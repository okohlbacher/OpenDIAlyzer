#pragma once
// One mmap for the whole .mzpeak, shared by every worker thread.
//
// WHY THIS EXISTS. mzPeak's own archive opens the file again for every member it hands out, and
// MzPeak::Index::spectra() gives each thread a Spectra holding its own Data::Signals. At 224
// workers that measured 908 open descriptors ON THE SAME FILE, blew past the default 1024 limit,
// and aborted:
//
//     what():  failed to open zip archive Can't open file: Too many open files
//
// Raising ulimit let it finish and showed the real cost: peak RSS 127.6 GB against 74.5 GB for the
// SAME analysis reading mzML. The columnar input -- the one that should be cheap -- was 71% more
// expensive than XML purely because of how it was opened.
//
// Centralising the Index alone was NOT enough and it is worth recording why, because the first
// attempt did exactly that and the descriptor count barely moved (1,060 -> 911). The Index is not
// where the handles come from; the per-thread Spectra are. What has to be shared is the BYTES.
//
// So: map the file once, parse the ZIP central directory once, and hand out Files that are plain
// spans of that mapping. Threads then share one read-only mapping with no locking (the page cache
// is doing the concurrency), the kernel pages in only what is touched, and identical members read
// by different workers cost one page each rather than one copy each.
//
// The descriptor is CLOSED after mmap. A mapping keeps the file alive on its own, so steady-state
// cost is zero descriptors, not one.
//
// FORMAT ASSUMPTION, checked rather than trusted: entries must be STORED (compression method 0).
// mzPeak writes parquet members uncompressed precisely so they can be read in place, but a
// DEFLATE-compressed member cannot be a span of the mapping, and silently returning its compressed
// bytes would surface far away as a parquet parse error. Such a member is rejected by name.

#include <mzpeak/io/archive.h>
#include <mzpeak/io/file.h>

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

namespace odia
{

/// The mapping itself, refcounted so every Archive and File keeps it alive.
class MappedFile
{
public:
  explicit MappedFile(const std::string& path)
  {
    const int fd = ::open(path.c_str(), O_RDONLY);
    if (fd < 0) { throw std::runtime_error("odia::MappedFile: cannot open " + path); }
    struct stat st {};
    if (::fstat(fd, &st) != 0) { ::close(fd); throw std::runtime_error("odia::MappedFile: fstat failed"); }
    size_ = static_cast<std::size_t>(st.st_size);
    void* p = ::mmap(nullptr, size_, PROT_READ, MAP_SHARED, fd, 0);
    // The mapping holds its own reference to the file, so the descriptor is released immediately.
    // This is the whole point: steady state is ZERO descriptors, however many threads read.
    ::close(fd);
    if (p == MAP_FAILED) { throw std::runtime_error("odia::MappedFile: mmap failed for " + path); }
    base_ = static_cast<const std::uint8_t*>(p);
  }
  ~MappedFile() { if (base_) { ::munmap(const_cast<std::uint8_t*>(base_), size_); } }
  MappedFile(const MappedFile&) = delete;
  MappedFile& operator=(const MappedFile&) = delete;

  const std::uint8_t* data() const { return base_; }
  std::size_t size() const { return size_; }

private:
  const std::uint8_t* base_ = nullptr;
  std::size_t size_ = 0;
};

/// A member of the archive, as a window on the shared mapping. Holds no descriptor and copies
/// nothing; read() memcpy's out of the page cache.
class SpanFile final : public MzPeak::IO::File
{
public:
  SpanFile(std::shared_ptr<const MappedFile> map, std::string name,
           std::size_t offset, std::size_t length)
    : map_(std::move(map)), name_(std::move(name)), off_(offset), len_(length) {}

  std::string name() const override { return name_; }
  std::size_t size() const override { return len_; }
  std::optional<std::size_t> tell() const override
  {
    if (!open_) { return std::nullopt; }
    return pos_;
  }
  /// EOF is an EMPTY optional, not 0.
  ///
  /// This mirrors ZipBuffer::read (mzpeak/src/zip_buffer.cpp:39), which returns {} whenever
  /// zip_fread yields <= 0, and the consumer relies on it: index.cpp:94 loops while the optional is
  /// engaged. Returning optional{0} at EOF -- the "natural" contract -- is always truthy and hangs
  /// the reader forever. It did: the access test went from 1.03 s to a timeout, with no output,
  /// because the very first metadata read never terminated.
  std::optional<std::size_t> read(std::uint8_t* buffer, std::size_t n) override
  {
    if (!open_ || buffer == nullptr || n == 0) { return std::nullopt; }
    const std::size_t avail = len_ - std::min(pos_, len_);
    const std::size_t take = std::min(n, avail);
    if (take == 0) { return std::nullopt; }
    std::memcpy(buffer, map_->data() + off_ + pos_, take);
    pos_ += take;
    return take;
  }
  bool seek(std::size_t pos) override
  {
    if (!open_ || pos > len_) { return false; }
    pos_ = pos;
    return true;
  }
  void close() override { open_ = false; }
  bool is_open() const override { return open_; }

private:
  std::shared_ptr<const MappedFile> map_;
  std::string name_;
  std::size_t off_ = 0, len_ = 0, pos_ = 0;
  bool open_ = true;
};

/// Archive over one shared mapping. Cheap to construct -- the directory is parsed once and shared,
/// so a per-thread Archive costs a shared_ptr copy, not a file open and not a directory scan.
class MmapArchive final : public MzPeak::IO::Archive
{
public:
  struct Entry { std::size_t offset, length; };
  using Directory = std::unordered_map<std::string, Entry>;

  /// Map @p path and parse its central directory. Use the shared_ptr overload for extra Archives
  /// over the same file so the mapping and the directory are not duplicated.
  explicit MmapArchive(const std::string& path)
    : map_(std::make_shared<MappedFile>(path)), dir_(std::make_shared<Directory>(parse_(*map_))) {}

  MmapArchive(std::shared_ptr<const MappedFile> map, std::shared_ptr<const Directory> dir)
    : map_(std::move(map)), dir_(std::move(dir)) {}

  std::vector<MzPeak::IO::fs::path> list() override
  {
    std::vector<MzPeak::IO::fs::path> out;
    out.reserve(dir_->size());
    for (const auto& kv : *dir_) { out.emplace_back(kv.first); }
    std::sort(out.begin(), out.end());          // deterministic order, whatever the hash layout
    return out;
  }

  std::unique_ptr<MzPeak::IO::File> read_file(const MzPeak::IO::fs::path& p) override
  {
    const auto it = dir_->find(p.string());
    if (it == dir_->end()) { return nullptr; }
    return std::make_unique<SpanFile>(map_, p.string(), it->second.offset, it->second.length);
  }

  std::shared_ptr<const MappedFile> mapping() const { return map_; }
  std::shared_ptr<const Directory> directory() const { return dir_; }

private:
  static std::uint32_t u32(const std::uint8_t* p) { return std::uint32_t(p[0]) | (std::uint32_t(p[1]) << 8) | (std::uint32_t(p[2]) << 16) | (std::uint32_t(p[3]) << 24); }
  static std::uint16_t u16(const std::uint8_t* p) { return std::uint16_t(std::uint16_t(p[0]) | (std::uint16_t(p[1]) << 8)); }

  /// Locate the end-of-central-directory record, walk the central directory, and record each
  /// member's DATA offset (past its local header, whose extra field may differ from the central
  /// one -- reading the local header is the only correct way to find where the bytes start).
  static Directory parse_(const MappedFile& m)
  {
    const std::uint8_t* b = m.data();
    const std::size_t n = m.size();
    if (n < 22) { throw std::runtime_error("odia::MmapArchive: file too small to be a zip"); }
    std::size_t eocd = 0;
    bool found = false;
    const std::size_t back = std::min<std::size_t>(n, 66000);      // EOCD + max comment
    for (std::size_t i = n - 22; i + 1 > n - back; --i)
    {
      if (u32(b + i) == 0x06054b50) { eocd = i; found = true; break; }
      if (i == 0) { break; }
    }
    if (!found) { throw std::runtime_error("odia::MmapArchive: no end-of-central-directory record"); }

    std::uint64_t count = u16(b + eocd + 10);
    std::uint64_t cd_off = u32(b + eocd + 16);
    // ZIP64: a >4 GB archive stores the real values in the ZIP64 EOCD, and the 32-bit fields are
    // sentinels. mzPeak files exceed 4 GB routinely, so this is the normal path, not an edge case.
    if (cd_off == 0xFFFFFFFFu || count == 0xFFFFu)
    {
      if (eocd < 20) { throw std::runtime_error("odia::MmapArchive: truncated zip64 locator"); }
      const std::size_t loc = eocd - 20;
      if (u32(b + loc) != 0x07064b50) { throw std::runtime_error("odia::MmapArchive: missing zip64 locator"); }
      const std::uint64_t z64 = *reinterpret_cast<const std::uint64_t*>(b + loc + 8);
      if (z64 + 56 > n || u32(b + z64) != 0x06064b50) { throw std::runtime_error("odia::MmapArchive: bad zip64 EOCD"); }
      count = *reinterpret_cast<const std::uint64_t*>(b + z64 + 32);
      cd_off = *reinterpret_cast<const std::uint64_t*>(b + z64 + 48);
    }

    Directory dir;
    dir.reserve(static_cast<std::size_t>(count) * 2);
    std::size_t p = static_cast<std::size_t>(cd_off);
    for (std::uint64_t k = 0; k < count && p + 46 <= n; ++k)
    {
      if (u32(b + p) != 0x02014b50) { break; }
      const std::uint16_t method = u16(b + p + 10);
      const std::uint16_t nlen = u16(b + p + 28);
      const std::uint16_t elen = u16(b + p + 30);
      const std::uint16_t clen = u16(b + p + 32);
      std::uint64_t usize = u32(b + p + 24);
      std::uint64_t lho = u32(b + p + 42);
      const std::string name(reinterpret_cast<const char*>(b + p + 46), nlen);

      // ZIP64 extra field: sizes and the local-header offset move here once they exceed 32 bits.
      if (usize == 0xFFFFFFFFu || lho == 0xFFFFFFFFu)
      {
        const std::uint8_t* ex = b + p + 46 + nlen;
        const std::uint8_t* ex_end = ex + elen;
        while (ex + 4 <= ex_end)
        {
          const std::uint16_t tag = u16(ex), sz = u16(ex + 2);
          if (tag == 0x0001)
          {
            const std::uint8_t* q = ex + 4;
            if (usize == 0xFFFFFFFFu && q + 8 <= ex_end) { usize = *reinterpret_cast<const std::uint64_t*>(q); q += 8; }
            if (u32(b + p + 20) == 0xFFFFFFFFu && q + 8 <= ex_end) { q += 8; }   // compressed size
            if (lho == 0xFFFFFFFFu && q + 8 <= ex_end) { lho = *reinterpret_cast<const std::uint64_t*>(q); }
            break;
          }
          ex += 4 + sz;
        }
      }

      if (method != 0)
      {
        // Not a defect in this class: a DEFLATE member simply cannot be a span of the mapping.
        // Named, because returning its compressed bytes would surface much later as a parquet
        // parse error with nothing pointing back here.
        throw std::runtime_error("odia::MmapArchive: member '" + name +
                                 "' is compressed (method " + std::to_string(method) +
                                 "); mzPeak members must be STORED to be read in place");
      }
      if (lho + 30 <= n)
      {
        // The LOCAL header's name/extra lengths are authoritative for where the data begins; they
        // are permitted to differ from the central directory's.
        const std::uint16_t lnlen = u16(b + lho + 26);
        const std::uint16_t lelen = u16(b + lho + 28);
        const std::size_t data = static_cast<std::size_t>(lho) + 30 + lnlen + lelen;
        if (data + usize <= n) { dir.emplace(name, Entry{data, static_cast<std::size_t>(usize)}); }
      }
      p += 46 + nlen + elen + clen;
    }
    if (dir.empty()) { throw std::runtime_error("odia::MmapArchive: central directory listed no readable members"); }
    return dir;
  }

  std::shared_ptr<const MappedFile> map_;
  std::shared_ptr<const Directory> dir_;
};

} // namespace odia
