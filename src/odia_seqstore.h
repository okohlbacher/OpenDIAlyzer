// odia_seqstore.h — interned string storage: one blob, 4-byte handles, no duplication.
//
// WHY. Measured on the Astral benchmark library (2026-08-01, phase-resolved probe):
//
//     7,149,966 compounds + 78,569,077 transitions  ->  38.79 GB RSS for library_load
//     of which  in_use 7.92 GB | RETAINED 11.29 GB | mmapped 18.94 GB
//
// 11.29 GB — 29% — is memory glibc holds after freeing, because ~471 million small string
// allocations fragment the arena around the live ones. The allocations themselves come from a
// layout that stores every string separately:
//
//   * `peptide_ref` is 78.6M copies of ~7.1M DISTINCT values — an 11x duplication. The parquet
//     source stored it dictionary-encoded, i.e. once; materialising into std::string exploded it.
//   * protein accessions: ~20k distinct, referenced from millions of precursors.
//   * every std::string over 15 chars (libstdc++ SSO limit) is its own malloc. "DECOY_PEPTIDEK_2"
//     is 16. So essentially all of them are.
//
// A std::string member costs 32 B inline PLUS a heap block PLUS allocator bookkeeping PLUS the
// fragmentation it leaves behind. An Id here costs 4 B and no allocation.
//
// DESIGN — why chunks rather than one growing buffer.
// The obvious implementation is a single std::string blob with string_views into it. That is a
// dangling-pointer trap: any growth reallocates the blob and every previously handed-out view
// becomes invalid, silently, long after the fact. Chunks of fixed size are never reallocated, so a
// view stays valid for the store's lifetime no matter how much is appended afterwards. A string
// that would straddle a chunk boundary starts a new chunk instead of being split.
//
// NOT a general string class. It is append-only and never erases: a proteomics library is built
// once and read many times, and supporting erase would cost either reference counts or compaction,
// both of which buy nothing here.

#ifndef ODIA_SEQSTORE_H
#define ODIA_SEQSTORE_H

#include <cstddef>
#include <cstdint>
#include <cstring>
#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>
#include <algorithm>
#include <unordered_map>
#include <vector>

namespace odia
{

class SequenceStore
{
public:
  using Id = std::uint32_t;
  /// Returned by find() for an absent string, and usable as "no value" in callers' records.
  static constexpr Id npos = ~Id(0);

  explicit SequenceStore(std::size_t chunk_bytes = 8u << 20) : chunk_bytes_(chunk_bytes)
  {
    if (chunk_bytes_ < 1024) { throw std::invalid_argument("SequenceStore: chunk too small"); }
  }

  /// Store `s` once. Repeated calls with equal content return the SAME Id and store nothing new —
  /// which is the entire point for peptide_ref and protein accessions.
  Id intern(std::string_view s)
  {
    const auto it = index_.find(s);
    if (it != index_.end()) { return it->second; }
    const Id id = appendRaw_(s);
    // Key the map with the STORED view, not the caller's: the caller's buffer may be a temporary,
    // and the stored one is stable for the store's lifetime by construction.
    index_.emplace(view(id), id);
    return id;
  }

  /// Store without deduplicating. For values unique by construction (e.g. transition names), where
  /// the hash lookup would cost more than it saves.
  Id append(std::string_view s) { return appendRaw_(s); }

  /// Id of `s` if present, else npos. Does not store.
  Id find(std::string_view s) const
  {
    const auto it = index_.find(s);
    return it == index_.end() ? npos : it->second;
  }

  /// Stable for the lifetime of the store — chunks are never reallocated.
  std::string_view view(Id id) const
  {
    if (id >= entries_.size()) { throw std::out_of_range("SequenceStore: bad Id"); }
    const Entry& e = entries_[id];
    return std::string_view(chunks_[e.chunk].get() + e.offset, e.length);
  }

  std::size_t size() const { return entries_.size(); }
  bool empty() const { return entries_.empty(); }

  /// Actual resident bytes: chunk capacity (not bytes used), the entry table, and the dedup index.
  /// Capacity is what the process holds, which is the number that matters for a memory comparison.
  std::size_t bytes() const
  {
    return sizeof(*this)
         + chunks_.size() * chunk_bytes_
         + chunks_.capacity() * sizeof(std::unique_ptr<char[]>)
         + entries_.capacity() * sizeof(Entry)
         + index_.size() * (sizeof(std::string_view) + sizeof(Id) + 16);   // node + bucket estimate
  }

  /// What the same strings cost as individual std::string objects, at the given reference count.
  /// `n_refs` is the number of REFERENCES, not distinct values -- that difference is the saving.
  std::size_t bytesAsStdString(std::size_t n_refs) const
  {
    std::size_t heap = 0;
    for (const Entry& e : entries_)
    {
      // libstdc++: <=15 chars live inline; longer ones take a heap block, rounded to a 16-byte
      // malloc bucket with a 16-byte header.
      if (e.length > 15) { heap += ((e.length + 1 + 15) / 16) * 16 + 16; }
    }
    // Each REFERENCE costs a 32-byte std::string; each DISTINCT value costs its heap block, and a
    // duplicated value costs its heap block once per reference.
    const double dup = entries_.empty() ? 1.0 : double(n_refs) / double(entries_.size());
    return n_refs * sizeof(std::string) + std::size_t(heap * dup);
  }


  /// A reference to a SPAN of a stored sequence: 12 bytes, no characters of its own.
  ///
  /// The point, for a proteomics library: a tryptic peptide is a SUBSTRING of its protein. Storing
  /// the ~20k FASTA protein sequences once (~11 MB for the human proteome) lets every one of the
  /// millions of peptides be (protein, offset, length) instead of its own std::string with its own
  /// heap block. Peptides that are NOT substrings -- variants, semi-tryptic products crossing a
  /// modification, anything synthesised -- are appended to the store and get a span covering their
  /// own entry, so callers handle one uniform type either way.
  struct Span
  {
    Id seq = npos;              ///< which stored sequence
    std::uint32_t offset = 0;   ///< start within it
    std::uint32_t length = 0;   ///< characters
    bool valid() const { return seq != npos; }
  };

  std::string_view view(Span sp) const
  {
    if (!sp.valid()) { return {}; }
    const std::string_view whole = view(sp.seq);
    if (std::size_t(sp.offset) + sp.length > whole.size())
    {
      throw std::out_of_range("SequenceStore: span leaves its sequence");
    }
    return whole.substr(sp.offset, sp.length);
  }

  /// Span for `needle` inside stored sequence `seq`, or an invalid span if it does not occur.
  /// The caller decides what to do about a miss -- this never guesses.
  Span locate(Id seq, std::string_view needle) const
  {
    if (needle.empty()) { return Span{}; }
    const std::string_view hay = view(seq);
    const std::size_t at = hay.find(needle);
    if (at == std::string_view::npos) { return Span{}; }
    return Span{seq, std::uint32_t(at), std::uint32_t(needle.size())};
  }


  /// Index every k-mer position of every stored sequence, so a peptide can be located anywhere in
  /// the stored set without knowing which sequence it came from.
  ///
  /// STRUCTURE, and why not a hash map. The first version used
  /// `unordered_multimap<uint64_t, Pos>`, which makes ONE NODE ALLOCATION PER RESIDUE -- ~11.4M for
  /// the human proteome -- and turns every probe into a pointer chase. That is the exact pathology
  /// this class exists to remove, reproduced inside it. This is a single sorted vector: one
  /// allocation, binary search, and the candidates for a key land in contiguous memory.
  ///
  /// k = 5 is below the length of essentially every tryptic peptide, so the short-needle fallback
  /// is a safety net rather than a code path. With 20 residues that is 3.2M possible keys over
  /// ~11.4M positions, about 3-4 candidates per probe, each rejected or confirmed by one memcmp.
  ///
  /// Call it again after adding more sequences (e.g. interned decoys) -- it rebuilds from scratch,
  /// and `indexAppended()` below extends it incrementally instead.
  void buildIndex(std::size_t k = 5)
  {
    k_ = k;
    kmers_.clear();
    indexed_upto_ = 0;
    indexAppended();
  }

  /// Index sequences added since the last call. Lets decoys -- which are only discovered as the
  /// library is read -- become searchable without rebuilding the whole index.
  void indexAppended()
  {
    if (k_ == 0) { k_ = 5; }
    std::size_t add = 0;
    for (Id id = Id(indexed_upto_); id < Id(entries_.size()); ++id)
    {
      const std::size_t n = entries_[id].length;
      add += n >= k_ ? n - k_ + 1 : 0;
    }
    if (add == 0 && indexed_upto_ == entries_.size()) { return; }
    kmers_.reserve(kmers_.size() + add);
    const std::size_t first_new = kmers_.size();
    for (Id id = Id(indexed_upto_); id < Id(entries_.size()); ++id)
    {
      const std::string_view sv = view(id);
      if (sv.size() < k_) { continue; }
      for (std::size_t off = 0; off + k_ <= sv.size(); ++off)
      {
        kmers_.push_back(KMer{kmerHash_(sv.substr(off, k_)), id, std::uint32_t(off)});
      }
    }
    // Sort only the new tail, then merge -- rebuilding from scratch on every incremental call
    // would be quadratic across a library read.
    std::sort(kmers_.begin() + std::ptrdiff_t(first_new), kmers_.end(),
              [](const KMer& a, const KMer& b) { return a.key < b.key; });
    if (first_new > 0)
    {
      std::inplace_merge(kmers_.begin(), kmers_.begin() + std::ptrdiff_t(first_new), kmers_.end(),
                         [](const KMer& a, const KMer& b) { return a.key < b.key; });
    }
    indexed_upto_ = entries_.size();
  }

  bool indexed() const { return !kmers_.empty(); }
  std::size_t indexBytes() const { return kmers_.capacity() * sizeof(KMer); }

  /// Find `needle` anywhere in the indexed sequences, or an invalid span.
  /// Candidates are verified by full comparison, so a hash collision costs a failed memcmp and
  /// never a wrong answer.
  Span findAnywhere(std::string_view needle) const
  {
    if (needle.empty() || kmers_.empty()) { return {}; }
    if (needle.size() < k_)
    {
      // Safety net only: with k=5 this is essentially unreachable for tryptic peptides. Returning
      // "not found" here instead would be a SILENT false negative -- measured, an earlier k=8 sent
      // 127,853 real peptides to standalone storage purely for being short.
      for (Id id = 0; id < Id(entries_.size()); ++id)
      {
        const std::string_view hay = view(id);
        const std::size_t at = hay.find(needle);
        if (at != std::string_view::npos) { return Span{id, std::uint32_t(at), std::uint32_t(needle.size())}; }
      }
      return {};
    }
    const std::uint64_t key = kmerHash_(needle.substr(0, k_));
    auto lo = std::lower_bound(kmers_.begin(), kmers_.end(), key,
                               [](const KMer& a, std::uint64_t v) { return a.key < v; });
    for (; lo != kmers_.end() && lo->key == key; ++lo)
    {
      const std::string_view hay = view(lo->seq);
      if (std::size_t(lo->offset) + needle.size() > hay.size()) { continue; }
      if (hay.compare(lo->offset, needle.size(), needle) == 0)
      {
        return Span{lo->seq, lo->offset, std::uint32_t(needle.size())};
      }
    }
    return {};
  }

  /// Store `s` in its own right and return a span covering it. For sequences that are not a
  /// substring of anything already stored.
  Span internAsSpan(std::string_view s)
  {
    const Id id = intern(s);
    return Span{id, 0u, std::uint32_t(s.size())};
  }

  /// The whole rule in one call: prefer a substring of `parent`, fall back to storing it.
  /// `out_was_substring` reports which happened, so a caller can log how much of a library is
  /// genuinely derived from its FASTA -- a number worth knowing rather than assuming.
  Span spanOfOrIntern(Id parent, std::string_view s, bool* out_was_substring = nullptr)
  {
    if (parent != npos)
    {
      const Span sp = locate(parent, s);
      if (sp.valid())
      {
        if (out_was_substring) { *out_was_substring = true; }
        return sp;
      }
    }
    if (out_was_substring) { *out_was_substring = false; }
    return internAsSpan(s);
  }

  /// Bytes a caller pays per peptide with spans, versus as std::string.
  static constexpr std::size_t bytesPerSpan() { return sizeof(Span); }

  /// Bytes a caller pays per reference with this store: one Id.
  static constexpr std::size_t bytesPerRef() { return sizeof(Id); }

  /// Pre-size for a known workload. Parquet's footer carries both numbers exactly (row count and
  /// total uncompressed column bytes), so the caller never has to guess.
  void reserve(std::size_t n_strings, std::size_t total_chars)
  {
    entries_.reserve(n_strings);
    index_.reserve(n_strings);
    chunks_.reserve(total_chars / chunk_bytes_ + 2);
  }

  /// Release the dedup index once loading is done. It exists only to answer intern(); reads use
  /// view(). On a 7.1M-entry store that is hundreds of MB returned for nothing given up.
  void freezeIndex()
  {
    std::unordered_map<std::string_view, Id> empty;
    index_.swap(empty);
  }
  bool frozen() const { return index_.empty() && !entries_.empty(); }

private:
  struct KMer { std::uint64_t key; Id seq; std::uint32_t offset; };   // 16 B, contiguous
  static std::uint64_t kmerHash_(std::string_view s)
  {
    std::uint64_t h = 1469598103934665603ull;                  // FNV-1a
    for (char c : s) { h = (h ^ std::uint8_t(c)) * 1099511628211ull; }
    return h;
  }

  struct Entry
  {
    std::uint32_t chunk;    ///< which chunk holds it
    std::uint32_t offset;   ///< byte offset within that chunk
    std::uint32_t length;   ///< bytes (not NUL-terminated)
  };

  Id appendRaw_(std::string_view s)
  {
    if (s.size() > chunk_bytes_)
    {
      // A single string larger than a chunk gets its own oversized chunk rather than being split:
      // splitting would break the contiguity that makes view() a plain string_view.
      chunks_.push_back(std::unique_ptr<char[]>(new char[s.size()]));
      std::memcpy(chunks_.back().get(), s.data(), s.size());
      entries_.push_back(Entry{std::uint32_t(chunks_.size() - 1), 0u, std::uint32_t(s.size())});
      used_ = chunk_bytes_;                       // force a fresh chunk for the next append
      return Id(entries_.size() - 1);
    }
    if (chunks_.empty() || used_ + s.size() > chunk_bytes_)
    {
      chunks_.push_back(std::unique_ptr<char[]>(new char[chunk_bytes_]));
      used_ = 0;
    }
    std::memcpy(chunks_.back().get() + used_, s.data(), s.size());
    entries_.push_back(Entry{std::uint32_t(chunks_.size() - 1), std::uint32_t(used_), std::uint32_t(s.size())});
    used_ += s.size();
    if (entries_.size() - 1 >= npos) { throw std::length_error("SequenceStore: more than 2^32-1 strings"); }
    return Id(entries_.size() - 1);
  }

  std::size_t chunk_bytes_;
  std::size_t used_ = 0;
  std::vector<std::unique_ptr<char[]>> chunks_;
  std::vector<Entry> entries_;
  std::unordered_map<std::string_view, Id> index_;
  std::size_t k_ = 5;
  std::vector<KMer> kmers_;          ///< ONE allocation, sorted by key
  std::size_t indexed_upto_ = 0;
};

} // namespace odia

#endif // ODIA_SEQSTORE_H
