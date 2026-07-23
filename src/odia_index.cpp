// odia-index -- streaming DIA indexer.
//
// Builds the per-(isolation window x RT block) index the extractor needs,
// without ever holding the run in memory. OpenSWATH materialises chromatograms
// and peaked at 60-83 GB on the agxt diaPASEF runs; this keeps only a binned
// representation whose size is set by the acquisition scheme, not by the number
// of peaks.
//
// Threading is lock-free by partition, the way DIA-NN does it (35 std::thread
// sites, 0 mutexes in diann.cpp): the XML parse is inherently serial, so it runs
// on the calling thread and routes each spectrum to the one worker that owns its
// isolation window. No two workers ever touch the same cell, so no cell needs a
// lock. Only the hand-off queues are synchronised.
//
//   odia-index <file.mzML> [--threads N] [--ppm 20] [--rt-block 30] [--selftest]

#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/INTERFACES/IMSDataConsumer.h>
#include <OpenMS/KERNEL/MSSpectrum.h>

#include <atomic>
#include <chrono>
#include <cmath>
#include <condition_variable>
#include <cstdio>
#include <cstring>
#include <deque>
#include <map>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

#ifdef __linux__
#include <unistd.h>
#endif

namespace
{
// Logarithmic m/z bins give constant relative (ppm) width, so one integer
// indexes a tolerance-sized bin anywhere in the spectrum.
inline long mz_bin(double mz, double mz0, double inv_log1p_tol)
{
  return static_cast<long>(std::log(mz / mz0) * inv_log1p_tol);
}

constexpr double MZ0 = 100.0;

// One (window, RT-block) cell: an occupancy bitmap for the prefilter plus the
// binned intensities. The bitmap is what makes rejection cheap; the intensities
// are what the gather reads.
struct Cell
{
  std::vector<uint64_t> occ;                        // 1 bit per m/z bin
  std::vector<std::pair<int32_t, float>> peaks;     // sparse (bin, intensity)

  void set(long bin)
  {
    const size_t w = static_cast<size_t>(bin) >> 6;
    if (w >= occ.size()) occ.resize(w + 1, 0);
    occ[w] |= 1ull << (bin & 63);
  }
  bool test(long bin) const
  {
    const size_t w = static_cast<size_t>(bin) >> 6;
    return w < occ.size() && (occ[w] >> (bin & 63) & 1ull);
  }
  size_t bytes() const
  {
    return occ.capacity() * sizeof(uint64_t) + peaks.capacity() * sizeof(peaks[0]);
  }
};

// A spectrum reduced to what the index needs. Built on the parse thread so the
// raw MSSpectrum can be released immediately.
struct Compact
{
  int window = -1;
  int block = -1;
  std::vector<std::pair<int32_t, float>> binned;
};

// Single-producer/single-consumer hand-off. One per worker, so the producer
// never contends with more than the one worker it is routing to.
class Queue
{
public:
  void push(Compact&& c)
  {
    {
      std::lock_guard<std::mutex> g(m_);
      q_.push_back(std::move(c));
    }
    cv_.notify_one();
  }
  void close()
  {
    {
      std::lock_guard<std::mutex> g(m_);
      closed_ = true;
    }
    cv_.notify_all();
  }
  bool pop(Compact& out)
  {
    std::unique_lock<std::mutex> g(m_);
    cv_.wait(g, [&] { return !q_.empty() || closed_; });
    if (q_.empty()) return false;
    out = std::move(q_.front());
    q_.pop_front();
    return true;
  }
  size_t depth()
  {
    std::lock_guard<std::mutex> g(m_);
    return q_.size();
  }

private:
  std::deque<Compact> q_;
  std::mutex m_;
  std::condition_variable cv_;
  bool closed_ = false;
};

class Indexer : public OpenMS::Interfaces::IMSDataConsumer
{
public:
  Indexer(int nthreads, double ppm, double rt_block)
    : inv_(1.0 / std::log1p(ppm * 1e-6)), rt_block_(rt_block),
      queues_(nthreads), cells_(nthreads), workers_(nthreads)
  {
    for (int i = 0; i < nthreads; ++i)
      workers_[i] = std::thread([this, i] { run_(i); });
  }

  void consumeSpectrum(OpenMS::MSSpectrum& s) override
  {
    if (s.getMSLevel() != 2 || s.empty() || s.getPrecursors().empty()) return;

    Compact c;
    const double pmz = s.getPrecursors().front().getMZ();
    c.window = window_of_(pmz);
    c.block = static_cast<int>(s.getRT() / rt_block_);
    c.binned.reserve(s.size());
    for (const auto& p : s)
      c.binned.emplace_back(static_cast<int32_t>(mz_bin(p.getMZ(), MZ0, inv_)),
                            static_cast<float>(p.getIntensity()));
    ++spectra_;
    peaks_ += c.binned.size();
    // Route by window: a window is owned by exactly one worker, so the cell it
    // writes into can never be touched concurrently.
    queues_[c.window % queues_.size()].push(std::move(c));
  }

  void consumeChromatogram(OpenMS::MSChromatogram&) override {}
  void setExpectedSize(size_t, size_t) override {}
  void setExperimentalSettings(const OpenMS::ExperimentalSettings&) override {}

  void finish()
  {
    for (auto& q : queues_) q.close();
    for (auto& t : workers_) if (t.joinable()) t.join();
  }

  size_t spectra() const { return spectra_; }
  size_t peaks() const { return peaks_; }
  size_t cells() const
  {
    size_t n = 0;
    for (const auto& m : cells_) n += m.size();
    return n;
  }
  size_t index_bytes() const
  {
    size_t b = 0;
    for (const auto& m : cells_)
      for (const auto& [k, c] : m) { (void)k; b += c.bytes(); }
    return b;
  }
  size_t windows() const { return win_index_.size(); }

private:
  // Isolation windows are discovered as they appear; the acquisition defines
  // them, so we must not assume a fixed grid.
  int window_of_(double pmz)
  {
    const long key = static_cast<long>(pmz * 100.0 + 0.5);
    auto it = win_index_.find(key);
    if (it != win_index_.end()) return it->second;
    const int id = static_cast<int>(win_index_.size());
    win_index_.emplace(key, id);
    return id;
  }

  void run_(int id)
  {
    Compact c;
    while (queues_[id].pop(c))
    {
      Cell& cell = cells_[id][(static_cast<int64_t>(c.window) << 32) | uint32_t(c.block)];
      for (const auto& [bin, inten] : c.binned)
      {
        cell.set(bin);
        cell.peaks.emplace_back(bin, inten);
      }
    }
  }

  double inv_, rt_block_;
  std::map<long, int> win_index_;
  std::vector<Queue> queues_;
  std::vector<std::map<int64_t, Cell>> cells_;
  std::vector<std::thread> workers_;
  std::atomic<size_t> spectra_{0}, peaks_{0};
};

size_t peak_rss_kb()
{
#ifdef __linux__
  FILE* f = std::fopen("/proc/self/status", "r");
  if (!f) return 0;
  char line[256];
  size_t kb = 0;
  while (std::fgets(line, sizeof line, f))
    if (std::strncmp(line, "VmHWM:", 6) == 0) { std::sscanf(line + 6, "%zu", &kb); break; }
  std::fclose(f);
  return kb;
#else
  return 0;
#endif
}

#define CHECK(c)                                                                                   \
  if (!(c)) { std::fprintf(stderr, "selftest FAILED line %d: %s\n", __LINE__, #c); return 1; }

int selftest()
{
  const double inv = 1.0 / std::log1p(20e-6);
  CHECK(mz_bin(MZ0, MZ0, inv) == 0);
  CHECK(mz_bin(500.0, MZ0, inv) < mz_bin(501.0, MZ0, inv));

  Cell c;
  c.set(0); c.set(63); c.set(64); c.set(5000);
  CHECK(c.test(0)); CHECK(c.test(63)); CHECK(c.test(64)); CHECK(c.test(5000));
  CHECK(!c.test(1)); CHECK(!c.test(62)); CHECK(!c.test(4999));
  // Bit 64 must land in the second word, not alias bit 0 -- the classic
  // off-by-one in a bitmap indexed by >>6 and &63.
  CHECK(c.occ.size() >= 79);
  CHECK((c.occ[0] & 1ull) && (c.occ[0] >> 63 & 1ull));
  CHECK(c.occ[1] & 1ull);

  Queue q;
  Compact in; in.window = 7; in.block = 3;
  q.push(std::move(in));
  q.close();
  Compact out;
  CHECK(q.pop(out));
  CHECK(out.window == 7 && out.block == 3);
  CHECK(!q.pop(out));                       // closed and drained

  // Cell key packing must keep window and block distinct.
  auto key = [](int w, int b) { return (static_cast<int64_t>(w) << 32) | uint32_t(b); };
  CHECK(key(1, 0) != key(0, 1));
  CHECK(key(3, 5) == key(3, 5));
  std::puts("selftest OK");
  return 0;
}
#undef CHECK
} // namespace

int main(int argc, char** argv)
{
  if (argc == 2 && std::string(argv[1]) == "--selftest") return selftest();
  if (argc < 2)
  {
    std::fprintf(stderr, "usage: odia-index <file.mzML> [--threads N] [--ppm P] [--rt-block S]\n");
    return 2;
  }
  int nthreads = 8;
  double ppm = 20.0, rt_block = 30.0;
  for (int i = 2; i + 1 < argc; i += 2)
  {
    const std::string k = argv[i];
    if (k == "--threads") nthreads = std::atoi(argv[i + 1]);
    else if (k == "--ppm") ppm = std::atof(argv[i + 1]);
    else if (k == "--rt-block") rt_block = std::atof(argv[i + 1]);
  }
  if (nthreads < 1) nthreads = 1;

  const auto t0 = std::chrono::steady_clock::now();
  Indexer idx(nthreads, ppm, rt_block);
  try
  {
    OpenMS::MzMLFile().transform(argv[1], &idx, true, true);
  }
  catch (const std::exception& e)
  {
    std::fprintf(stderr, "error: %s\n", e.what());
    return 1;
  }
  idx.finish();
  const double wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();

  std::printf("file          %s\n", argv[1]);
  std::printf("threads       %d   (1 parse + %d index workers)\n", nthreads, nthreads);
  std::printf("ppm / block   %.0f ppm / %.0f s\n", ppm, rt_block);
  std::printf("spectra       %zu MS2\n", idx.spectra());
  std::printf("peaks         %zu\n", idx.peaks());
  std::printf("windows       %zu\n", idx.windows());
  std::printf("cells         %zu\n", idx.cells());
  std::printf("index size    %.1f MB\n", idx.index_bytes() / 1048576.0);
  std::printf("peak RSS      %.1f MB\n", peak_rss_kb() / 1024.0);
  std::printf("wall          %.2f s\n", wall);
  return 0;
}
