// Self-check for odia_chromstore.h.
//
// The claims worth testing are the ones whose failure is SILENT: a slice that quietly returns the
// wrong points, and a quantisation error that is fine at the apex and ruinous at the baseline --
// which is where S/N is estimated.

#include "odia_chromstore.h"

#include <cmath>
#include <cstdio>
#include <random>
#include <vector>

using namespace odia;

static int g_fail = 0;
#define CHECK(cond)                                                                    \
  do {                                                                                 \
    if (!(cond)) { std::printf("  FAIL %s:%d: %s\n", __FILE__, __LINE__, #cond); ++g_fail; } \
  } while (0)

int main()
{
  std::printf("odia_chromstore_test\n");

  // A realistic window axis: ~1945 cycles over a 2333 s gradient.
  const std::size_t kCycles = 1945;
  std::vector<double> axis(kCycles);
  for (std::size_t i = 0; i < kCycles; ++i) { axis[i] = 12.0 + 1.1994 * static_cast<double>(i); }

  ChromStore store;
  const std::uint32_t ax = store.addAxis(axis);
  CHECK(store.axisCount() == 1);

  // ---- 1. a chromatogram with the benchmark's real dynamic range ------------------------------
  // Apex 1e6 on a 1e2 baseline: 4 orders. This is the case the header flags as the risk, so it is
  // the case the test uses rather than a friendly one.
  const std::uint32_t start = 300, len = 1196;
  std::vector<float> sig(len, 100.0F);
  for (std::uint32_t k = 0; k < len; ++k)
  {
    const double d = (static_cast<double>(k) - 600.0) / 4.0;         // ~9 s FWHM at 1.2 s/cycle
    sig[k] = static_cast<float>(100.0 + 1.0e6 * std::exp(-0.5 * d * d));
  }
  const std::size_t id = store.add(ax, start, sig);

  std::vector<double> rt;
  std::vector<float> got;
  store.get(id, rt, got);
  CHECK(rt.size() == len && got.size() == len);
  // The times must come from the AXIS, not be reconstructed by arithmetic.
  CHECK(std::abs(rt[0] - axis[start]) < 1e-12);
  CHECK(std::abs(rt[len - 1] - axis[start + len - 1]) < 1e-12);

  const double apex_err = static_cast<double>(std::abs(got[600] - sig[600]) / sig[600]);
  double base_err = 0.0;
  for (std::uint32_t k = 0; k < 50; ++k)                              // far from the peak
  {
    base_err = std::max(base_err, static_cast<double>(std::abs(got[k] - sig[k]) / sig[k]));
  }
  std::printf("  quantisation: apex rel err %.2e, baseline (1e2 on a 1e6 apex) rel err %.3f\n",
              apex_err, base_err);
  CHECK(apex_err < 1e-4);
  // The baseline error is LARGE by construction -- 1e6/65535 = 15.3 counts, so a 100-count baseline
  // resolves to ~7 levels. Asserted loosely and PRINTED, because the number is the point: anything
  // estimating noise from the baseline inherits it.
  CHECK(base_err < 0.20);

  // ---- 1b. ALL FOUR ENCODINGS on the same realistic signal ------------------------------------
  // The metric that actually matters is CORRELATION fidelity: the dominant sub-scores
  // (xcorr_shape, library_corr, dotprod, manhattan) are scale-invariant and read SHAPE, not
  // absolute counts. So a per-point relative error spread uniformly across the range costs almost
  // nothing, while an error concentrated at the apex -- which is what LINEAR quantisation gives --
  // destroys the baseline that S/N is estimated from.
  {
    auto corr = [](const std::vector<float>& a, const std::vector<float>& b) {
      const std::size_t n = std::min(a.size(), b.size());
      double ma = 0, mb = 0;
      for (std::size_t k = 0; k < n; ++k) { ma += a[k]; mb += b[k]; }
      ma /= static_cast<double>(n); mb /= static_cast<double>(n);
      double num = 0, da = 0, db = 0;
      for (std::size_t k = 0; k < n; ++k)
      {
        const double x = a[k] - ma, y = b[k] - mb;
        num += x * y; da += x * x; db += y * y;
      }
      return (da > 0 && db > 0) ? num / std::sqrt(da * db) : 0.0;
    };
    struct E { const char* name; ChromEncoding e; int bytes; };
    const E encs[] = {{"Quantised8Log ", ChromEncoding::Quantised8Log, 1},
                      {"Quantised8Lin ", ChromEncoding::Quantised8Lin, 1},
                      {"Quantised16   ", ChromEncoding::Quantised16, 2},
                      {"Float32       ", ChromEncoding::Float32, 4}};
    std::printf("  %-15s %6s %11s %13s %13s\n", "encoding", "B/pt", "apex relerr", "baseline rel",
                "1 - corr");
    for (const E& e : encs)
    {
      ChromStore st(e.e);
      const std::uint32_t a = st.addAxis(axis);
      const std::size_t i2 = st.add(a, start, sig);
      std::vector<double> r2; std::vector<float> g2;
      st.get(i2, r2, g2);
      double ap = static_cast<double>(std::abs(g2[600] - sig[600]) / sig[600]);
      double bs = 0.0;
      for (std::uint32_t k = 0; k < 50; ++k)
      { bs = std::max(bs, static_cast<double>(std::abs(g2[k] - sig[k]) / sig[k])); }
      std::printf("  %-15s %6d %11.2e %13.4f %13.2e\n", e.name, e.bytes, ap, bs,
                  1.0 - corr(sig, g2));
    }
    // Linear 8-bit cannot represent a 100-count baseline under a 1e6 apex at all: the code step is
    // 3922 counts, so it rounds to 0. That is the claim, so it is asserted rather than described.
    {
      ChromStore lin(ChromEncoding::Quantised8Lin);
      const std::uint32_t a = lin.addAxis(axis);
      const std::size_t i3 = lin.add(a, start, sig);
      std::vector<double> r3; std::vector<float> g3;
      lin.get(i3, r3, g3);
      CHECK(g3[0] == 0.0F);                      // baseline annihilated
    }
    // Log 8-bit keeps the baseline to a few percent AND correlates essentially perfectly.
    {
      ChromStore lg(ChromEncoding::Quantised8Log);
      const std::uint32_t a = lg.addAxis(axis);
      const std::size_t i4 = lg.add(a, start, sig);
      std::vector<double> r4; std::vector<float> g4;
      lg.get(i4, r4, g4);
      CHECK(g4[0] > 90.0F && g4[0] < 111.0F);    // baseline preserved within ~4%
      CHECK(1.0 - corr(sig, g4) < 1e-4);         // shape intact, which is what the scores read
    }
  }

  // ---- 2. slicing must return exactly the points a narrower window would --------------------
  // This is the operation the store exists for, and getting it wrong is invisible: a slice off by
  // one still looks like a chromatogram.
  const double lo = axis[start + 400], hi = axis[start + 799];
  std::vector<double> srt;
  std::vector<float> sint;
  const std::uint32_t n = store.getRtRange(id, lo, hi, srt, sint);
  std::printf("  RT slice [%.1f, %.1f] -> %u points (expected 400)\n", lo, hi, n);
  CHECK(n == 400);
  CHECK(std::abs(srt.front() - lo) < 1e-9);
  CHECK(std::abs(srt.back() - hi) < 1e-9);
  for (std::uint32_t k = 0; k < n; ++k) { CHECK(std::abs(srt[k] - axis[start + 400 + k]) < 1e-9); }
  // A slice must equal the corresponding part of the whole -- decoding a sub-range must not shift.
  double slice_diff = 0.0;
  for (std::uint32_t k = 0; k < n; ++k)
  {
    slice_diff = std::max(slice_diff, std::abs(static_cast<double>(sint[k]) - got[400 + k]));
  }
  CHECK(slice_diff == 0.0);

  // ---- 3. a request beyond the covered range clamps, it does not read off the end -------------
  std::uint32_t n2 = store.getRtRange(id, axis[0], axis[kCycles - 1], srt, sint);
  CHECK(n2 == len);                                                   // clamped to what it covers
  n2 = store.getRtRange(id, axis[0], axis[1], srt, sint);
  CHECK(n2 == 0);                                                     // entirely before the slice
  CHECK(srt.empty() && sint.empty());

  // ---- 4. an all-zero chromatogram is legitimate and must not become NaN ----------------------
  // A transition with no signal in the window is common; dividing by a zero scale would poison it.
  const std::size_t zid = store.add(ax, 0, std::vector<float>(100, 0.0F));
  store.get(zid, rt, got);
  bool finite = true;
  for (const float v : got) { finite = finite && std::isfinite(v) && v == 0.0F; }
  CHECK(finite);

  // ---- 5. the size claim, measured rather than asserted ---------------------------------------
  {
    ChromStore big;
    const std::uint32_t a2 = big.addAxis(axis);
    std::mt19937 rng(4);
    std::uniform_real_distribution<float> u(0.0F, 1.0e5F);
    for (int c = 0; c < 2000; ++c)
    {
      std::vector<float> v(len);
      for (auto& x : v) { x = u(rng); }
      big.add(a2, start, v);
    }
    const double ratio = static_cast<double>(big.bytesAsPeaks()) / static_cast<double>(big.bytes());
    std::printf("  2000 chromatograms x %u points: %.1f MB stored vs %.1f MB as ChromatogramPeak "
                "-> %.1fx\n", len, big.bytes() / 1048576.0, big.bytesAsPeaks() / 1048576.0, ratio);
    CHECK(ratio > 7.0);                                               // header claims 7.9x
  }

  // ---- 6. Float32 encoding is exact, so a quantisation suspicion can be ruled out -------------
  {
    ChromStore f(ChromEncoding::Float32);
    const std::uint32_t a3 = f.addAxis(axis);
    const std::size_t fid = f.add(a3, start, sig);
    CHECK(f.quantisationError(fid, sig) == 0.0);
    std::printf("  Float32 encoding: exact (rel err 0)\n");
  }

  if (g_fail) { std::printf("odia_chromstore_test FAILED (%d)\n", g_fail); return 1; }
  std::printf("odia_chromstore_test OK\n");
  return 0;
}
