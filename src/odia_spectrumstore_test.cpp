// Self-check for odia_spectrumstore.h.
//
// The claims worth testing are the ones that would fail SILENTLY: an m/z quantisation that is fine
// in the abstract and coarse relative to the 10 ppm extraction window, a binary search that returns
// the wrong first peak because it searches widened values instead of stored ones, and a summed
// range that quietly drops the boundary peaks extraction depends on.

#include "odia_spectrumstore.h"

#include <cmath>
#include <cstdio>
#include <random>
#include <vector>

using namespace odia;

static int g_fail = 0;
#define CHECK(c)                                                                       \
  do {                                                                                 \
    if (!(c)) { std::printf("  FAIL %s:%d: %s\n", __FILE__, __LINE__, #c); ++g_fail; } \
  } while (0)

int main()
{
  std::printf("odia_spectrumstore_test\n");

  // ---- 1. the precision claim, across the real m/z range ---------------------------------------
  // float32 must be far finer than the narrowest extraction window (10 ppm) everywhere the
  // instrument reports, or the storage becomes a term in the matching decision.
  std::printf("  %-12s %14s %14s\n", "m/z", "quant (ppm)", "vs 10 ppm win");
  for (const double mz : {200.0, 500.0, 1000.0, 2000.0})
  {
    const double ppm = SpectrumStore::quantisationPpm(mz);
    std::printf("  %-12.0f %14.4f %13.0fx finer\n", mz, ppm, 10.0 / ppm);
    CHECK(ppm < 0.1);                          // two orders below the window, everywhere
  }

  // ---- 2. a realistic spectrum round-trips within that precision --------------------------------
  std::mt19937 rng(11);
  std::uniform_real_distribution<double> um(150.0, 1800.0);
  std::uniform_real_distribution<float> ui(10.0F, 1.0e6F);
  const std::uint32_t n = 4000;
  std::vector<double> mz(n);
  std::vector<float> in(n);
  for (std::uint32_t k = 0; k < n; ++k) { mz[k] = um(rng); in[k] = ui(rng); }
  std::sort(mz.begin(), mz.end());

  SpectrumStore st;
  SpectrumStore::Meta m;
  m.rt = 1234.5; m.ms_level = 2; m.precursor_mz = 700.0; m.iso_lower = 699.0; m.iso_upper = 701.0;
  const std::size_t id = st.add(mz.data(), in.data(), n, m);
  CHECK(st.size() == 1 && st.count(id) == n);
  CHECK(st.meta(id).rt == 1234.5 && st.meta(id).ms_level == 2);

  std::vector<double> wm, wi;
  st.widen(id, wm, wi);
  CHECK(wm.size() == n && wi.size() == n);
  double worst_ppm = 0.0;
  for (std::uint32_t k = 0; k < n; ++k)
  {
    worst_ppm = std::max(worst_ppm, 1e6 * std::abs(wm[k] - mz[k]) / mz[k]);
    CHECK(wi[k] == in[k]);                     // intensities are float32 both sides: exact
  }
  std::printf("  worst m/z error over %u peaks: %.4f ppm\n", n, worst_ppm);
  CHECK(worst_ppm < 0.1);

  // ---- 3. ascending order is REQUIRED, not hoped for -------------------------------------------
  // Everything downstream binary-searches these arrays; an unsorted one returns wrong peaks rather
  // than failing, which is the worst possible outcome.
  {
    SpectrumStore bad;
    const double dm[3] = {500.0, 499.0, 501.0};
    const float di[3] = {1.0F, 1.0F, 1.0F};
    bool threw = false;
    try { bad.add(dm, di, 3, {}); } catch (const std::invalid_argument&) { threw = true; }
    CHECK(threw);
  }

  // ---- 4. lowerBound searches the STORED values -------------------------------------------------
  // Searching widened doubles while storing float32 puts the boundary in a different place than
  // sumRange uses, so a peak can be counted by one and missed by the other.
  for (std::uint32_t probe = 0; probe < n; probe += 373)
  {
    const std::uint32_t lb = st.lowerBound(id, mz[probe]);
    CHECK(lb <= probe);                                     // never past the peak itself
    CHECK(st.mz(id)[lb] >= static_cast<float>(mz[probe]) - 1e-3F);
  }

  // ---- 5. sumRange is inclusive at both ends ----------------------------------------------------
  // Extraction sums a +/- window around a fragment m/z; dropping the boundary peak is a silent
  // intensity loss that looks like a weaker peptide, not like a bug.
  {
    const std::uint32_t a = 1000, b = 1009;
    double want = 0.0;
    for (std::uint32_t k = a; k <= b; ++k) { want += in[k]; }
    const double got = st.sumRange(id, mz[a], mz[b]);
    // float32 m/z boundaries can admit an immediate neighbour; require the target peaks to be
    // present, not that no neighbour sneaks in.
    std::printf("  sumRange over peaks [%u,%u]: got %.1f, want %.1f\n", a, b, got, want);
    CHECK(got >= want * 0.999);
  }

  // ---- 6. many spectra: one arena, and the size claim measured ----------------------------------
  {
    SpectrumStore big;
    big.reserve(2000, 2000ull * 1200);
    std::vector<double> m2(1200);
    std::vector<float> i2(1200);
    for (int s = 0; s < 2000; ++s)
    {
      for (std::uint32_t k = 0; k < 1200; ++k) { m2[k] = um(rng); i2[k] = ui(rng); }
      std::sort(m2.begin(), m2.end());
      SpectrumStore::Meta mm;
      mm.rt = 0.6 * s; mm.ms_level = 2;
      big.add(m2.data(), i2.data(), 1200, mm);
    }
    const double ratio = double(big.bytesAsDoublePairs()) / double(big.bytes());
    std::printf("  2000 spectra x 1200 peaks: %.1f MB stored vs %.1f MB as double pairs -> %.2fx\n",
                big.bytes() / 1048576.0, big.bytesAsDoublePairs() / 1048576.0, ratio);
    CHECK(big.size() == 2000);
    CHECK(big.peakCount() == 2000ull * 1200);
    CHECK(ratio > 1.9);                       // 8 B/peak against 16 B/peak
    // Spectra stay independently addressable after all that appending.
    CHECK(big.count(1999) == 1200);
    CHECK(std::abs(big.meta(1999).rt - 0.6 * 1999) < 1e-9);
  }

  // ---- 7. an empty spectrum is legal and must not corrupt its neighbours ------------------------
  // Empty MS2 spectra occur in real runs; if one shifted the arena the whole store would be wrong
  // from that point on, which is exactly the failure that produces "0 peaks decoded" much later.
  {
    SpectrumStore s2;
    const double a[2] = {100.0, 200.0};
    const float b[2] = {5.0F, 6.0F};
    s2.add(a, b, 2, {});
    s2.add(nullptr, nullptr, 0, {});
    s2.add(a, b, 2, {});
    CHECK(s2.size() == 3);
    CHECK(s2.count(1) == 0);
    CHECK(s2.count(2) == 2);
    CHECK(s2.mz(2)[0] == 100.0F && s2.intensity(2)[1] == 6.0F);
    std::vector<double> e1, e2;
    s2.widen(1, e1, e2);
    CHECK(e1.empty() && e2.empty());
    CHECK(s2.sumRange(1, 0.0, 1e6) == 0.0);
  }

  if (g_fail) { std::printf("odia_spectrumstore_test FAILED (%d)\n", g_fail); return 1; }
  std::printf("odia_spectrumstore_test OK\n");
  return 0;
}
