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

  // ---- 8. a drift-less spectrum must NOT inherit the arena's back-filled -1 -------------------
  // The arena `dt_` is shared and gets back-filled with -1 as soon as ANY spectrum supplies
  // mobility. If drift() keyed on "the arena is non-empty", every drift-less spectrum in the same
  // store would report a full array of -1 -- which passes a null check and then fails every real
  // 1/K0 band (0.6-1.4), excluding all of its peaks from every window. Silently: the run finishes
  // and those scans simply contribute nothing. So drift() must key on the SPECTRUM.
  {
    SpectrumStore s3;
    const double a[3] = {100.0, 200.0, 300.0};
    const float b[3] = {1.0F, 2.0F, 3.0F};
    const double d[3] = {0.9, 1.0, 1.1};

    SpectrumStore::Meta m0;
    m0.drift = 0.95;                                  // per-slice layout: a SCALAR mobility
    const std::size_t no_im = s3.add(a, b, 3, m0);    // no per-peak array
    const std::size_t with  = s3.add(a, b, 3, {}, d); // per-peak array -> back-fills the arena
    const std::size_t after = s3.add(a, b, 3, {});    // drift-less again, arena now non-empty

    CHECK(s3.hasDrift());                             // the arena does hold mobility...
    CHECK(s3.drift(no_im) == nullptr);                // ...but these two spectra do not
    CHECK(s3.drift(after) == nullptr);
    CHECK(s3.drift(with) != nullptr);
    CHECK(s3.drift(with)[0] == 0.9F && s3.drift(with)[2] == 1.1F);
    CHECK(s3.meta(no_im).per_peak_drift == false && s3.meta(with).per_peak_drift == true);

    // A banded consumer keeps every peak of a drift-less spectrum, exactly as the raw path does.
    int kept = 0;
    const float* pd = s3.drift(after);
    for (std::uint32_t k = 0; k < s3.count(after); ++k)
    {
      if (pd && !(pd[k] >= 0.6F && pd[k] < 1.4F)) { continue; }
      ++kept;
    }
    CHECK(kept == 3);

    // widen() falls back to the SCALAR, not to a hard -1: that is what the raw path emits.
    std::vector<double> wm2, wi2, wd2;
    s3.widen(no_im, wm2, wi2, &wd2);
    CHECK(wd2.size() == 3 && wd2[0] == 0.95 && wd2[2] == 0.95);
    s3.widen(with, wm2, wi2, &wd2);
    CHECK(std::abs(wd2[1] - 1.0) < 1e-6);
  }

  // ---- 9. the window edge is not narrowed to float ---------------------------------------------
  // Narrowing the edge rounds it to nearest, which can move it INWARD by half a float ulp and drop
  // a peak that is genuinely inside the window. Widening the stored value instead is exact.
  {
    SpectrumStore s4;
    const double one = 1000.0009999;                  // sits between two float32 neighbours
    const float iv = 42.0F;
    s4.add(&one, &iv, 1, {});
    const float stored = s4.mz(0)[0];
    const double lo = std::nextafter(double(stored), 0.0);          // just below the stored value
    const double hi = std::nextafter(double(stored), 1e9);          // just above
    CHECK(s4.lowerBound(0, lo) == 0);                 // the peak is found...
    CHECK(s4.sumRange(0, lo, hi) == 42.0);            // ...and summed
    // An edge exactly AT the stored value is inclusive at both ends.
    CHECK(s4.sumRange(0, double(stored), double(stored)) == 42.0);
  }

  if (g_fail) { std::printf("odia_spectrumstore_test FAILED (%d)\n", g_fail); return 1; }
  std::printf("odia_spectrumstore_test OK\n");
  return 0;
}
