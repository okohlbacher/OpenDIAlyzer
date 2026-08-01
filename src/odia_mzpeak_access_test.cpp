// Test the mzPeak -> OpenSWATH streaming adapter.
//
// Validates the two things that make it usable as a drop-in input layer:
//   1. STRUCTURE: it produces a sane SwathMap vector (one MS1 + ~24 diaPASEF windows with
//      isolation bounds and an IM range), and every window's groups are RT-ordered.
//   2. CONTRACT: the spectra it hands OpenSWATH satisfy what the extractor requires --
//      m/z-SORTED peaks (extractChromatograms throws otherwise), parallel mz/intensity/drift
//      arrays of equal length, and getSpectraByRT returning the right window.
//   3. BOUNDED MEMORY: decoding a handful of spectra must not pull in the whole run.
//
// Usage: odia-mzpeak-access-test <file.mzpeak>
#include "odia_mzpeak_access.h"

#include <sys/resource.h>
#include <cmath>
#include <cstdio>
#include <string>

static long rss_mb() { struct rusage r; getrusage(RUSAGE_SELF, &r); return r.ru_maxrss / 1024; }

int main(int argc, char** argv)
{
  if (argc < 2) { std::fprintf(stderr, "usage: odia-mzpeak-access-test <file.mzpeak>\n"); return 2; }
  int fails = 0;
  auto fail = [&](const std::string& m) { std::fprintf(stderr, "FAIL: %s\n", m.c_str()); ++fails; };

  const long rss_before = rss_mb();
  odia::MzPeakIndexPtr idx;
  std::vector<OpenSwath::SwathMap> maps = odia::loadMzPeakSwathMaps(argv[1], idx);
  const long rss_index = rss_mb();
  std::printf("index built: %llu mzPeak spectra -> %zu SwathMaps (RSS %ld -> %ld MB)\n",
              (unsigned long long) idx->nSpectraTotal(), maps.size(), rss_before, rss_index);

  // ---- 1. structure -------------------------------------------------------
  std::size_t n_ms1 = 0, n_ms2 = 0;
  for (const auto& m : maps) { (m.ms1 ? n_ms1 : n_ms2)++; }
  std::printf("  MS1 maps %zu, MS2 (isolation windows) %zu\n", n_ms1, n_ms2);
  // NOTE: window count and IM presence are properties of the FILE (diaPASEF vs plain DIA),
  // not of the adapter -- report them, do not fail on them. Only structural invariants below
  // are adapter bugs.
  if (n_ms2 < 1) { fail("no MS2 isolation windows found at all"); }
  if (n_ms1 > 1) { fail("more than one MS1 map, got " + std::to_string(n_ms1)); }

  std::size_t total_groups = 0;
  for (const auto& m : maps)
  {
    const std::size_t n = m.sptr->getNrSpectra();
    total_groups += n;
    if (n == 0) { fail("a SwathMap has no spectra"); continue; }
    if (!m.ms1)
    {
      if (!(m.lower <= m.center && m.center <= m.upper)) { fail("isolation bounds not ordered"); }
      // IM range is only meaningful for ion-mobility data; absence is legal for plain DIA.
    }
    // RT ordering: OpenSWATH assumes spectra are RT-ordered
    double prev = -1e30;
    for (std::size_t i = 0; i < n; ++i)
    {
      const double rt = m.sptr->getSpectrumMetaById(static_cast<int>(i)).RT;
      if (rt < prev) { fail("spectra not RT-ordered"); break; }
      prev = rt;
    }
  }
  std::printf("  total (frame x window) groups: %zu\n", total_groups);

  // ---- 2. contract on real decoded spectra --------------------------------
  const OpenSwath::SwathMap* ms2 = nullptr;
  for (const auto& m : maps) { if (!m.ms1) { ms2 = &m; break; } }
  if (!ms2) { fail("no MS2 map to test"); std::fprintf(stderr, "%d failure(s)\n", fails); return fails ? 1 : 0; }

  // Sample from the MIDDLE of the gradient: the first groups of a run are usually
  // pre-elution and legitimately empty, which would make an "are peaks decoded?" check
  // vacuously pass (or look like a decode failure).
  const std::size_t n_tot = ms2->sptr->getNrSpectra();
  const std::size_t n_probe = std::min<std::size_t>(5, n_tot);
  const std::size_t first = (n_tot > n_probe) ? (n_tot / 2) : 0;
  std::size_t peaks_seen = 0;
  for (std::size_t i = first; i < first + n_probe; ++i)
  {
    OpenSwath::SpectrumPtr sp = ms2->sptr->getSpectrumById(static_cast<int>(i));
    auto mz = sp->getMZArray(); auto in = sp->getIntensityArray(); auto dt = sp->getDriftTimeArray();
    if (!mz || !in) { fail("spectrum missing mz/intensity array"); continue; }
    if (mz->data.size() != in->data.size()) { fail("mz/intensity length mismatch"); }
    if (dt && dt->data.size() != mz->data.size()) { fail("drift/mz length mismatch"); }
    for (std::size_t k = 1; k < mz->data.size(); ++k)
    {
      if (mz->data[k] < mz->data[k - 1]) { fail("peaks NOT m/z-sorted -- extractChromatograms would throw"); break; }
    }
    peaks_seen += mz->data.size();
  }
  std::printf("  decoded %zu spectra, %zu peaks, all m/z-sorted with parallel drift array\n",
              n_probe, peaks_seen);
  if (peaks_seen == 0) { fail("no peaks decoded from mid-gradient spectra -- decoding is broken"); }

  // getSpectraByRT must land in the requested window
  {
    const double rt_mid = ms2->sptr->getSpectrumMetaById(static_cast<int>(ms2->sptr->getNrSpectra() / 2)).RT;
    auto sel = ms2->sptr->getSpectraByRT(rt_mid, 5.0);
    if (sel.empty()) { fail("getSpectraByRT(+-5s) returned nothing at a real RT"); }
    for (std::size_t s : sel)
    {
      const double rt = ms2->sptr->getSpectrumMetaById(static_cast<int>(s)).RT;
      if (std::abs(rt - rt_mid) > 5.0 + 1e-6) { fail("getSpectraByRT returned an out-of-window spectrum"); break; }
    }
    auto one = ms2->sptr->getSpectraByRT(rt_mid, -1.0);
    if (one.size() != 1) { fail("getSpectraByRT(deltaRT<0) should return exactly the nearest spectrum"); }
    std::printf("  getSpectraByRT ok (%zu spectra within +-5s of %.2fs)\n", sel.size(), rt_mid);
  }

  // ---- 3. bounded memory --------------------------------------------------
  const long rss_after = rss_mb();
  std::printf("  RSS after decoding: %ld MB (index %ld MB)\n", rss_after, rss_index);
  // the whole run is ~10 GB on disk; a metadata index + a few spectra must stay far below that
  if (rss_after > 8000) { fail("RSS > 8 GB -- adapter is not streaming (materialising the run?)"); }

  // lightClone must be cheap and independent
  {
    auto c = ms2->sptr->lightClone();
    if (!c || c->getNrSpectra() != ms2->sptr->getNrSpectra()) { fail("lightClone lost state"); }
  }

  std::fprintf(stderr, fails ? "mzpeak access test FAILED (%d)\n" : "mzpeak access test OK\n", fails);
  return fails ? 1 : 0;
}
