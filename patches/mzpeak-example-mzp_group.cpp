// Validate the OpenDIAlyzer adapter's ALGORITHM against real mzPeak data, without OpenMS:
// group spectra by (RT, isolation-window centre), then decode a group and assemble the
// concatenated mz/intensity/drift arrays exactly as MzPeakIndex::decode() does.
#include <mzpeak/open.h>
#include <mzpeak/index.h>
#include <mzpeak/spectra.h>
#include <mzpeak/spectrum.h>
#include <mzpeak/spectrum_metadata.h>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <map>
#include <vector>
int main(int argc, char** argv)
{
  auto index = MzPeak::open(argv[1]);
  auto spectra = index.spectra();
  struct Grp { double rt; std::vector<uint64_t> slices; };
  std::map<long, std::map<long, Grp>> by_window;      // wkey -> rtkey -> group
  uint64_t i = 0; std::size_t n_ms1 = 0, n_ms2 = 0, n_im = 0;
  for (auto sp : spectra) {
    const double rt = sp.retention_time() ? *sp.retention_time() : 0.0;
    long wkey = -1; double centre = -1;
    if (sp.ms_level() != 1) {
      for (const auto& p : sp.precursors())
        if (p.isolation_window.target_mz) { centre = *p.isolation_window.target_mz; break; }
      if (centre < 0) { ++i; continue; }
      wkey = std::lround(centre * 1000.0); ++n_ms2;
    } else ++n_ms1;
    if (sp.ion_mobility()) ++n_im;
    auto& g = by_window[wkey][std::lround(rt * 10000.0)];
    g.rt = rt; g.slices.push_back(i); ++i;
  }
  std::size_t groups = 0, multi = 0;
  for (auto& w : by_window) for (auto& g : w.second) { ++groups; if (g.second.slices.size() > 1) ++multi; }
  std::printf("spectra=%llu  MS1=%zu MS2=%zu with-IM=%zu\n", (unsigned long long)i, n_ms1, n_ms2, n_im);
  std::printf("windows=%zu  groups=%zu  groups-with->1-slice=%zu (IM aggregation exercised: %s)\n",
              by_window.size(), groups, multi, multi ? "YES" : "no (1 spectrum per frame/window)");
  // assemble the middle group of the first MS2 window, exactly like the adapter
  for (auto& w : by_window) {
    if (w.first == -1) continue;
    auto it = w.second.begin(); std::advance(it, w.second.size() / 2);
    std::vector<double> mz, in, dt;
    for (uint64_t si : it->second.slices) {
      auto sp = spectra[static_cast<std::size_t>(si)];
      const auto& m = sp.mz(); const auto& y = sp.intensity();
      const double im = sp.ion_mobility() ? *sp.ion_mobility() : -1.0;
      const std::size_t n = std::min(m.size(), y.size());
      mz.insert(mz.end(), m.begin(), m.begin() + n);
      for (std::size_t k = 0; k < n; ++k) in.push_back((double)y[k]);
      dt.insert(dt.end(), n, im);
    }
    std::vector<std::size_t> ord(mz.size());
    for (std::size_t k = 0; k < ord.size(); ++k) ord[k] = k;
    std::sort(ord.begin(), ord.end(), [&](std::size_t a, std::size_t b){ return mz[a] < mz[b]; });
    bool sorted = true; double prev = -1;
    for (std::size_t k : ord) { if (mz[k] < prev) { sorted = false; break; } prev = mz[k]; }
    std::printf("assembled group (window %ld, rt %.2f): %zu peaks, arrays parallel=%s, m/z-sorted=%s\n",
                w.first, it->second.rt, mz.size(),
                (mz.size()==in.size() && mz.size()==dt.size()) ? "yes":"NO", sorted ? "yes":"NO");
    return mz.empty() ? 1 : 0;
  }
  std::printf("no MS2 window found\n"); return 1;
}
