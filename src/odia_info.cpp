// odia-info — report the DIA acquisition scheme of a run.
//
// Exists because the extraction cost model depends on the *actual* isolation
// window layout and cycle structure, not on assumptions about them. Every
// candidate-count estimate in the design must be derived from real numbers
// this tool prints.

#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <algorithm>
#include <cstdio>
#include <map>
#include <string>
#include <vector>

namespace
{
struct Window
{
  double lower, upper;
  bool operator<(const Window& o) const
  {
    return lower != o.lower ? lower < o.lower : upper < o.upper;
  }
};

// Median of an unsorted copy. Median not mean: a stalled scan or a gap between
// LC segments skews the mean badly, and cycle time drives the RT tolerance.
double median(std::vector<double> v)
{
  if (v.empty()) return 0.0;
  const size_t mid = v.size() / 2;
  std::nth_element(v.begin(), v.begin() + mid, v.end());
  return v[mid];
}

// Logarithmic m/z binning gives constant relative (ppm) width, so one integer
// index addresses a tolerance-sized bin anywhere in the spectrum. This is the
// addressing scheme the prefilter depends on.
inline long mz_bin(double mz, double mz0, double inv_log1p_tol)
{
  return static_cast<long>(std::log(mz / mz0) * inv_log1p_tol);
}

#define CHECK(cond)                                                                                \
  if (!(cond))                                                                                     \
  {                                                                                                \
    std::fprintf(stderr, "selftest FAILED at line %d: %s\n", __LINE__, #cond);                     \
    return 1;                                                                                      \
  }

int selftest()
{
  // ponytail: one runnable check on the only non-trivial logic here.
  CHECK(median({}) == 0.0);
  CHECK(median({5.0}) == 5.0);
  CHECK(median({3.0, 1.0, 2.0}) == 2.0);
  Window a{100.0, 110.0}, b{100.0, 120.0}, c{200.0, 210.0};
  CHECK(a < b); CHECK(b < c); CHECK(!(c < a));

  const double tol = 20e-6, mz0 = 200.0, inv = 1.0 / std::log1p(tol);
  CHECK(mz_bin(mz0, mz0, inv) == 0);

  // Binning is NOT tolerance-exact: two m/z within tolerance may straddle a bin
  // boundary and land in ADJACENT bins. The guaranteed invariant is only
  // "same or adjacent" — which is why the prefilter must probe neighbours.
  // An earlier version of this test asserted "same bin" and failed; the test
  // was wrong, and the design consequence is real (see MZ_BIN_NEIGHBOURS).
  for (double mz : {250.0, 500.0, 1234.5, 1999.0})
  {
    CHECK(std::abs(mz_bin(mz, mz0, inv) - mz_bin(mz * (1 + tol / 2), mz0, inv)) <= 1);
    CHECK(std::abs(mz_bin(mz, mz0, inv) - mz_bin(mz * (1 - tol / 2), mz0, inv)) <= 1);
  }
  // Well beyond tolerance must land strictly further away than a neighbour.
  CHECK(std::abs(mz_bin(500.0, mz0, inv) - mz_bin(500.0 * (1 + 10 * tol), mz0, inv)) > 1);
  CHECK(mz_bin(500.0, mz0, inv) < mz_bin(501.0, mz0, inv));

  std::puts("selftest OK");
  return 0;
}
#undef CHECK

// Measure how selective an occupancy-bitmap prefilter would be.
//
// The whole architecture rests on rejecting candidates by testing their
// fragment m/z values against a bitmap of occupied bins, before doing any
// gather. If observed spectra are dense enough that most random fragment sets
// pass, the prefilter is worthless and the design is wrong. This measures it
// rather than assuming it.
int occupancy(const OpenMS::MSExperiment& exp, double ppm)
{
  const double tol = ppm * 1e-6;
  const double inv = 1.0 / std::log1p(tol);
  const double mz0 = 100.0;

  std::vector<double> occ_frac, npeaks;
  double lo = 1e12, hi = 0.0;
  std::vector<long> bins;

  for (const auto& s : exp.getSpectra())
  {
    if (s.getMSLevel() != 2 || s.size() < 2) continue;
    bins.clear();
    bins.reserve(s.size());
    for (const auto& p : s)
    {
      bins.push_back(mz_bin(p.getMZ(), mz0, inv));
      lo = std::min(lo, p.getMZ());
      hi = std::max(hi, p.getMZ());
    }
    std::sort(bins.begin(), bins.end());
    const size_t distinct = std::unique(bins.begin(), bins.end()) - bins.begin();
    const long span = bins.back() - bins.front() + 1;
    if (span > 0) occ_frac.push_back(double(distinct) / double(span));
    npeaks.push_back(double(distinct));
  }

  if (occ_frac.empty()) { std::fprintf(stderr, "no MS2 spectra\n"); return 1; }

  // A probe must check the bin and both neighbours (see selftest): values
  // within tolerance can straddle a boundary. That triples the chance a random
  // fragment registers a hit, so the honest selectivity uses 3p, not p.
  constexpr int MZ_BIN_NEIGHBOURS = 3;
  const double p_raw = median(occ_frac);
  const double p = std::min(1.0, MZ_BIN_NEIGHBOURS * p_raw);
  const double med_peaks = median(npeaks);
  const long total_bins = mz_bin(hi, mz0, inv) - mz_bin(lo, mz0, inv) + 1;

  std::printf("tolerance        %.0f ppm\n", ppm);
  std::printf("MS2 spectra      %zu\n", occ_frac.size());
  std::printf("fragment m/z     %.1f - %.1f Th  (%ld bins)\n", lo, hi, total_bins);
  std::printf("distinct bins    %.0f per spectrum (median)\n", med_peaks);
  std::printf("occupancy        %.5f  (%.3f%% of bins in span)\n", p_raw, 100 * p_raw);
  std::printf("hit prob (x%d nb) %.5f\n", MZ_BIN_NEIGHBOURS, p);
  std::puts("\nP(random candidate passes prefilter), F=6 fragments:");
  std::puts("  need k of 6    probability      rejection factor");

  // Binomial tail: P(>= k hits of F) for independent random fragments.
  const int F = 6;
  for (int k = 6; k >= 2; --k)
  {
    double tail = 0.0;
    for (int i = k; i <= F; ++i)
    {
      double c = 1.0;
      for (int j = 0; j < i; ++j) c = c * (F - j) / (j + 1);
      tail += c * std::pow(p, i) * std::pow(1 - p, F - i);
    }
    std::printf("  %d            %12.3e     %12.3e x\n", k, tail, tail > 0 ? 1.0 / tail : 0.0);
  }
  return 0;
}
} // namespace

int main(int argc, char** argv)
{
  if (argc == 2 && std::string(argv[1]) == "--selftest") return selftest();
  if (argc < 2)
  {
    std::fprintf(stderr,
                 "usage: odia-info <file.mzML> [--occupancy <ppm>]\n"
                 "       odia-info --selftest\n");
    return 2;
  }
  double occ_ppm = 0.0;
  if (argc == 4 && std::string(argv[2]) == "--occupancy") occ_ppm = std::atof(argv[3]);

  OpenMS::MSExperiment exp;
  try
  {
    OpenMS::MzMLFile().load(argv[1], exp);
  }
  catch (const std::exception& e)
  {
    std::fprintf(stderr, "error: cannot read %s: %s\n", argv[1], e.what());
    return 1;
  }

  if (occ_ppm > 0.0) return occupancy(exp, occ_ppm);

  std::map<Window, size_t> windows;
  std::vector<double> ms1_rts;
  size_t ms1 = 0, ms2 = 0, other = 0, peaks = 0;
  bool has_im = false;

  for (const auto& s : exp.getSpectra())
  {
    peaks += s.size();
    if (s.getDriftTime() != -1.0) has_im = true;

    switch (s.getMSLevel())
    {
      case 1:
        ++ms1;
        ms1_rts.push_back(s.getRT());
        break;
      case 2:
      {
        ++ms2;
        if (s.getPrecursors().empty()) break;
        const auto& p = s.getPrecursors().front();
        windows[{p.getMZ() - p.getIsolationWindowLowerOffset(),
                 p.getMZ() + p.getIsolationWindowUpperOffset()}]++;
        break;
      }
      default: ++other;
    }
  }

  // Cycle time from consecutive MS1 scans — the interval over which one full
  // set of DIA windows is acquired, i.e. the sampling period of every
  // chromatographic trace we will later extract.
  std::vector<double> deltas;
  for (size_t i = 1; i < ms1_rts.size(); ++i) deltas.push_back(ms1_rts[i] - ms1_rts[i - 1]);

  std::printf("file            %s\n", argv[1]);
  std::printf("spectra         %zu  (MS1 %zu, MS2 %zu, other %zu)\n", exp.size(), ms1, ms2, other);
  std::printf("peaks           %zu\n", peaks);
  std::printf("ion mobility    %s\n", has_im ? "yes" : "no");
  if (!ms1_rts.empty())
    std::printf("RT range        %.1f - %.1f s\n", ms1_rts.front(), ms1_rts.back());
  std::printf("cycle time      %.3f s (median)\n", median(deltas));
  std::printf("distinct windows %zu\n", windows.size());

  if (!windows.empty())
  {
    std::vector<double> widths;
    for (const auto& [w, n] : windows) widths.push_back(w.upper - w.lower);
    std::printf("window width    %.2f Th (median)\n", median(widths));
    std::printf("window range    %.2f - %.2f Th\n", windows.begin()->first.lower,
                windows.rbegin()->first.upper);
    std::puts("\n  lower     upper    width   scans");
    for (const auto& [w, n] : windows)
      std::printf("  %8.2f %8.2f %7.2f %7zu\n", w.lower, w.upper, w.upper - w.lower, n);
  }
  return 0;
}
