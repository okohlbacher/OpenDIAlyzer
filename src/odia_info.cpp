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

int selftest()
{
  // ponytail: one runnable check on the only non-trivial logic here.
  if (median({}) != 0.0) return 1;
  if (median({5.0}) != 5.0) return 1;
  if (median({3.0, 1.0, 2.0}) != 2.0) return 1;
  Window a{100.0, 110.0}, b{100.0, 120.0}, c{200.0, 210.0};
  if (!(a < b) || !(b < c) || (c < a)) return 1;
  std::puts("selftest OK");
  return 0;
}
} // namespace

int main(int argc, char** argv)
{
  if (argc == 2 && std::string(argv[1]) == "--selftest") return selftest();
  if (argc != 2)
  {
    std::fprintf(stderr, "usage: odia-info <file.mzML> | --selftest\n");
    return 2;
  }

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
