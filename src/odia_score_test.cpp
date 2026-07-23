// Validate odia::score against OpenSWATH's own reference values.
//
// Every expected number below is copied from
// ext/OpenMS/src/tests/class_tests/openswathalgo/Scoring_test.cpp -- so this is
// differential testing against a decade-old, production-validated oracle, not
// against our own expectations. If a future CUDA port diverges from these, it
// is wrong.

#include "odia_score.h"

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

namespace
{
int failures = 0;

void near(const std::string& what, double got, double want, double tol = 1e-5)
{
  if (std::abs(got - want) > tol)
  {
    std::fprintf(stderr, "FAIL %-28s got %.9f want %.9f (|d|=%.2e)\n",
                 what.c_str(), got, want, std::abs(got - want));
    ++failures;
  }
}
void eq(const std::string& what, int got, int want)
{
  if (got != want) { std::fprintf(stderr, "FAIL %-28s got %d want %d\n", what.c_str(), got, want); ++failures; }
}
} // namespace

int main()
{
  using namespace odia;

  // --- standardize_data (Scoring_test.cpp:standardize_data_test) ---
  {
    std::vector<double> d1{0, 1, 3, 5, 2, 0}, d2{1, 3, 5, 2, 0, 0};
    standardize(d1);
    standardize(d2);
    near("standardize d1[0]", d1[0], -1.03479296);
    near("standardize d1[1]", d1[1], -0.47036043);
    near("standardize d1[2]", d1[2], 0.65850461);
    near("standardize d1[3]", d1[3], 1.78736965);
    near("standardize d1[4]", d1[4], 0.09407209);
    near("standardize d1[5]", d1[5], -1.03479296);
    near("standardize d2[2]", d2[2], 1.78736965);
  }

  // --- cross-correlation (Scoring_test.cpp:calcxcorr_test), maxdelay=2 ---
  // Reference: result.data[k].second at delay = k-2, on standardized inputs.
  {
    auto xc = xcorr({0, 1, 3, 5, 2, 0}, {1, 3, 5, 2, 0, 0}, 2);
    // slot index = delay + maxdelay
    near("xcorr delay=+2", xc[2 + 2], -0.7374631);
    near("xcorr delay=+1", xc[1 + 2], -0.567846);
    near("xcorr delay= 0", xc[0 + 2], 0.4159292);
    near("xcorr delay=-1", xc[-1 + 2], 0.8215339);
    near("xcorr delay=-2", xc[-2 + 2], 0.15634218);

    auto pk = xcorr_max(xc, 2);
    eq("xcorr_max delay", pk.delay, -1);          // 0.8215339 is the largest
    near("xcorr_max value", pk.value, 0.8215339);
  }

  // --- NormalizedManhattanDist (Scoring_test.cpp:74) ---
  {
    std::vector<double> a{0, 1, 3, 5, 2, 0}, b{1, 3, 5, 2, 0, 0};
    near("norm_manhattan", normalized_manhattan(a, b), 0.15151515);
  }

  // --- RootMeanSquareDeviation (Scoring_test.cpp:91) ---
  {
    std::vector<double> a{0, 1, 3, 5, 2, 0}, b{1, 3, 5, 2, 0, 0};
    near("rmsd", rmsd(a, b), 1.91485421551);
  }

  // --- SpectralAngle (Scoring_test.cpp:198-212) ---
  {
    near("spectral_angle 1",
         spectral_angle({0.03174064, 0.11582065, 0.63258941},
                        {0.71882213, 0.00087569, 0.36516896}), 0.7699453419277419);
    near("spectral_angle 2",
         spectral_angle({0.6608937, 0.0726909, 0.40912141},
                        {0.52081914, 0.71088, 0.0175557}), 0.9449782659258582);
    near("spectral_angle 3",
         spectral_angle({0.58858475, 0.08963515, 0.08578046},
                        {0.76180969, 0.72763536, 0.50090751}), 0.6547156284689354);
  }

  // --- all-pairs consistency: diagonal is autocorrelation, peak at delay 0 ---
  {
    std::vector<std::vector<double>> tr{{0, 1, 3, 5, 2, 0}, {1, 3, 5, 2, 0, 0}};
    auto ps = allpairs_xcorr(tr, 2);
    eq("allpairs count", int(ps.size()), 3);        // (0,0),(0,1),(1,1)
    // A standardized trace's autocorrelation at delay 0 is 1.0 by construction.
    for (const auto& p : ps)
      if (p.i == p.j) { near("autocorr peak@0", p.value, 1.0); eq("autocorr delay", p.delay, 0); }
    // The (0,1) off-diagonal peak must equal the pairwise xcorr max above.
    for (const auto& p : ps)
      if (p.i == 0 && p.j == 1) { near("pair(0,1) value", p.value, 0.8215339); eq("pair(0,1) delay", p.delay, -1); }
  }

  if (failures == 0) { std::puts("odia-score: all reference checks OK"); return 0; }
  std::fprintf(stderr, "odia-score: %d FAILURES\n", failures);
  return 1;
}
