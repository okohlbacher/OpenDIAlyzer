// Self-check for odia_features.h. Standalone: no OpenMS.
//
// The point of the class is a size and allocation-count claim, so the test asserts those against
// the MEASURED shape of the Astral benchmark rather than against round numbers.

#include "odia_features.h"

#include <cassert>
#include <cstdio>
#include <vector>

using odia::FeatureTable;

int main()
{
  std::printf("odia_features_test\n");

  // Measured on the benchmark: 423,079 precursors -> 2,070,081 features, 30,173,963 subordinates,
  // 113,855,060 meta values, i.e. 4.89 features/precursor, 14.58 subordinates/feature, 55 columns.
  const std::size_t n_prec = 423079;
  const std::size_t n_feat = 2070081;
  const std::size_t n_sub  = 30173963;
  const std::size_t n_cols = 55;

  std::vector<std::string> cols;
  for (std::size_t i = 0; i < n_cols; ++i) { cols.push_back("VAR_SCORE_" + std::to_string(i)); }

  FeatureTable t;
  t.setColumns(cols);
  assert(t.columnCount() == n_cols);

  // Populate a 1/1000 scale model, then extrapolate -- allocating the real thing would need the
  // memory this class exists to avoid.
  const std::size_t scale = 1000;
  const std::size_t feats = n_feat / scale, subs = n_sub / scale;
  std::vector<float> row(n_cols, 1.0f);
  t.reserve(feats);
  for (std::size_t i = 0; i < feats; ++i)
  {
    t.append(std::int64_t(i), std::int64_t(i / 5), 100.0f, 1.0f, 0.5f, 90.0f, 110.0f, 0.0f,
             row.data());
  }
  for (std::size_t i = 0; i < subs; ++i) { t.appendTransition(std::int64_t(i / 15), std::int64_t(i), 1.0f, 2.0f); }

  assert(t.size() == feats);
  assert(t.transitionCount() == subs);
  assert(t.row(7)[3] == 1.0f);
  assert(t.featureId(7) == 7);

  const double flat   = double(t.bytes()) * scale;
  const double nested = double(t.bytesAsFeatureMap()) * scale;
  const double flat_kb_per_prec   = flat   / n_prec / 1024.0;
  const double nested_kb_per_prec = nested / n_prec / 1024.0;

  std::printf("  per precursor: flat %.2f KB vs Feature/MetaInfo %.1f KB  (%.1fx)\n",
              flat_kb_per_prec, nested_kb_per_prec, nested / flat);
  std::printf("  allocations:   flat %zu total vs %.0f (one MetaInfo per feature and subordinate)\n",
              t.allocationCount(), double(t.allocationCountAsFeatureMap()) * scale);

  // The nested figure must land near the 30.7 KB/precursor the component walk measured. If this
  // drifts, either the model here or the constants in bytesAsFeatureMap() are wrong.
  assert(nested_kb_per_prec > 20.0 && nested_kb_per_prec < 45.0);

  // The whole justification is an order-of-magnitude reduction. Anything less and the change is
  // not worth replacing the bundle writer for.
  assert(nested / flat > 10.0);

  // Allocation count must be a CONSTANT, not proportional to the data -- that is what stops the
  // arena fragmenting.
  assert(t.allocationCount() < 100);
  assert(double(t.allocationCountAsFeatureMap()) * scale > 30e6);

  std::printf("odia_features_test OK\n");
  return 0;
}
