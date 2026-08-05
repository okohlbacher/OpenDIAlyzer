// Self-check for odia_scored.h, at the real table shapes.
//
// The claims worth testing are the ones that would fail SILENTLY: a CSR grouping that misplaces a
// feature's subordinates (the rows are still there, just attributed to the wrong parent), an
// absent column that quietly takes a storage slot, and a size claim that is asserted rather than
// measured.

#include "odia_scored.h"

#include <cstdio>
#include <random>
#include <string>
#include <vector>

using namespace odia;

static int g_fail = 0;
#define CHECK(c)                                                                       \
  do {                                                                                 \
    if (!(c)) { std::printf("  FAIL %s:%d: %s\n", __FILE__, __LINE__, #c); ++g_fail; } \
  } while (0)

int main()
{
  std::printf("odia_scored_test\n");

  // ---- 1. absent columns take no slot and no bytes ---------------------------------------------
  // features declares 65 columns; 16 are all-null and 2 constant on the measured run. A column
  // that is declared but absent must not consume a storage slot, or the block is 65 wide again.
  ColumnSet cs;
  std::vector<std::string> declared;
  for (int i = 0; i < 65; ++i) { declared.push_back("col" + std::to_string(i)); }
  cs.declare(declared);
  for (int i = 0; i < 47; ++i) { cs.present("col" + std::to_string(i)); }   // 47 of 65 carry data
  CHECK(cs.declared() == 65);
  CHECK(cs.stored() == 47);
  CHECK(cs.isPresent(0) && cs.slot(0) == 0);
  CHECK(!cs.isPresent(64));
  CHECK(cs.slot(64) == ColumnSet::kAbsent);          // absent: no slot, writer emits null
  bool threw = false;
  try { cs.present("not_a_column"); } catch (const std::invalid_argument&) { threw = true; }
  CHECK(threw);                                       // typo must fail, not silently add a column
  cs.present("col0");                                 // idempotent
  CHECK(cs.stored() == 47);

  // ---- 2. features at the real shape ------------------------------------------------------------
  const std::size_t kFeat = 200000;
  FeatureRows fr;
  fr.setColumns(cs);
  fr.setRunId(1168718576198231433LL);                 // constant: once, not 2.07M times
  fr.reserve(kFeat);
  std::mt19937 rng(3);
  std::uniform_real_distribution<double> us(0.0, 1.0);
  std::vector<double> row(cs.stored());
  for (std::size_t i = 0; i < kFeat; ++i)
  {
    for (auto& v : row) { v = us(rng); }
    row[0] = double(i);                               // a marker to verify the block indexing
    fr.append(std::int64_t(i) * 7919, std::int64_t(i), 100.0F + float(i) * 0.01F,
              1.5F, 0.3F, 99.0F, 101.0F, row.data());
  }
  fr.compact();
  CHECK(fr.size() == kFeat);
  CHECK(fr.runId() == 1168718576198231433LL);
  CHECK(fr.featureId(1234) == 1234 * 7919);
  CHECK(fr.scores(1234)[0] == 1234.0);                // the block is indexed per row, not per col
  CHECK(fr.scores(kFeat - 1)[0] == double(kFeat - 1));
  std::printf("  features           %6.1f B/row  (naive 8 B x 65 = 520)\n", fr.bytesPerRow());
  CHECK(fr.bytesPerRow() < 420.0);

  // ---- 3. CSR grouping: the parent key is never stored, and must still be recoverable ----------
  // If the offsets are off by one feature, every subordinate is attributed to its neighbour --
  // the rows are all present and the totals all match, which is exactly why it would not be
  // noticed downstream.
  TransitionRows tr;
  const std::size_t kSub = 11;                        // ~10.6 subordinates per feature, measured
  tr.reserve(kFeat, kFeat * kSub);
  for (std::size_t f = 0; f < kFeat; ++f)
  {
    const std::uint32_t n = std::uint32_t(f % 4 == 0 ? kSub + 2 : kSub);   // ragged on purpose
    for (std::uint32_t k = 0; k < n; ++k)
    {
      tr.append(std::uint32_t(f * 100 + k), float(f), double(k), double(f) * 2.0,
                float(k) * 0.5F, 3.25, -1.5);
    }
    tr.endFeature();
  }
  tr.compact();
  CHECK(tr.features() == kFeat);
  for (std::size_t f : {std::size_t(0), std::size_t(1), std::size_t(4), kFeat - 1})
  {
    const std::uint32_t want = std::uint32_t(f % 4 == 0 ? kSub + 2 : kSub);
    CHECK(tr.count(f) == want);
    // every row of feature f must carry f's marker -- catches an off-by-one in the offsets
    for (std::uint64_t r = tr.begin(f); r < tr.end(f); ++r)
    {
      CHECK(tr.area(r) == float(f));
      CHECK(tr.transitionId(r) == std::uint32_t(f * 100 + (r - tr.begin(f))));
    }
  }
  CHECK(tr.begin(0) == 0);
  CHECK(tr.end(kFeat - 1) == tr.size());              // no gap, no overlap, nothing past the end
  std::printf("  feature_transition %6.1f B/row  (naive 8 B x 44 = 352)  %zu rows\n",
              tr.bytesPerRow(), tr.size());
  CHECK(tr.bytesPerRow() < 60.0);

  // ---- 4. a feature with NO subordinates must not swallow the next one's -----------------------
  {
    TransitionRows t2;
    t2.reserve(3, 4);
    t2.append(1, 1.0F, 1.0, 1.0, 1.0F, 1.0, 1.0);
    t2.endFeature();                                   // feature 0: one row
    t2.endFeature();                                   // feature 1: EMPTY
    t2.append(2, 2.0F, 2.0, 2.0, 2.0F, 2.0, 2.0);
    t2.endFeature();                                   // feature 2: one row
    CHECK(t2.features() == 3);
    CHECK(t2.count(0) == 1);
    CHECK(t2.count(1) == 0);
    CHECK(t2.count(2) == 1);
    CHECK(t2.transitionId(t2.begin(2)) == 2);          // feature 2 gets ITS row, not feature 1's
  }

  // ---- 5. precursor isotopes, same grouping ----------------------------------------------------
  PrecursorRows pr;
  pr.reserve(kFeat, kFeat * 3);
  for (std::size_t f = 0; f < kFeat; ++f)
  {
    for (std::uint8_t k = 0; k < 3; ++k) { pr.append(k, float(f), double(f)); }
    pr.endFeature();
  }
  pr.compact();
  CHECK(pr.features() == kFeat);
  CHECK(pr.end(kFeat - 1) == pr.size());
  CHECK(pr.isotope(pr.begin(7)) == 0 && pr.isotope(pr.begin(7) + 2) == 2);
  std::printf("  feature_precursor  %6.1f B/row  (naive 8 B x 5 = 40)\n", pr.bytesPerRow());
  CHECK(pr.bytesPerRow() < 25.0);

  // ---- 6. the whole-run projection, measured rather than claimed -------------------------------
  {
    const double f_per = fr.bytesPerRow(), t_per = tr.bytesPerRow(), p_per = pr.bytesPerRow();
    const double gb = (2070355.0 * f_per + 21896391.0 * t_per + 1226331.0 * p_per) / 1073741824.0;
    std::printf("  -> full run projection: %.2f GB against 20.55 GB live as OpenMS Features"
                " (%.1fx)\n", gb, 20.55 / gb);
    CHECK(gb < 3.0);
  }

  if (g_fail) { std::printf("odia_scored_test FAILED (%d)\n", g_fail); return 1; }
  std::printf("odia_scored_test OK\n");
  return 0;
}
