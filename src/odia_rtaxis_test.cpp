// Self-check for odia_rtaxis.h.
//
// The properties worth testing are the ones that fail SILENTLY: a rounding rule that biases every
// stored RT in one direction, a missing value that reads as the start of the run, and a
// recalibration that reorders the run instead of re-timing it.

#include "odia_rtaxis.h"

#include <cmath>
#include <cstdio>
#include <vector>

using namespace odia;

static int g_fail = 0;
#define CHECK(c)                                                                        \
  do {                                                                                  \
    if (!(c)) { std::printf("  FAIL %s:%d: %s\n", __FILE__, __LINE__, #c); ++g_fail; }  \
  } while (0)

int main()
{
  std::printf("odia_rtaxis_test\n");

  // A realistic grid: 1945 cycles, ~1.1994 s apart, starting at 12 s.
  std::vector<double> t(1945);
  for (std::size_t i = 0; i < t.size(); ++i) { t[i] = 12.0 + 1.1994 * static_cast<double>(i); }
  RtAxis ax(t);
  CHECK(ax.size() == 1945);

  // ---- 1. round to the NEARER cycle, not always down -----------------------------------------
  // Always-down would bias every stored RT half a cycle early: a systematic shift over the whole
  // run, not noise, and invisible in any single value.
  const double mid = 0.5 * (t[100] + t[101]);
  CHECK(ax.index(t[100] + 0.1) == 100);
  CHECK(ax.index(t[101] - 0.1) == 101);
  CHECK(ax.index(mid - 1e-9) == 100);
  CHECK(ax.index(mid + 1e-9) == 101);
  // Over many random times the mean signed error must be ~0, which is what "nearest" buys.
  double signed_err = 0.0;
  const int kN = 20000;
  for (int k = 0; k < kN; ++k)
  {
    const double sec = t.front() + (t.back() - t.front()) * (static_cast<double>(k) / kN);
    signed_err += ax.seconds(ax.index(sec)) - sec;
  }
  signed_err /= kN;
  std::printf("  mean signed quantisation error over %d times: %+.4f s (cycle %.4f s)\n",
              kN, signed_err, t[1] - t[0]);
  CHECK(std::abs(signed_err) < 0.02);            // ~0; always-down would give about -0.6

  // ---- 2. clamping and the missing value ------------------------------------------------------
  CHECK(ax.index(t.front() - 1000.0) == 0);      // before the run
  CHECK(ax.index(t.back() + 1000.0) == 1944);    // after it
  CHECK(ax.index(std::nan("")) == kNoRt);
  // kNoRt must NOT read as the first cycle -- that is the whole reason it is not 0.
  CHECK(std::isnan(ax.seconds(kNoRt)));
  CHECK(!std::isnan(ax.seconds(0)));

  // ---- 3. quantisation is REPORTED, not assumed ----------------------------------------------
  const double q = ax.quantisationSeconds(500);
  std::printf("  worst quantisation at index 500: %.4f s "
              "(measured peak width 8.9 s, RT error 21.7 s -> %.1f%% of the error)\n",
              q, 100.0 * q / 21.7);
  CHECK(q > 0.5 && q < 0.7);

  // ---- 4. RECALIBRATION MUST NOT MOVE THE ACQUISITION TIMES -----------------------------------
  // The property that matters, and the one an earlier version of this class got wrong by rewriting
  // the seconds. When a cycle was acquired is a physical fact; a calibration cannot change it. What
  // a calibration changes is where the LIBRARY thinks that cycle sits, which is the iRT half.
  const std::uint32_t held_a = ax.index(600.0), held_b = ax.index(1800.0);
  const double sec_a = ax.seconds(held_a), sec_b = ax.seconds(held_b);
  CHECK(!ax.calibrated());                       // no calibration installed yet
  CHECK(std::isnan(ax.irt(held_a)));             // and iRT is absent, not 0 -- 0 is a legal iRT

  ax.setCalibration([](double s) { return 0.001 * s - 0.4; });
  CHECK(ax.calibrated());
  CHECK(ax.seconds(held_a) == sec_a);            // seconds UNTOUCHED
  CHECK(ax.seconds(held_b) == sec_b);
  CHECK(std::abs(ax.irt(held_a) - (0.001 * sec_a - 0.4)) < 1e-6);
  std::printf("  calibrate: index %u stays at %.2f s, iRT now %.4f\n",
              held_a, ax.seconds(held_a), ax.irt(held_a));

  // A SECOND calibration replaces the first and still leaves the seconds alone.
  ax.recalibrate([](double s) { return 0.0012 * s - 0.5; });
  CHECK(ax.seconds(held_a) == sec_a);
  CHECK(std::abs(ax.irt(held_a) - (0.0012 * sec_a - 0.5)) < 1e-6);

  // ---- 4b. the two halves are inverses over the same index ------------------------------------
  for (std::uint32_t i = 0; i < 1945; i += 53)
  {
    CHECK(ax.irtToIndex(ax.irt(i)) == i);        // iRT -> index -> iRT is a fixed point
  }
  CHECK(ax.irtToIndex(std::nan("")) == kNoRt);
  std::printf("  iRT<->index round trip exact over the grid\n");

  // ---- 5. a non-monotone calibration is a bug, not a recalibration ----------------------------
  // Elution order is physics. A map that reverses it leaves irtToIndex searching an unsorted array,
  // i.e. every lookup silently wrong.
  bool threw = false;
  try { ax.recalibrate([](double s) { return -s; }); }
  catch (const std::invalid_argument&) { threw = true; }
  CHECK(threw);
  CHECK(ax.seconds(held_a) == sec_a);                                  // seconds still intact
  CHECK(std::abs(ax.irt(held_a) - (0.0012 * sec_a - 0.5)) < 1e-6);     // and the good calibration kept
  std::printf("  non-monotone calibration rejected, both halves intact\n");

  // ---- 5b. re-installing the grid drops a calibration fitted to the old one --------------------
  {
    RtAxis fresh(std::vector<double>{0.0, 1.0, 2.0});
    fresh.setCalibration([](double s) { return s; });
    CHECK(fresh.calibrated());
    fresh.reset(std::vector<double>{0.0, 1.0, 2.0, 3.0});
    CHECK(!fresh.calibrated());                  // stale calibration must not survive a new grid
  }

  // ---- 6. a non-ascending axis is rejected at construction ------------------------------------
  threw = false;
  try { RtAxis bad(std::vector<double>{5.0, 3.0, 9.0}); }
  catch (const std::invalid_argument&) { threw = true; }
  CHECK(threw);

  // ---- 7. round trip: index -> seconds -> index is a fixed point ------------------------------
  // Still exact after calibration, because calibration never touched the seconds.
  for (std::uint32_t i = 0; i < 1945; i += 37)
  {
    CHECK(ax.index(ax.seconds(i)) == i);
  }

  if (g_fail) { std::printf("odia_rtaxis_test FAILED (%d)\n", g_fail); return 1; }
  std::printf("odia_rtaxis_test OK\n");
  return 0;
}
