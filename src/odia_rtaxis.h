#pragma once
// The run's retention-time axis, and the rule that RT is stored as an INDEX into it.
//
// THE RULE. Inside this tool a retention time is a uint32 index, not a double. Everything that
// needs seconds asks the axis. There is exactly one place where an RT becomes a number, and it is
// here.
//
// WHY, and it is not primarily about bytes. Two reasons, in order of importance:
//
//   1. RECALIBRATION BECOMES AN AXIS OPERATION. A run is recalibrated by fitting library RT to
//      observed RT and re-extracting on the corrected axis. When every structure stores its own
//      double, recalibration means rewriting all of them -- 14.9M scored rows, every feature, every
//      chromatogram -- and any structure missed silently keeps the old calibration. When they store
//      indices, recalibration rewrites ONE array of ~1,945 doubles and every index follows for
//      free. Nothing can be left behind, because nothing else holds an RT.
//
//   2. It is smaller. 14.9M scored rows x 2 RTs x 8 B = 238 MB as doubles, 119 MB as uint32, and
//      the chromatogram store drops 42.7 GB of duplicated timestamps (odia_chromstore.h). But the
//      correctness argument above is the one that matters -- the bytes are a side effect.
//
// PRECISION. Indices snap to the acquisition grid, so an RT is exact for anything OBSERVED (an
// apex, a boundary, a cycle time -- these ARE grid points) and quantised to one cycle for anything
// PREDICTED. On this instrument a cycle is ~1.2 s against a measured peak width of 8.9 s (median)
// and an RT prediction error of 21.7 s (median), so the quantisation is ~1.4% of the error it is
// added to. That is measured, not assumed: see docs/OpenDIAlyzer-classifier-backlog.md.
//
// Predictions outside the acquisition range clamp to the first or last index rather than wrapping
// or going out of bounds, and `kNoRt` marks "no retention time" so a missing value cannot be
// confused with the start of the run.

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

namespace odia
{

/// A retention time that is absent. Distinct from index 0, which is the run's first cycle.
inline constexpr std::uint32_t kNoRt = std::numeric_limits<std::uint32_t>::max();

/// The run's time grid. Ascending, one entry per acquisition cycle.
class RtAxis
{
public:
  RtAxis() = default;

  explicit RtAxis(std::vector<double> seconds) { reset(std::move(seconds)); }

  /// Install the grid. Must be ascending: the index<->time mapping is a binary search, and a
  /// non-monotone axis would silently return wrong indices rather than fail.
  void reset(std::vector<double> seconds)
  {
    if (!std::is_sorted(seconds.begin(), seconds.end()))
    {
      throw std::invalid_argument("odia::RtAxis: retention times must be ascending");
    }
    t_ = std::move(seconds);
  }

  bool empty() const { return t_.empty(); }
  std::size_t size() const { return t_.size(); }

  /// Seconds at an index. kNoRt yields NaN, so a missing RT propagates as a missing value rather
  /// than as the start of the run.
  double seconds(std::uint32_t i) const
  {
    if (i == kNoRt || t_.empty()) { return std::numeric_limits<double>::quiet_NaN(); }
    return t_[std::min<std::size_t>(i, t_.size() - 1)];
  }

  /// Nearest index to a time in seconds. Clamps at both ends; NaN maps to kNoRt.
  std::uint32_t index(double sec) const
  {
    if (t_.empty() || !std::isfinite(sec)) { return kNoRt; }
    const auto it = std::lower_bound(t_.begin(), t_.end(), sec);
    if (it == t_.begin()) { return 0; }
    if (it == t_.end()) { return static_cast<std::uint32_t>(t_.size() - 1); }
    // Round to the NEARER neighbour rather than always down: always-down biases every stored RT
    // half a cycle early, which over a whole run is a systematic shift, not noise.
    const std::size_t hi = static_cast<std::size_t>(it - t_.begin());
    const std::size_t lo = hi - 1;
    return static_cast<std::uint32_t>((sec - t_[lo] <= t_[hi] - sec) ? lo : hi);
  }

  /// Half the local cycle spacing at @p i -- the worst error the index representation can carry.
  /// Exposed so callers can report the quantisation rather than assume it is negligible.
  double quantisationSeconds(std::uint32_t i) const
  {
    if (t_.size() < 2 || i == kNoRt) { return 0.0; }
    const std::size_t k = std::min<std::size_t>(i, t_.size() - 1);
    const double a = (k > 0) ? t_[k] - t_[k - 1] : t_[1] - t_[0];
    const double b = (k + 1 < t_.size()) ? t_[k + 1] - t_[k] : a;
    return 0.5 * std::max(a, b);
  }

  /// RECALIBRATE: replace the grid's times, keeping every stored index valid.
  ///
  /// This is the operation the whole representation exists for. `map` transforms a time in seconds;
  /// applying it here re-times the entire run at once, and every index held anywhere -- scored
  /// rows, features, chromatograms -- follows without being touched. Nothing can be missed, because
  /// nothing else stores an RT.
  ///
  /// The result must stay ascending. A transform that reorders the run is not a recalibration, it
  /// is a bug, and it is rejected rather than silently producing an axis whose binary search is
  /// meaningless.
  template <typename Fn>
  void recalibrate(Fn&& map)
  {
    std::vector<double> out;
    out.reserve(t_.size());
    for (const double v : t_) { out.push_back(map(v)); }
    if (!std::is_sorted(out.begin(), out.end()))
    {
      throw std::invalid_argument("odia::RtAxis::recalibrate: transform is not monotone; that "
                                  "reorders the run rather than re-timing it");
    }
    t_.swap(out);
  }

  const std::vector<double>& seconds() const { return t_; }

private:
  std::vector<double> t_;
};

} // namespace odia
