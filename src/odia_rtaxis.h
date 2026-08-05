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
    irt_.clear();               // a new grid invalidates any calibration fitted to the old one
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

  // ---- the iRT half of the axis ---------------------------------------------------------------
  //
  // TWO PARALLEL ARRAYS OVER THE SAME INDEX. `t_` is when a cycle was acquired, in seconds. `irt_`
  // is where that cycle sits in the library's iRT space under the CURRENT calibration. One index
  // addresses both, so a caller asks for whichever it actually needs and never converts.
  //
  // An earlier version of this class re-timed the run by rewriting `t_`, which is wrong: the
  // acquisition times are physical facts and no calibration changes them. What a calibration
  // changes is the MAP into iRT space. Keeping them apart means a recalibration cannot corrupt the
  // observed times, and the two questions -- "when was this seen" and "where does the library think
  // it belongs" -- stop being answerable only by whichever one happened to be stored.
  //
  // f32 for iRT, f64 for seconds, deliberately. Seconds span 0..2333 and are compared against
  // windows a few seconds wide, so the ulp must stay far below that. iRT is a normalised prediction
  // whose own error is percent-scale, so f32's ~1e-7 relative precision is six orders below the
  // signal -- and it halves the array. That asymmetry is the point: precision follows what the
  // number is used for, not what it is called.

  /// Install the iRT for every cycle from a transform seconds -> iRT. Rejects a non-monotone
  /// result: elution order is physics, and a map that reorders it makes irtToIndex meaningless.
  template <typename Fn>
  void setCalibration(Fn&& to_irt)
  {
    std::vector<float> out;
    out.reserve(t_.size());
    for (const double v : t_) { out.push_back(static_cast<float>(to_irt(v))); }
    if (!std::is_sorted(out.begin(), out.end()))
    {
      throw std::invalid_argument("odia::RtAxis::setCalibration: seconds -> iRT is not monotone; "
                                  "that reorders the run rather than recalibrating it");
    }
    irt_.swap(out);
  }

  bool calibrated() const { return irt_.size() == t_.size() && !irt_.empty(); }

  /// iRT at a cycle. NaN when uncalibrated or the index is missing -- never 0, which is a legal iRT.
  double irt(std::uint32_t i) const
  {
    if (i == kNoRt || !calibrated()) { return std::numeric_limits<double>::quiet_NaN(); }
    return irt_[std::min<std::size_t>(i, irt_.size() - 1)];
  }

  /// Nearest cycle to an iRT, the inverse of irt(). Clamps at both ends; NaN yields kNoRt.
  std::uint32_t irtToIndex(double v) const
  {
    if (!calibrated() || !std::isfinite(v)) { return kNoRt; }
    const auto it = std::lower_bound(irt_.begin(), irt_.end(), static_cast<float>(v));
    if (it == irt_.begin()) { return 0; }
    if (it == irt_.end()) { return static_cast<std::uint32_t>(irt_.size() - 1); }
    const std::size_t hi = static_cast<std::size_t>(it - irt_.begin());
    const std::size_t lo = hi - 1;
    return static_cast<std::uint32_t>((v - irt_[lo] <= irt_[hi] - v) ? lo : hi);
  }

  /// RECALIBRATE: install a new seconds -> iRT map. Every index held anywhere stays valid and
  /// keeps meaning the same acquisition cycle; only where the library thinks that cycle belongs
  /// changes. Nothing else in the tool stores an iRT, so nothing can be left on the old one.
  template <typename Fn>
  void recalibrate(Fn&& to_irt) { setCalibration(std::forward<Fn>(to_irt)); }

  const std::vector<double>& seconds() const { return t_; }
  const std::vector<float>& irt() const { return irt_; }

private:
  std::vector<double> t_;      ///< acquisition seconds, immutable once installed
  std::vector<float> irt_;     ///< iRT of each cycle under the current calibration
};

} // namespace odia
