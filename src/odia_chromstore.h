#pragma once
// A compact chromatogram store: 1 byte per point instead of 16.
//
// THE PROBLEM. OpenMS holds an extracted chromatogram as a vector of ChromatogramPeak, which is a
// double RT plus a double intensity -- 16 bytes per point. On the Astral benchmark that is
// 4,463,919 transitions x ~1,196 cycles x 16 B = 85 GB, which is larger than the 6.0 GB raw file
// the data was extracted FROM. That is the wrong way round: extraction densifies a sparse signal,
// so it should cost more than the peaks but nothing like fourteen times the file.
//
// WHERE THE 85 GB GOES, and it is not the signal:
//
//   * 42.7 GB is the RETENTION TIME, stored once per transition. Every transition extracted from
//     the same SWATH window shares an IDENTICAL time axis -- the cycle times of that window -- so
//     those 1,196 timestamps are stored 4.46 million times over.
//   * 42.7 GB is intensity at double precision, from an instrument whose counts carry perhaps five
//     significant digits.
//
// THE LAYOUT. One global axis per SWATH window, and each chromatogram is a contiguous SLICE of it:
//
//   axis:   per window, the ascending MS2 cycle times           (150 x ~1945 x 8 B = 2.3 MB total)
//   chrom:  window id, start index, length, scale, uint8[length]
//
// A chromatogram therefore stores NO times at all: point i is at axis[window][start + i]. That is
// exactly the "index to the start RT and length, infer the rest" arrangement, and it also makes the
// slice operation this store exists to support -- narrowing an RT window -- a pair of integers,
// with no data movement.
//
//   per point:        1 B          (uint8 intensity, LOG-spaced against a per-chromatogram max)
//   per chromatogram: 16 B         (window, start, length, scale)
//   benchmark total:  4.46M x (16 + 1196) = 5.4 GB   vs 85 GB, a 15.5x reduction (measured)
//
// WHY LOG SPACING, AND WHY IT IS NOT A COMPROMISE. What the scores actually read is SHAPE:
// xcorr_shape, library_corr, dotprod and manhattan are all scale-invariant, so the quantity to
// preserve is RELATIVE accuracy per point, uniformly. Linear quantisation does the opposite -- it
// spends its resolution near the apex, where the data least needs it, and starves the baseline that
// S/N is estimated from. Log spacing gives (1e4)^(1/254) - 1 = 3.7% relative precision at EVERY
// intensity, so 255 codes cover four decades.
//
// Measured on a realistic 1e6 apex over a 1e2 baseline (odia_chromstore_test):
//
//     encoding         B/pt   apex relerr   baseline relerr   1 - corr
//     Quantised8Log       1      0.00e+00            0.0001   5.09e-05
//     Quantised8Lin       1      0.00e+00            1.0000   2.78e-06
//     Quantised16         2      0.00e+00            0.0682   1.28e-10
//     Float32             4      0.00e+00            0.0000   0.00e+00
//
// One byte LOG-spaced preserves the baseline 680x better than two bytes linear, at half the size.
// Linear 8-bit annihilates it outright (the code step is 3922 counts under a 1e6 apex, so a
// 100-count baseline rounds to zero) -- the failure is the SPACING, not the bit depth.
//
// Note also that correlation alone would have chosen wrongly: Quantised8Lin scores BETTER on
// 1 - corr than Quantised8Log, because correlation is dominated by the apex it preserved, while
// the baseline it destroyed does not show up there. That is why the test asserts both.
//
// ChromEncoding::Float32 remains, so a suspected quantisation effect can be ruled out by changing
// one value rather than by reasoning about it.

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

namespace odia
{

/// How intensities are stored. Quantised is the point of this file; Float exists so a suspected
/// quantisation effect can be tested by changing one value rather than by reasoning about it.
enum class ChromEncoding
{
  Quantised8Log,  ///< uint8 + f32 scale, LOG-spaced: 1 B/point. Default -- see the note below.
  Quantised8Lin,  ///< uint8 + f32 scale, linear: 1 B/point. Provided to demonstrate why not.
  Quantised16,    ///< uint16 + f32 scale, linear: 2 B/point
  Float32         ///< float: 4 B/point, no quantisation error
};

/// Dynamic range covered by the log encodings, as a ratio to the chromatogram's own maximum.
/// 1e4 spans the benchmark's observed 1e2 baseline under a 1e6 apex with room to spare; 255 codes
/// across it give (1e4)^(1/254) - 1 = 3.7% relative precision at EVERY intensity.
inline constexpr double kChromLogDynamicRange = 1.0e4;

/// The shared time axis of one SWATH window. Chromatograms reference it by index.
struct ChromAxis
{
  std::vector<double> rt;      ///< ascending cycle times, seconds

  /// First index with rt >= t (lower_bound). Used to turn an RT window into a slice.
  std::size_t lower(double t) const
  {
    return static_cast<std::size_t>(std::lower_bound(rt.begin(), rt.end(), t) - rt.begin());
  }
  /// First index with rt > t (upper_bound).
  std::size_t upper(double t) const
  {
    return static_cast<std::size_t>(std::upper_bound(rt.begin(), rt.end(), t) - rt.begin());
  }
};

/// A stored chromatogram: no times, just where it starts on its window's axis and how long it is.
struct ChromRef
{
  std::uint32_t axis = 0;      ///< which window's axis
  std::uint32_t start = 0;     ///< first index on that axis
  std::uint32_t length = 0;    ///< number of points
  std::uint64_t offset = 0;    ///< where the payload begins in the store's blob
  float scale = 0.0F;          ///< the chromatogram's maximum: the top code decodes to this
};

/// Compact store for many chromatograms over a few shared axes.
///
/// Deliberately append-only and index-addressed: the extractor produces chromatograms once, in an
/// order it chooses, and everything downstream refers to them by index. There is no erase, because
/// a store that could shrink would need the indices to be stable across it, which is a different
/// and much more expensive object.
class ChromStore
{
public:
  explicit ChromStore(ChromEncoding enc = ChromEncoding::Quantised8Log) : enc_(enc) {}

  /// Register a window's time axis; returns its id. Axes are few (one per SWATH window) and small.
  std::uint32_t addAxis(std::vector<double> rt)
  {
    if (!std::is_sorted(rt.begin(), rt.end()))
    {
      throw std::invalid_argument("odia::ChromStore: axis retention times must be ascending");
    }
    axes_.push_back(ChromAxis{std::move(rt)});
    return static_cast<std::uint32_t>(axes_.size() - 1);
  }

  const ChromAxis& axis(std::uint32_t id) const { return axes_.at(id); }
  std::size_t axisCount() const { return axes_.size(); }

  /// Append a chromatogram covering axis points [start, start + intensity.size()).
  ///
  /// The caller passes the intensities only -- the times are already on the axis, which is the
  /// entire saving. A mismatch between `start + n` and the axis length is a programming error and
  /// throws rather than silently storing a chromatogram that indexes off the end.
  std::size_t add(std::uint32_t axis_id, std::uint32_t start, const std::vector<float>& intensity)
  {
    const ChromAxis& a = axes_.at(axis_id);
    if (static_cast<std::size_t>(start) + intensity.size() > a.rt.size())
    {
      throw std::out_of_range("odia::ChromStore::add: slice extends past the axis");
    }
    ChromRef r;
    r.axis = axis_id;
    r.start = start;
    r.length = static_cast<std::uint32_t>(intensity.size());
    r.offset = blob_.size();

    if (enc_ == ChromEncoding::Float32)
    {
      const auto* p = reinterpret_cast<const std::uint8_t*>(intensity.data());
      blob_.insert(blob_.end(), p, p + intensity.size() * sizeof(float));
      refs_.push_back(r);
      return refs_.size() - 1;
    }

    float mx = 0.0F;
    for (const float v : intensity) { if (std::isfinite(v) && v > mx) { mx = v; } }
    r.scale = mx;
    // A flat-zero chromatogram keeps scale 0 and stores zero codes; decode then returns zeros,
    // which is the truth. Dividing by a zero scale would produce NaN for a case that is common
    // (a transition with no signal in the window) and entirely legitimate.
    const bool have = (mx > 0.0F);

    if (enc_ == ChromEncoding::Quantised16)
    {
      const double inv = have ? (65535.0 / static_cast<double>(mx)) : 0.0;
      blob_.reserve(blob_.size() + intensity.size() * 2);
      for (const float v : intensity)
      {
        const double c = (std::isfinite(v) && v > 0.0F) ? std::lround(static_cast<double>(v) * inv) : 0.0;
        const std::uint16_t q = static_cast<std::uint16_t>(std::min(65535.0, std::max(0.0, c)));
        blob_.push_back(static_cast<std::uint8_t>(q & 0xFF));
        blob_.push_back(static_cast<std::uint8_t>(q >> 8));
      }
    }
    else if (enc_ == ChromEncoding::Quantised8Lin)
    {
      const double inv = have ? (255.0 / static_cast<double>(mx)) : 0.0;
      blob_.reserve(blob_.size() + intensity.size());
      for (const float v : intensity)
      {
        const double c = (std::isfinite(v) && v > 0.0F) ? std::lround(static_cast<double>(v) * inv) : 0.0;
        blob_.push_back(static_cast<std::uint8_t>(std::min(255.0, std::max(0.0, c))));
      }
    }
    else   // Quantised8Log
    {
      // Code 0 is reserved for EXACT zero -- a transition with no signal in a cycle is common and
      // must decode back to zero, not to the bottom of the range. Codes 1..255 are log-spaced from
      // scale/kChromLogDynamicRange up to scale, so relative precision is constant everywhere
      // rather than concentrated at the apex where it is least needed.
      const double floor_v = have ? static_cast<double>(mx) / kChromLogDynamicRange : 0.0;
      const double lr = have ? std::log(static_cast<double>(mx) / floor_v) : 1.0;
      blob_.reserve(blob_.size() + intensity.size());
      for (const float v : intensity)
      {
        std::uint8_t code = 0;
        if (have && std::isfinite(v) && v > 0.0F)
        {
          const double vv = std::max(floor_v, std::min(static_cast<double>(mx), static_cast<double>(v)));
          const double t = std::log(vv / floor_v) / lr;                     // 0..1
          code = static_cast<std::uint8_t>(1 + std::lround(t * 254.0));
        }
        blob_.push_back(code);
      }
    }
    refs_.push_back(r);
    return refs_.size() - 1;
  }

  std::size_t size() const { return refs_.size(); }
  const ChromRef& ref(std::size_t i) const { return refs_.at(i); }

  /// Decode chromatogram @p i into @p rt / @p intensity.
  void get(std::size_t i, std::vector<double>& rt, std::vector<float>& intensity) const
  {
    getSlice(i, 0, refs_.at(i).length, rt, intensity);
  }

  /// Decode a SUB-RANGE, [from, from + n) within the chromatogram.
  ///
  /// This is the operation the whole store exists for: narrowing a chromatogram's RT window costs
  /// two integers and copies only what is asked for. Nothing is re-read and nothing is re-decoded
  /// outside the requested range.
  void getSlice(std::size_t i, std::uint32_t from, std::uint32_t n,
                std::vector<double>& rt, std::vector<float>& intensity) const
  {
    const ChromRef& r = refs_.at(i);
    if (from > r.length) { from = r.length; }
    if (from + n > r.length) { n = r.length - from; }
    const ChromAxis& a = axes_.at(r.axis);
    rt.resize(n);
    intensity.resize(n);
    for (std::uint32_t k = 0; k < n; ++k) { rt[k] = a.rt[r.start + from + k]; }

    if (enc_ == ChromEncoding::Float32)
    {
      const auto* src = reinterpret_cast<const float*>(blob_.data() + r.offset);
      for (std::uint32_t k = 0; k < n; ++k) { intensity[k] = src[from + k]; }
    }
    else if (enc_ == ChromEncoding::Quantised16)
    {
      const double s = static_cast<double>(r.scale) / 65535.0;
      const std::uint8_t* src = blob_.data() + r.offset + 2ULL * from;
      for (std::uint32_t k = 0; k < n; ++k)
      {
        const std::uint16_t q = static_cast<std::uint16_t>(src[2 * k] | (src[2 * k + 1] << 8));
        intensity[k] = static_cast<float>(static_cast<double>(q) * s);
      }
    }
    else if (enc_ == ChromEncoding::Quantised8Lin)
    {
      const double s = static_cast<double>(r.scale) / 255.0;
      const std::uint8_t* src = blob_.data() + r.offset + from;
      for (std::uint32_t k = 0; k < n; ++k)
      {
        intensity[k] = static_cast<float>(static_cast<double>(src[k]) * s);
      }
    }
    else   // Quantised8Log
    {
      const std::uint8_t* src = blob_.data() + r.offset + from;
      const double mx = static_cast<double>(r.scale);
      const double floor_v = (mx > 0.0) ? mx / kChromLogDynamicRange : 0.0;
      const double lr = (mx > 0.0) ? std::log(mx / floor_v) : 1.0;
      for (std::uint32_t k = 0; k < n; ++k)
      {
        const std::uint8_t c = src[k];
        intensity[k] = (c == 0) ? 0.0F
                     : static_cast<float>(floor_v * std::exp(lr * (static_cast<double>(c) - 1.0) / 254.0));
      }
    }
  }

  /// Slice by RETENTION TIME, inclusive on both ends, clamped to what this chromatogram covers.
  /// Returns the number of points written.
  std::uint32_t getRtRange(std::size_t i, double rt_low, double rt_high,
                           std::vector<double>& rt, std::vector<float>& intensity) const
  {
    const ChromRef& r = refs_.at(i);
    const ChromAxis& a = axes_.at(r.axis);
    const std::size_t lo = std::max<std::size_t>(a.lower(rt_low), r.start);
    const std::size_t hi = std::min<std::size_t>(a.upper(rt_high), r.start + r.length);
    if (hi <= lo) { rt.clear(); intensity.clear(); return 0; }
    const std::uint32_t from = static_cast<std::uint32_t>(lo - r.start);
    const std::uint32_t n = static_cast<std::uint32_t>(hi - lo);
    getSlice(i, from, n, rt, intensity);
    return n;
  }

  /// Bytes held, so a caller can report the saving rather than claim it.
  std::size_t bytes() const
  {
    std::size_t b = blob_.capacity() + refs_.capacity() * sizeof(ChromRef);
    for (const auto& a : axes_) { b += a.rt.capacity() * sizeof(double); }
    return b;
  }

  /// What the same data would cost as OpenMS ChromatogramPeak (double RT + double intensity).
  std::size_t bytesAsPeaks() const
  {
    std::size_t pts = 0;
    for (const auto& r : refs_) { pts += r.length; }
    return pts * 16;
  }

  /// Worst RELATIVE intensity error introduced by quantisation, over chromatogram @p i, measured
  /// against the values that were stored. Reported rather than assumed -- see the header note on
  /// baseline resolution.
  double quantisationError(std::size_t i, const std::vector<float>& original) const
  {
    std::vector<double> rt;
    std::vector<float> got;
    get(i, rt, got);
    double worst = 0.0;
    const std::size_t n = std::min(original.size(), got.size());
    for (std::size_t k = 0; k < n; ++k)
    {
      const double o = original[k];
      if (!(std::abs(o) > 0.0)) { continue; }
      worst = std::max(worst, std::abs(static_cast<double>(got[k]) - o) / std::abs(o));
    }
    return worst;
  }

private:
  ChromEncoding enc_;
  std::vector<ChromAxis> axes_;
  std::vector<ChromRef> refs_;
  std::vector<std::uint8_t> blob_;
};

} // namespace odia
