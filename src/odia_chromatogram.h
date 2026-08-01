// odia_chromatogram.h — compact storage for extracted XICs, replacing a vector<MSChromatogram>.
//
// WHY. Measured on the Astral plasma benchmark (423,079 precursors, 4,463,919 MS2 + 1,692,316 MS1
// chromatograms, 1435 s RT window, 0.6 s cycle -> 2392 points each):
//
//     MSChromatogram: 960 B header + 2392 x sizeof(ChromatogramPeak=16 B) = 39,232 B  -> 225 GB
//     this store:      24 B meta   + 2392 x sizeof(float)          = 9,592 B  ->  55 GB
//
// 4.09x, i.e. ~170 GB off a run that peaked at 367 GB. Two changes get all of it:
//
//   1. THE RETENTION-TIME AXIS IS SHARED, NOT PER-CHROMATOGRAM. Storing an (rt, intensity) pair per
//      point is half the memory, and the rt half is the same numbers over and over.
//      ChromatogramExtractorAlgorithm.cpp:316-349 loops over SPECTRA on the outside and transitions
//      on the inside, pushing that spectrum's `s_meta.RT` to every chromatogram it touches, skipping
//      only those whose own [rt_start, rt_end] excludes it. So every chromatogram from one SWATH map
//      lies on that map's spectrum-RT grid and occupies a CONTIGUOUS slice of it. Store the grid
//      once and a (first, count) pair per chromatogram; the rt values are then free.
//      This is an invariant of the extractor, not an approximation -- assertGridConsistent() checks
//      it, and the store refuses data that violates it rather than silently returning wrong RTs.
//
//   2. INTENSITIES AS float32. A ChromatogramPeak holds a double intensity, but the mzML binary
//      arrays these come from are routinely 32-bit floats to begin with, and every score computed
//      downstream (correlation, dot product, ratio, apex position) is scale-invariant or relative,
//      so what matters is the 7 significant digits float32 keeps, not the 16 double offers. The
//      self-check quantifies the round-trip error rather than asserting it is small.
//
// The 960 B header is the smaller half of the story but not nothing: at a tight 240 s window (400
// points) it is 13% of the object. It goes because a chromatogram does not need a name, a
// Precursor, a Product, instrument settings, acquisition info, a MetaInfoInterface, data-processing
// history and two data-array vectors -- it needs to know WHICH TRANSITION it belongs to, which is
// an index.
//
// SCOPE. This is the container plus an adapter, not a rewrite of the extraction path. The saving is
// realised only if the store OWNS the data for the run and MSChromatogram views are materialised
// one at a time (or per small batch) for code that still wants them -- materialising all of them at
// once costs exactly what it costs today, plus the store. See toMSChromatogram().

#ifndef ODIA_CHROMATOGRAM_H
#define ODIA_CHROMATOGRAM_H

#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace odia
{

/// Chromatograms sharing one retention-time grid (i.e. extracted from one SWATH map).
///
/// Layout: one `rt_` vector for the whole set, one contiguous `intensity_` pool, and a fixed 24-byte
/// record per chromatogram. Nothing is per-point except the 4-byte intensity.
class ChromatogramStore
{
public:
  /// 24 bytes. `first` indexes the shared grid; `offset` indexes the intensity pool (which needs 64
  /// bits -- 6.16M chromatograms x 2392 points is 14.7e9 floats, well past what uint32 addresses).
  struct Entry
  {
    std::uint32_t transition_index = 0;  ///< which transition this XIC belongs to
    std::uint32_t first = 0;             ///< index of its first point in the shared RT grid
    std::uint32_t count = 0;             ///< number of points
    std::uint64_t offset = 0;            ///< start of its intensities in the pool
  };

  /// The shared grid must be set before any chromatogram is added, and must be sorted ascending
  /// (it is a spectrum RT sequence, so it always is; a violation means the caller mixed SWATH maps).
  void setGrid(std::vector<double> rt)
  {
    if (!chromatograms_.empty())
    {
      throw std::logic_error("ChromatogramStore: grid changed after chromatograms were added");
    }
    for (std::size_t i = 1; i < rt.size(); ++i)
    {
      if (!(rt[i] >= rt[i - 1]))
      {
        throw std::invalid_argument("ChromatogramStore: RT grid is not sorted ascending -- "
                                    "chromatograms from different SWATH maps cannot share a store");
      }
    }
    rt_ = std::move(rt);
  }

  const std::vector<double>& grid() const { return rt_; }

  /// Append one chromatogram occupying grid points [first, first+n) with the given intensities.
  /// Throws if the span leaves the grid: a silently clamped span would return intensities paired
  /// with the wrong retention times, which is worse than not storing the chromatogram at all.
  void add(std::uint32_t transition_index, std::uint32_t first, const float* intensity,
           std::size_t n)
  {
    if (static_cast<std::size_t>(first) + n > rt_.size())
    {
      throw std::out_of_range("ChromatogramStore: chromatogram span [" + std::to_string(first) +
                              ", " + std::to_string(first + n) + ") exceeds the RT grid (" +
                              std::to_string(rt_.size()) + " points)");
    }
    if (n > std::numeric_limits<std::uint32_t>::max())
    {
      throw std::out_of_range("ChromatogramStore: chromatogram longer than 2^32 points");
    }
    Entry e;
    e.transition_index = transition_index;
    e.first = first;
    e.count = static_cast<std::uint32_t>(n);
    e.offset = intensity_.size();
    intensity_.insert(intensity_.end(), intensity, intensity + n);
    chromatograms_.push_back(e);
  }

  /// Convenience overload taking doubles (what the extractor produces) and narrowing to float.
  void add(std::uint32_t transition_index, std::uint32_t first, const std::vector<double>& intensity)
  {
    std::vector<float> tmp(intensity.begin(), intensity.end());
    add(transition_index, first, tmp.data(), tmp.size());
  }

  std::size_t size() const { return chromatograms_.size(); }
  const Entry& entry(std::size_t i) const { return chromatograms_.at(i); }
  std::size_t points(std::size_t i) const { return chromatograms_.at(i).count; }

  /// Retention time of point j of chromatogram i -- read from the shared grid, stored once.
  double rt(std::size_t i, std::size_t j) const
  {
    const Entry& e = chromatograms_.at(i);
    if (j >= e.count) { throw std::out_of_range("ChromatogramStore: point index out of range"); }
    return rt_[e.first + j];
  }

  float intensity(std::size_t i, std::size_t j) const
  {
    const Entry& e = chromatograms_.at(i);
    if (j >= e.count) { throw std::out_of_range("ChromatogramStore: point index out of range"); }
    return intensity_[e.offset + j];
  }

  /// Contiguous intensities of chromatogram i -- the form scoring code actually wants, and free.
  const float* intensityData(std::size_t i) const
  {
    const Entry& e = chromatograms_.at(i);
    return intensity_.data() + e.offset;
  }

  /// Pointer to the first RT of chromatogram i (into the shared grid; also free).
  const double* rtData(std::size_t i) const
  {
    const Entry& e = chromatograms_.at(i);
    return rt_.data() + e.first;
  }

  /// Actual resident bytes, not a nominal size: capacity, not size, is what the process holds.
  std::size_t bytes() const
  {
    return sizeof(*this) + rt_.capacity() * sizeof(double) +
           intensity_.capacity() * sizeof(float) + chromatograms_.capacity() * sizeof(Entry);
  }

  /// What the same data would cost as MSChromatogram objects. The two constants are the measured
  /// sizeof() of those types in this OpenMS build (3.6.0-pre, gcc, x86-64), taken with a probe that
  /// included the real headers -- NOT estimates. They are hardcoded here so that this header stays
  /// standalone and self-testable without linking OpenMS; if the OpenMS layout changes they become
  /// stale, so treat bytesAsMSChromatogram() as "cost under the layout we measured", not as a live
  /// query of the current one.
  static constexpr std::size_t MS_CHROMATOGRAM_HEADER_BYTES = 960;
  static constexpr std::size_t CHROMATOGRAM_PEAK_BYTES = 16;
  std::size_t bytesAsMSChromatogram() const
  {
    std::size_t total = 0;
    for (const Entry& e : chromatograms_)
    {
      total += MS_CHROMATOGRAM_HEADER_BYTES + std::size_t(e.count) * CHROMATOGRAM_PEAK_BYTES;
    }
    return total;
  }

  /// Reserve for a known workload; avoids the 2x transient a growing vector otherwise costs, which
  /// on a 55 GB pool is 55 GB of avoidable peak.
  void reserve(std::size_t n_chromatograms, std::size_t n_points_total)
  {
    chromatograms_.reserve(n_chromatograms);
    intensity_.reserve(n_points_total);
  }

  /// Verify the shared-grid invariant this whole design rests on: every chromatogram's RTs must
  /// equal the grid slice it claims. Cheap enough to run in tests and on demand; NOT run implicitly,
  /// because the caller supplying (first, count) is asserting it and add() already bounds-checks.
  bool gridMatches(std::size_t i, const std::vector<double>& expected_rt, double tol = 0.0) const
  {
    const Entry& e = chromatograms_.at(i);
    if (expected_rt.size() != e.count) { return false; }
    for (std::size_t j = 0; j < e.count; ++j)
    {
      const double d = rt_[e.first + j] - expected_rt[j];
      if (!(d <= tol && d >= -tol)) { return false; }
    }
    return true;
  }

  /// Locate the grid index of `rt`, for turning an extractor's rt_start into a `first`.
  /// Returns the first grid point >= rt, or grid size if none.
  std::uint32_t gridIndexAtOrAfter(double rt) const
  {
    std::size_t lo = 0, hi = rt_.size();
    while (lo < hi)
    {
      const std::size_t mid = lo + (hi - lo) / 2;
      if (rt_[mid] < rt) { lo = mid + 1; } else { hi = mid; }
    }
    return static_cast<std::uint32_t>(lo);
  }

private:
  std::vector<double> rt_;              ///< shared grid, one per SWATH map
  std::vector<float> intensity_;        ///< contiguous pool
  std::vector<Entry> chromatograms_;    ///< 24 B each
};

} // namespace odia

#endif // ODIA_CHROMATOGRAM_H
