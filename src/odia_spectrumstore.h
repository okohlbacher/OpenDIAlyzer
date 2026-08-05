#pragma once
// Every spectrum of the run, decoded ONCE, in one contiguous arena.
//
// WHY. The spectra are walked five times per run, not once:
//
//   1. prefilter/scan_targets        every MS1 + MS2 spectrum
//   2. prefilter/scan_decoys         the same spectra again -- same swath_maps, different candidates
//   3. setup/cirt_calibration        anchor extraction over the SWATH windows
//   4. extract_pass1_wide            every spectrum
//   5. extract_pass2_narrow          every spectrum again
//
// mzML hid this: it parses once and leaves everything resident, so passes 2-5 are walks over RAM.
// mzPeak streams honestly, so each pass genuinely re-decodes parquet -- which is why the same
// analysis went from 20:16 on mzML to over 83 minutes still inside the prefilter on mzPeak. The
// format did not get slower; it stopped concealing four redundant passes.
//
// So: decode once into this, and let all five passes read from memory.
//
// THE REPRESENTATION, and why not MSSpectrum. OpenSwath::Spectrum holds two vector<double> --
// 16 B/peak, plus two heap allocations per spectrum, plus MSSpectrum's metadata and DataArrays on
// top. Measured on this run the resident spectra cost ~9.06 GB. Here a peak is a float32 m/z and a
// float32 intensity, 8 B, in ONE arena for the whole run: two allocations total rather than two per
// spectrum, which also keeps ~1.2M allocations out of an allocator whose retained pool is already
// the largest single item at peak.
//
// PRECISION, measured rather than assumed (see the test). float32 carries ~24 bits of mantissa, so
// a m/z near 2000 Th resolves to ~1.2e-4 Th = 0.06 ppm. The narrowest extraction window in use is
// 10 ppm, and the mass calibration reports residuals in whole ppm. 0.06 ppm is ~170x finer than the
// window it feeds and well under the instrument's own accuracy, so the quantisation is not a term
// in any decision this tool makes. Intensities are float32 for the same reason they are float32
// everywhere else in OpenMS -- the detector does not deliver more.
//
// WHAT IS NOT DONE HERE, deliberately: no delta or varint encoding of m/z. The arrays are sorted
// and would compress well, but extraction binary-searches them per transition, and a
// variable-length encoding cannot be searched without either decoding or a second index. Random
// access is the requirement; density is second.

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

namespace odia
{

/// One run's spectra, decoded once. Indexed exactly as the source: index i is spectrum i.
class SpectrumStore
{
public:
  struct Meta
  {
    double rt = 0.0;              ///< seconds
    double precursor_mz = 0.0;    ///< isolation window centre; 0 for MS1
    double iso_lower = 0.0;       ///< isolation window, absolute Th (not offsets)
    double iso_upper = 0.0;
    double drift = -1.0;          ///< ion mobility, or -1
    std::uint64_t offset = 0;     ///< first peak in the arena
    std::uint32_t count = 0;      ///< peaks
    std::uint8_t ms_level = 0;
  };

  void reserve(std::size_t n_spectra, std::size_t n_peaks)
  {
    meta_.reserve(n_spectra);
    mz_.reserve(n_peaks);
    in_.reserve(n_peaks);
  }

  /// Append a spectrum. `mz` must be ascending -- everything downstream binary-searches it, and an
  /// unsorted array silently returns wrong peaks rather than failing.
  ///
  /// `drift` is optional and PER PEAK: PASEF carries an ion mobility per peak, and dropping it
  /// would silently disable IM scoring on exactly the instruments that need it. It costs nothing
  /// when absent -- the array is only grown if a spectrum supplies one.
  std::size_t add(const double* mz, const float* intensity, std::uint32_t n, const Meta& m,
                  const double* drift = nullptr)
  {
    for (std::uint32_t k = 1; k < n; ++k)
    {
      if (mz[k] < mz[k - 1])
      {
        throw std::invalid_argument("odia::SpectrumStore: m/z array is not ascending");
      }
    }
    Meta e = m;
    e.offset = mz_.size();
    e.count = n;
    if (drift && dt_.empty() && !mz_.empty()) { dt_.resize(mz_.size(), -1.0F); }   // back-fill
    for (std::uint32_t k = 0; k < n; ++k)
    {
      mz_.push_back(static_cast<float>(mz[k]));
      in_.push_back(intensity[k]);
    }
    if (!dt_.empty() || drift)
    {
      for (std::uint32_t k = 0; k < n; ++k)
      { dt_.push_back(drift ? static_cast<float>(drift[k]) : -1.0F); }
    }
    meta_.push_back(e);
    return meta_.size() - 1;
  }

  std::size_t size() const { return meta_.size(); }
  std::size_t peakCount() const { return mz_.size(); }
  const Meta& meta(std::size_t i) const { return meta_.at(i); }

  /// Peaks of spectrum @p i, as spans into the arena. No copy, no allocation.
  const float* mz(std::size_t i) const { return mz_.data() + meta_.at(i).offset; }
  const float* intensity(std::size_t i) const { return in_.data() + meta_.at(i).offset; }
  bool hasDrift() const { return !dt_.empty(); }
  const float* drift(std::size_t i) const
  { return dt_.empty() ? nullptr : dt_.data() + meta_.at(i).offset; }
  std::uint32_t count(std::size_t i) const { return meta_.at(i).count; }

  /// Widen into the double arrays the OpenSwath interface requires.
  ///
  /// This is the cost the design accepts: a widening copy per access, in exchange for never
  /// decoding parquet more than once. A copy of contiguous float32 is memory-bandwidth work; a
  /// parquet decode is page reads, decompression and column assembly. Trading the second for the
  /// first four times over is the whole point.
  void widen(std::size_t i, std::vector<double>& mz_out, std::vector<double>& in_out,
             std::vector<double>* dt_out = nullptr) const
  {
    const Meta& e = meta_.at(i);
    mz_out.resize(e.count);
    in_out.resize(e.count);
    const float* m = mz_.data() + e.offset;
    const float* v = in_.data() + e.offset;
    for (std::uint32_t k = 0; k < e.count; ++k) { mz_out[k] = m[k]; in_out[k] = v[k]; }
    if (dt_out)
    {
      dt_out->resize(e.count);
      const float* d = dt_.empty() ? nullptr : dt_.data() + e.offset;
      for (std::uint32_t k = 0; k < e.count; ++k) { (*dt_out)[k] = d ? d[k] : -1.0; }
    }
  }

  /// First peak with m/z >= @p lo, by binary search over the stored float32 values.
  std::uint32_t lowerBound(std::size_t i, double lo) const
  {
    const Meta& e = meta_.at(i);
    const float* b = mz_.data() + e.offset;
    return static_cast<std::uint32_t>(std::lower_bound(b, b + e.count, static_cast<float>(lo)) - b);
  }

  /// Summed intensity in [lo, hi], the operation extraction actually performs.
  double sumRange(std::size_t i, double lo, double hi) const
  {
    const Meta& e = meta_.at(i);
    const float* m = mz_.data() + e.offset;
    const float* v = in_.data() + e.offset;
    double s = 0.0;
    for (std::uint32_t k = lowerBound(i, lo); k < e.count && m[k] <= static_cast<float>(hi); ++k)
    {
      s += v[k];
    }
    return s;
  }

  std::size_t bytes() const
  {
    return (mz_.capacity() + in_.capacity() + dt_.capacity()) * sizeof(float)
         + meta_.capacity() * sizeof(Meta);
  }
  /// What the same peaks would cost as OpenSwath::Spectrum (two vector<double>), ignoring the
  /// two heap allocations per spectrum that this store also removes.
  std::size_t bytesAsDoublePairs() const { return mz_.size() * 2 * sizeof(double); }

  /// Worst relative m/z error the float32 storage can introduce at @p mz, in ppm. Reported rather
  /// than assumed: it has to stay far below the extraction window to be irrelevant.
  static double quantisationPpm(double mz)
  {
    const float f = static_cast<float>(mz);
    const double next = static_cast<double>(std::nextafter(f, std::numeric_limits<float>::max()));
    return mz > 0.0 ? 1e6 * 0.5 * (next - static_cast<double>(f)) / mz : 0.0;
  }

private:
  std::vector<float> mz_, in_, dt_;   // dt_ empty unless a spectrum supplied drift
  std::vector<Meta> meta_;
};

} // namespace odia
