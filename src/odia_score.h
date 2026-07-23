// odia-score -- the DIA scoring kernel, CPU reference implementation.
//
// This is the hot path a future CUDA port targets: for a candidate's F fragment
// traces over T retention-time points, reduce the F x T matrix to a score
// vector. The expensive term is the all-pairs cross-correlation (O(F^2 * D * T),
// D = 2*maxdelay+1), which is what OpenSWATH's MRMScoring computes and what
// carries enough arithmetic intensity (36-110 flop/byte FP32) to suit an H100.
//
// Correctness is defined by OpenSWATH's own algorithms: every function here is
// checked against the exact reference values in
// ext/OpenMS/src/tests/class_tests/openswathalgo/Scoring_test.cpp, so the CUDA
// port has a bit-comparable oracle. Idea and reference values from OpenMS
// (BSD-3); implemented independently in a layout suited to GPU porting (SoA,
// dense, no per-pair heap allocation).
//
// Header-only so the same source compiles for host and, later, device.

#pragma once

#include <cmath>
#include <cstddef>
#include <vector>

namespace odia
{

// Standardize in place: subtract mean, divide by *population* standard deviation
// (÷n, matching OpenSWATH's standardize_data, not the ÷(n-1) sample form).
inline void standardize(std::vector<double>& x)
{
  const size_t n = x.size();
  if (n == 0) return;
  double mean = 0.0;
  for (double v : x) mean += v;
  mean /= n;
  double var = 0.0;
  for (double v : x) var += (v - mean) * (v - mean);
  var /= n;                                  // population variance
  const double sd = std::sqrt(var);
  if (sd == 0.0) { for (double& v : x) v -= mean; return; }
  for (double& v : x) v = (v - mean) / sd;
}

// Cross-correlation at each integer delay in [-maxdelay, +maxdelay], normalised
// by 1/n. Inputs are assumed already standardized (the "...Post" variant in
// OpenSWATH). Returns one value per delay; the delay for slot k is k-maxdelay.
//
// GPU-shaped: output is a dense array indexed by delay, not a vector of
// (delay,value) pairs. One thread per delay is the natural device mapping.
inline std::vector<double> xcorr_post(const std::vector<double>& a,
                                      const std::vector<double>& b, int maxdelay)
{
  const int n = static_cast<int>(a.size());
  std::vector<double> out(2 * maxdelay + 1, 0.0);
  for (int delay = -maxdelay; delay <= maxdelay; ++delay)
  {
    // overlap of a[i] with b[i+delay]
    const int lo = delay < 0 ? -delay : 0;
    const int hi = delay < 0 ? n : n - delay;
    double s = 0.0;
    for (int i = lo; i < hi; ++i) s += a[i] * b[i + delay];
    out[delay + maxdelay] = s / n;
  }
  return out;
}

// Convenience: standardize copies, then cross-correlate. Mirrors
// normalizedCrossCorrelation.
inline std::vector<double> xcorr(std::vector<double> a, std::vector<double> b, int maxdelay)
{
  standardize(a);
  standardize(b);
  return xcorr_post(a, b, maxdelay);
}

// Largest |value| in a delay array, and its delay. This is the co-elution
// score OpenSWATH extracts from the xcorr array (xcorrArrayGetMaxPeak).
struct XCorrPeak { int delay; double value; };
inline XCorrPeak xcorr_max(const std::vector<double>& xc, int maxdelay)
{
  XCorrPeak best{0, 0.0};
  double bestabs = -1.0;
  for (size_t k = 0; k < xc.size(); ++k)
    if (std::abs(xc[k]) > bestabs) { bestabs = std::abs(xc[k]); best = {int(k) - maxdelay, xc[k]}; }
  return best;
}

// Normalised Manhattan distance (mQuest delta_ratio_sum): each vector scaled to
// unit mean, mean absolute difference. Matches NormalizedManhattanDist.
inline double normalized_manhattan(std::vector<double> x, std::vector<double> y)
{
  const size_t n = x.size();
  if (n == 0) return 0.0;
  double mx = 0.0, my = 0.0;
  for (size_t i = 0; i < n; ++i) { mx += x[i]; my += y[i]; }
  mx /= n; my /= n;
  double s = 0.0;
  for (size_t i = 0; i < n; ++i) s += std::abs(x[i] / mx - y[i] / my);
  return s / n;
}

// Root-mean-square deviation. Matches RootMeanSquareDeviation.
inline double rmsd(const std::vector<double>& x, const std::vector<double>& y)
{
  const size_t n = x.size();
  if (n == 0) return 0.0;
  double s = 0.0;
  for (size_t i = 0; i < n; ++i) s += (x[i] - y[i]) * (x[i] - y[i]);
  return std::sqrt(s / n);
}

// Spectral angle = acos(normalised dot product). Matches SpectralAngle.
inline double spectral_angle(const std::vector<double>& x, const std::vector<double>& y)
{
  const size_t n = x.size();
  double dot = 0.0, nx = 0.0, ny = 0.0;
  for (size_t i = 0; i < n; ++i) { dot += x[i] * y[i]; nx += x[i] * x[i]; ny += y[i] * y[i]; }
  const double denom = std::sqrt(nx) * std::sqrt(ny);
  if (denom == 0.0) return 0.0;
  double c = dot / denom;
  if (c > 1.0) c = 1.0; else if (c < -1.0) c = -1.0;
  return std::acos(c);
}

// The all-pairs cross-correlation matrix over F fragment traces -- the term
// that dominates cost and that a CUDA port parallelises. Upper triangle only
// (the matrix is symmetric in |delay|), matching MRMScoring::initializeXCorrMatrix.
// traces: F vectors, each length T, standardized in place. Returns, for each
// unordered pair (i<=j), the max-peak co-elution value.
struct PairScore { int i, j, delay; double value; };
inline std::vector<PairScore> allpairs_xcorr(std::vector<std::vector<double>> traces, int maxdelay)
{
  const int F = static_cast<int>(traces.size());
  for (auto& t : traces) standardize(t);
  std::vector<PairScore> out;
  out.reserve(static_cast<size_t>(F) * (F + 1) / 2);
  for (int i = 0; i < F; ++i)
    for (int j = i; j < F; ++j)
    {
      const auto xc = xcorr_post(traces[i], traces[j], maxdelay);
      const auto pk = xcorr_max(xc, maxdelay);
      out.push_back({i, j, pk.delay, pk.value});
    }
  return out;
}

} // namespace odia
