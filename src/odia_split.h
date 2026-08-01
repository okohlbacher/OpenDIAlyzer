#pragma once
// Extraction thread/batch split.
//
// OpenSwathWorkflow's legacy path sizes its inner batch team
//   omp_set_num_threads(std::max(1, total_nr_threads / threads_outer_loop_))
// so passing -1 yields max(1, 224/-1) == 1 and the inner loop runs SERIAL. All parallelism then
// comes from the outer SWATH-window loop, and with the benchmark's measured max/mean work ratio of
// 3.95 over 150 windows that averages 150/3.95 = 38 concurrent -- against 38.2 measured.
//
// The two knobs only work together. Measured on the Astral benchmark (423079 precursors, 150 MS2
// windows, 224 threads):
//
//   outer=-1, batch auto (10000) : inner 1 thread, 1 batch/window  ->  38.2 cores, peak 187 GB
//   outer=-1, batch 2000         : inner 1 thread                  ->  43.0 cores, 26% slower
//   outer=16, batch auto         : inner 14 threads, 1 batch/window->  16.4 cores (worse)
//   outer=16, batch 200          : inner 14 threads, ~14 batches   -> 105.6 cores, peak 88 GB
//
// So the derivation must supply BOTH inner threads and inner iterations to occupy them.

#include <algorithm>
#include <cstddef>

namespace odia
{

struct ExtractionSplit
{
  int outer_threads;   ///< SWATH windows processed concurrently
  int inner_threads;   ///< threads per window, = threads / outer_threads
  int batch_size;      ///< compounds per inner batch; aims for ~inner_threads batches per window
  std::size_t batches_per_window;
};

/// Derive from the thread budget and the actual data shape. Aims for ~14 inner threads: enough to
/// fill a team, few enough that a normal per-window compound count still divides into that many
/// batches.
inline ExtractionSplit deriveExtractionSplit(int threads, std::size_t n_ms2_windows,
                                             std::size_t n_compounds)
{
  ExtractionSplit s{};
  if (threads < 1) { threads = 1; }
  if (n_ms2_windows == 0) { n_ms2_windows = 1; }

  // Clamp INSIDE n_ms2_windows. max(2, min(n, ...)) would return 2 for a single-window run,
  // requesting more outer threads than there are iterations to give them.
  s.outer_threads = std::min<int>(static_cast<int>(n_ms2_windows), std::max(1, threads / 14));
  s.inner_threads = std::max(1, threads / s.outer_threads);

  const std::size_t per_window = (n_compounds + n_ms2_windows - 1) / n_ms2_windows;
  const std::size_t inner = static_cast<std::size_t>(s.inner_threads);
  // Floor of 50: below that the per-batch setup cost dominates the work in the batch.
  s.batch_size = static_cast<int>(std::max<std::size_t>(50, (per_window + inner - 1) / inner));
  s.batches_per_window = per_window / static_cast<std::size_t>(std::max(1, s.batch_size));
  return s;
}

} // namespace odia
