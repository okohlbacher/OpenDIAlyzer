# ODIA parallel efficiency — what is actually limiting it

Measured on the e2e run: 423,079 precursors, 4,463,919 transitions, 224 threads,
**2,437 s wall, 70,522 s CPU → 28.9 cores average of 224 requested (13%)**.

## 1. The ceiling is the serial tail, not the extraction

| phase | wall | share |
|---|---:|---:|
| extraction pass 1 (1435 s window) | 887 s | 36% |
| extraction pass 2 (240 s window) | 266 s | 11% |
| **everything else** | **1,284 s** | **53%** |

Extraction is under half the run. So:

```
extraction  2x faster -> 1,860 s   speedup 1.31x
extraction  4x faster -> 1,572 s   speedup 1.55x
extraction  8x faster -> 1,428 s   speedup 1.71x
extraction -> ZERO    -> 1,284 s   speedup 1.90x   <-- HARD CEILING
```

**Making extraction infinitely fast buys 1.90x.** Any effort spent on extraction parallelism alone
is bounded by that number. The non-extraction 1,284 s is where the headroom is, and it is largely
single-threaded:

- **parquet library load** — the reader is a row-by-row loop over 78.5M transitions, one thread
- **`remapFeaturePrecursorIds_`** — a SQLite hash-join rebuild over 2.07M features, run **twice**
  (once per pass), single-threaded
- **the semi-supervised LDA** — `src/odia_lda.h` contains no threading whatsoever; measured 142 s
  standalone on this feature count, and its folds are embarrassingly parallel
- calibration, OSW writing, feature-map handling

## 2. Extraction itself has only 150-way parallelism on a 224-core machine

The parallel unit is the SWATH window. There are 150 of them, so even with perfect load balance the
occupancy ceiling is **150/224 = 67%**. Measured: 41.5 cores, 18.5%.

The gap between 150 and 41.5 is load imbalance. Binning the searched precursors by m/z into 150
equal-width bins:

```
min 128    median 1,278    max 13,304 precursors per bin
max/median = 10.4x ;  top 10 bins hold 27.8% of all precursors
```

*Caveat, stated because it matters:* equal-width bins are a **proxy**. Astral's real scheme uses
variable-width windows chosen specifically to equalise ion load, so the true imbalance is smaller
than 10.4x. What the number establishes is that precursor density varies by an order of magnitude
across the m/z range, which is why `schedule(dynamic,1)` is already used — and why the *last* window
still sets the tail.

### 2.1 The nested inner parallelism is switched off

`OpenSwathWorkflow.cpp` offers two schedulers:

```cpp
use_swath_range_scheduler = batchSize <= 0 && load_into_memory && !ms1_only && !nested_requested;
```

ODIA satisfies `batchSize <= 0` but **not `load_into_memory`** — `-readOptions` defaults to `normal`
(streaming), deliberately, to match OpenSwathWorkflow. So the wave scheduler is never used and the
legacy path runs: one `#pragma omp parallel for schedule(dynamic,1)` over the 150 windows.

Inside that, the batch loop *can* nest (line 1121), but the team size is

```cpp
omp_set_num_threads(std::max(1, total_nr_threads / threads_outer_loop_));
```

and ODIA passes `outer_loop_threads = -1`, intending the wave scheduler. `224 / -1 = -224`, so
`max(1, -224) = 1`: **the nested team is one thread**. (If `MT_ENABLE_NESTED_OPENMP` is not defined
in this build the block is compiled out entirely — same outcome by a different route.)

So ODIA's setting of `-1` selects a scheduler it cannot reach, and lands in a path where that same
`-1` disables the only other source of parallelism. Extraction gets 150 coarse, unequal work units
and nothing finer.

## 3. Serialisation points

| where | what it serialises |
|---|---|
| `OpenSwathWorkflow.cpp:1212` `#pragma omp critical (osw_write_out)` | **all** feature+chromatogram output, across every thread |
| `MRMFeatureFinderScoring.cpp:917` `critical (openswath_feature_output_list)` | per-feature appends to the output list |
| `:1138`, `:634` `critical (osw_write_stdout)` | logging (cheap) |
| `:1221` `critical (progress)` | progress counter (cheap) |

The first is the one that matters: every thread that finishes a batch must queue behind a single
lock to write its features. At 224 threads and 2.07M features this is a real convergence point,
and it sits at the end of the longest phase.

## 4. An alternative explanation not yet ruled out

Low core occupancy during extraction is consistent with load imbalance, but **also** with memory
bandwidth saturation: extraction streams spectra and writes chromatograms, and the run held 367 GB.
A memory-bound phase shows exactly this signature — many idle cores, high wall time — and no amount
of scheduling fixes it. Distinguishing the two needs a run at reduced thread count (if wall time is
flat from 64 to 224 threads, it is bandwidth) or `perf stat` on memory stalls. **Not done.** Do not
assume the imbalance explanation is the whole story.

Note that the chromatogram work in this session cuts memory traffic directly: freeing `chrom_list`
after conversion removes a full duplicate of every chromatogram, and `ChromatogramStore` cuts the
per-point footprint 4x. If the phase is bandwidth-bound, those help throughput as well as memory.

## 5. Ranked levers

| # | lever | est. effect | cost |
|---|---|---|---|
| 1 | **Parallelise the LDA's folds** | ~95 s of 142 s | small; folds are independent by construction |
| 2 | **Run `remapFeaturePrecursorIds_` once, not per pass** | unmeasured, serial SQL over 2.07M rows | needs the pass-1 remap to be reusable |
| 3 | **Parallelise the parquet library reader** | part of ~230-300 s | row-range partitioning; the row loop is independent |
| 4 | Set `outer_loop_threads` > 0 to enable nested batch parallelism | raises the 150 ceiling | changes scheduling; needs measurement |
| 5 | `-readOptions cacheWorkingInMemory` to reach the wave scheduler | unknown | costs full SWATH residency — memory ODIA does not have to spare |
| 6 | Batch the `osw_write_out` critical section | unmeasured | buffer per thread, flush in bulk |

Levers 1-3 attack the 1,284 s serial tail, which is where the 1.90x ceiling lives. Levers 4-6
attack extraction, and are bounded by that ceiling — worth doing only after the tail shrinks.

**The honest summary:** ODIA at 13% occupancy is not primarily an extraction-parallelism problem. It
is a run that spends 53% of its wall clock in phases that were never parallelised at all.
