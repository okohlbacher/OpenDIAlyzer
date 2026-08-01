# Chromatogram extraction and scoring — where the time goes, and what can be done about it

An adversarial read of the OpenMS code ODIA actually executes for mass-trace extraction and
chromatogram analysis. Every number is measured on the Astral plasma benchmark unless labelled.

## The code that runs

| stage | file | lines |
|---|---|---:|
| m/z-window integration per spectrum | `ChromatogramExtractorAlgorithm.cpp` | 415 |
| batch orchestration | `OpenSwathWorkflow.cpp` | ~1700 |
| peak picking + smoothing per chromatogram | `PeakPickerChromatogram.cpp` | 409 |
| co-elution / shape sub-scores | `MRMScoring.cpp` | 836 |
| cross-correlation primitives | `openswathalgo/ALGO/Scoring.cpp` | 229 |
| DIA-specific (isotope, massdev) | `DIAScoring.cpp` | 590 |

## F1 — Invert the extraction loop: 12.6x less inner-loop work

**Measured inputs.** MS2 spectra on this instrument carry a median of **2,465 peaks**
(mean 2,360, min 220, max 4,237). One SWATH window's share of the library is
4,463,919 / 150 = **29,759 transitions**.

The current extractor is *already* good: `extract_value_tophat` advances a single monotonic
iterator (`mz_it`) across the spectrum's m/z array while transitions are visited in sorted m/z
order, so the sweep costs O(P) **in total**, not per transition. Per spectrum:

```
cost = O(P + T) = 2,360 + 29,759 = 32,119
```

But look at what those 29,759 per-transition visits find. Transition density in a SWATH window is
29,759 / ~1800 Da ≈ 16.5 per Da; a 10 ppm window at m/z 500 is 5.0 mDa wide. So a peak lands inside
about 0.08 transition windows, and:

```
expected matches per spectrum = P x density x width = 195
=> 99.3% of the per-transition visits find NOTHING
```

Indexing the transitions by m/z (uniform bins → O(1) lookup; the "hash" idea) and iterating over
**peaks** instead of transitions gives

```
cost = O(P + matches) = 2,360 + 195 = 2,555      12.6x less
whole run: 18.7e9 -> 1.5e9 inner iterations
```

**The catch, and it is the interesting part.** A chromatogram must have one point per spectrum in
its RT range, *including zeros*. Today that is `push_back(rt); push_back(0.0)` per transition — an
O(T) cost that inversion cannot remove. It only disappears if the destination is a **dense,
preallocated, zero-initialised array**, where "no peak" means "write nothing".

That is exactly `src/odia_chromatogram.h`. So the loop inversion and the compact store are not two
independent optimisations — **the store is the prerequisite that makes the inversion pay.** Neither
is worth much alone; together they remove 12.6x of the iteration and 4.09x of the memory.

## F2 — The cross-correlation matrix is materialised and thrown away

`MRMScoring::initializeXCorrMatrix` builds, for each of the N(N+1)/2 transition pairs, a full
`XCorrArrayType` (a `vector<pair<int,double>>` of 2·min(10, len)+1 ≈ 21 entries), then reads it
**once, on the very next line**, to take its maximum:

```cpp
xcorr_matrix_(i, j) = Scoring::normalizedCrossCorrelationPost(...);
auto x = Scoring::xcorrArrayGetMaxPeak(xcorr_matrix_(i, j));
xcorr_matrix_max_peak_(i, j)     = std::abs(x->first);
xcorr_matrix_max_peak_sec_(i, j) = x->second;
```

Every score that consumes this (`MRMScoring.cpp:394-466`) uses **only the two scalar matrices**. At
~10.5 transitions per precursor that is ~55 heap allocations per scored peak group, and 2.07M peak
groups were scored — order **10^8 allocations** for data that is discarded.

`Scoring::normalizedCrossCorrelationMaxPost()` is implemented (`Scoring.h/.cpp`) and computes the
same maximum in the same loop with no allocation, preserving the first-wins tie rule so the reported
delay is identical on ties.

**It is NOT enabled, and a correction is owed.** I first reported that `getXCorrMatrix()` "has no
consumer anywhere in the codebase". That was wrong, and the error was mine: the grep that
established it was

```
grep -rn "getXCorrMatrix" src/ ... | grep -v "MRMScoring"
```

whose `grep -v` filtered out **`MRMScoring_test.cpp`** — the one consumer. It asserts on
`getXCorrMatrix()(0,0).data.size() == 21` and on the array's values. Switching to the max-only path
made `MRMScoring_test` **SEGFAULT** (empty matrix, indexed) and took `IonMobilityScoring_test` with
it. Both changes are reverted and both tests pass again.

So the accurate statement is: no *production* consumer, but a tested public contract. Enabling the
fast path is an API decision, not a local optimisation — and the wall-time gain is **unmeasured**;
what is quantified is the allocation count, not the seconds.

## F3 — Window bounds recomputed per (spectrum, transition)

`extract_value_tophat` computes

```cpp
left  = mz - mz * mz_extraction_window / 2.0 * 1.0e-6;
right = mz + mz * mz_extraction_window / 2.0 * 1.0e-6;
```

on every call. Both depend only on the transition's m/z and the window setting — invariant across
spectra. At 3,889 spectra × ~29,759 transitions × 150 windows that is ~1.7e10 redundant evaluations
per pass, removable by precomputing `(left, right)` into the coordinates array once. Small per
call; free to fix.

## F4 — Two things NOT worth doing

Recording these because they are the obvious suggestions and both are wrong here.

**FFT for the cross-correlation is slower.** With ~40 points per peak region and
`XCORR_MAX_DELAY = 10`, direct evaluation is 21 × 40 ≈ 840 multiply-adds. An FFT-based correlation
is ~3·n·log₂n ≈ 1,150 plus transform setup and two zero-padded buffers. Direct wins, and the Eigen
`Map`-based dot product already vectorises it.

**Do not "fix" the extraction sweep.** It looks like a naive nested loop and is not: the monotonic
iterator makes it O(P + T) total, which is why `extractChromatograms` throws if the coordinates are
not sorted by m/z. That precondition is load-bearing.

## Ranked

| # | change | effect | state |
|---|---|---|---|
| 1 | Invert extraction loop + dense store | **12.6x** fewer inner iterations | store built & tested; inversion not written |
| 2 | Exact `reserve()` in the extractor | 42% capacity slack ≈ 156 GB; ~2.7M reallocs/batch | **done**, 15/15 extractor tests pass |
| 3 | Free `chrom_list` after conversion | removes a full duplicate of every chromatogram | **done**, under A/B |
| 4 | Max-only cross-correlation | ~1e8 allocations | primitive **done**, not enabled (tested contract) |
| 5 | Hoist `left`/`right` | ~1.7e10 redundant evaluations | not done |

Items 2 and 3 are already in the binary. Item 1 is the largest remaining, and it is blocked on
integrating the store into `ChromatogramExtractorAlgorithm` rather than on any new insight.
