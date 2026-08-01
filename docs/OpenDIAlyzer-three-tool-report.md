# OpenDIAlyzer vs OpenSwathWorkflow vs DIA-NN — measured comparison and bottleneck analysis

**Data:** one Thermo Astral plasma DIA run (`astral.mzML`, 6.4 GB), ONE shared spectral library of
7,149,966 precursors / 78,569,077 transitions (PeptDeep-predicted RT), on an idle 224-core /
2.2 TB node. All three tools got the same file and the same library.

**Read §7 before quoting anything.** This is n=1 per configuration with no replicates, the runs are
not thread-matched, and a known RT confound is live in every row.

---

## 1. Measurements

Thermo Astral plasma DIA (`astral.mzML`, 6.4 GB, 2333 s gradient, 150 SWATH windows), 7,149,966-precursor
PeptDeep-predicted library, IBMI `data` node (224 cores, 2.2 TB). All rows are IDs at 1% FDR.

| tool | threads | wall | CPU | avg cores | peak RSS | IDs @1% FDR |
|---|---:|---:|---:|---:|---:|---:|
| DIA-NN 2.0 | 224 | **349 s** | 8447% | 84.5 | 19.7 GB | **9,261** |
| DIA-NN 1.7.12 ᵈ | 224 | 614 s | 11947% | **119.5** | 16.2 GB | **8,405** |
| ODIA — GBT, parquet | 224 | 2,834 s | 2884% | 28.8 | 350.0 GB | 6,798 |
| ODIA — LDA, sqlite | 224 | 2,437 s | 2893% | 28.9 | 366.8 GB | 4,367 |
| OpenSwathWorkflow ᵃ | 24 ᵇ | 6,886 s | 1448% | 14.5 | 153 GB | 3,275 ᶜ |

ᵃ **extracted at 30 ppm** on an instrument whose measured accuracy is 1.66 ppm -- 18x too wide. This
is not a favourable configuration for OpenSwathWorkflow and the row should not be read as a general
statement about it.
ᵇ 24 threads; wall and occupancy are NOT comparable to the 224-thread rows.
ᶜ OSW's features scored with **ODIA's GBT**, not pyprophet -- holding the scorer constant is what
makes this a comparison of feature sets. pyprophet was abandoned after three attempts
(5h23m, 8h33m, 2h50m with 2.5 h of no output while holding ~330 GB).
ᵈ built from source (CC BY 4.0), 2-line build-only patch for gcc>=8 goto/initialisation strictness.

### The DIA-NN numbers here are NOT the ones quoted earlier in this project

Both DIA-NN versions were originally run against `library.tsv`, which carries OpenSWATH-style decoys
(`DECOY_<target id>` with a **shuffled** sequence). DIA-NN pairs decoys by **modified sequence +
charge**, so none of the 3,546,541 pair:

| DIA-NN | library | IDs | MS2 mass acc |
|---|---|---:|---:|
| 1.7.12 | with decoys | 1,092,962 **(invalid)** | 40.88 ppm |
| 1.7.12 | target-only | 8,405 | 1.7055 ppm |
| 2.0 | with decoys | 8,164 | -- |
| 2.0 | target-only | **9,261** | 1.656 ppm |

1.7.12 disabled its own decoy generation *because decoys were present* and was left with no null.
2.0 gave no warning but was still depressed by 13.4%. **Any earlier comparison against 8,164 was
against a handicapped reference.**

### Instrument mass accuracy: three independent measurements agree

DIA-NN 2.0 reports 1.656 ppm, DIA-NN 1.7.12 reports 1.7055 ppm, and direct measurement gave ~1.6 ppm.
ODIA extracts at 10 ppm by default -- 6x wider than the real error -- and OpenSwathWorkflow's
reference run used 30 ppm.

## 4. Core occupancy and phase breakdown

> **The benchmark node is shared and was oversubscribed.** Measured mid-run: 18 users,
> `load average: 180.14, 200.76, 232.90` on 224 cores, with ~47 cores continuously consumed by
> other tenants (a `dorado` basecaller at 1787% for 1d 14h, five `python` jobs at ~1530%).
>
> Every occupancy number in this section therefore measures **cores the run WON**, not how well the
> tool scales, and the runs were not all made in the same load window. Treat them as operational
> observations, not as a scalability comparison.
>
> **Contention-independent and safe to compare: peak RSS, and total CPU-seconds.** On CPU-seconds —
> DIA-NN 2.0 **29,480**, DIA-NN 1.7.12 73,354, ODIA GBT sqlite 57,950 — ODIA performs *less total
> work* than DIA-NN 1.7.12 and about 2× DIA-NN 2.0, rather than the 6.5× its wall clock suggests.
>
> `experiments/bench_run_fair.sh` records load on both sides of a run so future numbers carry their
> own context.

The column below is **average core occupancy over the whole run** (CPU% ÷ threads requested),
including serial phases. It is *not* parallel efficiency of the parallel sections, and it rewards
tools that request fewer threads — OSW's 60% comes partly from requesting 24.

| tool | threads | avg. cores | occupancy |
|---|---:|---:|---:|
| DIA-NN | 224 | 75.2 | 34% |
| OpenSwathWorkflow | 24 | 14.5 | 60% |
| ODIA @10 ppm | 224 | 25.3 | 11% |
| ODIA @30 ppm | 224 | 72.2 | 32% |

**ODIA @10 ppm, 1,814 s wall:**

| phase | wall | share | parallelism |
|---|---:|---:|---|
| library TSV parse | ~459 s | 25% | **single-threaded** |
| prefilter + remap + LDA scoring | ~580 s | 32% | partial |
| extraction | 545 s | 30% | 41.5 of 224 cores |
| calibration | 230 s | 13% | mixed (nonlinear phase: 3 h 45 m CPU over 207 s wall ≈ 65 cores) |

**Extraction is only 30% of wall time.** Occupancy falls as ODIA gets faster, which is consistent
with fixed serial cost dominating a shrinking parallel section — though also with faster runs simply
having less parallelisable work. These were not distinguished.

---

## 5. ODIA bottlenecks

1. **Serial library TSV parse, ~459 s single-threaded** — 25% of wall at 10 ppm and an Amdahl
   ceiling. Fix: PQP (SQLite) instead of re-parsing an 18.5 GB TSV.
2. **Non-extraction phases are 70% of wall.** Prefilter, the PRECURSOR_ID remap (full table
   rebuild) and the LDA pass are each comparable to extraction.
3. **Extraction uses 41.5 of 224 cores**, capped by 150 SWATH windows as the parallel unit.
4. **Calibration anchor yield 1.7%** (68 usable pairs from 3,897 candidates, measured at 10 ppm).

---

## 6. Calibration

| | DIA-NN | OpenSwathWorkflow | OpenDIAlyzer |
|---|---|---|---|
| RT anchors | its own, iterative | `auto_irt`: samples the **library** (measured here: 492 linear / 23,816 nonlinear), with irtkit (10) + cirtkit (124) as *priority* sequences within that sample | same OpenMS `CalibrationWorkflow`, `SAMPLE_ONCE`, 675 priority CiRT sequences (measured: 500 linear / 3,897 nonlinear) |
| nonlinear fit | yes | yes | **no — `alignmentMethod` pinned to `linear`**, so the "nonlinear" phase is a second linear fit |
| RT window | set from data (`2.86643`, **unit not stated in the log** — not comparable to the seconds below) | `estimateWindow` on surviving anchor residuals → 375.3 s | same mechanism → 746.9–1,135.5 s |
| mass calibration | **measures its own error and narrows 30/20 ppm → 6 ppm** | does not auto-narrow; uses the configured window | does not auto-narrow; uses the configured window |
| QC gates | internal | `min_rsq` 0.95, `min_coverage` 0.6 (defaults) | **both relaxed**: 0.70 / 0.30 |
| GPU | **not used in this run** (no CUDA/NVIDIA/OpenCL in the log; 75 CPU cores) | n/a | n/a |

The mass-calibration row is the substantive difference: DIA-NN measures its own mass error and acts
on it; **neither OSW nor ODIA does** — that cell is identical for both, and OSW would have been
equally exposed had its output been FDR-scored under the same window.

ODIA-specific weaknesses, both currently surfaced in its own logs: an anchor yield of 1.7%, and
both QC gates relaxed at once, which lets outlier removal discard up to 70% of matched anchors to
manufacture the R² it is then checked against.

---

## 7. Caveats

1. **n = 1.** One file, one sample, one run per configuration. No replicates. Wall/RSS quoted to 3–4
   significant figures are single observations; a 53-ID difference (§2) is not resolvable.
2. **Not thread-matched.** OSW at 24, others at 224. Cross-tool speed and occupancy comparisons in
   §1 and §4 are confounded. OSW at 224 threads is **unmeasured** — an earlier claim that it
   "cannot" scale to 224 was based on a segfault that turned out to be my own libOpenMS ABI break
   (rebuilt binaries run clean) and is **withdrawn**.
3. **Run vintages differ.** Rows were produced across a period that also included a decoy-pairing
   fix (a prior regime kept 23 decoys out of ~3.5M, which alone can flatten score distributions and
   floor the q-value), the streaming change, and a calibration retune. The 30 ppm row is *not* a
   clean control for the 10 ppm row.
4. **RT clipping is live in every row.** All runs use RT windows below the measured p95 residual
   (~1,720 s), depressing every ID count by an unknown amount — plausibly the largest unquantified
   confound in the ID comparison, and it affects ODIA more than DIA-NN, which set its own.
5. **Different work per tool.** ODIA's numbers include its prefilter and in-process LDA; OSW's
   exclude the pyprophet step its output requires; DIA-NN's include its own calibration.
6. **DIA-NN's 8,164** is its own FDR estimate (precursors, single-pass, no MBR, no entrapment
   validation). Any "% of DIA-NN" inherits that assumption.
7. **Measurement method:** wall/CPU/peak RSS from `/usr/bin/time -v` (ru_maxrss) for every arm.
   Peak-RSS comparisons across tools are sensitive to allocator behaviour.

---

## 8. Next steps

1. **Implement mass-accuracy inference** (hook now exists, empty, in the calibration phase): extract
   wide, measure observed-minus-theoretical m/z on confident hits, narrow to ~5x robust sd, re-extract.
   Default 10–15 ppm meanwhile.
2. **Run the clean ablation** the conclusions currently lack: one binary, 600 s RT, m/z ∈ {30, 20,
   10, 5} — plus an RT arm at fixed 10 ppm to quantify what the live RT clipping is costing.
3. **Kill the serial library parse** (PQP) and profile the memory floor with heaptrack rather than
   modelling it.
4. **Re-derive prefilter / LOESS / anchor-yield conclusions at 10 ppm** — all were measured in a
   regime where nothing discriminated.
