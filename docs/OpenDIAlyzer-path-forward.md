# Path forward — sequenced by what can be attributed, not by what looks promising

Written after adversarial reviews of a night's work by two external models (kimi, codex). Both
independently predicted the same failure for the plan as originally stated ("calibration → memory →
occupancy"):

> *optimize a proxy, change the system, then instrument too late to attribute the result.*

Every conclusion reached by inference in that night was wrong or confounded; every conclusion
reached by instrumentation held. The ordering below follows from that, not from expected payoff.

## Standing rules, derived from specific failures

| Rule | The failure it prevents |
|---|---|
| **Measure before intervening, in the same build.** | The prefilter rewrite and its timers landed together, so it has no before-number, permanently. |
| **One variable per run.** | The sqlite-vs-parquet comparison used different extraction runs; the variable was never isolated. |
| **Count the null, not just the benefit.** | 20 ppm screening predicted +800 IDs and delivered −327: reachability was counted, the FDR burden of 1M extra candidates was not. |
| **A gate that passes is not a result.** | Peakedness, inferred ppm and anchor count are internal diagnostics. IDs at controlled FDP are the result. |
| **Nothing under the noise floor is signal.** | Three equivalent configurations gave 6,522 / 6,607 / 6,622. Differences below ~100 IDs have been over-read repeatedly. |
| **Never reset a pre-existing reference checkout.** | `git checkout -- .` in the gitignored OpenMS tree destroyed load-bearing uncommitted work. `setup_node.sh` automated the same act until it was gated. |

## Phase 0 — establish what the numbers mean (blocking, cheap)

Nothing below is interpretable without these.

1. **Noise floor.** Three repeats of one fixed configuration on the same binary. Report mean and
   spread of IDs. Every later claim is measured against this, not against a single run.
2. **Phase-resolved memory.** `PhaseTimer` now reports RSS at phase entry/exit from
   `/proc/self/statm`. One run identifies which phase actually holds peak RSS.
   **This is required before any memory work**: the §2 decomposition of
   `OpenDIAlyzer-memory-review.md` was refuted by its own §6 A/B (a 62% cut of the chromatogram term
   moved peak by 4.6%), and the dominant term is still unidentified.
3. **An exclusive allocation, or `taskset`.** The benchmark node has 18 users and a load average of
   180–233. Until this exists, wall clock and core occupancy measure cores *won*, and no scalability
   claim can be made.

## Phase 1 — quality (IDs), in dependency order

**1.1 Validate the mass calibration that now fires.** The RT-direction fix made it fit for the first
time: MS2 5,709 anchors at peakedness 4.9, MS1 560 at 5.5, MS1 sigma 1.5 ppm against an instrument
independently measured at 1.66–1.71 ppm. That corroborates the *measurement*. It does not establish
that the identifications gained are real: anchors are the same run's top-2,000 targets ranked by the
score the narrowing then improves — selection on the optimised quantity. Required before claiming a
win: **entrapment FDP** (`-entrapment_tag`, already implemented) and cross-fitted anchors.

**1.2 `-ms1_scores` A/B.** 12 MS1 sub-scores are withheld from the classifier. The measurement that
justified excluding them ("35-36 sub-scores gave 4,913 vs 29 giving 6,607") is now explained by the
`library_rt` defect — the same path reaches 6,522 with no MS1 scores at all. The exclusion rests on
a confounded comparison and must be re-measured, not assumed either way.

**1.3 MS2 anchor quality.** MS2 sigma came out at 4.5 ppm on a 1.66 ppm instrument, so the MS2
anchor picking still admits interference, and its 13.4 ppm estimate was correctly rejected by the
never-widen guard. Note the trap: widening the anchor **search** is not the same as accepting a
wider **estimate**. Conflating them produces "fixes that change inputs but not the symptom" — which
is what three earlier mass-calibration fixes did.

**1.4 Named RT transform types.** Two directional defects in one subsystem in one session
(`library_rt`; `native_trafo` at *two* call sites). `LibraryToRun` / `RunToLibrary` as distinct
types makes the class unrepresentable rather than fixed case by case.

**1.5 Library fragment intensities may be predicted for the wrong collision energy.** Four of the
24 sub-scores are library-correlation terms (`var_library_corr`, `_dotprod`, `_manhattan`,
`_sangle`); they can only be as good as the library's predicted intensities. The AlphaPeptDeep MS2
model is conditioned on (charge, **NCE**, **instrument**), and:

| | value | status |
|---|---|---|
| benchmark run | Orbitrap Astral (`MS:1003378`), beam-type CID, **collision energy 25.0** | **verified from the mzML** |
| OpenDIALibGen default | **NCE 30.0**, instrument `QE` | verified in `docs/OpenDIALibGen.md` |
| the 7.1M library actually used | unknown | **NOT verified** — the `.oswpq` metadata records the OpenMS converter, not the prediction parameters |

So the premise is half-established and must not be built on. **The decisive test needs no search
and no library rebuild**: take the confident IDs, and compare their observed fragment intensities
against PeptDeep predictions at NCE 25 vs NCE 30 (Astral vs QE). Whichever correlates better is the
answer, measured directly. Only if 25/Astral wins is re-predicting the library worth its cost.

Two things this also exposes, independent of the outcome:
- **The library format does not record how it was predicted.** A predicted library whose NCE and
  instrument are unrecoverable cannot be validated against the run it is used on. Whatever is
  decided here, `OpenDIALibGen` should write NCE/instrument/model-version into the bundle metadata.
- `PeptDeepMS2Inference::predictMS2` is callable at search time, so per-run re-prediction at the
  observed NCE needs no new dependency — but that is an optimisation to consider *after* the
  correlation test says the mismatch matters.

## Decided against: hand-rolling a replacement for `TransitionParquetFile`

The compact reader loads the same library in 15.0 s against the OpenMS reader's 144.8 s, so
replacing the production loader looks obvious. It was attempted and abandoned, deliberately.

Producing a `LightTargetedExperiment` requires reimplementing that reader's whole field mapping:
modification parsing (`(UniMod:4)` -> `compound.modifications`), the TraML-id fallback that decoy
pairing depends on, `/`-separated accession splitting, decoy flags, fragment types, and the
detecting/identifying/quantifying flags. A subtle mismatch anywhere produces a library that is
WRONG BUT STILL RUNS, and still emits q-values.

The arithmetic does not support that risk:

| | |
|---|---|
| gain | 145 s of a 1,903 s run (**4%**); ~6 GB transient of a 186 GB peak (**3%**) |
| risk | silent corruption of the most load-bearing structure in the tool |
| meanwhile | **extraction is 65% of wall and 90% of CPU** |

The compact representation's value is real but it is realised by the pipeline NOT materialising
`LightTargetedExperiment` for all 7.1M precursors -- only for the ~5.7% that survive prefiltering.
That is a pipeline change, not a loader swap, and it should be judged after extraction is fixed,
because extraction dominates everything else.

`-compact_probe` stays: it is the measurement that quantified the representation (1.94 GB vs
16.70 GB; 0.07 GB retained vs 11.29 GB) and it costs nothing to keep.

## Phase 2 — memory, gated on Phase 0.2

Do not start until the phase-resolved measurement names the dominant term. If it is chromatograms,
`ChromatogramStore` applies — but it **cannot** be integrated by editing OpenSwathWorkflow (the
OpenMS core is off-limits); ODIA must own the extract→score loop locally on public APIs. Both
reviewers warn the integration will otherwise either recreate the forbidden core edit or hold both
representations at once and *raise* peak RSS.

## Phase 3 — speed and occupancy, last

Serial phases measured so far: `library_load` 155 s at 1.1 cores (~215M string allocations and
~1.1e9 atomic refcount operations in the parquet→`LightTargetedExperiment` conversion),
`prefilter/build_decoy_view` 62.5 s, `pair_and_rebuild` 53.5 s, `precursor_index`, `score_load`,
`write_scores`, `context_fdr` — all at 1.0 cores. The structural fix is caching the *converted*
library, so the conversion is paid once rather than per run.

Not to be attempted before Phase 0.3 exists: any before/after on a shared node measures contention.

## Explicitly not on this list

- **Widening the prefilter screen.** Measured dead: 10 → 20 ppm cost 327 IDs. 30 ppm would be worse.
  The "30 ppm recovers 1,015 of 1,371 unreachable" line in earlier planning is superseded; it counts
  benefit without the null.
- **`ChromatogramStore` integration ahead of Phase 0.2.** Its headline saving is computed from the
  refuted decomposition.

---

## Results, 2026-08-01 (all against a measured noise floor of 6,487 +/- 83)

| # | change | IDs | verdict |
|---|---|---:|---|
| — | baseline (loader fix only, MS2 scope) | 6,506 | control |
| **1.2** | `-ms1_scores` (24 -> 36 features) | **6,574** | **+68: inside the noise. No effect.** |
| 1.1 | mass-cal RT direction fixed, window applied | 6,468 | -154 vs 6,622: outside noise, a LOSS |
| — | 20 ppm prefilter screen | 6,195 | -327: measured dead |

### 1.2 settled, and the original reasoning was wrong

The exclusion of the 12 MS1 sub-scores was justified in code by: *"sqlite (29 sub-scores) gave
6,607 identifications, this path (35-36 sub-scores) gave 4,913 -- 26% fewer for having MORE
features"*, concluding they were sparse noise the fit wasted capacity on.

That 4,913 was the `library_rt`-defective path. With that fixed the same path reaches 6,522 **with no
MS1 scores at all**, and a controlled A/B now shows the MS1 features are worth **+68 IDs against a
+/-83 noise floor** -- i.e. nothing.

(Correction: these were first reported as 6,437 vs 6,506. That compared the baseline's PASS-1
number against the treatment's final -- a two-pass run emits one `in-process LDA FDR` line per
pass. The corrected pair is 6,506 vs 6,574. The conclusion is unchanged, but it was reached by
comparing the wrong two numbers.) They never caused the deficit and they do not fix it.

Default stays `false`: 12 more features cost compute for no return. Right conclusion, wrong reason,
now corrected in the source comment.

### Performance results

| change | effect |
|---|---|
| loader: `resolveRaw` + exact `reserve` | library_load 144.8 s -> 133.9 s (**7.5%**); IDs unchanged |
| `-readOptions cacheWorkingInMemory` | **+10% wall, +13% peak RSS.** Extraction occupancy 38.9 -> 97.8 cores, but 2.4x the CPU for 3% less extraction wall |
| `MALLOC_ARENA_MAX=4` | **>2x slower** for -4.5% memory |
| parallelising the compact reader | **1.14x and no more** -- memory-latency bound |

### The pattern worth keeping

Everything that worked came from reading code or instrumenting: `library_rt` (+33% IDs), the
mass-calibration RT direction (both call sites), the fragmentation attribution (83-88% of peak), the
hardcoded `kNCE`/`kInstrument`, `OPENMS_TMPDIR`, and the accessor-level `shared_ptr`/`std::string`
churn.

Everything that tried to add parallelism or capacity to existing work regressed or gained nothing:
20 ppm screening, `MALLOC_ARENA_MAX`, `cacheWorkingInMemory`, three parallel-reader attempts, and
`reserve()` (predicted 24 GB, measured 0).

---

## 2026-08-01: determinism, the un-phased 27%, and what the allocator is not

### The ID "noise floor" was a bug, not noise

Three runs on byte-identical input: **6487 / 6565 / 6433**. This had been treated as a +/-1% noise
floor and used to size A/B tests. It was neither classifier -- GBT does no subsampling and LDA is
seeded. A precursor group's index is its FIRST-OCCURRENCE POSITION in row order, and row order is
whatever the parallel extraction produced, so a seeded shuffle of *indices* still assigned the same
precursor to a different fold on every run.

Fixed by sorting the group lists by precursor id before the shuffle (`odia_lda.h`). Fold sizes stay
exactly balanced. `odia_lda_test` now permutes the rows and requires identical per-row d-scores:
**8.290e-01** max delta without the fix, **1.954e-14** with it.

Consequence: every A/B measured before this was compared against a band that was partly
self-inflicted, and several "no effect" verdicts sit inside it.

### 27% of the run was outside every phase

| | wall | CPU-s |
|---|---:|---:|
| sum of all phases | 1449 s | 43857 |
| actual run | 1987 s | 60409 |
| **un-phased** | **537 s (27%)** | **16552 (27%)** |

The global RSS peak lands in that gap. MemProbe sampled the phase stack at 5 Hz and discarded it; it
now charges un-phased wall time to the phase it follows, so gaps have names.

One occupant found: the `mem/component` walk ran unguarded on every run -- serial, a `vector<string>`
per subordinate and a keyed lookup per meta value across 2.07M features and 113.9M meta values,
~150M allocations against a fragmented 180 GB heap. Timed live at **>=41 s** (joined mid-walk, so a
lower bound). Now behind `-mem_components`, default off.

### The instrumentation that measured nothing

Extraction-region instrumentation was "published and verified" three times and benchmarked for ~75
minutes while `strings(libOpenMS.so)` contained none of it, on the run node and on Ceph. The publish
target builds ODIA only, so an edit under `src/OpenMS` is compiled by nothing and the stale `.so`
ships as if fresh. `VERIFY_SYMBOL` passed throughout because the symbol given to it already existed.

A missing log line reads exactly like a region that cost nothing. `build_and_publish.sh` now refuses
when OpenMS sources are newer than the installed library, and names the file. The OpenMS-side
instrumentation was reverted: it cannot deploy under a deliberately read-only install, and the
question is answerable ODIA-side.

### Allocator contention is NOT what caps extraction

`LD_PRELOAD` tcmalloc, byte-identical run otherwise:

| phase | glibc | tcmalloc |
|---|---:|---:|
| `library_load` | 106.5 s | 86.5 s |
| `dia_run_load` | 130.6 s | 106.1 s |
| `prefilter` | 240.6 s | 173.5 s |
| `extract_pass1_wide` | 570.8 s / **41.9 cores** | 922.3 s / **40.8 cores** |
| `extract_pass2_narrow` | 266.5 s | 482.4 s |
| RSS after extraction | 186.01 GB | **80.15 GB** |

Extraction occupancy is **unchanged** (41.9 vs 40.8), so malloc is not the limiter. (A 30 s
steady-state sample showed 123/224 and was not representative; the phase average supersedes it.)
`MALLOC_ARENA_MAX=4` being >2x slower was a real effect read as evidence for the wrong cause.

tcmalloc's actual effect is memory. Complete run, same input, same binary:

| | glibc | tcmalloc | delta |
|---|---:|---:|---:|
| peak RSS | 189.4 GB | **105.8 GB** | **-44%** |
| wall | 32:07.65 | 38:58.74 | +21% |
| user CPU | 53195 s | 71298 s | +34% |
| sys CPU | 5466 s | 6347 s | +16% |
| IDs @ q<0.01 | 6433 | 6417 | none |

44% off the peak for 21% more wall time and identical IDs. Document as an `LD_PRELOAD` option for
memory-constrained nodes -- 189 GB needs a large node, 106 GB fits most of the cluster -- not a
default, since wall time is what the tool optimises.

Caveat not yet separated: glibc's peak lands in the un-phased gap, where the ~150M-allocation
`mem/component` walk also lived. Part of that 189 GB may have been the diagnostic itself. The `gap`
run (walk gated off, gap accounting on) settles it.

### The live candidate: batch granularity

Per-batch compound counts from the extraction log, 224 threads:

| batches | min | p50 | p90 | max | mean | max/mean |
|---:|---:|---:|---:|---:|---:|---:|
| 324 | 130 | 1416 | 10730 | 13418 | 3398 | **3.95x** |

324 units for 224 threads, and the largest is ~2.7x a perfectly-balanced share -- the makespan is
capped before anything else applies. `calculateInnerBatchSize` sizes batches as 5% of available RAM
/ 2 KB per compound, clamped to [2000, 10000], and **never consults the thread count**, so a bigger
machine gets coarser work units. On a 2.2 TB node it pins at the 10000 ceiling.

`-innerBatchSize` is already an ODIA option, so this is testable with no code change.
