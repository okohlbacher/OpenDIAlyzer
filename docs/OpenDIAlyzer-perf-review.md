# Quality & performance — adversarial review and priorities (2026-07-29)

Reference point, measured on ONE Thermo Astral plasma DIA run (~4 GB mzML) against the SAME
complete 78.5M-transition human library, on an idle 224-core / 2.2 TB node:

| tool | wall | IDs | CPU | peak RSS |
|---|---|---|---|---|
| DIA-NN 2.0 | **509 s** | 8,164 precursors @1% q | 7519% (~75 cores) | **22.7 GB** |
| OpenSwathWorkflow | >71 min, running | — | 936% (~9 cores) | 133 GB |
| OpenDIAlyzer | >85 min, running | — | 7044% (~70 cores) | **861 GB** |

So on the only complete three-way datapoint we have: **>10× slower and ~38× more memory than
DIA-NN.** Everything below is ranked against closing that gap.

---

## P0 — Memory: 861 GB for one run

**Symptom.** RSS grows with work completed, superlinearly at first, apparently plateauing near
861 GB. Earlier on a 1/2000 subset it went 1.7 → 9.8 GB over 17 min.

**What has been ruled out.** Recycling the per-thread mzPeak reader every 256 decodes changed
nothing (tested, reverted). So it is NOT reader decode state.

**Leading hypotheses, untested:**
1. `OpenSwathWorkflow` accumulates every extracted chromatogram and feature in memory until the
   run ends. 78.5M transitions × float arrays is easily hundreds of GB. `batchSize` exists
   precisely to bound this and we pass **0 (no batching)**.
2. `load_into_memory=true` (our `auto` default for mzML) materialises the working SWATH map on
   top of that.
3. The 7.1M-compound library itself is resident twice — once in `LightTargetedExperiment`, once
   in the OSW writer's seeded tables.

**Why this is P0.** A tool needing 861 GB per run cannot be used on the "16 cores, 128 GB, no
GPU" workstation the README promises. It fails the project's own design principle #2.

**Cheapest decisive experiments (run in parallel, they are independent):**
- `-batchSize 5000` and `-batchSize 50000` — if RSS drops proportionally, hypothesis 1 is right
  and the fix is a sane default, not new code.
- `-readOptions normal` (load_into_memory=false) — isolates hypothesis 2. NOTE this also
  disables the wave scheduler (`OpenSwathWorkflow.cpp:486`), so read the CPU number too.
- `heaptrack`/`massif` on a small subset — the only thing that will actually name the allocator.

---

## P1 — The serial 18.5 GB TSV library parse (7:39, single-threaded)

Measured directly: `Load TSV file -- done [took 07:50 m (CPU), 07:39 m (Wall)]`, then
`Loaded 40371 proteins, 7149966 compounds with 78569077 transitions.` One core throughout, and
it is where much of the resident library memory is born.

**DIA-NN completes its ENTIRE run in 509 s** — less than OpenMS spends parsing the library. This
is an architectural difference (indexed binary vs text re-parse), not an extraction problem.

**Fix, cheap and benefits BOTH OpenMS tools:** convert the TSV to **PQP** (SQLite) once with
`TargetedFileConverter`, then pass `.pqp`. `TransitionPQPFile` already exists; OpenDIAlyzer
already advertises `-tr` as accepting PQP. Expected: 7:39 → seconds, plus a smaller resident
footprint because PQP can be read incrementally.

**Risk to check:** whether PQP loading is genuinely lazy or just a faster full load. If the
latter it fixes the time, not the memory.

---

## P2 — Only ~1% of the library is addressable, and the two tools disagree on how much

`TransitionListEvidenceFilter retained 34165 of 3603425 target precursors` (OpenSWATH), while
OpenDIAlyzer reported `67429 peptides`. Same file, same library, ~2× discrepancy.

Two consequences:
1. **Any ID comparison is suspect** until this is understood. If ODIA extracts precursors OSW
   correctly excludes, extra "IDs" are noise, not sensitivity.
2. **We are doing ~100× more work than necessary.** Filtering the library to the acquisition's
   m/z windows BEFORE extraction would shrink every cost above. DIA-NN effectively does this.

This is a correctness question first and a performance opportunity second — which is why it
outranks raw speed work.

---

## P3 — Parallel efficiency: ~70 of 224 cores (31%)

Better than it looked earlier today (the "1.5 cores" and "24-core ceiling" findings were
contention and library-size artifacts, now retracted), but still 69% idle.

Known contributors, in order:
- the serial library parse (P1) — pure serial time in Amdahl terms
- `min_per_chunk = 2000` disables within-window chunking on small per-window coordinate counts;
  correct in principle, but it means parallelism silently tracks library size
- OSW is at ~9 cores in the same run, so some of this is inherent to the workflow, not us

**Do not optimise this before P0/P1.** Amdahl: if 7:39 of a 90-minute run is serial, perfect
extraction parallelism still leaves the serial parse.

---

## P4 — Correctness debt that could invalidate results

- `Read chromatogram while reading SWATH files, did not expect that!` on BOTH mzML inputs.
  Something in our load path does not expect chromatograms alongside spectra.
- **PASS00779 needs re-verification.** The 495,651 vs 491,377 result predates the IM fix,
  the readOptions rework, the progress logging, and an OpenMS rebuild.
- **Feature counts are not IDs.** Until pyprophet runs, OSW/ODIA numbers are raw features and
  are NOT comparable with DIA-NN's q-value-filtered 8,164. Any table mixing them is misleading.

---

## P5 — Originally requested, still not started

`Calibration:*` registration, CiRT/LOESS recalibration, and online library finetuning. Deferred
deliberately: layering finetuning onto a tool that needs 861 GB and 90 minutes optimises the
wrong thing. Revisit once P0/P1 land.

---

## Method notes for the overnight iteration

Today produced six wrong claims, every one *plausible-but-unverified*. Guards now in force:

1. `features == 0` is a FAILURE, never a timing datapoint. Three silent aborts (rsq threshold,
   missing CiRT anchors in sampled subsets, overlapping windows) each looked like fast results.
2. **Never benchmark concurrently.** Concurrency made OpenSWATH read 84 cores vs 156 idle and
   produced two retracted conclusions.
3. **Never extrapolate from one datapoint.** 13.6 days → 3.3 days once a second point existed.
4. **Conservation-style invariants can be tautological.** A wrong split conserves too.
5. **Match parameters across arms.** Stock OSW defaults differ from ours on mz window, min_rsq
   and IM; unmatched arms measure defaults, not tools.
6. Do not wait on `pgrep -f <pattern>` your own command contains — it self-matches.
