# Handoff: open questions in ODIA recalibration and tuning

**Written 2026-08-02.** Self-contained: everything needed to discuss these questions in a fresh
context is stated here, including the measurements that motivate them and the ones that refute
earlier reasoning.

Benchmark throughout: one Thermo Astral plasma DIA run (`astral.mzML`, 6.4 GB, 2333 s gradient,
150 SWATH windows) against a 7,149,966-precursor PeptDeep-predicted library (78,569,077
transitions), prefiltered to 423,079 precursors / 4,463,919 transitions before extraction.

**Scoring is now deterministic.** Five runs returned exactly 6430 IDs across different schedulers,
thread topologies, node loads and a 48% timing perturbation. Any ID difference below is real, not
noise. This was not true before 2026-08-01; two order-dependency bugs (LDA fold assignment by row
position, GBT histogram chunking by row index) produced a +/-1% spread that earlier work treated as
a noise floor.

---

## 1. What is measured and settled

### 1.1 The RT transform's shape is correct; its WIDTH is what costs IDs

`experiments/calmap.py` builds a reference map `library_rt -> exp_rt` from a run's own q<0.01 rank-1
target features, then reports the residual around it. Result: the residual median sits within
+/-10 s of zero across the whole gradient, and the tool's own `delta_rt` tracks the reference
residual almost exactly (median 31.9 vs 30.3 s, p99 113.0 vs 114.8 s). **Global RT calibration is
not where IDs are lost.**

### 1.2 Window widths cost 500 IDs, and the defaults are wrong for this data

| pass 1 | pass 2 | IDs | note |
|---:|---:|---:|---|
| 600 s | 240 s | 6430 | current defaults |
| 1435 s | 240 s | 6552 | +122 -- pass-1 window explicit, bypassing the yield gate |
| 1435 s | 400 s | 6695 | +265 |
| 1435 s | 600 s | **6930** | **+500 (+7.8%)** |

Cost of the widest: peak RSS 89.9 -> 92.9 GB, extraction 961 -> 993 s. Negligible. The
counter-hypothesis (extra candidates cost more in FDR than they return) is refuted over this range.

`recal900` and `recal1435` were queued to find the turnover and had not finished at handoff.

### 1.3 Pass 2 was clipping the residual distribution

Widening pass 2 decompresses the residuals, which is the direct signature of clipping:

| pass-2 window | p99 residual | p99 delta_rt |
|---:|---:|---:|
| 240 s (+/-120) | 113.4 s | 110.9 s |
| 400 s (+/-200) | 155.5 s | 168.3 s |
| 600 s (+/-300) | 163.2 s | 198.2 s |

At 240 s the p99 of 113-115 s sits pressed against the +/-120 s boundary. At 600 s it is 163 s,
comfortably inside +/-300 for the first time.

### 1.4 Calibration quality varies 14x along the gradient

Median `|delta_rt|` by RT octile, once clipping is removed (600 s pass 2):

| octile | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| median &#124;delta_rt&#124; (s) | **6.8** | 20.0 | 29.9 | 33.1 | 35.1 | 40.6 | 46.5 | **95.9** |

Removing the clip made early RT *better* (16.0 -> 6.8 s: the narrow window had been forcing a
compromise apex) and revealed late-RT deviations that were previously cut away (56.6 -> 95.9 s).

The mechanism is visible in the reference map: it steepens sharply at the end. Library RT
0.83 -> 0.87 (a 4% span) covers 1810 -> 2041 s, so the mapping compresses and a fixed library-RT
error becomes a much larger run-RT error.

### 1.5 A real bug fixed: RT rescale composed across passes

`rescaleLibraryRT_` read the CURRENT `c.rt`, so pass 1 gave `rt := t1(rt_orig)` and pass 2 gave
`rt := t2(t1(rt_orig))`. But `t2` comes from `recalibrate_`, which fits `library_rt -> exp_rt` using
`library_rt` from `precursor_index_` -- built once, before any mutation, so in ORIGINAL units. Pass 2
evaluated the transform far outside its fitted domain and every window landed nowhere.

Effect: `-rt_calibration bootstrap` returned **0 IDs**; after the fix, **6039**. Pass 1 had been
healthy throughout (bootstrap's pass-1 p95 anchor residual was 112.86 s, *better* than cirt's
117.24 s). `cirt` was never affected because `use_native_trafo` skips the function entirely -- which
is why the defect survived: the default path does not touch it. Fixed in `367f25c`.

`-rt_calibration none` still returns 0, correctly: it uses normalized library RT verbatim as
seconds, so windows sit at 0-1 s on a 2333 s gradient.

---

## 2. Open question A: the hard gates

Three gates fire on every run and discard usable estimates. `setup/mass_calibration` costs 1.2 s and
its output is thrown away every time -- the machinery is inert.

```
MS1: residuals are FLAT (peakedness 1.47059 < 3). That is what a mostly-noise anchor set looks like
MS2 mass calibration returned 18.7327 ppm, WIDER than the configured 10 ppm -- rejecting it
MS1: only 141 anchor errors (need 200) -- not inferred
yield 1.796254% is below -rt_calib_min_yield_pct. Keeping 600 s. Calibration may only narrow.
```

The last one alone costs 122 IDs.

**Questions.**

- A1. The MS2 estimate of 18.7 ppm against a configured 10 ppm is rejected as "wider". But is the
  data saying the extraction window is genuinely too narrow? Nothing has tested a wider m/z window.
  Note MS1 sigma was measured at 1.5 ppm, matching the instrument spec of 1.66-1.71, so MS1 and MS2
  disagree by ~10x and pooling them is wrong.
- A2. Is peakedness a good proxy at all? It asks "does this look peaked", where the quantity of
  interest is "would this window have retained held-out anchors".
- A3. `rt_calib_min_yield_pct` rejects a 1.8% anchor yield. Is 70 anchor pairs from 3897 candidates
  actually too few for a *window* estimate? mzRefinery needs >=2-3 calibrants for identifiability
  and >=100 per time segment for its LC model.

## 3. Open question B: replacing the gates without hard thresholds

Two designs were written and both were killed by adversarial review before implementation. The
failures are recorded because they are easy to reinvent. Full text:
`docs/OpenDIAlyzer-adaptive-window-plan.md`.

**v1 -- monotone contraction with hold-out validation.** Killed by three findings:
- *Truncation feedback.* Anchors collected AT the current window are truncated at it, so the scale
  is underestimated, proposing a smaller window, truncating harder. The hold-out check cannot detect
  this because it is computed on the same truncated sample: it measures capture *conditional on
  having survived the cut*.
- *Contraction cannot express "too narrow"*, which is exactly what the MS2 18.7 ppm figure says.
- *The Wilson lower bound relocates the minimum-anchor threshold rather than removing it*:
  `n/(n+z^2) >= 0.99` needs `n >= 380`, worse than the explicit 200 it replaced, and un-auditable.
- Also: hashing *precursor* id splits a peptide's charge states across estimation and validation
  with correlated residuals, breaking the iid assumption Wilson needs.

**v2 -- distribution-free quantile ladder on a fixed sample.** Not yet re-reviewed in full, but the
partial review found the truncation bias survives: with `G` the truncated CDF and `F` the true one,
`Quantile_p(G) = Quantile_{p*F(W0)}(F)`, and validation measures `G(w) = F(w)/F(W0) >= F(w)`. At 1%
clipped mass, true coverage is 0.9801 while measured capture reads 0.99 and passes.

**The general obstacle: you cannot estimate a width from a sample that was selected by that width.**

**Questions.**

- B1. Does the escape route work -- collect anchor residuals at a deliberately WIDER tolerance than
  any window under consideration, so the censoring point sits far out in the tail? Affordable in
  principle (a few hundred stratified anchors, not 423,079 precursors), but it needs an extraction
  at that wider tolerance, which is the thing being avoided.
- B2. Empirically, the censoring bias measured *small* for RT: a 2.4x wider pass-1 window moved the
  final p99 by 1.4 s. Does that generalise to m/z, or is RT unusually benign here?
- B3. Is hold-out capture even the right acceptance criterion? It is a proxy for identification
  yield. The one direct test available (widening pass 2) moved IDs and capture in the same
  direction, but n=1.

## 4. Open question C: RT-segmented windows

§1.4 says a scalar window is sized for the worst region and is ~40x too wide where the calibration
is good. Both reference implementations are RT-resolved: OpenMS `InternalCalibration` defaults
`RT_chunking` to 300 s one-sided (-1 means global, documented as the exception), and mzRefinery bins
in 75 s windows merged until >=100 identifications contribute.

**Questions.**

- C1. `OpenSwathWorkflow::performExtraction` takes a SCALAR `rt_extraction_window`. How does ODIA get
  per-precursor widths without modifying OpenMS (CLEAN-ROOM: OpenMS stays unvendored)?
  - Partitioning the library by predicted RT and running pass 2 once per segment is ODIA-side
    orchestration. Whether it costs K extra passes over the spectra, or approximately nothing
    because an RT-segmented library also needs only an RT-segmented slice of spectra, has not been
    worked out.
- C2. `applyEmpiricalRT_` already replaces pass-1-identified precursors' RT with the measured apex.
  Those precursors could take a much tighter window than unidentified ones. Same scalar obstacle.
- C3. Does the late-gradient degradation (95.9 s median in the last octile) reflect genuine LC drift,
  or the library's RT prediction being worse for late-eluting (more hydrophobic) peptides? These
  imply different fixes.

## 5. Open question D: is CiRT worth its cost now that bootstrap works?

`setup/cirt_calibration` is a real extraction pass: **424.9 s wall, 15,277 CPU-s** for 500 linear
anchors / 5293 transitions plus 3897 nonlinear candidates -- 69% of a full pass-1 extraction's CPU
for 0.12% of the compounds, returning 70 anchor pairs (1.8% yield). Its CPU cost is invariant to
thread configuration (15,277.5 vs 15,270.9 CPU-s across two very different splits).

Before the rescale fix, bootstrap scored 0 and CiRT looked indispensable. After the fix bootstrap
scores **6039** against CiRT's 6552 at matched pass-1 width.

**Questions.**

- D1. Is 513 IDs (7.8%) worth 424.9 s? That is roughly the same magnitude as the window-tuning gain,
  obtained for free.
- D2. Bootstrap has not been tested with the tuned pass-2 window (600 s). CiRT's advantage may shrink
  or vanish.
- D3. Could bootstrap seed a cheaper CiRT (fewer anchors, narrower nonlinear window) that recovers
  the difference for a fraction of 424.9 s?

## 6. Open question E: library provenance is unrecoverable

`library/metadata.json` records generator, counts and fragment-type breakdowns but **no NCE, no
instrument, no model version**. For every library shipped so far, "was this predicted for the
instrument we are searching?" cannot be answered from the artifact.

`OpenDIALibGen` had `kNCE = 35.0f` and `kInstrument = 2` (timsTOF) as compile-time constants,
measured from a Bruker run's `analysis.tdf`; they are now `-nce` / `-instrument` options with those
values as defaults, and a provenance sidecar is written. But the *benchmark library's* settings
remain unknown.

**Questions.**

- E1. What was this library predicted at? If timsTOF/NCE 35 for an Astral run, the fragment
  intensities are systematically wrong and every downstream score is degraded.
- E2. peptdeep has no Astral embedding. `QE` and `Lumos` are the defensible choices (HCD in the ion
  routing multipole -- fragmentation is what the model conditions on; the analyser difference is not
  in the embedding). Which, and at what NCE? A grid over `-nce`/`-instrument` is now a flag sweep.

## 7. Methodological caveats on everything above

- **Node contention.** `spock` acquired other users mid-session (load 160-240 against 224 cores).
  Wall-clock comparisons spanning that window are unreliable and several were withdrawn. ID counts
  are unaffected. `experiments/bench_odia.sh` now records `uptime` before and after every run;
  `scripts/ibmi-nodes.sh` classifies node state ("BUSY - IDs ok, timings not").
- **The reference map is built from a run's own IDs**, so it is censored at that run's extraction
  window. §1.3 exploits this deliberately (comparing maps at different windows measures the
  censoring); it also means a single map's p99 must not be read as the underlying RT accuracy.
- **n=1 per configuration.** Determinism makes each number exact, but exactness is not
  generalisation: every figure here is one run, one instrument, one library.

## 8. Files

| path | what |
|---|---|
| `docs/OpenDIAlyzer-mz-rt-calibration-review.md` | the research handoff this work draws on: HSM, anchor budgets, OpenSwathWorkflow's Q0.99x1.3 defaults, the Q0.70 misuse |
| `docs/OpenDIAlyzer-adaptive-window-plan.md` | v1 and v2 designs, with §7 recording v1's refuted claims |
| `docs/OpenDIAlyzer-path-forward.md` | the wider optimisation log this sits inside |
| `experiments/calmap.py` | reference-map tool; `calmap.py run1.oswpq run2.oswpq` prints residual stats and a censoring comparison |
| `experiments/bench_odia.sh` | one benchmark invocation, records node load |
| `src/opendialyzer.cpp` | `rescaleLibraryRT_`, `recalibrate_`, `calibrateCiRT_`, `calibrateMassFromPass_`, `buildPrecursorIndex_` |
