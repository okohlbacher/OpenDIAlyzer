# RT/predictor calibration: OpenDIAlyzer vs DIA-NN — how we differ, and how we close the gap

Status: 2026-07-27. Companion to the DIA-NN calibration research (agent report, same date)
and `docs/OpenDIAlyzer-merger-plan.md`. Sources tagged **[DOC]** documented, **[CODE]** read
from source, **[INF]** inferred; full citations in the research report.

## 1. How OpenDIAlyzer calibrates today (`src/opendialyzer.cpp`)

Two-pass, in-process, no external tools:

1. **Linear bootstrap** (`main_`): library RT is normalized ~[0,1]; a 3-point
   interpolated-linear map sends it to the run's actual RT range `[rt_min, rt_max]`
   (read from any SWATH map). This just gives pass 1 a *bounded, centred* window
   (`rt_extraction_window`, default 600 s). Whole-range extraction (`rt_win=-1`) exploded
   the candidate count ~100× and choked the writer, so the window is always bounded.
2. **Extract pass 1** over that window → raw `.osw`.
3. **Recalibrate** (`recalibrate_`): in-process semi-supervised LDA scores every peak group;
   the top target peak groups become RT anchors; `fitTrafo_` fits a monotone
   `library_RT → observed_RT` map.
4. **Extract pass 2** over a narrow window (`rt_extraction_window_recal`, default 240 s).
5. **Final FDR** (`finalScore_`): LDA + q-values, count targets at q<0.01.

### What `fitTrafo_` does now (after the 2026-07-27 change)
- 60 quantile bins over the anchors; per-bin **median** x and y (outlier-robust).
- **PAVA isotonic regression** on the bin medians → provably monotone fit.
  *(This replaced a forward-max clamp that pinned y to a running max and created flat
  plateaus wherever a noisy bin dipped.)*
- Piecewise-linear interpolation between the isotonic knots.
- Reports the 95th-pct **in-sample** anchor residual (diagnostic only — see §3).

### Anchor policy
Top-2000 targets by d-score (lenient), NOT the honest q<0.01 set — because the wide pass-1
often has ~0 targets at honest 1% FDR. Fine for a robust RT map; explicitly *not* safe as
training truth (matters for any future fine-tuning). RT-derived sub-scores
(`VAR_NORM_RT_SCORE`) are excluded when selecting RT anchors (anti-circularity).

## 2. How DIA-NN / alphaDIA / Spectronaut differ

All three implement the same three-part idea; we were missing pieces 1 and (partly) 3.

| | OpenDIAlyzer (before) | DIA-NN | alphaDIA |
|---|---|---|---|
| RT transform form | binned median + **forward-clamp** | **isotonic + spline, monotone** [DOC] | LOESS (6 local regressors) [CODE] |
| Iterations | **1** recalibration, 2 passes total | calibrate→ID→recal on the growing 1%-FDR set→narrow, **3–5×** [DOC] | iterative, per-dimension [CODE] |
| Window sizing | **fixed** 600→240 s | **from the observed residual distribution**, auto-narrowed [DOC] | `k × Pxx(residual)`: RT P99.9×1.5, MS2 P99×1.1 [CODE] |
| Unseen-peptide residual | predicted RT only (PeptDeep sd≈**810 s** here) | **empirical/observed-value 2nd pass** collapses seen-bulk residual to the measurement floor [DOC] | same, **plus retrains** the predictor (RT err 317→11 s) [CODE] |
| Predictor adaptation | none | per-run RT/mass/IM **calibration** + a per-run NN FDR classifier; weights **not** retrained by default [DOC] | genuine transfer-learning fine-tune of RT/MS2/Charge [CODE] |
| Mass window | fixed 30 ppm | **dynamic, per-precursor**, auto-tightened [DOC] | 99th-pct residual, floor 5/10 ppm [CODE] |

### The one insight that explains our gap
Our recalibration residual (~43 s on the anchors) is the **in-sample** fit residual. The
residual that actually sets the window is the **predictive** residual on the ~½ of the
library *not yet identified* — governed by PeptDeep's prediction error, measured here at
**sd≈810 s, P95≈1720 s** (`run_full_arm.sh:27-32`). DIA-NN/alphaDIA reach ~11–21 s windows
**not by fitting anchors better** but by (a) sizing the window from the real residual
distribution and (b) **eliminating the unseen-bulk residual** — via an observed-value
second pass (seen peptides) and, in alphaDIA, retraining the predictor (unseen peptides).

## 3. What we changed on 2026-07-27, and what we deliberately did NOT

**Implemented (deterministic, unit-tested locally; cluster build pending maintenance):**
1. **Isotonic PAVA transform** — the documented functional-form difference (Calib-RT: DIA-NN
   = isotonic+spline). Robust median binning kept; forward-clamp replaced with true isotonic
   regression. Helps the *unseen* peptides (better global map). *(pava self-test: 5 cases
   incl. weighted pooling + textbook; engine `--selftest` now also asserts monotonicity.)*
2. **Empirical library second pass** (`-empirical_rt`, **default OFF pending FDR
   validation**) — the **#1 lever** from the research. For every pass-1 target scoring above
   the best decoy (model-free high-confidence gate), its **measured apex RT** replaces the
   predicted RT in pass 2. Seen precursors' residual → measurement floor, so the same-width
   window is now *centred on the measured pass-1 apex* instead of on a mis-prediction. `recalibrate_`
   emits the map keyed by `TRAML_ID` (== `LightCompound::id`); `applyEmpiricalRT_` applies it
   after the global map.
   - **FDR safety — the symmetric-procedure design (2026-07-27, post-benchmark):** the first
     attempt (target-only + paired-decoy at the *target's* RT) was asymmetric — codex's
     winner's-curse critique: a false target locked onto its own fluke keeps the advantage; its
     decoy, given the *target's* RT, was never selected for its own extremeness. **Fix:** apply
     empirical RT to **every precursor — target AND decoy — using its OWN pass-1 apex** (the
     best-scoring pass-1 peak group's observed RT). The procedure is now *identical* for both
     classes, so a false target locked onto its own fluke is mirrored by a decoy locked onto
     *its* own fluke → the decoy null captures the inflation → FDR preserved. No pairing /
     `"DECOY_"` assumption. And since every seen precursor is centred on the peak it actually
     found, the residual → ~0 and the narrow window stops collapsing. Still default OFF pending
     the entrapment/decoy-diagnostic confirmation on the EMP8 benchmark, but FDR-safe by design.
3. **p95 anchor-residual diagnostic** — logged each pass to watch calibration quality
   improve across iterations.

**Deliberately NOT done (and why):**
- **Auto window-sizing from the anchor residual** — the anchor residual is *in-sample*
  (~43 s), not the *predictive* residual (~810 s). Sizing the window from it would collapse
  the window and lose every unseen peptide — the exact trap the research warns about. Safe
  auto-sizing needs a held-out predictive-residual estimate (backlog).
- **Predictor fine-tuning (neural weight updates)** — infeasible in this build:
  `ONNXPredictorBase` is inference-only (no optimizer/gradient/checkpoint symbols, verified).
  "Fine tuning" here means per-run *calibration*, not weight training. Neural fine-tune is a
  separate research spike (merger-plan Phase 3a, NO-GO today).

## 4. Ranked backlog to close the rest of the gap (gain-per-effort, all deterministic)

1. **Iterate calibrate→search→narrow 3–5×** (not just once). Low effort; the loop DIA-NN
   credits Spectronaut for. Needs the extraction to be fast enough (task #22).
2. **Empirical MS2 + IM** replacement (extend the observed-value pass beyond RT). Needs the
   `FEATURE_TRANSITION` intensity reader + observed 1/K0 (merger-plan Phase 2b).
3. **Held-out predictive-residual estimate** → then residual-driven window sizing becomes
   safe. Medium effort.
4. **Per-run mass recalibration** (fit ppm-vs-m/z from confident IDs; dynamic ppm window).
   OpenMS `InternalCalibration` covers most of it.
5. **Per-run IM (1/K0) calibration** — affine-align predicted CCS/IM to observed.
6. **Predictor transfer-learning** (alphaDIA-style) — highest effort, needs training-enabled
   ONNX + the merger's live-predictor ownership. Biggest gains for non-standard mods /
   instruments; for standard tryptic on a covered instrument, steps 1–2 capture most of it.

**Do 1+2 first.** They target the dominant lever (collapse the unseen/seen residual and
iterate) and are deterministic. 6 (neural) is last and gated on a feasibility spike.
