# Answers from DIA-NN 1.7.12 to the calibration handoff questions

**Written 2026-08-02.** Source: `DiaNN/reconstructed/diann-1.7.12/src/diann.cpp` (11,260 lines, the
whole engine in one file). Line numbers below refer to that file. Cross-checked against
`DiaNN/docs/adaptive_mz_calibration_supplement.md`, an existing line-verified account of the m/z
path.

Answers the questions in `OpenDIAlyzer-calibration-handoff.md` §2–§6.

---

## Summary: DIA-NN has none of the three gates

The single largest finding is structural. DIA-NN does not decide whether a calibration estimate is
*trustworthy* and then discard it. It has no peakedness test, no "wider than configured, reject"
rule, and no yield-percentage threshold. Every question in §2 is asking how to make a better
accept/reject rule for a quantity DIA-NN never accepts or rejects.

What DIA-NN does instead, in order:

1. **Extract at a deliberately wide tolerance** so the anchor sample is not censored by the width
   being estimated. MS1 and MS2 calibration extraction is at **100 ppm** (`CalibrationMassAccuracy`,
   :217, applied at :6671), against a 20 ppm default operating width. RT calibration runs with
   `RT_windowed_search = false` (:6583) — **the full gradient, no RT window at all**.
2. **Make that affordable by subsampling precursors, not by narrowing tolerance.** `BatchMode`
   (:137) splits the library into `entries/2000` random batches with a fixed seed (:9000–9014).
   Calibration processes batches until 1000 IDs at q≤0.10 accumulate (`MinCal`, :243; loop
   :10371–10392). For your 423,079 prefiltered precursors that is 211 batches, of which typically
   the first few suffice.
3. **Take the width as a residual quantile, with no acceptance test** — `Q0.70(|r|)` (:9890).
4. **Then grow it against identification count** until IDs stop improving (:10450–10464).

Steps 1–2 together are the escape route from question B1, and they answer the cost objection that
made it look unaffordable. Step 4 answers B3.

---

## §2 Open question A: the hard gates

### A1 — "MS2 returned 18.7 ppm, wider than the configured 10 ppm, rejecting it"

**DIA-NN would accept 18.7 ppm and then try widening past it.** There is no comparison against a
configured width anywhere in the calibration path. The only comparison made is
corrected-vs-uncorrected on the same anchors (see A2). Three specific points:

- **18.7 ppm is a legitimate measurement here, not an outlier**, because DIA-NN's calibration
  extraction is at 100 ppm. Your 18.7 ppm figure is presumably measured at a 10 ppm extraction, in
  which case it is *itself* censored and the true value is larger still. The measurement is telling
  you the window is too narrow, and it is understating by how much.
- **After fitting, DIA-NN explicitly widens.** The recalibration path multiplies the width by 5,
  capped at the 100 ppm calibration ceiling (`MassAccuracy = Min(CalibrationMassAccuracy,
  MassAccuracy * 5.0)`, :10399), resets all precursors and weights, and re-searches. This is the
  operation your v1 design could not express ("contraction cannot express *too narrow*"). DIA-NN
  expresses it as an unconditional ×5 blowup whenever the first pass under-delivers.
- **MS1 and MS2 are never pooled.** Separate correction models (`MassCorrection` /
  `MassCorrectionMs1`), separate RT bin grids (`MassCalBins` / `MassCalBinsMs1`), separate widths,
  separate anchor counts, fitted in separate blocks (:9814–10041). MS2 is fitted first; MS1 failure
  does not block MS2. Your observation that MS1 (1.5 ppm) and MS2 (18.7 ppm) disagree by 10× is not
  a problem to reconcile — DIA-NN assumes they will differ and additionally inflates the MS1
  width by **5×** over its own residual quantile (:10004, :10035) because MS1 anchors are a
  quality-filtered subset and their residual spread understates the population.

Note the scale convention when comparing numbers: DIA-NN's width is a **half-width in relative
units**, the extraction interval being `μ(1−w) < m/z < μ(1+w)`, and `w = Q0.70(|r|)` — so the
starting width deliberately **clips ~30% of the calibration anchors**. It is not a coverage
statistic and was never meant to be one.

### A2 — is peakedness a good proxy?

**No, and DIA-NN's replacement is directly decision-relevant.** DIA-NN's only fit-quality test is:

```
MassAccuracy  = Q0.70(|residual after correction|)     // :9890
acc_no_cal    = Q0.70(|residual with no correction|)   // :9917
if (acc_no_cal < MassAccuracy) { MassAccuracy = acc_no_cal; MassCorrection = 0; }  // :9920–9924
```

That is, "does applying this transform make the residuals tighter than doing nothing?" If not, the
transform is zeroed and the *no-correction* width is used — the machinery is not discarded, its
output is replaced by the identity. The same test is repeated for MS1 (:10029–10036). Your setup
currently answers a question ("does this look peaked") whose answer does not determine an action;
DIA-NN's test's two outcomes are two different corrections, both usable.

Your peakedness statement is also framed against the wrong alternative: 1.47 vs. a threshold of 3.
DIA-NN would ask instead whether the fitted MS1 correction beats zero on the MS1 anchors, and if
the residuals really are flat noise it will beat nothing and get zeroed, at zero cost and with no
threshold.

### A3 — is 70 anchor pairs too few?

**For a mass-calibration fit, no — DIA-NN's minimum is 20.** `MinMassDeltaCal = 20` (:239), and
the outlier-refinement pass needs half that, 10 (`mass_cnt < min_mass_cnt / 2`, :9877). A reference
run drops it to 8 (`MinMassDeltaCalRef`, :240). 70 anchors comfortably passes.

**For an RT window, it is marginal — DIA-NN's minimum is 100** (`RTWinSearchMinCal`, :189), applied
at :9752. But the crucial difference is **what happens on failure**:

```
RT_windowed_search = (rt_delta.size() >= RTWinSearchMinCal);   // :9752
```

Below 100 anchors DIA-NN does not keep a default window — it **turns RT-windowed search off
entirely and searches the whole gradient** (the `RT_windowed_search` guards at :7179, :7201, :8146
all become no-ops). "Not enough evidence to narrow" resolves to *do not narrow*, which is the safe
direction. Your `rt_calib_min_yield_pct` gate resolves it to "keep 600 s", which is a narrowing
decision made on no evidence, and it costs 122 IDs.

There is a deeper point behind A3. DIA-NN never has a yield problem because it does not use a
special anchor set. Its anchors are whatever came out of the ordinary search at q≤0.10
(`MassCalQvalue = 0.1`, :160; `iRTMaxQvalue = 0.1`, :159) over random batches of the real library.
The 1.8% yield you measure is a property of running a 500-peptide CiRT set against a plasma sample,
not of the calibration problem. See §5 (question D).

**Concrete change:** the yield gate should be deleted, not tuned. Replace it with DIA-NN's failure
semantics — if the anchor count is too low to size a window, do not apply a window.

---

## §3 Open question B: replacing the gates

### B1 — does the "collect at a deliberately wider tolerance" escape route work?

**Yes. It is exactly what DIA-NN does, on both axes, and the cost objection is answered by
subsampling precursors rather than by avoiding the wide extraction.**

| axis | calibration-phase tolerance | operating tolerance | ratio |
|---|---|---|---|
| m/z (MS1 and MS2) | 100 ppm (:217, :6671) | ~10–20 ppm | 5–10× |
| RT | none — full gradient (:6583) | `RT_window` (§4) | ∞ |

The censoring point sits so far out in the tail that `Quantile_p(G) ≈ Quantile_p(F)` to within
nothing that matters, which is precisely the condition your v2 analysis identified as necessary.

The cost is bounded structurally rather than by tolerance:

- Library split into `entries / MinBatch` random batches, `MinBatch = 2000` (:138, :9000).
  Shuffle seed is fixed at 1 (:8999) — the batching is deterministic, so this does not reintroduce
  the order-dependency class of bug you just fixed.
- The calibration loop (:10371) processes batch 0, then 1, … stopping as soon as `Ids10 >= MinCal`
  (1000) **and** at least `MinCal/2` (500) of those carry a usable MS2 mass error (:10385–10390).
- So the wide-tolerance extraction touches a few thousand precursors, not 423,079. Order 1% of a
  full pass, at 5–10× the m/z tolerance and unlimited RT.

There is a second-level fallback for the case where even the wide pass under-delivers: if IDs stay
below `MinCal`, DIA-NN takes whatever fit it can get, blows the width up ×5, resets precursors and
weights, and re-runs the whole calibration loop (`goto calibrate`, :10393–10406), with
`MinCalRec = 100` (:242) as the lower bar for that second attempt. It never terminates in a state
where an under-wide window was inferred from an under-wide sample.

**For ODIA:** the pass-1 extraction already exists and is already wide (1435 s). The missing piece
is the m/z side — a subsampled, deliberately-wide-tolerance extraction to source mass anchors. The
"it needs an extraction at that wider tolerance, which is the thing being avoided" objection
dissolves once the extraction is over ~2,000 precursors instead of 423,079.

### B2 — does the small RT censoring bias generalise to m/z?

**DIA-NN does not measure it because it does not incur it**, and its design contains the reason your
v1/v2 kept tripping: it **truncates for fitting but never for width estimation**.

The outlier-refinement pass (`RemoveMassAccOutliers`, :10411–10417) *does* re-truncate anchors at
the current width and refit the bias model on the survivors (:9838–9841) — truncation used
deliberately, for robustness of the *shape* fit. But when it then computes the width, the residual
loop runs over **every eligible anchor with no outlier condition** (:9895–9901 for MS2,
:10022–10028 for MS1), including the ones excluded from the refit. The two uses are separated by
construction.

That separation is the missing rule in both your killed designs. v1 and v2 both estimated the width
from the same truncated sample they fitted on. Keep the trimming for the transform, and always
compute the quantile on the full anchor set collected at the wide tolerance.

Your empirical result (2.4× wider pass 1 moved the final p99 by 1.4 s) is consistent with the
censoring point already being far enough out that it does not bind — which is what DIA-NN engineers
on purpose. It is not evidence that RT is unusually benign; it is evidence that 1435 s was already
wide enough. Expect the same behaviour on m/z **at 100 ppm**, and not at 10 ppm.

### B3 — is hold-out capture the right acceptance criterion?

**No. DIA-NN uses identification count directly, and it is cheap enough to do so.**

After calibration and refinement, DIA-NN runs a one-dimensional search on the MS2 width
(:10446–10466):

```
w0 = MassAccuracy                      // Q0.70 of the refined residuals
loop:  w *= 1.2
       if (w > 10*w0) break
       reset_precursors(); reset_weights();
       re-extract, re-score, recompute q-values
       ids = Ids10                     // targets at q <= 0.10
       if (ids > best_ids) best = w, fail = 0;  else fail++
       if (fail >= 3) break
MassAccuracy = best_w
```

Three details worth copying:

- **The objective is IDs at q≤0.10, not q≤0.01** — ~5–10× more events, so the objective is much less
  noisy per unit of compute, and it is evaluated on a single batch (`curr_batch`), not the full
  library.
- **It only searches upward.** The residual quantile is treated as a lower bound on a sensible
  width, which is the correct prior given that the anchors came from a wide pass and the quantile
  is Q0.70.
- **Stop rule is 3 consecutive non-improvements**, capped at 10×w0. Cheap, and it tolerates a
  non-monotone curve — which is exactly the shape your §1.2 table has
  (6552 → 6695 → **6930** → 6835).

That last point matters for your open question about `recal1435`. DIA-NN would have stopped at 600 s
only after 900 s **and** two more candidates failed to beat it. With a ×1.2 ladder from 600 s the
candidates are 720, 864, 1037 — so the turnover you found at 900 s is one failure, not a stop.

Note also that DIA-NN's chosen width satisfies **no** coverage property. It starts from a quantile
that clips 30% of anchors and then moves in whatever direction IDs prefer. Hold-out capture would
have rejected the starting point outright. Your one direct test (widening pass 2 moved IDs and
capture the same way) had n=1; DIA-NN's design says do not build on it — measure the thing you
care about.

**Concrete change:** replace the whole adaptive-window plan with a ×1.2 upward ladder on
`-rt_calib_*` width and on the m/z tolerance, objective = IDs at q≤0.10 on one library subsample,
stop after 3 failures. That is ~40 lines and it retires the truncation-bias problem entirely,
because the criterion no longer references the anchor distribution.

---

## §4 Open question C: RT-segmented windows

### C1 — how to get per-precursor widths against a scalar API

**DIA-NN also has a scalar RT half-width, and declines to segment it.** It solves §1.4 a different
way. Three separate mechanisms:

**(a) The window *center* is RT-resolved; only the *width* is scalar.**

```
pred_RT = calc_spline(RT_coeff, RT_points, pep->iRT);           // :7179
RT_min = pred_RT - RT_window;  RT_max = pred_RT + RT_window;    // :7181
```

The transform is a **monotone cubic Hermite spline** (`spline()`, :1092–1151; `calc_spline()`,
:1153–1179) with up to 20 segments, `segments = Min(20, Max(1, 2*sqrt(n/20)))` (:9732, with
`RTSegments = 20` and `MinRTPredBin = 20`, :168–169). Knots are placed at **equal anchor counts**,
not equal RT spacing (:1153 in `spline`), so the knot density automatically follows the anchor
density — which is the direct structural answer to your §1.4 observation that the map steepens at
the end. The fit is run **twice**: fit, drop the worst 20% of residuals, refit (`map_RT`,
:1184–1203). Same trim-for-shape / full-sample-for-width discipline as B2.

**(b) The width is generous and does no discrimination.**

```
RT_window = Max( 2.0 * Q0.80(|residual|),                        // RTWindowMargin=2.0, RTWindowLoss=0.20
                 1.0 * (RT_max - RT_min) / 40 );                 // MinRTWinFactor=1.0, RTWindowFactor=40
                                                                 // :9751, constants :185–188
```

Two× a robust quantile, floored at 1/40 of the observed RT span. For your 2333 s gradient the floor
alone is ±58 s. The window is a prefilter, not a discriminator.

**(c) The discrimination is a learned score that sees RT position and RT deviation separately.**

```
sc[pRT]  = pos / span;                              // where in the gradient, normalised   :7756
sc[pdRT] = sqrt(Min(Abs(delta), span) / span);      // how far off, normalised, sqrt-compressed
```
(`score_RT`, :7742–7758 and :7760–7776.)

`pdRT` feeds the linear classifier from iteration `CalibrationIter+2` (:6621). `pRT` is withheld
from the linear classifier (:6644) but the **neural network is trained on the full score vector**
(`training[i][...] = it->target.info[0].scores`, :9353), so the NN sees both. A network with both
inputs can learn "a 90 s deviation is unremarkable at the end of the gradient and damning at the
start" — which is an RT-dependent tolerance, learned from the data, obtained without segmenting
anything and without a second extraction pass.

**For ODIA, in order of cost:**

1. **Free, no OpenMS change, no extra pass.** DIA-NN enforces its window *twice*: as a scan-range
   prefilter (:7201–7204) and again as a post-hoc rejection of the winning peak
   (`if (RT_windowed_search && |RT_apex − spline(iRT)| > RT_window) continue;`, :8146). The second
   enforcement is entirely downstream of extraction. Run `performExtraction` with a generous scalar
   window and apply an **RT-segmented or per-precursor acceptance rule on the candidate peaks**
   afterwards. You get per-precursor widths with zero changes to OpenMS and zero extra passes.
2. **Better.** Add `delta_rt` and normalised RT *position* as two scoring features and let the LDA
   or GBT learn the interaction, instead of hard-rejecting on `delta_rt` at all. Note the sqrt
   compression on the deviation — a raw `|delta|` in seconds is a poor feature for a linear model.
3. Only if 1 and 2 fail: K segmented passes.

Option 1 is also what makes the §1.4 table actionable without a redesign: the 14× spread across
octiles is a spread in *acceptance thresholds*, and thresholds are applied after extraction.

### C2 — tighter windows for pass-1-identified precursors

DIA-NN does not do this. Every iteration re-searches everything at the same window and lets the
classifier separate them. The nearest analogue is `--gen-ref` (:1971), which persists
high-confidence precursor indices to a reference library (`lib->save(gen_ref_file, &rt_ref, ...)`,
:9760) for use as calibrants in *future* runs.

The reason not to tighten is visible in §1.4 of your own handoff: removing the pass-2 clip made
early-RT calibration *better* (16.0 → 6.8 s) because the narrow window had been forcing a
compromise apex. A tight window on an already-identified precursor risks reproducing exactly that
on the precursors you are most confident about. Under C1 option 1, the same effect is obtained
safely: extract wide, accept tightly, and a bad acceptance decision costs one precursor rather than
distorting the apex.

### C3 — LC drift or a worse RT prediction for late-eluting peptides?

DIA-NN cannot distinguish these either, and is designed not to need to: the equal-count-knot spline
absorbs both as long as the bias is monotone in library RT. It does report a diagnostic —
`iRT prediction: median error` (:6508) — comparing its own in-silico predictions against measured
values.

**The test you can run with what you already have.** The reference map is
`library_rt → exp_rt`. Build a second map from a *different run on the same instrument and
gradient*, using **measured** RT from that run as the x-axis instead of library-predicted RT.

- If the last-octile residual degrades in both → LC drift (or genuine late-gradient peak-shape
  degradation). Fix is RT-resolved acceptance (C1).
- If it degrades only against library RT → the library's predictor is worse for hydrophobic
  peptides. Fix is on the library side, and it connects directly to question E.

`--gen-ref` shows how DIA-NN would source that second map cheaply.

---

## §5 Open question D: is CiRT worth its cost?

### D1 — 513 IDs for 424.9 s

**DIA-NN's default is bootstrap. The reference-peptide path exists and is off by default.**

`RefCal = false` (:181). `reference_run()` (:10228) is only invoked when `--ref <file>` is given
(:11042, :11125). Its thresholds are all *looser* than the bootstrap path's:

| | reference run | normal calibration |
|---|---:|---:|
| min anchors for RT window | 10 (`RTWinSearchMinRef`, :188) | 100 (`RTWinSearchMinCal`, :189) |
| min anchors for mass cal | 8 (`MinMassDeltaCalRef`, :240) | 20 (`MinMassDeltaCal`, :239) |
| RT window floor | span/10 (`RTRefWindowFactor`, :183) | span/40 (`RTWindowFactor`, :184) |

A reference run produces a *wider*, *less confident* calibration from *fewer* anchors. It is the
fallback for when bootstrapping cannot get off the ground, not a better estimator. Nothing in the
code treats it as an accuracy improvement.

So the DIA-NN-shaped answer to D1 is: **CiRT should not be on the default path.** Your own numbers
support this independently — 424.9 s and 15,277 CPU-s for +513 IDs, against +500 IDs from a window
change that costs 32 s.

But the more useful finding is *why* CiRT is expensive in ODIA and cheap in DIA-NN. DIA-NN's
bootstrap is not "search everything, then calibrate". It is "search a **random 2,000-precursor
subsample** of the real library over the full RT range at 100 ppm, repeat with the next subsample
until 1,000 IDs at q≤0.10". Your CiRT pass spends 69% of a full extraction's CPU on 500 curated
anchors plus 3,897 nonlinear candidates and returns 70 usable pairs — 1.8% yield — because a
500-peptide CiRT set has poor coverage in plasma. A random subsample of your own library would have
had a yield in the tens of percent for the same or less compute.

### D2 — bootstrap at the tuned 600 s pass-2 window

Untested, and DIA-NN offers no prediction. But note the mechanism: CiRT's advantage over bootstrap
is a better *transform*, and §1.1 already established that the transform's shape is not where IDs
are lost. If that holds, most of the 513-ID gap is a width effect and should shrink at 600 s. Worth
running — it is one flag.

### D3 — could bootstrap seed a cheaper CiRT?

DIA-NN does the reverse and it is the better trade. It does not make the anchor set smarter; it
makes the anchor set *ordinary* (a random subsample) and stops as soon as enough anchors accumulate.
A curated anchor set has to be searched in full before you know the yield; a random subsample lets
you stop at the first batch that clears the bar.

**Concrete change:** replace `setup/cirt_calibration` with a batched bootstrap — random subsample of
the prefiltered library, full RT range, wide m/z tolerance, accumulate until N anchors, stop. Keep
`-rt_calibration cirt` as an explicit opt-in for runs where bootstrap fails, with the DIA-NN-style
looser thresholds. This subsumes D1–D3 and also supplies the wide-tolerance anchor sample that
question B1 needs.

---

## §6 Open question E: library provenance

### E1/E2 — what was the library predicted at?

**DIA-NN 1.7.12 cannot help with this and does not have the problem.** The open source has no
spectrum predictor: the `predictor::Predictor` class is declared (:756) but the call site is behind
`#ifdef PREDICTOR` (:5195–5203) and the implementation is not in the distribution. Its only
in-silico stage is **RT**, via `predict_irt` (:1726–1748) — a linear additive model over amino-acid,
modification and terminal-position features. There is no NCE parameter, no instrument parameter, and
no fragment-intensity prediction anywhere in the file. (The DL predictor with instrument/NCE
settings arrives in 1.8; not verifiable from this source.)

That said, the code says something useful about the *shape* of the answer. DIA-NN never treats
library fragment intensities as fixed truth:

- Intensity agreement enters as **scores** (`pTimeCorr`, `pCos`, `pCorr`, the `pSig` family) that a
  **run-specific** classifier re-weights every iteration. A systematically wrong intensity pattern
  is partly absorbed by re-weighting.
- Each library entry carries a `lib_qvalue` (:3082, read at :4830) which is propagated into the
  protein-level error: `err = Max(dratio, lib_qvalue) + lib_qvalue` (:10634 region). Library
  quality is a per-entry input to FDR, not an assumption.
- `--gen-spec-lib` / `--out-lib` writes an **empirical** library from the run's own identifications,
  with measured rather than predicted RT by default (`iRTOutputLibrary = true`, :266;
  `--out-measured-rt` disables, :1901).

So the DIA-NN answer to "was this library predicted for the right instrument" is **stop asking and
replace it**: run once with the predicted library, export an empirical library from the confident
IDs, re-search. That is a direct answer to E1 that does not require recovering the provenance.

For E2, nothing in 1.7.12 constrains the NCE/instrument choice. Your reasoning (HCD in the ion
routing multipole, so `QE` or `Lumos`, the analyser difference not being in the embedding) is sound
and unaddressed by this code. Since it is now a flag sweep, run the grid — but score it the way
DIA-NN scores its width sweep: **IDs at q≤0.10 on one subsample**, not on a proxy, and expect the
empirical-library route to beat any grid point.

---

## What this changes, ranked

1. **Delete all three gates** (§2). Replace the yield gate's action with DIA-NN's: too few anchors
   ⇒ do not apply a window, rather than keep the default. Free, +122 IDs measured.
2. **Adopt the ×1.2 upward ID-yield ladder** for both the RT window and the m/z tolerance (§3, B3),
   objective = IDs at q≤0.10 on one library subsample, stop after 3 consecutive failures, cap at
   10× the starting value. This retires the entire adaptive-window design problem — no coverage
   statistic, no truncation bias, no hold-out split. ~40 lines.
3. **Replace CiRT with a batched bootstrap** over random subsamples of the real library, full RT
   range, wide m/z tolerance (§5). Saves ~425 s/run and simultaneously produces the uncensored
   anchor sample B1 needs.
4. **Separate trim-for-fit from full-sample-for-width** everywhere (§3, B2). Trim the worst residuals
   before fitting the transform; never trim before computing the width quantile.
5. **RT-segmented acceptance, not RT-segmented extraction** (§4, C1 option 1). Extract with a
   generous scalar window, apply the RT-resolved threshold on candidate peaks downstream.
   Clean-room-safe, no extra pass.
6. **Add normalised RT position and sqrt-compressed |delta_rt| as scoring features** (§4, C1 option
   2) so the classifier can learn the RT-dependent tolerance instead of it being a hard threshold.
7. **Export an empirical library and re-search** (§6) rather than trying to recover the benchmark
   library's NCE.

## Key line references

| what | file:line |
|---|---|
| 100 ppm calibration extraction tolerance | `diann.cpp:217`, applied `:6671` |
| RT-windowed search off during calibration | `:6583`, guards at `:7179`, `:7201`, `:8146` |
| Batching: `entries/2000`, fixed seed | `:137–138`, `:8999–9014` |
| Calibration batch loop, stop at 1000 IDs | `:10371–10392` |
| ×5 width blowup + full restart on under-delivery | `:10393–10406` |
| Width = `Q0.70(|r|)`; corrected-vs-uncorrected test | `:9890`, `:9917–9924` |
| Width quantile computed on **untrimmed** anchor set | `:9895–9901`, `:10022–10028` |
| Outlier-refinement trims for the **fit** only | `:9838–9841`, `:10411–10417` |
| MS1 width inflated ×5 over its quantile | `:10004`, `:10035` |
| ID-yield width ladder (×1.2, cap 10×, 3 fails) | `:10446–10466` |
| `RT_window = Max(2·Q0.80(|r|), span/40)` | `:9751`, constants `:184–188` |
| RT-window anchor minimum ⇒ disable windowing | `:9752`, `:189` |
| Monotone Hermite spline, equal-count knots | `:1092–1179` |
| `map_RT`: fit, trim worst 20%, refit | `:1184–1203` |
| `pRT` / `pdRT` scoring features | `:7742–7776` |
| NN trained on full score vector | `:9353` |
| Reference-peptide run, opt-in, looser thresholds | `:10228`, `:11042`, `:183`, `:188`, `:240` |
| No spectrum predictor in open source | `:756`, `:5195` (`#ifdef PREDICTOR`) |
| Empirical library export with measured RT | `:266`, `:1892`, `:1901` |
