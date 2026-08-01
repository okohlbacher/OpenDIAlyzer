# Adversarial review: adaptive m/z calibration supplement — and a fast plan for ODIA

**Reviewed document:** `~/Downloads/adaptive_mz_calibration_supplement.md` ("Supplementary Methods:
Adaptive Calibration of Mass-to-Charge Bias and Extraction Tolerances").

**Provenance of this report.** Two workflow runs were needed. The first hit the account's monthly
spend limit (71/107 agents); the resume completed 102/106, so the 3-vote adversarial verification
DID run this time — but the synthesis agent died again, on a session limit. This report is
synthesised by hand.

**Read the verification verdicts carefully.** The verifiers were harsh: most claims drew ≥2 refutes
and were killed. Inspecting the evidence, the majority are *split verdicts* — the factual core is
confirmed verbatim, but the load-bearing inference drawn from it overreaches, and the verifier
refutes on that. So "refuted" here usually means "the fact is right, the conclusion was too strong",
not "the fact is wrong".

**The figures this report depends on were verified by me directly against the vendored BSD OpenMS
3.6.0 source**, not taken from an agent claim:

| figure | source | status |
|---|---|---|
| `mz_estimation_percentile` = **99.0**, `mz_estimation_padding_factor` = **1.3** | `SwathMapMassCorrection.cpp:102-108` | **verified** |
| `windows:rt_percentile` = **95.0**, `rt_estimation_padding_factor` = **1.3** | `CalibrationWorkflow.cpp:97-103` | **verified** |
| `goodness:median` = **4.0** ppm, `goodness:MAD` = **2.0** | `InternalCalibration.cpp:195-196` | **verified** |
| `RT_chunking` default = **300 s** one-sided; −1 means a global model | `InternalCalibration.cpp:180` | **verified** |
| RANSAC outlier removal present, **disabled by default**, threshold 10.0 ppm² | `InternalCalibration.cpp:182-188` | **verified** |

---

## 1. The finding that matters most: Q0.70 is ~1σ, not a coverage window

The document sets the extraction half-width to the empirical **70th percentile** of absolute
relative residuals, and acknowledges it is "an adaptive extraction parameter rather than a
conventional estimate of instrumental standard deviation."

The robust-statistics literature makes that acknowledgement sharper than the document does. In
robust Chauvenet rejection, **the 68.3rd percentile of absolute deviations is defined as a robust
estimate of one standard deviation** — the Gaussian-equivalent scale. So Q0.70 ≈ 1σ, and an
extraction interval of ±Q0.70 retains ~68% of true signal and **clips ~32% by construction**. A
coverage-based window needs a *multiple* of that quantity (≈2.5–3× for ~99% retention), not the
quantity itself.

The counter-example is in our own dependency. **OpenSwathWorkflow sizes windows from HIGH
percentiles with padding > 1**: m/z at the **99th percentile × 1.3**
(`Calibration:MassIMCorrection:mz_estimation_percentile = 99.0`,
`mz_estimation_padding_factor = 1.3`) and RT at the **95th × 1.3**
(`Calibration:windows:rt_percentile = 95.0`). The BSD reference implementation deliberately
*over-covers* the residual distribution; the reviewed document deliberately under-covers it and
then recovers the loss with a downstream identification-yield search.

**Verdict: the Q0.70 choice is not wrong so much as mislabelled.** It is a scale estimate being
used as a coverage window, and the 1.2^k yield search exists to undo that. Removing the search
without also raising the multiplier would clip a third of the signal.

## 2. Where the document is out of step with published practice

| point | document | published practice | verdict |
|---|---|---|---|
| anchor confidence (C) | q ≤ 0.10, **no confidence weighting** | mzRefinery: **q < 0.01**; TRIC: strict 1% to *anchor*, loose 5% to *report* | **document is an outlier** — 10× the contamination, weighted equally with true anchors |
| half-width (D) | Q0.70 (≈1σ) | OpenSwathWorkflow Q0.99 × 1.3 (m/z); OpenMS acceptance: median ppm < 4.0, **MAD < 2.0** | **document ~3× too narrow** before its yield search |
| bootstrap | ±100 ppm → fit → 5× expand | mzRefinery: two-pass wide-then-narrow, first search ~±50 ppm with restricted search space | same shape, document starts 2× wider |

Two further notes. mzRefinery's model needs only **≥2 calibrants (linear) or ≥3 (quadratic)** to be
identifiable, so the document's minimum of 20 is conservative rather than risky. And robust
outlier rejection has a published **operating ceiling of 40–70% contamination** — beyond that the
authors state the technique stops working. That is decisive for ODIA: with ~99.8% of library
targets genuinely absent, anchors **must** be drawn from an identification-filtered set, never from
raw m/z matches, or every robust estimator is past its breakdown point before it starts.

## 3. Where I was wrong, and where the document is defensible

**B — "most intense peak within ±100 ppm biases the residuals."** I asserted this. The Astral
literature complicates it: **64% of centroids contain ≤5 ions while contributing only 6.5% of total
ion current.** The interferent population is numerous but individually *weak*, so selecting the
most intense peak is **safer than selecting the nearest** whenever true signal sits above the
single-ion noise floor — "nearest" would preferentially grab single-ion noise sitting closer to the
query centre. My criticism was wrong in sign for this instrument class. It would still bite for a
low-abundance precursor whose true peak is at the noise floor.

**F — "B_max = 1 makes the RT machinery inert." I now think the document is out of step here, and I
withdraw my earlier softening.** I had defended `B_max = 1` on the strength of a claim that drift is
predominantly *between-run* (σ 0.93 ppm single-run → 1.54 ppm over 60 runs / 136 h). That claim was
**refuted for scope overreach** in verification, and two reference implementations default the other
way:

- **OpenMS `InternalCalibration`: `RT_chunking` defaults to 300 s** (one-sided), with −1 documented
  as the opt-in for "ALL calibrants for all scans, i.e. a global model" — so RT-resolved is the
  default and global is the exception (verified, `InternalCalibration.cpp:180`).
- **mzRefinery** ships an LC-dependent model binning errors in 75-second windows, merged until
  ≥100 identifications contribute.

The reviewed document defaults to the opposite. That does not make it wrong for a given instrument,
but it is the minority position and the document offers no measurement defending it.

*Falsifiable test, cheap, and now worth running:* fit the ppm offset separately in early/mid/late RT
thirds of one ODIA run. If the between-third spread is small relative to the within-third σ
(~0.9 ppm), RT segmentation is genuinely inert on our data and `B_max = 1` is fine here.

**New — the model form contradicts the canonical lock-mass assumption.** Cox, Michalski & Mann
(2011), *Software Lock Mass by Two-Dimensional Minimization of Peptide Mass Errors* (JASMS), model
the error as an **additive** sum of two functions, one per variable, stating that "the
non-linearities in the mass scale should be independent of elution time", and fit each as piecewise
linear. The reviewed document instead uses `a·m² + b(t) + c(t)·m`, whose `c(t)·m` term is an
**m/z × time interaction** — precisely the dependence the canonical method argues is absent. Either
the document has evidence for an interaction that the lock-mass literature says should not exist, or
the term is unnecessary and is spending degrees of freedom on noise.

**A — the m⁻⁴ weighting.** Not resolved by the surviving agents. My objection stands on
first principles (constant-ppm error ⇒ variance ∝ m² ⇒ inverse-variance weight m⁻²; the document
uses m⁻⁴, over-weighting low m/z by a further m²), but I have no source confirming or refuting it
and the verification agents that would have tested it were killed by the spend limit. **Treat as
open.**

**E — yield-driven width search.** Also unverified. The concern (selecting the extraction width by
maximising IDs at q ≤ 0.10, then reporting at a stricter threshold, is selection on the test
statistic) is standard, but no source was retrieved that settles whether it measurably inflates
FDR. **Treat as open.** It is moot for ODIA regardless — see §4.

## 4. Plan for ODIA

The document's own procedure is **not affordable here**. Its 1.2^k search re-extracts, re-scores
and re-estimates confidence per candidate; at ~10–30 min per ODIA extraction pass that is hours
against a 1–2 minute budget. The parts worth taking are the *model* and the *bootstrap*, not the
search.

### 4.1 Estimator: half-sample mode, not KDE

For the m/z offset and width, the literature points away from KDE mode-finding:

- **Silverman's critical-bandwidth bootstrap is not correctly calibrated** for general modality
  testing — the calibration depends on unknown tail quantities (Hall & York corrections exist).
- The **half-sample mode (HSM)** (Bickel & Frühwirth 2006) recursively takes the shortest interval
  containing half the remaining points and returns **both** a mode location (→ the offset δ) **and**
  a scale from the retained half-range (→ σ̂). It is claimed more outlier-resistant than other
  robust location estimators *and* than low-bias mode estimators including KDE-based ones.

HSM is O(n log n) after a sort, needs no bandwidth choice, and yields exactly the two quantities the
hook must return. **Use HSM; skip KDE.** MassTraceFinder-style density estimation is the right tool
for finding traces in m/z–RT space, but it is the wrong instrument for estimating a 1-D calibration
mode from an already-paired residual sample.

### 4.2 Anchors: reuse what OpenSwathWorkflow already does

OpenSwathWorkflow already implements stratified subsampling with published budgets: **RT-binned
sampling (100 bins × 5 peptides linear ≈ 500 anchors; 2000 × 50 nonlinear), intensity-ranked top
fraction (0.4 linear / 0.7 nonlinear), then an evidence prefilter against raw data before fitting.**
That is the speed answer and the bootstrap answer in one, and ODIA already calls this machinery.

Anchor budget: mzRefinery needs ≥2–3 calibrants for identifiability and ≥100 per time segment for
its LC model. **Several hundred stratified anchors is ample** for a scalar offset + width — orders
of magnitude below the 67,429 candidates ODIA's calibration was extracting.

### 4.3 Concrete recipe for `inferMassAccuracyPpm_`

1. Take the existing auto-iRT stratified anchor sample (RT-binned, intensity-ranked). **No new
   extraction pass.**
2. Filter to identifications at **q < 0.01**, not 0.10 — matching mzRefinery, and mandatory in the
   99.8%-absent regime where robust estimators are otherwise past their breakdown point.
3. Pair observed vs theoretical m/z; compute signed ppm error. **MS1 and MS2 separately** (measured
   0.60 vs 1.63 ppm — a pooled number is 3× wrong for MS1).
4. **HSM** → offset δ and scale σ̂.
5. **Guard before trusting it:** peakedness. A true error distribution is centrally peaked; a
   noise-dominated one is flat. Require `density_center ≥ 3 × density_edge`. Fail → return ≤ 0
   ("not inferred") and keep the configured window. This is what stops the estimator widening
   itself to death.
6. Window = **k · σ̂**, with **k ≈ 3** (≈99% coverage for a Gaussian) — *not* Q0.70. Cross-check
   against OpenSwathWorkflow's Q0.99 × 1.3. Floor at ~2 ppm so a high-S/N anchor set cannot
   collapse the window below instrument capability.
7. Let **`SwathMapMassCorrection` own the offset** (it already fits `MZTrafoModel` and rewrites
   spectra); the hook returns only the residual *width*. Do not double-correct.

### 4.4 Validation — the part that decides whether it worked

"Narrower" is the failure mode, not the success criterion. Estimate on anchors only, then measure
**hold-out capture** on the 2,957 identifications from the 10 ppm run — a set never used in
estimation. Require **≥99.5%** of their fragments inside the inferred window. Accept the *smallest*
window that still passes. Secondary: OpenMS's own acceptance criteria (median ppm error < 4.0,
MAD < 2.0) and restored target/decoy separation.

### 4.5 Ranking by (expected gain)/(implementation risk)

1. **HSM + peakedness guard + k·σ̂ window in `inferMassAccuracyPpm_`** — highest gain (this is the
   0 → 2,957 ID axis), low risk, reuses existing anchors and `SwathMapMassCorrection`. **Build first.**
   Sanity-check the output against OpenMS's own acceptance thresholds, which are verified and free:
   median ppm error < 4.0 and MAD < 2.0.
2. **RT-thirds diagnostic** (§3, F) — one cheap run; decides whether RT-dependent mass correction is
   worth building at all. Do it before writing any segmentation code. Now higher priority than I
   first ranked it, because both reference implementations default to RT-resolved.
3. **Anchor threshold q < 0.01 and drop the yield search** — a parameter change plus a deletion.
4. **Consider enabling RANSAC** — OpenMS `InternalCalibration` already implements it for calibration
   outlier removal but ships it **disabled** (threshold 10.0 ppm²). Free to try, and directly aimed
   at the contaminated-anchor problem; measure against HSM rather than assuming it helps.
5. **RT window from a held-out residual** in `inferRtWindowSeconds_` — deferred until the 810 s
   figure underpinning the current RT reasoning is re-measured; it traces to a script comment and
   has never been verified.

---

## 5. What would falsify this plan

- HSM and median/MAD agree to within ~10% on real anchors → the estimator choice does not matter;
  take the simpler one.
- The peakedness ratio exceeds 3 even on a deliberately noise-only anchor set → the guard does not
  discriminate and needs replacing.
- Hold-out capture at k = 3 falls below 99.5% → σ̂ is underestimated; raise k rather than changing
  the estimator.
- RT-thirds offsets differ by ≫ 0.9 ppm → `B_max = 1` is wrong for this instrument and RT-dependent
  correction is needed after all.
