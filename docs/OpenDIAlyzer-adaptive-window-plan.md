# Plan v2: distribution-free window selection by validated quantile ladder

**Status: PLAN. Not implemented.** v1 was reviewed adversarially (kimi) and three findings killed
it. This is the rewrite. §7 records what v1 got wrong, because the failures are instructive and two
of them are easy to reintroduce.

## 1. The problem, from this run's own logs

Every hard gate in the mass-calibration path fires and rejects, on every run:

```
MS1: residuals are FLAT (peakedness 1.47059 < 3). That is what a mostly-noise anchor set looks like
MS2 mass calibration returned 18.7327 ppm, WIDER than the configured 10 ppm -- rejecting it
MS1: only 141 anchor errors (need 200) -- not inferred
mass accuracy still not inferable from 2000 pass-1 anchors; keeping the configured window
```

`setup/mass_calibration` costs 1.2 s and its result is discarded every time. The machinery is inert.

## 2. Design

Anchors are collected **once**, at the window already in force. The search then happens entirely on
that fixed sample -- nothing is re-extracted, so the sample's truncation point never moves.

```
R        <- residuals of anchors at q<0.01, collected ONCE at W0    # MS1 and MS2 separately
est, val <- split R by hash of PEPTIDE SEQUENCE                     # not precursor id: see 3.4
delta    <- HSM location on est

# Is the sample even informative about width, or is it clipped by W0?
if quantile(|r - delta| for r in est, 0.99) >= W0 * 0.95:
    return "not inferable"                                          # keep W0. See 3.2.

# Iteratively tighten: walk the ladder from loose to tight, keep the tightest that validates.
best <- W0
for p in [0.9995, 0.999, 0.995, 0.99, 0.98, 0.95]:                  # decreasing coverage
    w   <- quantile(|r - delta| for r in est, p)                    # DISTRIBUTION-FREE. No k*sigma.
    if w >= best: continue
    cap <- WilsonLower(#{r in val : |r - delta| <= w}, |val|, 0.95)
    if cap >= target: best <- w                                     # tighter AND validated
    else: break                                                     # stop at first failure
return max(best, floor)
```

The iteration the request asks for is the ladder: progressively tighter candidate bounds, each one
earning its place on data that did not propose it, stopping at the first that cannot.

## 3. Why each piece, and what it replaces

### 3.1 One fixed anchor sample, not re-collection per iteration

v1 re-collected anchors at the contracted window each pass. That is a truncation feedback loop: a
sample truncated at W has an underestimated scale, which proposes a smaller W, which truncates
harder. Worse, the validation is computed on the same truncated sample, so it measures capture
*conditional on having survived the cut* -- it can never see the signal it clipped.

Collecting once removes the loop entirely and is strictly cheaper.

### 3.2 A clipping test, so "the window is too NARROW" is representable

The logs say the MS2 residual spread wants 18.7 ppm against a configured 10 ppm. A contraction-only
design cannot express that; it would tighten into an already-miscalibrated run.

The guard is structural, not a tuned constant: if the 99th percentile of `|r - delta|` reaches the
collection boundary, the sample is censored and carries no usable width information -- return "not
inferable" and keep `W0`. That is today's behaviour, reached for a stated reason rather than by a
peakedness proxy.

**This is the honest answer, not a workaround.** Widening the extraction window is a decision about
the *extraction*, and it costs a re-extraction to validate. It is out of scope here, and §6 records
it as the follow-up it is.

### 3.3 Empirical quantiles, not `k * sigma`

v1 used `k = 3` for "~99% coverage", which is a Gaussian statement. Mass-error distributions are
routinely heavy-tailed, and bimodal when two populations mix (charge states, a lock-mass step
mid-run). Under a heavy tail, `3*sigma` under-covers; under bimodality, `sigma` describes neither
mode.

Taking the quantile directly makes the coverage claim mean what it says, whatever the shape, and
removes both `k` and the Gaussian assumption. HSM is kept for *location* only, where its
outlier-resistance is the point.

### 3.4 Split by peptide sequence, not precursor id

The same peptide at 2+ and 3+ hashes to different precursor ids but shares its chemistry and its
local calibration state; the residuals are correlated. Splitting on those puts correlated
observations on both sides, and the Wilson interval assumes iid Bernoulli -- correlation makes it
overconfident. Hashing the stripped sequence keeps every charge state of a peptide on one side.

### 3.5 The implied minimum anchor count, stated rather than hidden

The Wilson lower bound does **not** remove a minimum-n threshold; it relocates one. At `x = n`
(perfect observed capture), `LB = n / (n + z^2)`, so `LB >= 0.99` at 95% confidence requires
`n >= z^2 * 0.99 / 0.01 ~= 380`.

Replacing an explicit, auditable "200 anchors" with a larger implicit one buried in two other
constants is worse, not better. So: compute it, log it, and let the ladder's failure be attributed.

```
mass window: 141 val anchors, target 0.99 @95% needs >=380 for any step; keeping 10.0 ppm
```

### 3.6 Floor

Floor at the instrument's capability (~2 ppm) so a high-S/N anchor set cannot propose a window below
what the hardware can deliver. If the ladder bottoms out at the floor on every dataset, the floor is
doing the work and the estimator is uninformative -- that is falsification #3.

## 4. Unchanged from the research handoff

- Anchors from the existing stratified auto-iRT sample. **No new extraction pass.**
- `SwathMapMassCorrection` owns the offset; this returns only the width (handoff §4.3 step 7).
- MS1 and MS2 estimated separately (0.60 vs 1.63 ppm measured; pooling is 3x wrong for MS1).
- Anchors at q < 0.01 (handoff §2: at ~99.8% library absence, looser anchors put robust estimators
  past their 40-70% breakdown ceiling before they start).

## 5. Falsification

1. The ladder never takes a step -> our anchor sets carry no width information. Report it; do not
   ship a no-op that looks like a feature.
2. It takes steps and IDs drop -> hold-out capture is not a sufficient proxy for identification
   yield, and acceptance must become ID-based (expensive; handoff rejects on cost).
3. It bottoms out at the floor on every dataset -> the floor is doing the work.
4. The clipping test fires on every dataset -> W0 is genuinely too narrow and the real finding is
   about the extraction window, not the estimator.
5. Validation capture and identification yield disagree in sign on any dataset -> the whole
   acceptance criterion is wrong.

## 6. Out of scope, deliberately

- **Widening** the extraction window when the clipping test fires. Needs a re-extraction to
  validate. This is the natural follow-up and the MS2 18.7 ppm figure says it may matter.
- RT-segmented mass correction. Handoff §3/§5 proposes a cheap prior test: fit the ppm offset in
  early/mid/late RT thirds, compare between-third spread against within-third sigma (~0.9 ppm).
- RANSAC (handoff §4.5 item 4): measure against HSM rather than assume.

## 7. What v1 got wrong

| v1 claim | why it was wrong |
|---|---|
| "monotone contraction is safe by construction" | safe against *widening*, but it cannot express "your window is too narrow", which is what the MS2 18.7 ppm figure says. It would contract into a miscalibrated run. |
| "re-collect anchors at the current window each iteration" | truncation feedback: truncated sample -> underestimated scale -> tighter window -> harder truncation. Validation computed on the same truncated sample cannot detect it, because everything it observes survived the cut. |
| "the Wilson bound removes the need for a minimum-anchor constant" | it relocates it. ~380 at target 0.99 / 95% conf, versus the explicit 200 it replaced -- and un-auditable. |
| "split by hash of precursor id" | a peptide's charge states land on opposite sides with correlated residuals; Wilson then assumes an independence that does not hold. |
| "window = k * sigma, k ~ 3" | a Gaussian coverage claim on a distribution that is routinely heavy-tailed and sometimes bimodal. |
