# Why ODIA's LDA disagreed with pyprophet — line-by-line comparison

Both implementations read side by side: `src/odia_lda.h` against pyprophet 3.0.15
(`pyprophet/scoring/semi_supervised.py`, `classifiers.py`, `data_handling.py`, `_config.py`),
installed at `/ceph/ibmi/abi/oliver/envs/pyprophet`.

Reference measurement being explained: on the identical feature table
(`bench/narrowmz/mz10.osw`, 2.07M features, 10 ppm), pyprophet reported **4,302** target
precursors at q<0.01 and ODIA's in-process LDA **2,957**. Same features, same labels — so the
features carry the signal and the difference is entirely in the scoring algorithm.

## The differences, in order of measured impact

### 1. The semi-supervised loop never ignites — cold start (ROOT CAUSE)

**ODIA (before):** the loop is seeded with a *single* feature — the one with the largest absolute
Welch t between target and decoy training rows, given weight ±1. Iteration 0 then ranks precursors
by that one column, computes q-values, and selects training positives at q ≤ 0.15.

**pyprophet:** never seeds from a single feature. OpenSWATH hands it a composite `main_score`, and
`learn_randomized` does `train.rank_by("main_score")` before anything else. Its first *model* is a
full LDA (`LDALearner.learn`, `sklearn.LinearDiscriminantAnalysis`) on all top-decoy peaks against
the targets that main_score put above the 15% cutoff.

**Why the single feature is not enough.** Measured on `testdata/lda_fixture.txt` (44,568 real
OpenSWATH rows, 2,500 precursors, 29 sub-scores):

```
[fold 0 SEED] feature=4 |t|=13.69 finite=29482/29482 (100.0%) distinct_z=23920
[fold 0 it 0] m=29 need=31 | train groups=1666 (T=590) | minq=0.9654 medq=0.9654
              | q<=0.15: 0  q<=0.50: 0 | got=0
```

The seed feature is not degenerate — fully populated, 23,920 distinct values, |t| = 13.7. But
|t| = 13.7 over 29,482 rows is a **0.16 sd** per-row effect, and the fixture has ~18 candidate peak
groups per precursor. The statistic that actually ranks precursors is the *maximum* over those ~18
rows, and the max of 18 draws is governed by extreme-value spread, which swamps a 0.16 sd shift.
Result: `minq == medq == 0.9654` — a **flat q-value for every target group** — so zero positives are
selected, the fit is skipped, `w` is unchanged, and every subsequent iteration skips identically.

`trained=0, skipped=9` on all three folds. The "semi-supervised LDA" was ranking by one sub-score.

**Fix:** seed with a real multivariate direction. Target/decoy labels are known outright — no FDR
estimate is needed — so an LDA of all top-*target* rows against all top-*decoy* rows is available for
free and is far stronger than any one column. The target class is contaminated (most target
precursors are false), which shrinks the fitted direction but does not rotate it: the contaminant
*is* the decoy distribution, so it biases μ⁺ toward μ⁻ and costs magnitude, not orientation. It only
has to be good enough to ignite the loop.

Effect on the fixture:

| | target mean | decoy mean | distinct d-scores |
|---|---|---|---|
| single-feature seed | 1.011 | 1.011 | 35,964 |
| LDA seed | **0.119** | -0.000 | **44,347** |

Separation goes from nil to 0.119 sd and the tie blocks disappear.

### 2. Asymmetric class construction — top peaks vs all peaks

**ODIA (before):** positives are one row per precursor (`candidate.best_row`, the top-scoring peak
group) but negatives were **every** decoy row.

**pyprophet:** both classes are top-peaks-only — `td_peaks = train.get_top_decoy_peaks()`,
`bt_peaks = train.get_top_target_peaks().filter_(score >= cutoff)`.

With ~4–18 candidate peak groups per precursor, the old negative class was dominated by runner-ups,
so the two classes differed in *two* ways at once: target vs decoy (wanted) and rank-1 vs runner-up
(not wanted). The fitted direction partly separates "best peak in its group" from "not the best
peak" — a real, learnable axis carrying no target/decoy information — and it pulls both the negative
mean and the within-class covariance.

**Fix:** `LDAParams::top_decoys_only` (default true), matching pyprophet.

### 3. Cross-fold scores were pooled without being commensurable

**ODIA (before):** each of the F folds fits its own `w` and scores its own held-out groups; the raw
`dot(w, z)` values are then pooled into a single ranking for the final q-values. An LDA direction is
defined only up to scale, and its offset depends on that fold's training set, so fold A's 4.0 and
fold B's 4.0 mean different things. The pooled ranking was therefore partly sorted by *which fold a
precursor landed in*.

**pyprophet:** sidesteps this by averaging the fold weight vectors into one model
(`learner.averaged_learner(ws)`) and rescoring everything with it — at the cost of the leakage-free
property ODIA keeps (its averaged model has seen every group).

**Fix:** `LDAParams::normalize_folds` (default true) rescales each fold's held-out scores to that
fold's own decoy null (mean 0, sd 1). That is the one scale target-decoy FDR cares about — "how many
decoy sds above the null" — and it makes pooling valid while keeping folds leakage-free.

### 4. Parameter differences (not yet changed)

| | pyprophet default | ODIA |
|---|---|---|
| CV repeats | `ss_num_iter` = 10 | `n_folds` = 3 |
| semi-supervised iterations | `xeval_num_iter` = 10 | `n_iter` = 3 |
| train fraction | `xeval_fraction` = 0.5 | (F−1)/F = 0.67 |
| initial training FDR | `ss_initial_fdr` = 0.15 | 0.15 ✓ |
| iteration training FDR | `ss_iteration_fdr` = 0.05 | 0.05 ✓ |
| π₀ in the *training* cutoff | bootstrap, λ ∈ (0.1, 0.5, 0.05) | off |
| first fit uses main score | **no** (`use_main_score=False`) | n/a — ODIA has no main score |
| solver | sklearn SVD, rank-truncated at tol=1e-4 | Cholesky + ridge 1e-6 |

Two of these are worth noting even though they are not yet changed:

- **π₀ off makes ODIA's training cutoff stricter than pyprophet's at the same nominal 0.15**, so it
  selects fewer training positives — which directly feeds the cold-start failure in §1.
- **Solver regularisation differs in kind, not degree.** sklearn's SVD solver *discards* directions
  whose singular value falls below tol; ODIA's Cholesky with ridge 1e-6 *inverts* them, producing
  large weights along near-collinear directions. OpenSWATH sub-scores are strongly collinear
  (`xcorr_shape` vs `xcorr_shape_weighted`, `library_corr` vs `library_dotprod`, …), so this is a
  real difference in noise amplification, not a numerical detail.

## What the test fixtures can and cannot show

The synthetic oracles (`odia_lda_test`, `odia_lda_adversarial_test`) are **insensitive to all three
fixes** — 1094→1082 IDs on one, 1052→1056 on the other, FDR control identical at 0.008/0.009. Their
decoys are drawn i.i.d. with no rank structure, so the confound in §2 does not exist in them and the
seed in §1 always ignites.

`testdata/lda_fixture.txt` shows the cold-start failure clearly (§1) but **cannot measure the gap**:
a faithful reimplementation of pyprophet's loop — same sklearn `LinearDiscriminantAnalysis`, same
defaults (10×10, 0.15/0.05, xeval_fraction 0.5) — run on that same fixture **also never trains**
(0/10 folds). The fixture has no bootstrappable signal, so 0 IDs is the correct answer for it and it
cannot discriminate the two algorithms.

The gap measurement therefore has to be run on the real table (`-score_osw bench/narrowmz/mz10.osw`),
against the recorded references of 2,957 (ODIA) and 4,302 (pyprophet).

## Result on the real table

`-score_osw bench/narrowmz/mz10.osw -threads 32`, same 2.07M features as both references:

```
OpenDIAlyzer: LDA fitted 9 iteration(s), skipped 0.
OpenDIAlyzer: in-process LDA FDR -> 4367 target precursors at q<0.01.
2:18.58 wall, 1.4 GB peak RSS
```

| | IDs @ q<0.01 | decoys @ q<0.01 | empirical FDR |
|---|---|---|---|
| ODIA, single-feature seed | 2,957 | — | — |
| pyprophet 3.0.15 | 4,302 | — | — |
| **ODIA, LDA seed** | **4,367** | 42 | **0.96%** |

`skipped 0` — the loop now trains on every fold and iteration, where before it was starving. IDs
+47.7% over the old ODIA, and 1.5% ahead of pyprophet on identical input.

The FDR is calibrated: 42 decoys against 4,367 targets is 0.96% empirical at a 1% nominal cutoff.
Worth noting that ODIA reaches this while being the *more conservative* of the two on the q-value
side — π₀ correction is off here (`use_pi0=false`) and on in pyprophet — so the extra identifications
come from a better discriminant, not from a looser error model.

The §4 parameter differences (3×3 vs 10×10 folds/iterations, π₀ in the training cutoff, SVD
rank-truncation vs ridge) remain unaddressed and are the obvious next lever.
