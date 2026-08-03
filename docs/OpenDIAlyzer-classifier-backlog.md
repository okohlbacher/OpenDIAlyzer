# Classifier / anti-circularity backlog

Open items from the night of 2026-08-02, after two adversarial review passes (kimi, codex) over
commits `3956706` → `f9307b3` → `8955640`. Everything listed here was **confirmed against the
source** and deliberately *not* fixed yet, with the reason given.

Fixed already, for reference: M5 threshold mismatch, fake no-CV arm, dead iteration cap, silent
mask discard (`f9307b3`); mechanism-1 seed-row leak, name/column desync, unlinked OpenMP in the
NN test, failed-fit accounting, dead sort in `NNEnsemble::fit`, backwards FDR-floor message,
overstated `medianGap`, fixture precision, unknown-arm silence (`8955640`).

---

## 1. M5 compares peak ROWS, not precursor identities — MEDIUM

`src/odia_lda.h`, the `stop_on_composition` block.

The positive set is compared as a set of **best-row indices**. When a precursor's best peak group
changes between iterations — which is a normal, healthy thing for the model to do — the row index
changes, so one precursor leaving and re-entering looks like a removal *plus* an addition. Jaccard
is depressed and the shrink test can fire on a set that did not actually lose a single precursor.

**Fix:** compare `group[best_row]`, not `best_row`. The rule is *about* precursors; it currently
measures peaks.

**Why not tonight:** it changes what the stop rule measures, so it needs its own test (a fixture
where the best peak swaps but the precursor set is constant must report Jaccard 1.0) and a re-run
of the M5 arms. Both M5 arms are in the running ablation under the current semantics; changing it
mid-flight would make the table incomparable.

## 2. M5 commits the fit before testing for collapse, and skips can bypass the detector — MEDIUM

Same block. Two structural problems:

* The new model is installed by `fit_learner()` *before* the composition test runs, so on the
  iteration where a collapse is detected the collapsed model is already the live one. The rule
  stops further damage but does not roll back.
* The `continue` paths above it (`positive_rows.size() < m + 2`, `negative_rows.size() < 2`) exit
  the iteration before the composition test, so the most severe shrinkage — the kind that starves
  the fit entirely — never reaches the detector.

**Fix:** evaluate composition *before* committing, keep the last non-collapsed model, and run the
test on the skip paths too.

## 3. Fold normalisation is a label-dependent transform of the test set — MEDIUM, deepest item here

`src/odia_lda.h`, `normalize_folds`.

Each fold's held-out scores are standardised to *that fold's own decoy* mean and SD, and those same
decoys are then pooled to compute q-values. The standardising statistics are therefore computed
from the labels of the very rows the null is built on. With few decoys in a fold the effect is
stark: two decoys are forced to ±1/√2 while an exchangeable null target is unbounded — so target
and decoy are no longer exchangeable, which is the assumption target-decoy FDR rests on.

This is not new tonight; it predates the NN work and affects LDA and GBT equally. It is also the
item most likely to be quietly costing calibration accuracy.

**Fix directions:** label-blind calibration (standardise on all held-out rows, not decoys only), or
cross-fitted statistics, or pool fold-local valid statistics rather than rescaling. Needs an
entrapment measurement to choose between them — see item 5.

**Why not tonight:** it changes every arm's q-values, so it must not land in the middle of an
ablation, and the choice between fixes is an empirical question this project cannot answer without
item 5.

## 4. Harness reports 6,421 where the pipeline reports 6,433 — LOW, but unexplained

Same fixture, same classifier, same parameters; 12 precursors (0.19%) differ.

Most likely best-per-group tie-breaking: `odia_ablate.cpp:idsAt()` keeps the first maximum in row
order, and `bestPerGroup_()` in the tool may resolve equal d-scores differently. Harmless for the
ablation (it cancels across arms) but it is an unexplained difference between two things that
should be identical, and those have a history in this project of turning out to matter.

**Fix:** make `idsAt()` use the same tie-break as `bestPerGroup_()`, then assert equality.

## 5. External null validation — the one that gates items 2 and 3

Every internal diagnostic here is necessary and none is sufficient: `medianGap()` rises for genuine
discrimination *and* for label leakage, and it cannot distinguish them. The scheme's own header
says so — a uniformly biased seed is invisible to all five mechanisms.

The entrapment library is **already built** (177,763 promoted decoys). Running it gives an external
answer to "is the reported 1% actually 1%", which is what decides whether the fold-normalisation
change in item 3 helps or merely moves numbers.

**Do this first tomorrow.** It is the only measurement that can adjudicate the rest.

## 6. Bagging preserves the class ratio in expectation only — LOW (wording, possibly design)

`src/odia_anchor_training.h:trainBaggedOnAnchors`.

Positives and negatives are independent threshold samples at the same rate, so the retained ratio
matches only in expectation, not exactly; and XORing a constant into the seed is *domain
separation*, not probabilistic independence. Neither is wrong in effect — the per-bag class
weighting in `NNEnsemble::fit` compensates — but the comment claims more than the code delivers.

**Fix:** either say "equal expected sampling rate", or sample a fixed stratified count per class.

## 7. NN scoring is still the serial tail — DONE (`8f6a44c`, `5c03c12`)

All three per-group scans now share one parallel helper with a pre-sized output written by
position; the duplicated decoy re-scan is gone; and the nested thread budget (`max_threads/folds`)
is applied by every inner region rather than by training alone, which was leaving 3 folds x 180
threads on 224 cores.

Still open underneath it: measured CPU went 435% -> 699% with the scans parallelised but before the
oversubscription fix. **The post-fix number has not been measured yet** — check it on the next run
before assuming the fix worked.

## 9. `assert()` in tests — AUDITED, `odia_nn_test` fully converted

Swept every `src/*_test.cpp` for side-effecting calls inside `assert`. **`odia_nn_test.cpp` was the
only file affected**, and it had two more beyond the one `de38a1b` fixed — the XOR-capacity and
refit-determinism blocks, whose `fit()` calls also vanished under `-DNDEBUG`. The file now contains
no `assert` at all; every check is an always-compiled `CHECK` that reports its line and sets the
exit code. Proved by injecting a failing threshold and confirming a Release build exits 1.

No other test file puts a side-effecting call inside an assert. They do still lose their checks
under `-DNDEBUG`, which is the ordinary cost of assert-based tests and not the same defect — but
`ctest` runs the Release binaries, so **their assertions are not actually being evaluated on the
cluster either**. Converting them is mechanical and worth doing.

## 10. Batch-size tuning came from a synthetic — NEW

`batch_size = 256` was chosen from a measured curve on a 7k/600k Gaussian synthetic that is
LINEARLY SEPARABLE. It establishes the direction (updates matter, ~3x cheaper via smaller batches
than via more epochs) but not the magnitude for real sub-scores. Re-tune on the benchmark fixture
once the ablation has a working baseline.

## 11a. THE SEED-STEP PROXY DID NOT TRANSFER — new, and it invalidates a tuning claim

`f68ad14` tuned the network on a seed-step proxy (targets clearing the 99th percentile of decoys,
held-out groups of one subsample) and reported it beating the GBT 2.09% to 1.99%. End-to-end on the
full fixture:

| arm | IDs@1% |
|---|---:|
| gbt | 6,421 |
| nn (tuned) | **5,236** |

**18% behind**, where the proxy predicted ahead. The architecture fix is real and necessary — the
network went from 0 IDs and 0/9 iterations trained to 5,236 and 9/0 — but "beats the GBT" was an
artefact of the proxy, the tuning set, or both.

Two distinguishable causes, cheapest test first:
1. **The proxy is wrong.** It measures ONE fit on seed rows; production runs three semi-supervised
   iterations per fold, and a model that starts better can end worse by selecting a narrower
   positive set to retrain on. Test: instrument the per-iteration positive-set size and composition
   for both learners on the same fixture. This also directly exercises mechanism 5.
2. **Overfitting to the tuning subsample.** ~20 configurations, best kept, one subsample, margin
   0.10 points. Test: re-evaluate the top 3 configurations on a DIFFERENT subsample (`id%13==1`)
   without retuning.

Until one of these is settled, `-classifier gbt` remains the honest default and the NN defaults in
`odia_nn.h` are "best known", not "validated".

## 11. THE NETWORK DOES NOT WORK ON REAL SUB-SCORES — CAUSE FOUND (`f68ad14`), see 11a for what remains

`nn` returns 0 identifications where `gbt` returns 6,421, and the cause is upstream of everything
the five mechanisms do: the seed model never produces 26 confident positives, so all 9
fold-iterations skip. Measured on held-out groups with `odia-nn-diag`:

| model | held-out AUC | targets above 99th pct of decoys |
|---|---:|---:|
| one-feature bootstrap | 0.5118 | 1.53% |
| GBT | 0.5250 | 1.99% |
| NN (defaults) | 0.5159 | **1.38%** |

**At the top of the ranking the network is worse than the single feature it was seeded from**, and
its logits sit in [-0.15, 0.12]. AUCs near 0.5 are expected — most target candidates are genuinely
false — so the whole contest is decided in the last percentile, which is what the third column
measures and what an AUC hides.

Refuted already (do not re-try): heavy-tail saturation (clip3/tanh/asinh move AUC <0.002);
underfitting from a low learning rate (lr 0.5 and 2.0 are WORSE); histogram binning (see below).

Hypotheses still open, cheapest first:
1. **Twelve nets averaged from near-identical starts.** Each member barely moves from Xavier init in
   one epoch, and averaging 12 such nets shrinks the logit range further. Test: report the AUC of a
   SINGLE member against the ensemble.
2. **Depth.** 5 tanh layers on 24 features, one epoch. Test 1-2 hidden layers.
3. **The loss is wrong for the objective.** Logistic loss optimises the whole distribution; the FDR
   only cares about the top percentile. A ranking loss on the top-k may be the real answer.
4. **Class weighting.** `w_neg = n_pos/n_neg` is ~1.0 at the seed step (all targets vs all decoys),
   so it is not the imbalance — but the positives are ~98% noise, which is a different problem
   from imbalance and may need a robust loss.

## 12. THE SCORE FIXTURE IS WRITTEN CLASS-ORDERED — hazard, not yet a known bug

`-score_fixture` emits all decoy rows before all target rows (measured: first 200k rows 100% decoy,
last 200k 100% target), because it walks the FeatureMap in its natural order. Anything downstream
that is order-dependent therefore has the label available as a positional signal.

This already caused one false finding tonight — a rank transform that broke ties by sort position
reported a held-out GBT AUC of 0.9811 instead of 0.5250, purely by reading the file layout, and it
survived a group-wise held-out split because the leak was inside the feature rather than across the
split.

`OswRows::canonicalize()` exists to defend against exactly this class of problem in the production
path. **Check whether the fixture dump happens before or after it**, and either canonicalise before
dumping or interleave the classes.

## 8. Carried over, unrelated to tonight

* Report upstream: needless deep copy of every `Feature` inside `omp critical (osw_write_out)`.
* Report upstream: native string IDs used as join keys where indices would do (see also
  `docs/OSWPQ-parquet-library-format-spec.md` §7.2).
* `-compact_library true` still yields zero prefilter support, undiagnosed.
* `library_intensity` float32 change not contributed upstream.
