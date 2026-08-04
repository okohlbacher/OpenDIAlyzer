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

## 13. THE .oswpq SCORE TABLE CANNOT BE JOINED TO ITS FEATURE TABLE — FIXED

**Cause: two id conventions in one bundle.** `OpenSwathOSWParquetWriter` writes
`clearSignBit(feature.getUniqueId())` (`OpenSwathOSWParquetWriter.cpp:603`); our score writer cast
the raw `uint64` to `int64`. `UniqueIdGenerator` is effectively random, so half the ids have the
high bit set and came out negative on one side and positive on the other — hence exactly 50%
overlap.

Both earlier hypotheses were wrong and are recorded so they are not re-tried: it is **not** two
different runs, and it **is** a sign issue — but the transform is `+2^63` (clear the sign bit), not
the `+2^64` I tested, which is why that test matched zero rows and wrongly exonerated the idea.

Verified after the fix: **2,070,355 of 2,070,355 score rows join, 100%.**

It cost more than downstream convenience: the entrapment validation could map only 2,328 of 4,588
identifications back to a precursor, halving the sample its FDP estimate rests on.

### original text, kept for the diagnosis trail

In `bench/nnfix/nnfix.oswpq`, `runs/run_id=*/score_ms2.parquet` and `runs/run_id=*/features.parquet`
each hold **14,943,922 unique `feature_id`s and share only 7,472,175** — exactly half. Any consumer
joining scores to features silently loses half the rows.

Diagnosis so far, with two hypotheses already refuted:
* NOT a signed/unsigned wraparound. 7,471,747 score ids are negative and the feature table has none,
  which looks exactly like a uint64 written through `arrow::int64()` — but adding 2^64 to the
  negative ids matches **zero** feature ids, and reinterpreting the whole column as unsigned still
  gives 50.00%.
* NOT two different runs. Both tables have the same row count as the candidate set, and two
  independent 63-bit id sets of 15M elements would collide ~never, not 50% of the time.

Most likely: OpenMS assigns a **new** `UniqueId` when a Feature is copied, so the features written
to disk are not the objects that were scored. That is consistent with the existing upstream item
"native IDs are join keys, should be indices".

**This did not affect any identification count** — the pipeline computes q-values in memory and
never performs this join; its 6,433 is authoritative. It affects downstream consumers and any
external analysis, and it invalidated one gap number I computed tonight before I noticed.

**Workaround that works today:** `score_peptide.parquet` keys on `modified_sequence` and reproduces
the pipeline's own peptide count exactly (5,625). Use it for identity-level comparisons.

## 14. THE DIA-NN REFERENCE USED FOR ALL GAP ANALYSIS WAS WRONG — corrected

`bench/diann_ids.txt` (7,787 entries, 30 July) was used for every gap analysis in this project.
It is not DIA-NN 2.0's real output: against a target-only library DIA-NN 2.0 reports 9,261
precursor rows. Rebuilt from `bench/dn_fair/d20.parquet` as `bench/diann_ids_correct.txt`:

| key space | count |
|---|---:|
| precursor rows (DIA-NN's own headline) | 9,261 |
| unique stripped-sequence + charge | 8,812 |
| unique stripped sequence (peptide) | 7,831 |

Against the old list: 7,582 shared, **1,230 real IDs missing from it**, and 205 entries that were
never DIA-NN 2.0 IDs at q<0.01. So every earlier gap was understated by ~14% and slightly polluted.

**The real, trustworthy comparison** (peptide level, both sides keyed on sequence, using the tables
that reproduce each tool's own reported counts):

| | peptides @1% |
|---|---:|
| ODIA | 5,625 |
| DIA-NN 2.0 | 7,831 |
| shared | 4,766 |
| **DIA-NN only — the gap** | **3,065** |
| ODIA only | 859 |
| **recall of DIA-NN** | **60.9%** |

ODIA is not a subset: it finds 859 peptides DIA-NN does not. Re-do the prefilter/loss attribution
against `diann_ids_correct.txt`; the earlier attribution used the wrong target.

## 15. PASS-2 RE-EXTRACTION: measured, and NOT worth eliminating — CLOSED

The proposal (user's): pass 2 re-extracts the whole library after recalibration, but the underlying
signal is unchanged — only the RT prediction moved — so slice the pass-1 chromatograms instead of
re-reading the run.

**The premise is correct.** Adversarial review verified that a slice is *bit-identical* to a
re-extraction in the vendored extractor. The idea is sound.

**The payoff is not.** Measurement 0 — the extraction/scoring split of pass 2 — from three existing
runs, no new measurement required:

| run | rt_win | pass-2 wall |
|---|---:|---:|
| ms1_on | 600 s | 300.8 s |
| rtfeat_on | 600 s | 303.7 s |
| bestknown | 864 s | 337.4 s |

A 44% wider window costs 11.6% more time. Fitting `t = a + b·W` gives b ≈ 0.133 s per second of
window and a ≈ 222 s, so at 864 s: **extraction ≈ 115 s (34%), scoring ≈ 222 s (66%)**. Removing
re-extraction entirely saves ~115 s of 1,847 s = **6.2% of wall**, against a plan that claimed 337 s
/ 18%. Overstated 3x, because the plan attributed the whole phase to extraction.

Caveat, stated because the fit is thin: two window widths, so `b` is poorly determined. The two
600 s runs agree to 1% (300.8 / 303.7) which is reassuring, and the direction is not in doubt even
if the coefficient is.

**Also refuted along the way** (all from the review, all confirmed):
* Option B (rescore without narrowing) — dead. The ID ladder 240/400/600/720/864/900 →
  6552/6695/6930/6941/6980/6835 is NOT confounded, so the narrower window does real work beyond
  centring.
* Option D (reorder the mass calibration before pass 1) — dead, and it rested on **my error**: I
  claimed a pre-pass mass calibration on CiRT anchors already exists. It does not;
  `calibrateMassFromPass_` is called only between passes. A CiRT-anchored version would have tens of
  anchors where the estimator is already rejected as unsupported at 2,000.
* Option A (retain + slice) — incomplete alone: MS1's window genuinely narrows between passes, so
  MS1 must be re-extracted regardless.

**Decision: not building it.** 6.2% is not worth a retained chromatogram store (which does not
exist), a containment fallback, and an MS1 exception — in a project whose worst metric is memory.
Revisit only if scoring gets much cheaper, which would change the ratio.

## 16. ENTRAPMENT RAN, AND 5% IS TOO SMALL A FRACTION

First external FDR check. Library: 2,223,453 real targets + 86,257 entrapment (5% of decoys
promoted). At q<0.01, rank 1:

| | |
|---|---:|
| real target precursors | 2,301 |
| **entrapment precursors** | **2** |
| observed FDP | 0.09% |
| scaled x20 (entrapment was 5% of the decoy pool) | **1.7%** vs nominal 1% |

**Consistent with calibration, and too weak to say more.** Two events: the Poisson 95% interval on 2
is [0.24, 7.2], so the scaled FDP lies somewhere in ~0.2%-6%. That rules out gross miscalibration --
we are not claiming 1% and delivering 20% -- and cannot separate "well calibrated" from "twice as
bad as claimed".

**Next run: raise the fraction to 25-50%.** The uncertainty here is entirely event-count driven. The
cost is a smaller decoy null, which the script's own header already flags as the trade.

Note the sample was also halved by item 13 (the feature_id join defect), now fixed -- so a re-run
gets both a larger fraction and twice the joinable data.

## 17. -compact_library: SEVEN DEFECTS, the last of them mine and the worst of them silent

The option's help said "not diagnosed; do not enable" for weeks. It is now diagnosed four times over
and still must not be enabled, for a fifth reason that is not a loader defect.

**Fixed, each verified by a load-time assertion that now exists:**

| # | commit | never populated | how it surfaced |
|---|---|---|---|
| 1 | `480b7a3` | `PRECURSOR_MZ`, `LIBRARY_RT`, `DECOY` on precursors | prefilter: "no supported precursors" |
| 2 | `48fd9a8` | transition `DETECTING` flags | prefilter: "0 targets, 0 decoys" |
| 3 | `2b4da1a` | transition decoy flag (wrong source) | prefilter: "118,902 target / **0 decoy**" |
| 4 | `4642e64` | decoy↔target pairing (design, not omission) | prefilter: still 0 decoy |

The first three are the same shape — a setter that exists, is correct, and is never called — which
is why each fix only revealed the next. The fourth is a design error: `syntheticId()` encodes a
peptide's OWN row, but pairing downstream is by string (`"DECOY_" + target id`), so every decoy id
matched no target. Fixed by resolving the pairing at load from `TRAML_ID` and storing it as a 4-byte
index, then releasing the strings.

**Fixed since, and these are the ones worth reading:**

| # | commit | defect | how it surfaced |
|---|---|---|---|
| 5 | `4bc2aee` | `library_intensity` is **float32**, read as `arrow::DoubleArray` | 97% of targets at q<0.01 |
| 6 | this | the **unmodified** sequence stored as the peptide identity | nothing — silent wrong mass |
| 7 | this | null fragment charge defaulted to 1, ordinary uses 0 | nothing — inert today |

**#5 was mine, twice over.** I narrowed that column from float64 to save 654 MB and never revisited
the one path that reads it with a hard cast. An 8-byte stride over a 4-byte buffer:

    library_intensity : float     one chunk, 78,569,077 rows
    valid double reads end at row 39,284,538
    first decoy row               39,589,429

Targets are written before decoys, so the half-stride boundary falls INSIDE the target block:
targets read in-buffer (wrong but finite), every decoy read runs off the end of it. One cause, both
symptoms — in the prefilter intensity only ranks the top-6 fragments, so garbage picked a more
coincidence-prone decoy subset; in scoring it is used quantitatively, and NaN correlation maps to
the score floor (`MRMScoring.cpp:579`), which is a degenerate null.

| prefilter evidence | targets | decoys |
|---|---|---|
| ordinary | 120,513 | 109,284 |
| compact, before | 118,902 | **309,344** |
| compact, after | **120,525** | **109,298** |

Found by codex; kimi reached the same line independently but ranked it second behind a NULL
hypothesis I had already refuted by measurement. The lesson is in `arrowText()`, which exists in the
same file for exactly this reason and which I did not generalise to the numeric columns:
**never cast an Arrow numeric chunk to a fixed width.**

**#6 traded correctness for storage and documented the trade as if it were free.** The line stored
the unmodified sequence with a comment saying the identity for scoring "is still the modified form,
which is why a library needs both" — and then stored only one. `OpenSwathScoring.cpp:349` derives
the peptide formula, and therefore the isotope scores, from `compound.sequence`. On this library
44.1% of precursors carry a modification, and keying on the unmodified form merged 741,887 distinct
identities (3,218,138 modified vs 2,476,251 unmodified).

**What is still not at parity, deliberately:** `LightTransition::fragment_nr` is never set (ordinary
reads the ordinal). Storing it costs 2 B x 78.6M = 157 MB for a field only `MRMAssay`/`MRMIonSeries`
read, and neither is in this pipeline. Also unset: `peptide_group_label`, `rt_start`, `rt_end`,
`gene_name`, `sum_formula`, `compound_name`, and `LightCompound::modifications` (read only by
`MRMDecoy`, which ODIA does not invoke — the library arrives with decoys).

**Separate, still open:** `restoreRealIds_` maps positionally (`exp.compounds[i]` vs `Peptide(i)`)
while post-prefilter `exp.compounds` is a 786k subset of 7.1M. Latent only because `setOriginalId`
is never called in production, so the compact path ships synthetic ids in its output bundle.

**Two guards now stand where these got through:**

1. If more than 50% of retained targets pass q<0.01 the run fails with `UNEXPECTED_RESULT`. At a
   nominal 1% FDR that fraction is impossible with a valid null, and the broken run exited 0.
   Verified to fire on the known-broken build: `381457 of 393753 (96.8772%) ... Exit status: 13`.
2. The load check on library intensities is **per-label by construction**. An aggregate check would
   have passed on #5: it read one label correctly enough to look plausible and the other as garbage.
   Anything that can be wrong for one label only must be checked for each label separately.

## 8. Carried over, unrelated to tonight

* Report upstream: needless deep copy of every `Feature` inside `omp critical (osw_write_out)`.
* Report upstream: native string IDs used as join keys where indices would do (see also
  `docs/OSWPQ-parquet-library-format-spec.md` §7.2).
* `-compact_library true` still yields zero prefilter support, undiagnosed.
* `library_intensity` float32 change not contributed upstream.
