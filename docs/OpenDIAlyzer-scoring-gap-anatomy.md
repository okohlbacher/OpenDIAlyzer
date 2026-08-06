# Why we lose identifications at scoring

**2026-08-06.** Measured on the Astral benchmark, `-classifier gbt`, idle 224-core node, against the
corrected DIA-NN reference (`bench/diann_ids_correct.txt`, 8,812 keys). Adversarially reviewed by
Codex and Kimi; both reviews are cited where they changed the conclusion.

The prefilter is settled: an oracle that force-admits every deleted reference precursor yields
**+215 precursors and −45 proteins**. So the candidate set is not the constraint. This document is
about the stage that is.

## 1. Where the gap is

| | count |
|---|---:|
| reference precursors | 8,812 |
| removed before scoring (prefilter) | ~1,370–1,430 |
| reaching scoring | 7,442 |
| ODIA identifies | 5,819 of those |
| ODIA-only (DIA-NN misses) | 1,189 |
| **the gap** | **1,623** |

**A correction to an earlier reading of mine.** I first reported the 1,370 as "absent from ODIA's
library entirely — not searchable", and floated it as possibly the real story. It is not.
`experiments/id_gap.py` reads `library/precursors.parquet` **from the output bundle**, and ODIA
replaces `transition_exp` with the retained pairs before writing it (`opendialyzer.cpp:4127`,
`:5964`) — the bundle library is the **post-prefilter** set, 423,079 precursors, not the 7.1M
input. So that bucket is prefilter deletion wearing a library-absence costume: the same population
the oracle already bounded at +215. Kimi caught the contradiction against
`prefilter-loss-anatomy.md:38` ("all 8,812 present as candidates"); Codex independently confirmed
the resolution. **There is no separate 15.5% library-coverage problem.**

By best q, the 1,623 break down as 506 (31.2%) at q ≥ 0.50, 472 (29.1%) at 0.01–0.05, 414 (25.5%)
at 0.05–0.20, 231 (14.2%) at 0.20–0.50. **886 of 1,623 (54.6%) are near-misses.** No
"never-extracted" bucket exists — every searchable reference precursor was extracted and scored.

## 2. The missed precursors are not scoring like noise

Best score per precursor:

| stratum | n | p10 | median | p90 |
|---|---:|---:|---:|---:|
| decoys | 206,173 | −0.846 | −0.296 | 1.170 |
| all targets | 208,052 | −0.828 | −0.207 | 1.996 |
| reference FOUND | 6,426 | 8.143 | **9.730** | 10.597 |
| reference MISSED | 1,681 | −0.470 | **4.333** | 6.936 |

decoy IQR = [−0.641, 0.321], decoy p99 = 3.966.

| | inside decoy IQR | above decoy p99 |
|---|---:|---:|
| all targets | 46.5% | 5.4% |
| reference MISSED | **13.7%** | **53.7%** |
| reference FOUND | 0.0% | 100.0% |

**Over half the missed precursors score above the decoy p99.** They sit in a band between the
decoy tail and the found population. The classifier separates the *clean* population perfectly
well — every found reference is above decoy p99. The problem is the band, and what holds the
threshold above it.

## 3. What holds the threshold up

Operating threshold at 1% FDR: **score 7.248**, with 7,012 targets and **70 decoys** above it.

| threshold | missed refs recovered | decoys above |
|---:|---:|---:|
| 7.248 (current) | 4 | 70 |
| 6.523 | 330 | 202 |
| 5.436 | 609 | 626 |
| 3.966 | 903 | 2,063 |

**Codex's correction, which I had wrong:** 70 extreme decoys is not an anomaly — at 1% FDR with
~7,000 discoveries, roughly that many decoys *must* occupy the reported tail. It is arithmetic,
not pathology. And 7.248 is a fold-normalised discriminator, not seven Gaussian σ. The real
question is **whether those decoys faithfully represent false targets, or exploit decoy-specific
artefacts.** Codex also notes the counting SE on 70 events is ~12% (±23% at 95%) before
clustering, so no conclusion here survives a single decoy realisation — it must be repeated across
independent shuffles.

To reach 6.523 legitimately, ~130 of the 202 decoys there (64%) must disappear *while the
recovered targets survive*. That is the whole problem, stated exactly.

## 4. The mechanism: interference, and how decoys acquire real evidence

Both reviewers put this first, and the same mechanism explains **both** arms of the gap.

ODIA's decoys are shuffled sequences with **`product_mz_shift = +20.0`** and
`precursor_mz_shift = 0.0` (`src/opendialibgen.cpp:650-654`). So every decoy fragment is extracted
from an m/z region densely populated by real fragments of ~200k library precursors, in the decoy's
target's own SWATH and RT neighbourhood. If a co-eluting interferent lands in several of a decoy's
bins, that decoy acquires **genuine** chromatographic evidence — real peaks, real co-elution among
"its" transitions, real S/N, real mass accuracy. Only `library_corr` / `library_dotprod` stay low,
and with ~24 sub-scores a classifier can trade one weak feature against several strong ones.

ODIA computes every sub-score on **raw traces**: `background_subtraction="none"`
(`opendialyzer.cpp:2807`), and a repo-wide grep of `src/` for interference handling returns
nothing.

**Supporting measurement.** RT distance from each decoy to the nearest confident (q<0.01) target
in the same isolation window:

| group | median dRT | within 2 s |
|---|---:|---:|
| top 70 decoys (set the cutoff) | **0.21 s** | **95.7%** |
| decoys 71–500 | 0.47 s | 81.6% |
| decoys 501–2000 | 0.61 s | 78.5% |
| random bulk decoys | 0.78 s | 73.3% |

The threshold-setting decoys co-elute almost exactly with confidently identified peptides. This is
consistent with interference but does **not** prove it — co-elution with a real peptide is also
what a dense region looks like. §6 gives the cheap test that discriminates.

### Amplifiers, all verified in code

- **Winner's curse.** Up to 5 candidate peak groups per precursor (`opendialyzer.cpp:2796`) and the
  max is taken (`odia_lda.h:965`), so 206,173 decoys contribute ~1M peak-group opportunities. This
  is not FDR-invalid if opportunities are exchangeable across labels — but it fattens the decoy
  tail, which is exactly what sets the threshold.
- **GBT tiny leaves — a real defect.** `min_child_rows = 20` is documented as "minimum rows in a
  child" (`odia_gbt.h:61`) but **is not enforced per child**: `nL` is declared at `:383`, never
  incremented, and discarded with `(void)nL;` at `:390`. Only the *parent* is checked
  (`:407`). Children are constrained solely by `min_child_weight = 1.0` on the Hessian sum, which
  permits leaves of ~4 rows. That is an overfitting surface precisely in the tail.
- **Target–decoy pairs split across CV folds.** Folds are assigned per precursor group
  (`odia_lda.h:439`), independently for a target and its paired decoy. With 3 folds, ~2/3 of decoys
  are scored by a model that could have trained on their correlated parent.
- **Label-dependent fold normalisation.** Held-out scores are rescaled to *that fold's own decoy
  null* (`odia_lda.h:924-959`), and those decoys then also build the pooled null — already recorded
  as a known issue (`classifier-backlog.md:45`). Both reviewers rate it real but too small to
  produce the tail alone.

## 5. How DIA-NN differs (published literature only)

From the Nature Methods paper and its supplement. **No DIA-NN source was read** — its licence
forbids deriving from it and this project contains none of its code.

| | ODIA | DIA-NN |
|---|---|---|
| sub-scores per peak group | 24 (36 with `-ms1_scores`) | **73** |
| interference handling | **none** | identification-level arbitration + profile clipping |
| competing-peak context | none in pipeline | 2 dedicated scores |
| isotope scores | 3 (MS2) | 11, incl. testing whether "fragments" are heavy isotopologues of something lighter |
| co-elution | 1 tolerance | 3 mass tolerances, vs a designated best fragment |
| peak shape | EMG fit | 5 window segments |
| RT residual | `var_norm_rt_score` | apex RT + sqrt&#124;measured − predicted&#124; |
| classifier | LDA or 120-tree depth-4 GBT | linear stage + **12-DNN ensemble** |
| mass/RT windows | fixed | iteratively tightened from high-confidence IDs |

DIA-NN's interference correction has two published parts. **Identification-level:** precursors
matched to the *same retention time* that share interfering fragments are arbitrated, and only the
highest-scoring one is reported. **Quantification-level:** a best fragment is chosen as the one
maximising summed correlation with the others, smoothed into a reference profile, and every other
fragment's profile is clipped at 1.5·r·ref — explicitly independent of library intensity accuracy.

The authors attribute their advantage to the classifier plus interference correction: *"The use of
neural networks allows all 73 scores calculated for each elution peak to be used effectively"*,
with the linear classifier explicitly identified as the OpenSWATH/mProphet family ODIA belongs to.

**Which difference matters most here: interference correction, and it is not close.** Kimi's
argument, which the data supports: ODIA's classifier already separates the clean population fine
(found references are 100% above decoy p99), so the binding constraint is the null *tail*, and
interference is the only listed mechanism that manufactures real-looking decoy chromatograms. It
also uniquely explains the target side — missed references diluted to median 4.333 by interfered
fragments. One mechanism, both arms.

**The DNN ensemble is not the answer, and this was already measured.** ODIA built a 12-network
tanh ensemble of DIA-NN's published shape and it came in **5,236 against GBT's 6,421 — 18% behind**
(`classifier-backlog.md:136-144`). Copying the classifier without the evidence does not work.

Second choice is DIA-NN's iterative window tightening, which shares the "shrink the decoy null"
logic by reducing wrong candidate peaks per precursor.

## 6. What to do

### Phase 0 — falsify before building (cheap, existing data)

Codex's point: an O(P·F) pass over the **existing subordinate rows** can kill the hypothesis
before any trace scorer is written. Per peak group compute effective fragment count
`(ΣAᵢ)²/ΣAᵢ²`, largest-one/two-fragment area share, apex-RT MAD, FWHM CV, mass-error MAD, and
leave-one-fragment-out library-agreement stability.

Compare the **202 decoys above 6.523** against **reference misses in the same score band**, matched
on charge, precursor m/z, RT, transition count and library intensity. Matching matters: comparing
against all missed references, or against 200k easy decoys, mostly rediscovers the classifier.

- Tail driven by one or two transitions, discordant widths, shoulders → **interference confirmed**.
- Stable split-half support from several fragments with aligned apexes and mass errors → **refuted**,
  and the whole plan below is wrong.

### Phase 1 — one label-blind score, which is both fix and measurement

Per candidate peak group, using the retained `ChromStore` (1 B/point, shared RT axis,
`odia_chromstore.h`): correlate each fragment trace against the group's median trace and against
its expected library intensity ratio; mark failures as interfered; rescore on the survivors.

Predicted behaviour, which is the test: **missed references lose 1–3 fragments and gain score on a
coherent core**; **interference-lit decoys collapse**, because their signal is one interferent
whose intensity ratios match the interferent's peptide, not the decoy's shuffled library — so they
fail the ratio arm globally rather than fragment-wise.

This is why "lift targets vs suppress the tail" is a **false dichotomy**: the same label-blind
filter does both at once. Lifting targets without touching decoys would need ~3 score units of new
separation in exactly the band where false targets are most enriched — the hardest place to find
signal.

### Non-negotiable guards

- **The filter must be a pure function of the observed data, never of the label.** Anything trained
  on, tuned on, or *selected by* its effect on decoys breaks exchangeability: q stays 0.01 while
  true FDP climbs. Choosing between two label-blind filters by their decoy counts is still adaptive
  selection on the null.
- **Gate on entrapment FDP, not ID count** — `-entrapment_tag` already exists
  (`odia_fdr.h:219-233`). Run it on every variant *before* looking at identifications.
- **`-fdr_pi0` is out of bounds.** It offers the largest single-jump ID gain in the codebase (the
  code comment itself measures 22,959 → 37,539), and the uniformity assumption it needs is violated
  exactly by the mid-band pile-up in §2. It is the most seductive fake gain available.
- **Repeat across independent decoy shuffles.** 70 tail events carry ±23% at 95%. Decoy shuffles
  are seeded with `time(nullptr)` in release builds (`opendialibgen.cpp:616-618`), so replicates
  come from rebuilding the library.
- **Keep target–decoy pairs in the same CV fold**, and equalise peak-group opportunities across
  labels, before attributing anything to the new score.
- **Hold out a second dataset.** Tuning against this one Astral run until ODIA beats 9,261 is
  overfitting to one sample.

### Independent defects to fix regardless

1. **`min_child_rows` is not enforced** (`odia_gbt.h:383-390`). Real, verified, and it lives in the
   tail this analysis is about.
2. **`-ms1_scores` must be re-measured.** The recorded verdict (+68 against a ±83 "noise floor",
   called no effect) is **stale**: that noise floor was a fold-assignment determinism *bug*, since
   fixed. Several "no effect" verdicts sit inside it. MS1–MS2 coelution contrast and MS1 isotope
   agreement are computed and then discarded by default — and MS1 evidence is exactly what an
   interferent at a different precursor mass fails.
3. **`main_var_xx_swath_prelim_score` never reaches the classifier** — it starts with `MAIN_VAR_`
   and fails the `VAR_` prefix test in `scoreColumnsOf_` (`opendialyzer.cpp:1937`). OpenSWATH's
   composite main score is silently absent from ODIA's feature set.
4. **EMG peak-shape scoring is probably ON, not off.** `makeFeatureFinderParam_` does
   `ff.remove("Scores:use_elution_model_score")` (`:2828`) believing it disables it, but
   `DefaultParamHandler::setParameters` restores removed keys to their default, which is `true`.
   The TOPP tool removes *and separately sets* it; ODIA copied only the removal. Inferred from
   code, corroborated by the column arithmetic; worth confirming by measurement.
