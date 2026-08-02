# Where the missing IDs go: evidence, hypotheses, backlog

**2026-08-02.** Measured against a DIA-NN reference ID list (`/scratch/kohlbach/bench/diann_ids.txt`,
7,787 precursors, same `astral.mzML`, same library). ODIA run: `recal600` (6,930 IDs, tuned windows).

Tooling: `experiments/id_gap.py` (gap decomposition), `experiments/calmap.py` (RT reference map).

Scoring is deterministic (five runs at exactly 6430 across schedulers, topologies and loads), so
every count here is exact rather than sampled.

---

## 1. The measured gap

| | count | share of total gap |
|---|---:|---:|
| DIA-NN reference precursors | 7,787 | |
| — absent from ODIA's **prefiltered** library | **923** | **38.8%** |
| searchable by ODIA | 6,864 | |
| ODIA q<0.01 targets | 6,728 | |
| overlap | 5,406 | |
| **ODIA-only** (DIA-NN misses these) | **1,322** | -- |
| **reference-only** (the gap) | **1,458** | 61.2% |

**The ID sets are not nested.** Jaccard overlap is 66%. ODIA finds 1,322 precursors DIA-NN does not,
so this is not "ODIA is a weaker DIA-NN" -- the tools disagree in both directions, which constrains
any single-cause explanation.

### 1.1 Decomposition of the 1,458 searchable misses

| bucket | n | share |
|---|---:|---:|
| FDR: 0.01 <= q < 0.05 | 423 | 29.0% |
| DISCRIMINATION: q >= 0.50 | 341 | 23.4% |
| FDR: 0.05 <= q < 0.20 | 317 | 21.7% |
| DETECTION: never extracted | 197 | 13.5% |
| DISCRIMINATION: 0.20 <= q < 0.50 | 180 | 12.3% |

Grouped: **scoring-side 1,261 (86.5%)**, detection-side 197 (13.5%).

The never-extracted 197 show no systematic signature: m/z median 569.3 against 570.3 for identified
precursors, charge mix 2:108 / 3:80 / 4:9, library RT spanning the full 0.003-0.883 range. Nothing
points at a window or mass-range explanation.

### 1.2 Target and decoy scores barely separate

Rank-1 features from `recal600`:

| | target | decoy |
|---|---:|---:|
| n | 104,242 | 103,250 |
| median score | **-0.205** | **-0.295** |
| p90 | 2.019 | 1.217 |
| p99 | -- | 3.928 |

**46.8% of target rank-1 features lie inside the decoy interquartile range.** Only 5,568 targets
clear the decoy p99. With ~104k target candidates yielding ~6.7k IDs, the classifier -- not the
extractor -- is what sets the yield.

This is expected in principle (most library precursors are genuinely absent from any given run), but
it locates the binding constraint.

---

## 2. Hypotheses, ranked by the evidence behind them

### H1 -- Prefilter recall 88.1%. 923 precursors (38.8% of the gap) deleted before extraction ran
**CONFIRMED.** All 7,787 DIA-NN IDs are present in the FULL library (zero naming artefacts). The
prefilter keeps 201,022 of 2,730,206 target sequence_charge keys (7.4%) and removes 923 of DIA-NN's
IDs -- recall **88.1%**. This is a hard ceiling: scoring and calibration cannot recover a precursor
that was deleted before extraction.

**The criteria (defaults):**

| parameter | value | effect |
|---|---:|---|
| `ms2_top_transitions_per_precursor` | 6 | only the top-6 fragments BY PREDICTED LIBRARY INTENSITY are indexed |
| `ms2_min_fragment_hits` (`-prefilter_min_fragments`) | 4 | 4 of those 6 must hit |
| `ms2_top_peaks_per_spectrum` (`-prefilter_top_peaks`) | 1000 | ...within the top 1000 peaks... |
| `ms2_min_qualifying_spectra` | 1 | ...of a SINGLE spectrum |
| `-prefilter_mz_extraction_window` | 10 ppm | fragment match tolerance |

So: **4 of a precursor's top-6 predicted fragments must co-occur in one spectrum's top-1000 peaks.**
A low-abundance peptide whose fragments do not all reach the top 1000 of any single spectrum fails
and is deleted.

**H1 and H4 are coupled, not independent.** The top-6 are chosen by PREDICTED intensity. If the
library was predicted at the wrong NCE/instrument, the prefilter indexes the wrong six fragments and
hunts for peaks that were never going to be the intense ones. A provenance error would show up here
as prefilter recall loss, and separately as degraded sub-scores -- one root cause, two symptoms.

**Original assessment retained below for the record.**
**Evidence: strong, one test outstanding.** The prefilter cuts 7,149,966 precursors to 423,079
(94.1% removed) before extraction. 923 DIA-NN IDs are absent from the prefiltered library.
**Outstanding:** confirm they are present in the FULL library (i.e. the prefilter dropped them)
rather than absent for naming/modification reasons. That single query decides whether this is the
largest bucket in the whole analysis or an artefact of my sequence-key matching.
**Action:** measure prefilter recall directly -- what fraction of DIA-NN's IDs survive it, and at
what threshold. If recall is ~88%, the prefilter is capping achievable yield at 88% before any
scoring happens.

### H2 -- Scoring discrimination. 1,261 precursors (86.5% of searchable misses) extracted but not separated
**Evidence: strong.** 740 sit at 0.01 <= q < 0.20 -- found, ranked, just under the cutoff. Another
521 score like noise. Target/decoy medians differ by 0.09 on a scale where target p90 is 2.0.
**Actions:**
- H2a. Re-test `-ms1_scores`. Previously measured +68 IDs and dismissed against a +/-83 "noise
  floor" that was itself a bug. With determinism, +68 would be a real gain. **Cheapest test here.**
- H2b. Compare sub-score sets against what DIA-NN's classifier consumes; identify missing evidence
  (co-elution shape, isotope pattern agreement, MS1/MS2 consistency).
- H2c. Classifier capacity: GBT on ~23 sub-scores versus DIA-NN's neural-network ensemble with
  per-run retraining. Test whether more trees/depth moves anything before concluding the features
  are the limit rather than the model.

### H3 -- Interference. Sub-scores computed on un-deconvolved traces
**Evidence: indirect but mechanistically strong.** DIA-NN performs interference correction before
scoring; ODIA does not. Every sub-score is computed on raw extracted traces, so a co-eluting
interferent degrades the score of a genuinely present peptide. This is a candidate root cause *for*
H2 rather than an alternative to it.
**Action:** quantify on the 740 boundary cases -- do their traces show co-elution structure absent
from the confidently-identified set?

### H4 -- Library provenance. Fragment intensities possibly predicted for the wrong instrument
**Evidence: unresolved and cheap to settle.** `library/metadata.json` records no NCE, no instrument,
no model version. `OpenDIALibGen` had `kNCE = 35.0f` / `kInstrument = timsTOF` hardcoded, measured
from a Bruker run. If this library was predicted at those settings for an Astral run, **every**
intensity-dependent sub-score is systematically degraded -- which would show up exactly as H2.
**Action:** grid `-nce` x `-instrument` (now flags) on a library subset and correlate predicted
against observed fragment intensities for confidently-identified precursors.

### H5 -- Detection. 197 precursors never extracted
**Evidence: measured, and it is the smallest bucket.** No systematic m/z, charge or RT signature.
**Action:** low priority given size; revisit if H1 turns out to be a naming artefact and this bucket
absorbs those precursors.

### H6 -- Decoy model. Shuffled decoys may misestimate the null
**Evidence: weak, untested.** 104,242 target vs 103,250 decoy rank-1 features is a near-1:1 balance,
so the decoy library is well-sized. Whether shuffled decoys produce a null that matches the true
false-target distribution is a different question and is unmeasured.
**Action:** entrapment FDP. Requires a library rebuild with entrapment sequences -- the only method
here that yields ground truth rather than another relative comparison.

### H7 -- FDR estimator conservatism
**Evidence: none either way.** If ODIA's q-values are conservative relative to DIA-NN's at the same
true FDP, the 423 precursors at 0.01 <= q < 0.05 are a bookkeeping loss rather than a scoring one.
Cannot be settled by comparing two tools to each other -- both could be wrong.
**Action:** blocked on H6's entrapment set. Until then, treat "our q is too strict" as unfalsifiable
and do not act on it.

### H8 -- RT calibration
**Evidence: largely refuted as a primary cause.** The reference map shows the transform's shape is
correct (residual median within +/-10 s of zero) and the tool's own `delta_rt` tracks the reference
residual almost exactly. Window *widths* were worth +500 IDs and are now tuned. What remains is the
14x RT-dependent quality swing (6.8 s median in the first octile, 95.9 s in the last), which a
scalar window cannot express.
**Action:** RT-segmented windows -- see `docs/OpenDIAlyzer-calibration-handoff.md` §C. Bounded upside
now that the scalar optimum is found.

---

## 3. Adversarial reading

**Against "it is an FDR problem":** 86.5% of searchable misses being scoring-side is real, but 38.8%
of the *total* gap may never reach scoring at all (H1). Fixing the classifier cannot recover a
precursor the prefilter deleted. H1 must be settled first or effort goes to the wrong layer.

**Against "it is a scoring-power problem":** ODIA finds 1,322 precursors DIA-NN misses. A uniformly
weaker classifier does not do that. This looks more like two differently-biased scoring systems than
one strictly dominating the other -- which weakens any argument that simply copying DIA-NN's
approach recovers the difference, and raises the possibility that the union is the interesting
object.

**Against "it is interference":** plausible mechanism, zero direct measurement here. Should not be
prioritised over H1/H2a, both of which are cheap and decisive.

**Against the whole comparison:** the reference is DIA-NN's output, not ground truth. Every "miss"
assumes DIA-NN is right. Its own FDR could be anti-conservative, in which case some fraction of the
1,458 are false positives we are correctly rejecting. **Only entrapment (H6) breaks this circularity**,
and until it exists every number here is relative.

---

## 4. Ranked next actions

| # | action | cost | decisiveness |
|---|---|---|---|
| 1 | Confirm the 923 are in the FULL library (prefilter, not naming) | one query | settles the largest bucket |
| 2 | Measure prefilter recall against DIA-NN IDs | one query | bounds achievable yield |
| 3 | Re-test `-ms1_scores` with determinism | one run | previously +68, dismissed against a bogus floor |
| 4 | `-nce` x `-instrument` grid vs observed fragment intensities | a few runs | settles H4, which would masquerade as H2 |
| 5 | Interference structure in the 740 boundary cases | analysis only | tests H3 without new runs |
| 6 | Entrapment library | library rebuild | the only ground truth; unblocks H6 and H7 |
| 7 | RT-segmented windows | design + implementation | bounded; the scalar optimum is already taken |
