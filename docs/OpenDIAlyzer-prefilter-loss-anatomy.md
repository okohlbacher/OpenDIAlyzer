# Anatomy of the prefilter's losses

**2026-08-06.** Measured on the Astral benchmark against the DIA-NN reference list
(`/scratch/kohlbach/bench/diann_ids.txt`, 7,787 precursors, same run, same library). Evidence
dumped at `-prefilter_min_fragments 1` so that **sub-threshold** hit counts are recorded — at the
default 4 the filter short-circuits and every failure reports depth 0, which makes the losses
unattributable. `experiments/run_pfdiag1.sh`; dump is 382 MB, 3,603,425 candidate rows.

Analysis scripts are in the scratchpad (`depth_hist.py`, `depth_pop.py`, `ceiling2.py`,
`falsify.py`). Keys are stripped-sequence + charge; a precursor's depth is the **max** over its
modified forms, because a precursor survives if any of its forms passes.

## 1. Where the losses actually sit

`depth` = `ms2_best_fragment_hits`, the count of the top-6 predicted fragments matched within the
top-1000 peaks of a single spectrum, at ±5 ppm. The rule keeps `depth >= 4`.

| depth | reference precursors | | full target population |
|---:|---:|---|---:|
| 0 | 0 | | 847,154 |
| 1 | 3 | | 75,563 |
| 2 | 188 | | 1,645,520 |
| 3 | **766** | | 914,675 |
| 4 | 1,237 | | 107,686 |
| 5 | 1,612 | | 8,059 |
| 6 | 3,981 | | 4,768 |

**957 reference precursors are deleted by the default rule. 80.0% of them sit at depth 3 — one
fragment short. 19.6% at depth 2. Only 3 (0.3%) below depth 2.**

All 7,787 reference precursors are present in the dump as candidates, so none of the gap is a
naming or modification artefact.

This is the answer to the question that gated the whole line of work: **99.7% of the losses have
at least 2 matched fragments**, so evidence that needs matched fragments — a spectral dot product
against the library, a mass-error spread across matched fragments — is *defined* for essentially
all of them. Had they sat at depth 0–1 there would have been nothing to score.

## 2. Two hypotheses killed outright

- **`depth >= 4` but only ONE qualifying spectrum: 0.** The single-qualifying-spectrum criterion
  removes nothing on this run. `experiments/pf_lost.py` lists it as a distinct recovery axis; it
  is empty.
- **Losses with MS1 evidence despite failing MS2: 0.** There is no MS1 rescue route here.

## 3. The budget already contains the room

The prefilter keeps **120,513** evidence-bearing targets (`depth >= 4`). But a pair-union then
adds the decoy partner of every kept target *and* the target partner of every kept decoy, giving
**212,292 targets + 210,787 decoys = 423,079** retained.

So **~91,779 retained targets — 43% of the target budget — have no MS2 evidence of their own.**
They are retained only because their shuffled decoy partner passed. Meanwhile 914,675 depth-3
candidates, containing at least 766 real precursors, compete for zero slots.

Those 91,779 slots are reallocatable **without growing the search space**, which is the one thing
that must not change: a learned prefilter was measured end-to-end and identifications fell
*monotonically* as the budget grew (423k → 6,024 peptides; 707k → 5,855; 1.96M → 5,705), because
admitted candidates enlarge the target-decoy null faster than they add true positives.

Corollary: the retained set is 212,292 target vs 210,787 decoy, a **target/decoy ratio of 1.007**.
At population level the current criterion barely distinguishes a real peptide from a shuffled one.

## 4. Ranking by existing evidence is worse than random

Ranking the 786,187 target precursors whose best depth is exactly 3, taking the top 91,779 (the
reallocatable budget, 11.7% of the pool), and counting how many of the 766 reference losses are
captured. **Random selection captures 89.**

| ranked by | captured of 766 | vs random |
|---|---:|---:|
| qualifying spectra (recurrence) | 22 | **0.25×** |
| ms2 hit count | 23 | 0.26× |
| ms2 summed intensity | 42 | 0.47× |
| ms2 max intensity | 42 | 0.47× |
| ms1 max intensity | 76 | 0.85× |

Every one is *worse than random*. This directly kills the recurrence term (`run_len`) that was
proposed as a cheap co-elution surrogate: on this data it would actively harm.

## 5. The correction: it is pool purity, not "weak losses"

The first reading of §4 was **"the depth-3 pool is dominated by high-abundance spurious matches
and the real weak peptides sit underneath them."** That is half right, and the conclusion drawn
from it was wrong. The falsification test — do reference precursors that *pass* also look weak? —
settles it:

| depth | pool | reference | pool purity | ref/pool median summed intensity |
|---:|---:|---:|---:|---:|
| 3 | 786,187 | 766 | 0.1% | 0.53 |
| 4 | 102,723 | 1,237 | 1.2% | **0.18** |
| 5 | 7,861 | 1,612 | 21% | **0.17** |
| 6 | 4,482 | 3,981 | **89%** | 0.88 |

Reference precursors are weaker than the pool median at **every** depth, and most extremely at
depths 4 and 5 — where they *pass*. And the ratio tracks **pool purity**: at depth 6 the pool is
89% reference precursors and the ratio converges to 0.88.

So the effect is not a property of the losses. **The pool median is dominated by whatever noise it
contains, and that noise is high-intensity** — abundant peaks generate many chance 3-of-6 matches.
Intensity ranking promotes coincidences wherever noise dominates, which is every depth below 6.

The useful consequence: any ranking term that scales with **magnitude** is contaminated by this.
What survives is evidence that is **scale-free by construction** — the agreement between observed
and predicted fragment *ratios* (restricted to matched fragments, since the top-6 are chosen *by*
predicted intensity and an unrestricted dot product just re-encodes depth), and the *consistency*
of mass error across a candidate's matched fragments. Neither is a function of how intense the
peaks are.

## 6. What cannot be measured yet

**The evidence dump contains zero decoy rows** — only the target scan is dumped
(`prefilter_evidence_ = res.evidence` takes the target result). So target/decoy enrichment *at a
given depth* is not measurable, and neither is "target fraction in the top-2N", which was the
intended **reference-free** quality gate — the one that does not depend on DIA-NN being right.
Adding decoy rows to `-prefilter_out` is a prerequisite for evaluating any new score honestly.

## 7. What this does and does not license

It licenses building the scale-free terms: they are defined for 99.7% of the losses, and the
budget to spend on them already exists.

It does **not** license optimism about magnitude. The pool is ~1000:1 noise at depth 3 (786,187
candidates for 766 reference precursors), 86.5% of searchable misses are scoring-side rather than
detection-side, and 46.8% of target rank-1 features already lie inside the decoy interquartile
range. Recovering a precursor into the search space is necessary, not sufficient — it still has to
survive a classifier that is itself the binding constraint.
