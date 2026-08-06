# Anatomy of the prefilter's losses

> **CORRECTED 2026-08-06, after adversarial review (Codex).** The first version of this document
> used `bench/diann_ids.txt` (7,787 entries) as the reference. That list is **known-bad and this
> repository already documented it** — `docs/OpenDIAlyzer-classifier-backlog.md` §14: it omits
> **1,230 real DIA-NN 2.0 IDs** and contains **205 entries that were not IDs at q<0.01**. The
> corrected list is `bench/diann_ids_correct.txt` (8,812 stripped-sequence+charge keys). §1–§3
> below are re-run against it. **§4 and §5 have been withdrawn**: they rested on evidence columns
> that do not mean what I assumed — see §4.
>
> The headline conclusion is unchanged by the correction, which is the one reassuring part.

**2026-08-06.** Measured on the Astral benchmark against the corrected DIA-NN reference list
(`/scratch/kohlbach/bench/diann_ids_correct.txt`, 8,812 precursors, same run, same library). Evidence
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
| 1 | 4 | | 75,563 |
| 2 | 335 | | 1,645,520 |
| 3 | **1,091** | | 914,675 |
| 4 | 1,528 | | 107,686 |
| 5 | 1,757 | | 8,059 |
| 6 | 4,097 | | 4,768 |

**1,430 reference precursors are deleted by the default rule. 76.3% of them sit at depth 3 — one
fragment short. 23.4% at depth 2. Only 4 (0.3%) below depth 2.**

All 8,812 reference precursors are present in the dump as candidates, so none of the gap is a
naming or modification artefact.

This is the answer to the question that gated the whole line of work, and it is **robust to the
reference correction** — it was 99.7% against the bad list too, and the shape barely moved:
**99.7% of the losses have at least 2 matched fragments**, so evidence that needs matched fragments — a spectral dot product
against the library, a mass-error spread across matched fragments — is *defined* for essentially
all of them. Had they sat at depth 0–1 there would have been nothing to score.

## 2. Two hypotheses killed outright

- ~~**`depth >= 4` but only ONE qualifying spectrum: 0.**~~ **WITHDRAWN.** `ms2_qualifying_spectra`
  counts spectra meeting the *configured* threshold, which in this dump is 1 — so the measurement
  actually says "no precursor at depth >= 4 had only one spectrum with at least ONE hit", which is
  near-trivially true and not the intended test. This axis is **untested**, not empty. Testing it
  needs a dump at `min_fragments=3` or separate >=1/>=2/>=3 counters.
- **Losses with MS1 evidence despite failing MS2: 0.** This one stands — `ms1_hit_count` is
  accumulated independently of the MS2 threshold. There is no MS1 rescue route here.

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

## 4. WITHDRAWN — the ranking experiment and the "pool purity" reading

The first version reported that ranking the depth-3 pool by `ms2_qualifying_spectra`,
`ms2_hit_count`, `ms2_sum_intensity`, `ms2_max_intensity` and `ms1_max_intensity` captured
22/23/42/42/76 of the losses against a random expectation of 89 — i.e. **every ranking worse than
random** — and concluded from a follow-up test that this reflected *pool purity*: the noise is
high-intensity, so intensity ranking promotes coincidences.

**Both are withdrawn. The columns do not mean what I assumed.** The accumulator increments
`ms2_qualifying_spectra` only when a spectrum has already met the *configured* `min_fragment_hits`
— stated in the vendored patch itself:

> Called only when this spectrum already met min_fragment_hits, so this counts SPECTRA in which the
> precursor qualified -- i.e. how often the evidence RECURRED.

This dump was taken at `min_fragments=1`. So "qualifying spectra" counts spectra with **at least
one** fragment hit, and `ms2_hit_count` / the intensity accumulators aggregate that same ≥1-hit
background across the whole run — potentially hundreds of unrelated events. **None of it is the
evidence belonging to the spectrum that established depth 3.** I ranked candidates by a whole-run
abundance proxy and then interpreted the result as a statement about coincidence structure.

The same defect voids the depth-stratified median comparison that produced the "ratio tracks pool
purity" reading. Whether real weak peptides genuinely sit beneath abundant coincidences is, as of
now, **unmeasured**.

What it would take to measure it properly: a dump at `min_fragments=3` (so the accumulators count
the ≥3-hit event), or separate ≥1/≥2/≥3 counters, and evidence recorded from the *argmax-depth*
spectrum rather than summed over the run.

Two further confounds, raised in the same review and not yet controlled:

- **Depth < 4 is not the same as deletion.** The pair-union keeps a target if *either* it or its
  decoy passes, so some "losses" were searchable through their partner. Separating "own target
  passed" / "partner rescued it" / "actually extracted" / "identified" needs a join against a
  default-threshold dump.
- **Max-over-modified-forms is the wrong unit.** Forms differ in precursor m/z, fragments and RT;
  taking the max gives groups with more enumerated forms more chances, and presence of *some* form
  does not show the reference peptidoform was present.

## 6. Target/decoy enrichment by depth — reference-free, and the number to beat

The dump now carries both arms (`prefilter_evidence_decoy_`) at `min_fragments=3`, so the
enrichment can be computed **without any reference list**. That matters: our DIA-NN reference is
positive-unlabelled — it omits real IDs, contains non-IDs, and ODIA finds many precursors DIA-NN
does not — so every recall number here is suspect. Decoys are shuffled by construction, so target
excess at a given score IS the discriminating power, and it is the ratio the FDR is computed from.

7,149,966 candidates: 3,603,425 target + 3,546,541 decoy.

| depth | targets | decoys | T/D | target frac |
|---:|---:|---:|---:|---:|
| ≤3 (bulk) | 2,568,237 | 2,544,186 | 1.009 | 50.2% |
| 3 | 914,675 | 893,071 | 1.024 | 50.6% |
| 4 | 107,686 | 103,449 | **1.041** | 51.0% |
| 5 | 8,059 | 5,662 | 1.423 | 58.7% |
| 6 | 4,768 | **173** | **27.561** | **96.5%** |

Cumulative, which is how the rule actually operates:

| cut | targets | decoys | T/D | target frac | excess (T−D) |
|---|---:|---:|---:|---:|---:|
| depth ≥ 3 | 1,035,188 | 1,002,355 | 1.033 | 50.8% | 32,833 |
| **depth ≥ 4** (today) | 120,513 | 109,284 | **1.103** | **52.4%** | 11,229 |
| depth ≥ 5 | 12,827 | 5,835 | 2.198 | 68.7% | 6,992 |
| depth ≥ 6 | 4,768 | 173 | 27.561 | 96.5% | 4,595 |

Three things follow, and they reframe the design question:

1. **At the operating point the criterion is barely better than a coin flip** — 52.4% target. This
   corroborates, from a completely independent direction, the 212,292/210,787 = 1.007 ratio of the
   retained set.
2. **The discrimination is concentrated almost entirely at depth 6.** Depth 3 (1.024) and depth 4
   (1.041) are noise; the *marginal* admissions the rule makes at depth exactly 4 are 51.0%
   target / 49.0% decoy.
3. **Yield per slot is wildly uneven.** The rule spends 120,513 target slots for an excess of
   11,229 — 9%. Depth 6 alone returns 4,595 excess from 4,768 candidates — 96%. So 41% of all the
   discriminating signal sits in 4% of the admitted candidates.

**52.4% is the number any new score must beat at the same budget.** It is measurable offline from
this dump, before any extraction run, and it does not depend on DIA-NN being right about anything.

## 7. What this does and does not license

It licenses building the scale-free terms: they are defined for 99.7% of the losses, and the
budget to spend on them already exists.

It does **not** license optimism about magnitude. The pool is ~1000:1 noise at depth 3 (786,187
candidates for 766 reference precursors), 86.5% of searchable misses are scoring-side rather than
detection-side, and 46.8% of target rank-1 features already lie inside the decoy interquartile
range. Recovering a precursor into the search space is necessary, not sufficient — it still has to
survive a classifier that is itself the binding constraint.
