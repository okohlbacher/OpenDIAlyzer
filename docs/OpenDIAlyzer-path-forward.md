# Path forward — sequenced by what can be attributed, not by what looks promising

Written after adversarial reviews of a night's work by two external models (kimi, codex). Both
independently predicted the same failure for the plan as originally stated ("calibration → memory →
occupancy"):

> *optimize a proxy, change the system, then instrument too late to attribute the result.*

Every conclusion reached by inference in that night was wrong or confounded; every conclusion
reached by instrumentation held. The ordering below follows from that, not from expected payoff.

## Standing rules, derived from specific failures

| Rule | The failure it prevents |
|---|---|
| **Measure before intervening, in the same build.** | The prefilter rewrite and its timers landed together, so it has no before-number, permanently. |
| **One variable per run.** | The sqlite-vs-parquet comparison used different extraction runs; the variable was never isolated. |
| **Count the null, not just the benefit.** | 20 ppm screening predicted +800 IDs and delivered −327: reachability was counted, the FDR burden of 1M extra candidates was not. |
| **A gate that passes is not a result.** | Peakedness, inferred ppm and anchor count are internal diagnostics. IDs at controlled FDP are the result. |
| **Nothing under the noise floor is signal.** | Three equivalent configurations gave 6,522 / 6,607 / 6,622. Differences below ~100 IDs have been over-read repeatedly. |
| **Never reset a pre-existing reference checkout.** | `git checkout -- .` in the gitignored OpenMS tree destroyed load-bearing uncommitted work. `setup_node.sh` automated the same act until it was gated. |

## Phase 0 — establish what the numbers mean (blocking, cheap)

Nothing below is interpretable without these.

1. **Noise floor.** Three repeats of one fixed configuration on the same binary. Report mean and
   spread of IDs. Every later claim is measured against this, not against a single run.
2. **Phase-resolved memory.** `PhaseTimer` now reports RSS at phase entry/exit from
   `/proc/self/statm`. One run identifies which phase actually holds peak RSS.
   **This is required before any memory work**: the §2 decomposition of
   `OpenDIAlyzer-memory-review.md` was refuted by its own §6 A/B (a 62% cut of the chromatogram term
   moved peak by 4.6%), and the dominant term is still unidentified.
3. **An exclusive allocation, or `taskset`.** The benchmark node has 18 users and a load average of
   180–233. Until this exists, wall clock and core occupancy measure cores *won*, and no scalability
   claim can be made.

## Phase 1 — quality (IDs), in dependency order

**1.1 Validate the mass calibration that now fires.** The RT-direction fix made it fit for the first
time: MS2 5,709 anchors at peakedness 4.9, MS1 560 at 5.5, MS1 sigma 1.5 ppm against an instrument
independently measured at 1.66–1.71 ppm. That corroborates the *measurement*. It does not establish
that the identifications gained are real: anchors are the same run's top-2,000 targets ranked by the
score the narrowing then improves — selection on the optimised quantity. Required before claiming a
win: **entrapment FDP** (`-entrapment_tag`, already implemented) and cross-fitted anchors.

**1.2 `-ms1_scores` A/B.** 12 MS1 sub-scores are withheld from the classifier. The measurement that
justified excluding them ("35-36 sub-scores gave 4,913 vs 29 giving 6,607") is now explained by the
`library_rt` defect — the same path reaches 6,522 with no MS1 scores at all. The exclusion rests on
a confounded comparison and must be re-measured, not assumed either way.

**1.3 MS2 anchor quality.** MS2 sigma came out at 4.5 ppm on a 1.66 ppm instrument, so the MS2
anchor picking still admits interference, and its 13.4 ppm estimate was correctly rejected by the
never-widen guard. Note the trap: widening the anchor **search** is not the same as accepting a
wider **estimate**. Conflating them produces "fixes that change inputs but not the symptom" — which
is what three earlier mass-calibration fixes did.

**1.4 Named RT transform types.** Two directional defects in one subsystem in one session
(`library_rt`; `native_trafo` at *two* call sites). `LibraryToRun` / `RunToLibrary` as distinct
types makes the class unrepresentable rather than fixed case by case.

## Phase 2 — memory, gated on Phase 0.2

Do not start until the phase-resolved measurement names the dominant term. If it is chromatograms,
`ChromatogramStore` applies — but it **cannot** be integrated by editing OpenSwathWorkflow (the
OpenMS core is off-limits); ODIA must own the extract→score loop locally on public APIs. Both
reviewers warn the integration will otherwise either recreate the forbidden core edit or hold both
representations at once and *raise* peak RSS.

## Phase 3 — speed and occupancy, last

Serial phases measured so far: `library_load` 155 s at 1.1 cores (~215M string allocations and
~1.1e9 atomic refcount operations in the parquet→`LightTargetedExperiment` conversion),
`prefilter/build_decoy_view` 62.5 s, `pair_and_rebuild` 53.5 s, `precursor_index`, `score_load`,
`write_scores`, `context_fdr` — all at 1.0 cores. The structural fix is caching the *converted*
library, so the conversion is paid once rather than per run.

Not to be attempted before Phase 0.3 exists: any before/after on a shared node measures contention.

## Explicitly not on this list

- **Widening the prefilter screen.** Measured dead: 10 → 20 ppm cost 327 IDs. 30 ppm would be worse.
  The "30 ppm recovers 1,015 of 1,371 unreachable" line in earlier planning is superseded; it counts
  benefit without the null.
- **`ChromatogramStore` integration ahead of Phase 0.2.** Its headline saving is computed from the
  refuted decomposition.
