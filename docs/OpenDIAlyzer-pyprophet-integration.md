# PyProphet as OpenDIAlyzer's FDR model — analysis + requirements

> **SUPERSEDED framing (user directive 2026-07-26): OpenDIAlyzer must not exec
> external tools, so it does NOT call the `pyprophet` binary.** The semi-supervised
> model is **reimplemented in-process in C++** ("pyprophet-*like*"): start with LDA,
> keep boosted (XGBoost/HGB) in the backlog ([scoring-fdr-backlog](OpenDIAlyzer-scoring-fdr-backlog.md)).
> It is used in TWO places: (a) the recalibration phase, to select RT anchors with a
> real semi-supervised d-score instead of raw `VAR_XCORR_SHAPE`+`VAR_LIBRARY_CORR`
> thresholds; (b) the final FDR step. The §3 "subprocess wrapper" recommendation
> below is retained only as the analysis of what pyprophet does — the *mechanism*
> we replicate — NOT as the integration path. Everything else (the algorithm,
> cross-validation, context FDR, decoy-quality dependency, determinism, thread
> hygiene) carries over as requirements for the in-process implementation.

Scope: the FDR + semi-supervised scoring model, reimplemented **in-process** and
used both for recalibration-anchor confidence and final error control. Grounded in
the installed **pyprophet 3.0.15** as the reference algorithm.

---

## 1. What PyProphet is (analysis)

PyProphet is the semi-supervised FDR engine of the OpenSWATH stack — a Python
reimplementation of mProphet. It does **not** extract or score peaks; it consumes
the ~20 per-peak-group DIA sub-scores OpenSWATH already wrote (`FEATURE_MS2.VAR_*`)
and turns them into one discriminant score + calibrated error rates.

**Algorithm (per level):**
1. Take target and **decoy** peak groups with their `VAR_*` sub-scores.
2. Semi-supervised loop (`--ss_num_iter`): start from one strong "main" score,
   train a classifier (default **LDA**; also SVM / XGBoost / HistGradientBoosting)
   to separate targets from decoys, re-select a confident target set, repeat.
3. **Cross-validated** application (`--xeval_num_iter` folds) so a peak group is
   scored by a model it did not train on — this is what prevents the classifier
   from memorising the peaks it will judge.
4. Collapse sub-scores → one **d-score**; estimate the null from decoys, estimate
   π0, and convert to **p-value → q-value → PEP**. Writes `SCORE_MS2` (SCORE, RANK,
   PVALUE, QVALUE, PEP) back into the `.osw`.

**Levels:** `ms1 | ms2 | ms1ms2 | transition | alignment`. We use `ms2` today;
`ms1ms2` fuses MS1+MS2 evidence and usually helps diaPASEF.

**Context FDR (`infer` subcommand):** peak-group q-values are not the reporting
unit. `pyprophet infer peptide` / `infer protein` roll groups up to peptide/protein
with `--context [run-specific | experiment-wide | global]`. Single run ⇒
`run-specific`. (Top-level `peptide`/`protein`/`ipf`/`gene` are **deprecated** in
3.x → use `infer`; pin the version, the surface is moving.)

**Interfaces:** CLI-first (`pyprophet score`, `pyprophet infer …`, `export`,
`merge`, `reduce`, `subsample`). An importable package exists
(`pyprophet.{scoring,stats,infer,io,_config}`) but the public, tested,
version-stable surface is the CLI. Input/out is `.osw` (SQLite); 3.x also reads
`.parquet` / split-parquet.

---

## 2. How it's used today, and the friction

`OpenSwathWorkflow → features.osw` then, by hand:
`pyprophet score --in features.osw --level ms2 --threads 32` → `count_osw.py`.

Friction this integration removes:
- **Two tools, one result.** A run is "done" only after a separate manual step;
  easy to forget `--level`, seed, or the context rollup.
- **No feedback path.** OpenDIAlyzer is two-pass; today it cannot use FDR-calibrated
  confidence to pick pass-2 anchors (it uses raw `VAR_XCORR_SHAPE`+`VAR_LIBRARY_CORR`).
- **Silent thread oversubscription.** With HistGradientBoosting/XGBoost the learner
  spawns OpenMP threads *inside* each of `--threads` workers; pyprophet's own docs
  warn you must set `OMP_NUM_THREADS ≈ ceil(20/threads)` or it thrashes. A wrapper
  should own this, not the user.
- **No packaged provenance.** Version, classifier, seed, π0, #decoys are not
  recorded with the result.

---

## 3. Integration shape (recommendation)

**Subprocess, not library-import.** OpenDIAlyzer (C++) already shells to
OpenSwathWorkflow via `std::system`; it should likewise invoke a single Python
entry `experiments/odia_fdr.py <osw>` that drives the pyprophet **CLI**. Rationale
(ponytail rung 4 — reuse the supported surface):
- pyprophet's internal Python API is unstable across 3.x (active deprecations);
  the CLI is the contract that is tested and documented.
- Its heavy deps (numpy/sklearn/xgboost) stay isolated in the odia env; the C++
  tool needs no Python linkage.
- Matches the existing OpenSwathWorkflow-wrapper pattern (env fixed by the wrapper).

A `-fdr` flag on OpenDIAlyzer (default on) runs the postprocessor on the final
`.osw`; `-fdr false` reproduces today's raw-output behaviour for A/B.

```
OpenDIAlyzer main_():  … final pass writes out.osw
   if (getFlag_("fdr"))  std::system("python odia_fdr.py out.osw --level ms2 …")
```

`odia_fdr.py` (the Python postprocessing step) is the deliverable this doc scopes.

---

## 4. Functional requirements — `odia_fdr.py`

- **PP1 Input contract.** Accept an OpenSWATH-schema `.osw` containing `FEATURE`,
  `FEATURE_MS2` (all `VAR_*`), `PRECURSOR.DECOY`. Validate up front: tables exist,
  decoys present, sub-scores non-degenerate; fail with an actionable message, not a
  pyprophet stack trace.
- **PP2 Score.** Run `pyprophet score` at a configurable `--level` (default `ms2`;
  expose `ms1ms2`), `--classifier` (default `LDA` — deterministic, no OpenMP;
  `XGBoost`/`HistGradientBoosting` opt-in), fixed `--xeval_num_iter`/`--ss_num_iter`.
  Writes `SCORE_MS2` into the `.osw`.
- **PP3 Context FDR.** Run `pyprophet infer peptide` (and `infer protein`) with
  `--context run-specific` for a single run; expose `experiment-wide`/`global` for
  the future multi-run path (needs `merge` first — out of MVP scope).
- **PP4 Thread hygiene.** Own the `--threads` ↔ `OMP_NUM_THREADS` relationship so
  learner OpenMP threads do not oversubscribe; one knob in, correct env out.
- **PP5 Determinism.** Seed the cross-validation/subsampling; record the seed.
  Same input + config ⇒ same q-values (needed for the two-pass A/B to be honest).
- **PP6 Report.** Emit a compact summary next to the `.osw`: IDs@1% (distinct
  target peptides, via `count_osw.py` logic), the **d-score tail** (targets
  outscoring every decoy), #decoys, π0, classifier, seed, pyprophet version. This
  is the run's provenance record.
- **PP7 Failure modes.** Detect and report, not crash: too few decoys/targets for a
  stable model; a `VAR_*` column all-NaN; a degenerate null (π0→1, no separation).
  Return a distinct exit code so OpenDIAlyzer surfaces "extracted OK, FDR
  inconclusive" rather than a generic failure.
- **PP8 Idempotent + non-destructive.** Operate on a copy (or a documented in-place
  contract) so re-running does not corrupt a half-scored `.osw`; `score` refuses to
  re-score an already-scored file otherwise.

## 5. Non-functional

- **PPN1 Single-run MVP.** Run-specific FDR only; multi-run (merge/backpropagate,
  global context, MBR alignment) is explicitly deferred.
- **PPN2 Performance.** ms2 scoring on a full-proteome `.osw` is minutes; must not
  dominate the extraction passes. LDA default keeps it CPU-light and deterministic.
- **PPN3 Env pinned.** pyprophet **3.0.15** pinned in the odia env; the wrapper
  asserts the version and warns on drift (CLI is changing under us).
- **PPN4 Schema contract.** The `.osw` schema is a hard interface. If OpenDIAlyzer
  ever emits a native/columnar output, it must still produce a pyprophet-readable
  `.osw` (or `.parquet`) view — do not fork the schema silently.

## 6. Adversarial notes / open decisions

1. **Decoy quality feeds the null.** pyprophet's error rates are only as honest as
   the decoys. Our P2 independently-predicted decoys were built precisely so the
   null is well-matched — this integration inherits that dependency; a decoy
   regression silently miscalibrates q-values. Report #decoys + π0 so it is visible.
2. **Postprocessing vs in-loop.** The user asked for a *postprocessing* step —
   correct for the MVP (FDR after the final pass). But the two-pass engine could
   also use pyprophet-calibrated confidence to select pass-2 anchors instead of raw
   sub-score thresholds. That is a **separate** feedback requirement; note it,
   don't build it yet (adds a full FDR pass between extraction passes — cost).
3. **LDA vs boosted.** LDA is deterministic and dependency-light; XGBoost/HGB can
   lift IDs but add nondeterminism + the OpenMP thread trap. Default LDA; make the
   boosted learners an explicit opt-in measured against LDA.
4. **`ms2` vs `ms1ms2`.** diaPASEF has real MS1 signal; `ms1ms2` may recover IDs.
   A/B it once the two-pass baseline is fixed — but only after PP2 is stable.
5. **Version churn.** 3.x deprecated `peptide/protein/ipf/gene` for `infer`. Pinning
   is not optional; the wrapper must not assume a command exists — probe once.

Related: docs/OpenDIAlyzer-requirements.md · docs/OpenDIAlyzer-plan-synthesized.md ·
[[Library dissection and fix plan]] (P2 decoys → the null this FDR trusts).
