# OpenDIAlyzer — synthesized plan (after codex + vibe review)

> **Architecture directive (user, 2026-07-26): OpenDIAlyzer is ONE self-contained
> executable. It MUST NOT exec external tools — no `std::system`, no external
> `OpenSwathWorkflow`, no external `pyprophet`.** Extraction is in-process by
> inheriting `TOPPOpenSwathBase` (libOpenMS) and calling
> `OpenSwathWorkflow::performExtraction` directly (loaders `loadSwathFiles` /
> `loadTransitionList` are inherited; `.d` + library load ONCE, reused across
> passes). The confidence/FDR model is an **in-process LDA** (pyprophet-like),
> used both to select recalibration anchors and for final FDR — start with LDA,
> keep boosted (XGBoost/HGB) in the backlog. Literature scan on stronger scoring /
> FDR / recalibration options runs in the background.


Both reviewers independently reached the **same** verdict, and it overrides the
v1 requirements: **the requirements assume the cause and scope a solution to an
unproven problem.** Do not build the engine yet.

## The synthesis (codex ∩ vibe)

1. **Root-cause GATE first (1 day), before any tool code.** We have not isolated
   *why* OpenSWATH gets 5,594 vs DIA-NN's 43,640 on the same library. Run the
   cheap experiment that isolates it. Decision tree (vibe NR1, codex MVP-3/11):
   - **If OpenSWATH + accurate (DIA-NN/observed) RT ≥ ~35k** → OpenSWATH's
     *single-pass calibration* is the whole problem. The fix is a small
     two-pass / calibration change (or an OpenMS patch), **not a new engine.**
   - **If < ~20k** → calibration is not enough; the two-pass MVP is warranted.
   - **If both the gate and an optimized-OpenSWATH resource run hit target** →
     upstream a prefilter/calibration interface to OpenMS; kill the standalone
     engine.
2. **The MVP is a two-pass OpenSWATH orchestrator + calibration sidecar** — NOT a
   new extraction engine. Pass 1 (wide, OpenSWATH) → pyprophet → fit robust
   RT/mass/IM transforms from confident anchors → Pass 2 (narrow, OpenSWATH).
   Reuse OpenSWATH scoring and peak-picking **unchanged** (parity tests required
   before any new kernel).
3. **Also run the resource decomposition** (codex MVP-9): optimized OpenSWATH
   (PQP, local NVMe, physical cores, `maxConcurrentSwaths`, `innerBatchSize`,
   null-writer arm) — the 1.7 TB / 5.6 h may already be fixable with existing
   controls, killing the "new engine for performance" motivation.

## Cut from v1 (both reviewers)
C++ ONNX training (R3), MS2 fine-tuning, new scoring kernels, F4 local streaming
extraction kernel, interference subtraction/arbitration, in-process FDR, parallel
SQLite writer, columnar store, mzML/Bruker fallbacks, TSV loading, quantification,
NUMA scheduling, multi-pass>2.

**mzPeak streaming — the PREFERRED input mode (user directive 2026-07-26).** No
longer "just a deferred spike": it is the top input-layer backlog item, to be wired
in after the extraction MVP is validated. The reader is already built
(`libmzpeak.a/.so`); the integration is a single `MzPeakSpectrumAccess :
OpenSwath::ISpectrumAccess` adapter feeding `SwathMap`s to `performExtraction`
unchanged, with `SwathFile` (`.d`/mzML) demoted to fallback. Full plan:
docs/OpenDIAlyzer-mzpeak-backlog.md.

## Missing requirements to add (codex)
Canonical precursor/peptidoform identity + discovery unit; exact target/decoy
construction + paired filtering; cross-fitting to prevent calibration leakage;
peak-group discovery + multiple-peak retention; full-vs-half window conventions +
units for every RT/mass/IM param; min anchor counts + coverage + fallback +
failure exit; per-stratum recall curves; output schema versioning; external
validation beyond one S08 run; benchmark hardware/storage/cache/phases stated.

## Immediate action: THE GATE EXPERIMENT
`experiments/rt_gate.py` + benchmark: our P1+P2 library, but RT rewritten to the
**observed** run-RT (DIA-NN report.parquet, for the ~43k precursors DIA-NN found;
predicted-in-run-scale for the rest). Search S08 (1/8 fast). Metric: IDs@1% and
the d-score tail. This isolates whether **RT accuracy** is the dominant lever.
- ≥ ~20k (scaled) → calibration/RT is the story → build the two-pass MVP.
- ~5k → RT is not it → the cause is extraction locality / scoring / interference,
  and the requirements must be re-derived from that.

Related: [[Library dissection and fix plan]] (task-3: 43,640) ·
[[Bottleneck review - what survived]] · docs/OpenDIAlyzer-requirements.md (v1)
