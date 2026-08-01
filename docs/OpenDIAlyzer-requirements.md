# OpenDIAlyzer — requirements (v1, for adversarial review)

A new DIA extraction/identification tool, re-imagining `OpenSwathWorkflow`. It
USES OpenMS (links libOpenMS + libmzpeak) but does not live in the OpenMS repo;
it is the search companion to OpenDIALibGen.

**Why it exists — the measured mandate.** The same predicted library scores
**43,640 IDs in DIA-NN's engine but 5,594 in OpenSwathWorkflow** on agxt S08
(1% FDR). Our library is not the limiter (DIA-NN beats its own library with it).
The ~7× gap is entirely the extraction engine. OpenDIAlyzer targets that gap.

---

## 1. Deep adversarial analysis of OpenSwathWorkflow (the baseline)

Grounded in `ext/OpenMS/src/topp/OpenSwathWorkflow.cpp` (1855 lines) and this
session's measurements.

### 1.1 What it does (the actual control flow)
`loadSwathFiles()` → **`CalibrationWorkflow::performCalibration()` (ONCE)** →
one `TransformationDescription trafo_rtnorm` → **`OpenSwathWorkflow::
performExtraction(swath_maps, trafo_rtnorm, …)`** (all precursors, fixed
transform, full window) → `OpenSwathOSWWriter` → `.osw`. pyprophet does FDR
afterward, out of process.

### 1.2 Strengths (keep these)
- **Competitive on measured libraries**: OpenSWATH+empirical = 37,539 (≈85% of
  DIA-NN). The scoring (≈20 DIA sub-scores) and peak-picking are sound.
- **diaPASEF-capable**: IM windowing + IM sub-scores work (its own tests + our runs).
- **Mature, correct chemistry**: fragment m/z exact; the DIA scores are validated.
- **A real wave scheduler exists** (`maxConcurrentSwaths`, `innerBatchSize`,
  in-memory reads) — bounds memory *somewhat*.

### 1.3 Weaknesses (the mandate) — ranked by measured impact
1. **Single global RT/mass/IM calibration, no iteration.** `performCalibration`
   runs once on ~120 CiRT anchor peptides; the transform is frozen for the whole
   search. DIA-NN recalibrates iteratively on thousands of confident IDs and
   *fine-tunes its predictors*. Measured: our post-calibration RT residual is
   **42.9 s** vs DIA-NN's **21.1 s** — and a plain transform re-fit on all our
   confident IDs *does not beat* OpenSWATH's 42.9 s (it is PeptDeep's floor). So
   the gain requires **model-level recalibration**, not just a better transform.
2. **Exhaustive extraction, no candidate prefilter.** Every transition of every
   precursor is extracted over the full window and materialised. 3.6M precursors
   × ~12 transitions ⇒ ~43M chromatograms, **1.7 TB peak RSS**, ~5.6 h wall on
   S08. DIA-NN sends only ~244k of ~4M precursors to its expensive scorer.
3. **Single-threaded OSW writer.** `OpenSwathOSWWriter` drains a queue from 224
   producers on one thread; the 54–70 GB `.osw` write serialises (the 159 M
   voluntary context switches + system-CPU-dominance we measured).
4. **No interference/peak-group discovery beyond the picker.** DIA-NN subtracts
   fragment interference and arbitrates competing precursors at one RT; OpenSWATH
   scores each group in isolation.
5. **Out-of-process, decoupled FDR.** pyprophet is a separate tool; the engine
   cannot use FDR feedback to drive a second pass.
6. **Object-model + I/O overhead.** General-purpose `TargetedExperiment` /
   `MSExperiment`, TSV library parse (18 min), Bruker load (6 min) — all serial.

---

## 2. Functional requirements

**F1 Input.** Read DIA runs as **mzPeak (streaming, bounded memory)**, with
mzML and Bruker `.d` as fallbacks. diaPASEF (ion mobility) mandatory.
**F2 Library.** Read OpenSWATH TSV and PQP (targets + decoys), incl. IM (1/K0)
and predicted RT. Reuse `TransitionTSVFile`/`TransitionPQPFile`.
**F3 Candidate routing.** Map each library precursor to its SWATH isolation
window (m/z) and IM band; skip precursors with no compatible window/frames.
**F4 Extraction.** Extract fragment (and MS1) chromatograms **locally** around
each candidate's expected apex, streamed and discarded — bounded memory, not
full materialisation.
**F5 Scoring.** Compute the OpenSWATH DIA sub-scores (co-elution, library
dot-product/spectral-angle, mass accuracy, isotope, IM, MS1 contrast). Reuse
`MRMFeatureFinderScoring`/`DIAScoring`/`odia_score.h` where possible; do **not**
reinvent validated scores.
**F6 On-the-fly recalibration (the headline capability).** Iterative passes:
  (a) first pass with wide windows on a high-confidence subset;
  (b) fit RT + mass + IM transforms from confident hits (many, not just CiRT);
  (c) **fine-tune the RT/MS2 predictor** on those hits (to beat the 42.9 s floor);
  (d) narrow windows; (e) re-extract/re-score the full library.
  `recal_passes` controls it (1 = single-pass, OpenSWATH-like, for A/B).
**F7 FDR.** Emit an OpenSWATH-schema `.osw` that pyprophet scores unchanged
(v1); design toward an in-process, entrapment-validated null later.
**F8 Output.** `.osw` (features + sub-scores + q-values), compatible with the
existing `count_osw.py` / pyprophet flow.

---

## 3. mzPeak streaming requirements (from the built reader, M0)

- Use `MzPeak::open()` → `Index` → `Index::spectra()` (`fetch(i)`, **lazy peak
  decode**; metadata-only routing already proven: 128 MB RSS on a 9.7 GB file).
- **Ion mobility is at `selected_ion` level, not `scan`** — read
  `precursors()[].selected_ions[].ion_mobility_value` (the reader's `scan`-level
  `ion_mobility()` returns null on real diaPASEF; a one-line reader fix pending).
- Iterate spectra **in RT order, demultiplexed by SWATH isolation window**; hold
  only a bounded ring buffer per active SWATH.
- Caveat measured: the current per-spectrum `fetch()` iterate cost is ~40 ms
  (11 min for 16k spectra). M1 must batch-read the metadata table once, or the
  streaming source becomes the new bottleneck.

---

## 4. On-the-fly recalibration requirements (literature-grounded)

DIA-NN (Demichev 2020; docs): a first pass calibrates RT and mass, **retrains
the RT/spectral predictors on run-specific confident IDs**, then narrows the RT
scan window hard (our run: "RT window set to 1.4 min"). This is why DIA-NN gets
~43k on *either* predicted or empirical libraries. Requirements:
- **R1** Select a first-pass anchor set (confident, abundance/charge-stratified —
  not just CiRT), avoiding bias toward abundant tryptic peptides.
- **R2** Fit RT, MS1/MS2 mass, and IM transforms from the anchor set (robust to
  outliers; report residual quantiles, not sd).
- **R3** Optionally **fine-tune the PeptDeep RT (and MS2) ONNX model** on the
  anchor set — the measurement says a transform alone can't beat 42.9 s; DIA-NN's
  21.1 s comes from model adaptation. (Biggest, riskiest requirement; may be
  phased.)
- **R4** Narrow windows to `2·[p99(|resid|) + ½ peak-width + smoothing margin]`
  and re-extract. Convergence criterion + max `recal_passes`.
- **R5** Every recalibration decision auditable (logged transforms, residuals,
  anchor counts) — auditable error control is the project's differentiator.

---

## 5. Non-functional requirements
- **N1 Memory bounded**: peak RSS < 128 GB on full-proteome S08 (vs 1.7 TB).
- **N2 Wall time**: ≤ one gradient, ≤ 16 physical cores, for a full-proteome run
  (DIA-NN does ~5 min; target within a small factor).
- **N3 Parallelism**: NUMA-aware, physical-core pinning, a **non-single-threaded
  writer** (batch/parallel `.osw` insert or a columnar sink).
- **N4 Determinism** where feasible (seedable decoys/threads) for A/B.
- **N5 Instrumented**: per-phase timers + `/proc` counters built in (perf is
  blocked on the cluster).

## 6. Explicit non-goals (v1)
Not a new FDR method (pyprophet for now); not a library generator (OpenDIALibGen);
not GPU (CPU first — the gap is algorithmic, both tools are CPU-only); not
non-tryptic/immunopeptidomics-specific yet; not multi-run/MBR yet.

## 7. Open questions for the adversarial review
1. Is per-precursor **local** extraction actually safe, given our RT residual is
   42.9 s and true apices for *missed* peptides may be worse? (The streaming
   design must handle wide-residual candidates without unbounded memory.)
2. Is **model fine-tuning (R3)** worth the complexity, or do transform-only
   recalibration + a candidate prefilter capture most of the 7×?
3. What actually causes 5,594 vs 43,640 — RT, extraction locality, scoring, or
   interference handling? (We have not isolated it; the requirements assume all
   four contribute.)
4. Reuse `MRMFeatureFinderScoring` in a streaming loop, or extract the scoring
   kernels standalone? (The OpenMS class assumes materialised chromatograms.)
5. mzPeak streaming vs. pre-convert to a cached columnar store — which meets N1/N2?
6. Where is the writer replaced — parallel SQLite, `.oswpq`, or Parquet?

Related: [[Why DIA-NN is 100x faster than OpenSWATH]] ·
[[Bottleneck review - what survived]] · [[Streaming extractor design v1]] ·
[[Library dissection and fix plan]] (task-3 result).
