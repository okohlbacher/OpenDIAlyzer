# OpenDIAlyzer — manual

Targeted DIA analysis: takes a DIA run and a spectral library, extracts chromatograms for every
library precursor, scores the peak groups, and controls FDR against decoys.

This page documents the **executable**. For why the tool exists see [`../README.md`](../README.md);
for the architecture see [`DESIGN.md`](DESIGN.md).

## Synopsis

```bash
OpenDIAlyzer -in <run> -tr <library> -out <results> [options]
```

```bash
OpenDIAlyzer \
  -in astral.mzML \
  -tr library.oswpq \
  -out results.oswpq \
  -threads 224 \
  -mz_extraction_window 10 \
  -classifier gbt
```

Run `OpenDIAlyzer -selftest` to verify a build without any data: it exercises the RT transform,
isotonic fit, and the degenerate/small-anchor paths, and exits non-zero on failure.

## Inputs and outputs

| Flag | Formats | Notes |
|---|---|---|
| `-in` | `.mzML`, `.mzXML`, `.sqMass`, Bruker `.d` (with OpenTIMS), `.mzpeak` | The DIA run. mzPeak streams by (frame, window) group, so peak memory is O(one group). |
| `-tr` | `.oswpq`, `.pqp`, `.TraML`, `.tsv` | Assay library. A `.tsv` is parsed once and cached beside it as parquet; later runs skip the parse. |
| `-out` | `.oswpq`, `.osw` | `.oswpq` is a ZIP of parquet tables and scores **in memory** — no SQLite round-trip. `.osw` is the pyprophet-compatible SQLite schema. |

**Decoys must be paired by id.** A decoy's id is `<decoy_tag><target id>` (`DECOY_PEPTIDEK_2` for
`PEPTIDEK_2`); the sequence is shuffled, so sequence matching cannot pair them. The run aborts if
no decoy id carries the tag — that state produces q-values with no null distribution behind them,
which is worse than stopping.

> **Reading a `.pqp`**: the library must retain its original ids. OpenDIAlyzer forces
> `legacy_traml_id=true` for exactly this reason; the default renames every compound to its row
> number and silently destroys decoy pairing.

## How a run proceeds

1. **Library load** — parquet or cached TSV.
2. **Prefilter** — screen the library against MS1/MS2 evidence actually present in the run.
   Both label classes are screened by the identical criterion; a pair is kept when *either*
   member has evidence, so the selection is invariant under swapping the labels.
3. **RT calibration** — linear then non-linear iRT/CiRT, fitted from anchors found in the run.
4. **Pass 1** — extract and score at the wide RT window.
5. **Recalibration** — refit the RT map from the pass-1 confident identifications; optionally
   re-centre each identified precursor's window on its measured apex.
6. **Pass 2** — extract and score at the narrow window.
7. **FDR** — semi-supervised classifier, then q-values; optionally peptide/protein context FDR.

## Options

Defaults are the registered ones. Anything not listed is inherited from TOPP.

| Option | Default | Purpose |
|---|---|---|
| `-rt_extraction_window` | 600.0 | Pass-1 RT window in seconds (bounded; centered by the linear bootstrap calibration). |
| `-prefilter` | "true" | Screen the library against MS1/MS2 evidence in the run before  |
| `-decoy_tag` | "DECOY_" | Prefix that identifies a decoy's target partner in the library ids.  |
| `-prefilter_evidence` | "ms2" | Evidence required by -prefilter. ms2 = at least  |
| `-prefilter_mz_extraction_window` | -1.0 | m/z window used by -prefilter ONLY (ppm; <=0 = same as extraction). The  |
| `-prefilter_min_fragments` | 4 | MS2 evidence threshold: distinct library fragments that must match in one  |
| `-prefilter_top_peaks` | 1000 | Most intense peaks kept per spectrum when screening for -prefilter  |
| `-prefilter_min_spectra` | 1 | Number of spectra in which -prefilter_min_fragments must be met. 1 = best  |
| `-mz_calib_bootstrap_ppm` | 50.0 | Wide window used only to COLLECT calibration anchor errors, before any  |
| `-mz_calib_max_anchors` | 2000 | Precursors sampled for m/z calibration. A scalar offset+width needs  |
| `-mz_calib_min_anchors` | 200 | Minimum anchor ERRORS before a window may be inferred at all. Below this  |
| `-mz_calib_sigma_multiple` | 3.0 | Window = k x robust sigma. k=3 covers ~99% of a Gaussian error  |
| `-mz_calib_floor_ppm` | 2.0 | Never infer a window below this. Guards against a high-S/N anchor set  |
| `-mz_calib_min_peakedness` | 3.0 | Reject the inference unless residual density near zero exceeds density  |
| `-lda_folds` | 3 | Cross-validation folds (by precursor) for the in-process LDA. pyprophet  |
| `-lda_iterations` | 3 | Semi-supervised iterations per fold. pyprophet default is 10. |
| `-lda_train_fdr_initial` | 0.15 | FDR for selecting the FIRST training set, before any discriminant  |
| `-lda_train_fdr` | 0.05 | FDR for training-set selection in later iterations  |
| `-library_cache` | "true" | Cache a parsed TSV library as parquet ('<library>.oswpq') beside it and  |
| `-convert_library` | "" | Convert -tr to an OpenSwath parquet bundle at this path and EXIT without  |
| `-prefilter_out` | "" | Write the surviving precursors (id, sequence, decoy) here and EXIT  |
| `-mz_extraction_window_unit` | "ppm" | Unit for -mz_extraction_window / -mz_extraction_window_ms1. |
| `-extra_rt_extraction_window` | 0.0 | Extend the RT extraction window by this much on BOTH sides (seconds). |
| `-min_upper_edge_dist` | 0.0 | Minimum distance to the upper edge of a SWATH window for a precursor to be  |
| `-extraction_function` | "tophat" | Function used to extract a chromatogram point from the m/z window. |
| `-ms1_isotopes` | 3 | Number of MS1 isotopes to extract per precursor. |
| `-batchSize` | 0 | Precursors per batch (0 = no batching). Bounds peak memory on large libraries. |
| `-innerBatchSize` | -1 | Inner (per-SWATH) batch size; -1 = auto. |
| `-readOptions` | "normal" | How SWATH maps are loaded AND whether the workflow keeps the working map  |
| `-mz_extraction_window_ms1_unit` | "ppm" | Unit for -mz_extraction_window_ms1. SEPARATE from the MS2 unit, as in  |
| `-im_extraction_window_ms1` | -1.0 | Ion-mobility window for MS1 extraction (-1 = none). Used only when MS1 ion  |
| `-irt_mz_extraction_window` | 50.0 | m/z window for the iRT/CiRT calibration extraction. |
| `-irt_mz_extraction_window_unit` | "ppm" | Unit for -irt_mz_extraction_window. |
| `-irt_im_extraction_window` | -1.0 | Ion-mobility window for the iRT/CiRT calibration extraction (-1 = off). |
| `-rt_calibration` | "cirt" | Pass-1 RT handling: 'cirt' = OpenSWATH-equivalent iRT/CiRT calibration from the data (recommended); 'bootstrap' = linear map of normalized library RT  |
| `-use_estimated_rt_window` | "true" | Use the RT extraction window estimated by the calibration (as OpenSwathWorkflow does) instead of -rt_extraction_window for pass 1. |
| `-calibration_min_rsq` | 0.70 | Minimum R^2 for the iRT regression (OpenSWATH default 0.95 assumes spike-in iRT kits; predicted libraries are noisier). |
| `-calibration_min_coverage` | 0.30 | Minimum fraction of iRT anchor peptides that must survive outlier removal. |
| `-calibration_nonlinear_bins` | 400 | RT bins for nonlinear iRT sampling (OpenMS default 2000). Controls RT  |
| `-calibration_nonlinear_per_bin` | 10 | Nonlinear iRT anchors sampled per RT bin (OpenMS default 50). Anchor  |
| `-calibration_nonlinear_rt_window` | 2400.0 | FULL-width RT window (i.e. +/-1200 s) for finding nonlinear iRT anchors,  |
| `-rt_fit_loess_span` | 0.3 | LOESS neighbourhood as a fraction of anchors for the RT recalibration  |
| `-rt_calib_min_yield_pct` | 10.0 | Minimum calibration anchor yield (percent of candidates that produced a  |
| `-rt_extraction_window_min` | 0.0 | Lower bound (full width, seconds) on the pass-1 RT window when  |
| `-max_concurrent_swaths` | -1 | Cap SWATH windows extracted concurrently (-1 = auto). |
| `-rt_extraction_window_recal` | 240.0 | Pass-2 (narrow) RT window in seconds, after recalibration. |
| `-empirical_rt` | "false" | Pass-2+: replace predicted RT of every pass-1-seen precursor (target AND decoy) with its OWN measured apex RT (DIA-NN-style empirical library; symmetr |
| `-mz_extraction_window` | 30.0 | MS2 m/z window (ppm). |
| `-mz_extraction_window_ms1` | 30.0 | MS1 m/z window (ppm). |
| `-ion_mobility_window` | -1.0 | Ion-mobility extraction window. -1 = OFF, which is required for data  |
| `-recal_passes` | 2 | 1 = single wide pass (OpenSWATH-like); 2 = recalibrated two-pass. |
| `-threads` | 1 | OpenMP worker threads for extraction. |
| `-tempDirectory` | "" | Scratch dir for cached SWATH data. |
| `-pasef` | "auto" | diaPASEF IM windowing: auto (on if SWATHs carry IM) | true | false. |
| `-score_osw` | "" | Score an existing raw .osw in-process (LDA FDR) and exit — validation/standalone use. |
| `-fdr_context` | "global" | Also control FDR at the peptide and protein level and write  |
| `-picked_protein` | "true" | Protein-level FDR by picked target-decoy competition (Savitski 2015;  |
| `-entrapment_tag` | "" | Id prefix marking entrapment sequences in the library. When set, report  |
| `-classifier` | "lda" | Semi-supervised learner: lda (linear, default) or gbt (histogram  |

## Operational notes

**Set the m/z window to the instrument, not to a habit.** Chromatogram memory and interference both
scale with it. On the Astral benchmark the instrument's measured MS2 accuracy is 1.66–1.71 ppm;
extracting at 30 ppm yielded 27.3M features, 881 GB and *zero* identifications at any usable FDR,
while 10 ppm yielded 2.07M features and 189 GB.

**Screening width is not free.** `-prefilter_mz_extraction_window` sets the screen independently of
extraction, because the two want opposite things — a narrow screen discards real precursors before
any scorer sees them, a wide one admits interference. But widening it also inflates the search
space: measured on the Astral benchmark, going 10 → 20 ppm raised the searched library from 423,079
to 1,459,897 precursors and **lost 327 identifications**, because the larger null tightens the 1%
threshold faster than reachability improves.

**Memory scales with the RT window.** Chromatograms are ~55% of peak RSS and are linear in the
window: the same search at 600 s vs 1435 s peaked at 189.5 GB vs 366.8 GB.

**Calibration may only narrow.** A mass or RT calibration that returns a *wider* window than
configured is reporting a failed fit, not a worse instrument, and is rejected. An early version
without this guard inferred 72.6 ppm on a 1.66 ppm instrument and cost 2,302 identifications.

## Diagnostics

Every run logs a phase profile with wall **and** CPU time, so a serial phase is self-identifying
(`avg cores = cpu/wall`):

```
OpenDIAlyzer[phase] library_load: 155.7 s wall, 167.6 s cpu, 1.1 avg cores
OpenDIAlyzer[phase] dia_run_load: 138.4 s wall, 6673.5 s cpu, 48.2 avg cores
```

Other lines worth reading:

| Prefix | Tells you |
|---|---|
| `[prefilter]` | how many precursors survived, split by label — a large target/decoy asymmetry means the criterion is not discriminating |
| `[mzcal]` | whether a mass window was inferred, and if not, *why* (too few anchors, or flat residuals) |
| `[scoreload/*]` | rows loaded, rows skipped and the reason, columns before/after pruning, target count |
| `CiRT calibration OK` | anchor pairs **and the yield** — single-digit yield means the anchors were not found |

## Exit codes

`0` success. Non-zero on unreadable input, an unusable library (no decoy pairing, no supported
precursors), or a failed `-selftest`. The tool prefers to abort loudly over emitting q-values it
cannot stand behind.
