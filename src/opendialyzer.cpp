// OpenDIAlyzer — self-contained two-pass targeted DIA extraction with on-the-fly
// RT recalibration. ONE executable: it links libOpenMS and drives the OpenSWATH
// extraction machinery IN-PROCESS. It never execs an external tool.
//
// Architecture (user directives 2026-07-26): (a) no std::system / no external
// OpenSwathWorkflow / no external pyprophet; (b) derive from plain TOPPBase, NOT
// TOPPOpenSwathBase — the DIA-run and library loaders are inlined here (calling the
// libOpenMS classes SwathFile / TargetedDataFileLoader / TransitionTSVFile|PQP
// directly), so OpenDIAlyzer is its own tool that USES OpenMS rather than an
// OpenSWATH-derived one. The DIA run (.d/mzML) and library are loaded ONCE and
// reused across both passes.
//
// Two passes (codex+vibe synthesis): the same predicted library scores 43,640 in
// DIA-NN but 5,594 in single-pass OpenSWATH; the dominant suspect is OpenSWATH's
// single frozen RT calibration.
//   pass 1: extract over the WHOLE RT range (rt_extraction_window = -1, identity
//           transform) — measure each candidate apex empirically, no fragile
//           pre-calibration.
//   recal : from confident pass-1 peaks (high VAR_XCORR_SHAPE + VAR_LIBRARY_CORR;
//           an in-process LDA replaces these raw thresholds next) fit a piecewise
//           library_RT -> observed_RT transform.
//   pass 2: re-extract with that transform and a NARROW RT window.
// recal_passes=1 = single wide pass (OpenSWATH-like) for clean A/B.
//
// This OpenMS build migrated its public API to std::string (no OpenMS::String).

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathWorkflow.h>       // OpenSwathWorkflow, ChromExtractParams
#include <OpenMS/ANALYSIS/OPENSWATH/TransitionListEvidenceFilter.h>  // Library:prefilter
#ifdef WITH_MZPEAK
#include "odia_mzpeak_access.h"                              // odia::loadMzPeakSwathMaps
#endif
#include <OpenMS/ANALYSIS/OPENSWATH/CalibrationWorkflow.h>     // CiRT/iRT RT calibration (as in TOPP)
#include <OpenMS/ANALYSIS/OPENSWATH/MRMRTNormalizer.h>         // iRT outlier-detection defaults
#include <OpenMS/ANALYSIS/OPENSWATH/SwathMapMassCorrection.h>  // m/z + IM calibration defaults
#include <OpenMS/ANALYSIS/TARGETED/MRMMapping.h>               // chromatogram->assay mapping defaults
#include <OpenMS/SYSTEM/File.h>                                // File::find for the iRT kits
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathOSWWriter.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMFeatureFinderScoring.h>
#include <OpenMS/ANALYSIS/OPENSWATH/TransitionTSVFile.h>
#include <OpenMS/ANALYSIS/OPENSWATH/TransitionPQPFile.h>
#include <OpenMS/ANALYSIS/OPENSWATH/TransitionParquetFile.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathOSWParquetWriter.h>
#include <OpenMS/FORMAT/ParquetFile.h>
#include <OpenMS/FORMAT/ZipArchiveFile.h>
#include <arrow/api.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationDescription.h>
#include <OpenMS/CONCEPT/UniqueIdGenerator.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/FORMAT/DATAACCESS/MSDataWritingConsumer.h>    // NoopMSDataWritingConsumer
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/SwathFile.h>
#include <OpenMS/FORMAT/TargetedDataFileLoader.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/METADATA/ExperimentalSettings.h>

#include "odia_lda.h"                                          // in-process semi-supervised LDA
#include "odia_fdr.h"
#include "odia_library.h"                                   // compact library probe
#include <OpenMS/FORMAT/FASTAFile.h>
#include <unordered_set>
#include <arrow/array.h>
#include <OpenMS/FORMAT/ArrowSchemaRegistry.h>        // OSWPrecursorSchema / OSWTransitionSchema
#ifdef __GLIBC__
#include <malloc.h>                                            // malloc_trim: return the freed library to the OS
#endif

#include <sqlite3.h>
#include <sys/resource.h>
#include <malloc.h>
#include <thread>
#include <atomic>
#include <mutex>
#include <map>
#include <unistd.h>
#include <cstdio>
#include <chrono>
#include <memory>
#include <unordered_set>
#include <iomanip>
#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <limits>
#include <map>
#include <memory>
#include <filesystem>
#include <fstream>
#include <set>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace OpenMS;

class TOPPOpenDIAlyzer final : public TOPPBase
{
public:
  TOPPOpenDIAlyzer()
    : TOPPBase("OpenDIAlyzer",
               "Self-contained two-pass targeted DIA extraction with on-the-fly RT "
               "recalibration (in-process OpenSWATH). Companion to OpenDIALibGen.",
               /*official=*/false)
  {
  }

protected:
  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "DIA run (mzPeak .mzpeak [streaming, bounded memory] / Bruker .d / mzML / mzXML / sqMass).", false);
    registerInputFile_("tr", "<file>", "", "Spectral library (OpenSWATH TSV / PQP / OSWPQ).", false);
    registerOutputFile_("out", "<file>", "", "Output features (.osw). Pass-1 .osw lands beside it.", false);

    registerDoubleOption_("rt_extraction_window", "<s>", 600.0, "Pass-1 RT window in seconds (bounded; centered by the linear bootstrap calibration).", false);
    // "none" makes the engine use the library RT verbatim -- exactly what standalone
    // OpenSwathWorkflow does when no -tr_irt is supplied. Required for a like-for-like
    // extraction-fidelity comparison against that tool.
    // --- TOPP OpenSwathWorkflow parity -------------------------------------------
    // These change extraction RESULTS and were previously hardcoded here, so a run could not
    // be made to match stock OpenSwathWorkflow even in principle. Names and semantics follow
    // TOPP_OpenSwathWorkflow. NOT yet a drop-in INI: our defaults differ deliberately
    // (mz_extraction_window 30 vs TOPP 50, ion_mobility_window 0.047 vs -1) and TOPP options
    // outer_loop_threads / use_ms1_ion_mobility / sort_swath_maps / matching_window_only /
    // enable_ipf are still unregistered. Do not claim INI transferability until they are.
    // Library prefilter. WITHOUT this, extraction runs over EVERY library precursor: on a
    // whole-proteome library that is 7,149,966 candidates, and with the (TOPP-default)
    // permissive peak-picking settings nearly all of them emit features from noise --
    // measured 27,344,163 features, 1753 GB RSS and 2.9 h on one plasma run, where
    // OpenSwathWorkflow's own prefilter retained 34,165 precursors. This is a correctness
    // problem first and a performance catastrophe second.
    registerStringOption_("prefilter", "true|false", "true",
                          "Screen the library against MS1/MS2 evidence in the run before "
                          "extraction (OpenMS Library:prefilter). Off reproduces the unfiltered "
                          "behaviour, which is only sensible for a library already matched to "
                          "the acquisition.", false);
    setValidStrings_("prefilter", {"true", "false"});
    // The prefilter scores TARGETS only, so decoys are kept via their target partner, and the
    // partner link is the id prefix (MRMDecoy: decoy id = decoy_tag + target id). Must match the
    // tag the library was built with -- OpenSwathDecoyGenerator's own default is "DECOY_".
    registerStringOption_("decoy_tag", "<string>", "DECOY_",
                          "Prefix that identifies a decoy's target partner in the library ids. "
                          "Only used by -prefilter; must match the tag used to build the library.", false, true);
    // Default ms2, NOT the OpenMS default hybrid: hybrid is (MS1 OR MS2) and its MS1 arm keeps a
    // precursor on a single m/z coincidence anywhere in the run, which at proteome scale supports
    // 87% of the library and makes the prefilter a no-op. See prefilterLibrary_ for the numbers.
    registerStringOption_("prefilter_evidence", "ms1|ms2|hybrid", "ms2",
                          "Evidence required by -prefilter. ms2 = at least "
                          "-prefilter_min_fragments distinct library fragments in ONE spectrum "
                          "(selective). ms1 = precursor m/z seen at all (near-useless alone). "
                          "hybrid = MS1 OR MS2, i.e. as permissive as ms1.", false, true);
    setValidStrings_("prefilter_evidence", {"ms1", "ms2", "hybrid"});
    registerDoubleOption_("prefilter_mz_extraction_window", "<ppm>", -1.0,
                          "m/z window used by -prefilter ONLY (ppm; <=0 = same as extraction). The "
                          "screen and the scoring want opposite widths: narrow loses real "
                          "precursors before scoring, wide only costs the screen some specificity. "
                          "Try 2-3x -mz_extraction_window.", false, true);
    registerIntOption_("prefilter_min_fragments", "<n>", 4,
                       "MS2 evidence threshold: distinct library fragments that must match in one "
                       "spectrum for a precursor to survive -prefilter.", false, true);
    // The MS2 null scales as p^(min_fragments) in the per-fragment coincidence rate p, and p is
    // linear in the number of peaks retained per spectrum -- so this is the highest-leverage knob
    // for selectivity (1000 -> 300 cuts the 4-fragment null by ~(10/3)^4 ~ 120x).
    registerIntOption_("prefilter_top_peaks", "<n>", 1000,
                       "Most intense peaks kept per spectrum when screening for -prefilter "
                       "evidence. Lower = more selective.", false, true);
    // The one criterion that actually discriminates. m/z matching alone, maximised over every
    // spectrum in the run, gives a precursor ~3889 independent chances at a coincidence on this
    // gradient -- measured target:decoy enrichment 1.02x, i.e. none. Requiring the evidence to
    // RECUR costs a real precursor almost nothing (DIA-NN measures FWHM.Scans = 2.53 here) and
    // makes an isolated coincidence fail.
    registerIntOption_("prefilter_min_spectra", "<n>", 1,
                       "Number of spectra in which -prefilter_min_fragments must be met. 1 = best "
                       "single spectrum anywhere in the run (not discriminating). >1 requires the "
                       "evidence to recur across cycles.", false, true);
    // --- data-driven m/z window (see inferMassAccuracyPpm_) --------------------------------
    // Bootstrap width. Deliberately WIDE: the anchors are known-present iRT/CiRT peptides sought
    // with isotope, co-elution and SWATH-window constraints, and the noise-match budget scales with
    // the NUMBER OF QUERIES, so a window that is fatal across 7.15M library targets is safe across
    // a few hundred anchors. Do not raise much beyond this -- past ~50 ppm even constrained
    // matching admits too many coincidences.
    registerFlag_("mem_components", "PROBE: walk the feature map to attribute memory by component. "
                                   "Serial, and it constructs a vector<string> per subordinate plus a "
                                   "lookup per meta value -- ~150M allocations against a fragmented "
                                   "180 GB heap, measured at >=41 s. Off by default: it is a diagnostic, "
                                   "not part of the search.");
    registerStringOption_("compact_probe", "<lib.oswpq>", "", "PROBE: load this library into the "
                          "compact representation (src/odia_library.h) and report what it costs, "
                          "then exit. For comparing against the LightTargetedExperiment path on the "
                          "same file.", false, true);
    registerStringOption_("compact_probe_fasta", "<fasta>", "", "FASTA for -compact_probe. With it, "
                          "a peptide that occurs in its protein costs (protein, offset, length) and "
                          "no characters of its own; without it, sequences are interned instead.",
                          false, true);
    registerFlag_("ms1_scores", "Give the classifier the MS1-level sub-scores (var_ms1_*: isotope "
                  "correlation/overlap, mass deviation, MS1 xcorr shape/coelution) alongside the MS2 "
                  "ones. Off by default -- not because it is known worse, but because the evidence "
                  "that originally excluded them was confounded by the library_rt defect and the "
                  "alternative has not been measured since. In-memory (parquet) scoring only; the "
                  "sqlite path reads FEATURE_MS2 and cannot see them.", true);
    registerDoubleOption_("mz_calib_bootstrap_ppm", "<ppm>", 50.0,
                          "Wide window used only to COLLECT calibration anchor errors, before any "
                          "window is inferred. 0 disables m/z inference.", false, true);
    registerIntOption_("mz_calib_max_anchors", "<n>", 2000,
                       "Precursors sampled for m/z calibration. A scalar offset+width needs "
                       "hundreds, not thousands; this is what keeps calibration inside its time "
                       "budget.", false, true);
    registerIntOption_("mz_calib_min_anchors", "<n>", 200,
                       "Minimum anchor ERRORS before a window may be inferred at all. Below this "
                       "the scale estimate is too noisy to act on.", false, true);
    registerDoubleOption_("mz_calib_sigma_multiple", "<k>", 3.0,
                          "Window = k x robust sigma. k=3 covers ~99% of a Gaussian error "
                          "distribution. NB the reviewed supplement's 70th percentile is ~1 sigma, "
                          "i.e. it clips ~32% of true signal by construction; OpenMS's own "
                          "SwathMapMassCorrection instead uses the 99th percentile x 1.3.", false, true);
    registerDoubleOption_("mz_calib_floor_ppm", "<ppm>", 2.0,
                          "Never infer a window below this. Guards against a high-S/N anchor set "
                          "producing a window tighter than the instrument can actually deliver.", false, true);
    registerDoubleOption_("mz_calib_min_peakedness", "<ratio>", 3.0,
                          "Reject the inference unless residual density near zero exceeds density "
                          "at the window edge by this factor. Uniform (all-noise) residuals give "
                          "~1; a real error distribution gives >>1.", false, true);
    // --- semi-supervised LDA -----------------------------------------------------------------
    // Measured on the SAME .osw (2.07M peak groups, 10 ppm): pyprophet reports 4,302 target
    // precursors at 1% FDR where this LDA reported 2,957 -- a 45% gap on IDENTICAL features, so it
    // is the model, not the data. Line-by-line against pyprophet 3.0.15 the substantive differences
    // are iteration count and, more importantly, the FDR used to pick the FIRST training set:
    // pyprophet uses ss_initial_fdr=0.15 then ss_iteration_fdr=0.05, whereas this used 0.05
    // throughout. These are exposed so the gap can be closed empirically rather than guessed at.
    registerIntOption_("lda_folds", "<n>", 3,
                       "Cross-validation folds (by precursor) for the in-process LDA. pyprophet "
                       "instead resamples 50/50 for 10 iterations.", false, true);
    registerIntOption_("lda_iterations", "<n>", 3,
                       "Semi-supervised iterations per fold. pyprophet default is 10.", false, true);
    registerDoubleOption_("lda_train_fdr_initial", "<q>", 0.15,
                          "FDR for selecting the FIRST training set, before any discriminant "
                          "exists. Too strict here and too few positives are selected to fit one, "
                          "the iteration is skipped, and scoring silently falls back to a single "
                          "feature. pyprophet's ss_initial_fdr is 0.15.", false, true);
    registerDoubleOption_("lda_train_fdr", "<q>", 0.05,
                          "FDR for training-set selection in later iterations "
                          "(pyprophet ss_iteration_fdr = 0.05).", false, true);
    registerStringOption_("library_cache", "true|false", "true",
                          "Cache a parsed TSV library as parquet ('<library>.oswpq') beside it and "
                          "reuse it while it is not older than the TSV. The TSV parse is "
                          "single-threaded and took 7:39 on the whole-proteome library.", false);
    setValidStrings_("library_cache", {"true", "false"});
    registerStringOption_("convert_library", "<oswpq>", "",
                          "Convert -tr to an OpenSwath parquet bundle at this path and EXIT without "
                          "searching. Pass the result as -tr for later runs.", false);
    registerStringOption_("prefilter_out", "<tsv>", "",
                          "Write the surviving precursors (id, sequence, decoy) here and EXIT "
                          "before extraction. For tuning/validating -prefilter without paying for "
                          "a full run.", false, true);
    registerStringOption_("mz_extraction_window_unit", "ppm|Th", "ppm",
                          "Unit for -mz_extraction_window / -mz_extraction_window_ms1.", false);
    setValidStrings_("mz_extraction_window_unit", {"ppm", "Th"});
    registerDoubleOption_("extra_rt_extraction_window", "<s>", 0.0,
                          "Extend the RT extraction window by this much on BOTH sides (seconds).", false, true);
    registerDoubleOption_("min_upper_edge_dist", "<Th>", 0.0,
                          "Minimum distance to the upper edge of a SWATH window for a precursor to be "
                          "assigned to it.", false, true);
    registerStringOption_("extraction_function", "<name>", "tophat",
                          "Function used to extract a chromatogram point from the m/z window.", false, true);
    setValidStrings_("extraction_function", {"tophat", "bartlett"});
    registerIntOption_("ms1_isotopes", "<n>", 3,
                       "Number of MS1 isotopes to extract per precursor.", false, true);
    registerIntOption_("batchSize", "<n>", 0,
                       "Precursors per batch (0 = no batching). Bounds peak memory on large libraries.", false, true);
    registerIntOption_("innerBatchSize", "<n>", -1,
                       "Inner (per-SWATH) batch size; -1 = auto.", false, true);
    // Default 'normal' = STREAMING, chosen to MATCH OpenSwathWorkflow's own default so that ODIA is
    // not worse than the reference it wraps. It used to be 'auto', which resolved to
    // cacheWorkingInMemory for mzML to enable OpenSWATH's SWATH-range (wave) scheduler.
    //
    // Honest statement of the evidence: the comparison that motivated this was ODIA-wave at 224
    // threads (6d07h CPU / 2:03:46 wall / 1538 GB) against stock OSW-streaming at 24 threads
    // (1d05h CPU / 1:25:49 wall / 153 GB) for the same 27.34M features. That is CONFOUNDED -- it
    // cannot attribute the difference to the scheduler, because the thread counts differ 9x and
    // ODIA additionally ran a calibration and an LDA pass. ODIA's own utilisation (74 of 224 cores,
    // 33%) points at contention on the wave path, and the earlier 2.7x wave win was measured on a
    // smaller configuration -- consistent with "wave wins below some thread count and loses above".
    // The ODIA-vs-ODIA A/B at matched thread counts is the experiment that would settle it and has
    // NOT been run. Until it is, 'normal' is justified on reference-parity and bounded-memory
    // grounds, not on a proven scheduler verdict.
    registerStringOption_("readOptions", "auto|normal|cache|cacheWorkingInMemory|workingInMemory", "normal",
                          "How SWATH maps are loaded AND whether the workflow keeps the working map "
                          "in RAM. Default 'normal' STREAMS, bounding memory, and matches "
                          "OpenSwathWorkflow's default. 'auto' restores the old per-format behaviour "
                          "(Bruker .d: normal; mzML/mzXML: cacheWorkingInMemory). In-memory working "
                          "maps are a PRECONDITION for the SWATH-range (wave) scheduler "
                          "(OpenSwathWorkflow.cpp:486, needs batchSize<=0 AND load_into_memory), "
                          "which may still win on smaller inputs or lower thread counts -- try "
                          "cacheWorkingInMemory there and measure.", false, true);
    setValidStrings_("readOptions", {"auto", "normal", "cache", "cacheWorkingInMemory", "workingInMemory"});
    registerStringOption_("mz_extraction_window_ms1_unit", "ppm|Th", "ppm",
                          "Unit for -mz_extraction_window_ms1. SEPARATE from the MS2 unit, as in "
                          "TOPP: inheriting the MS2 unit would silently reinterpret a ppm MS1 window "
                          "as Th.", false, true);
    setValidStrings_("mz_extraction_window_ms1_unit", {"ppm", "Th"});
    registerDoubleOption_("im_extraction_window_ms1", "<1/K0>", -1.0,
                          "Ion-mobility window for MS1 extraction (-1 = none). Used only when MS1 ion "
                          "mobility is in play (diaPASEF); without it MS1-IM scores are computed from "
                          "chromatograms extracted with no IM windowing.", false, true);
    registerDoubleOption_("irt_mz_extraction_window", "<w>", 50.0,
                          "m/z window for the iRT/CiRT calibration extraction.", false, true);
    registerStringOption_("irt_mz_extraction_window_unit", "ppm|Th", "ppm",
                          "Unit for -irt_mz_extraction_window.", false, true);
    setValidStrings_("irt_mz_extraction_window_unit", {"ppm", "Th"});
    registerDoubleOption_("irt_im_extraction_window", "<1/K0>", -1.0,
                          "Ion-mobility window for the iRT/CiRT calibration extraction (-1 = off).", false, true);
    registerStringOption_("rt_calibration", "cirt|bootstrap|none", "cirt", "Pass-1 RT handling: 'cirt' = OpenSWATH-equivalent iRT/CiRT calibration from the data (recommended); 'bootstrap' = linear map of normalized library RT onto the run's RT range; 'none' = library RT verbatim (matches OpenSwathWorkflow without any iRT).", false);
    setValidStrings_("rt_calibration", {"cirt", "bootstrap", "none"});
    registerStringOption_("use_estimated_rt_window", "true|false", "true", "Use the RT extraction window estimated by the calibration (as OpenSwathWorkflow does) instead of -rt_extraction_window for pass 1.", false);
    setValidStrings_("use_estimated_rt_window", {"true", "false"});
    // OpenSWATH's defaults (rsq 0.95 / coverage 0.6) assume spike-in iRT kits, whose iRT is
    // near-perfectly linear in run RT. A PREDICTED library (PeptDeep) is much noisier -- we
    // measure rsq ~0.76 on this data -- so stock thresholds reject a calibration that is
    // still far better than no calibration at all. Relaxed here, and exposed for tuning.
    registerDoubleOption_("calibration_min_rsq", "<r2>", 0.70, "Minimum R^2 for the iRT regression (OpenSWATH default 0.95 assumes spike-in iRT kits; predicted libraries are noisier).", false);
    registerDoubleOption_("calibration_min_coverage", "<frac>", 0.30, "Minimum fraction of iRT anchor peptides that must survive outlier removal.", false);
    // Nonlinear iRT sampling. OpenMS defaults (2000 bins x 50/bin, 600 s window) were measured to
    // be a very bad trade here: 67,429 anchors / 968,800 transitions extracted for 32.7 CPU-HOURS
    // -- more than stock OpenSwathWorkflow's ENTIRE extraction -- yielding 136 usable anchor pairs
    // (0.2%). The yield is low for the same reason everything else here fails: at a 600 s window
    // against a PeptDeep RT residual of sd ~810 s, most anchors' peaks are not inside the window
    // being searched. Sampling MORE anchors cannot fix that; widening the search window can.
    // So: sample far fewer (a smooth curve needs hundreds of good points, not tens of thousands)
    // and look for them over a window wide enough to actually contain them.
    // Cost is ~ bins x per_bin x window_width. Cutting per_bin alone while widening the window
    // 4x is very nearly COST-NEUTRAL ((10/50) x (2400/600) = 0.8x) -- so the bin count has to come
    // down too or the "stop burning 32 CPU-hours" goal is not delivered.
    // (400/2000) x (10/50) x (2400/600) = 0.16x, i.e. ~5 CPU-hours instead of 32.7.
    // Bins are the RT-coverage knob and per_bin the density knob: 400 bins still places candidate
    // anchors every ~9 s of a 3900 s gradient, which is ample for a smooth transform, while
    // per_bin protects against over-sampling regions that were never the constraint.
    registerIntOption_("calibration_nonlinear_bins", "<n>", 400,
                       "RT bins for nonlinear iRT sampling (OpenMS default 2000). Controls RT "
                       "COVERAGE of candidate anchors; too few coarsens sparse regions.", false, true);
    // NB the "nonlinear" phase is a second LINEAR fit here: alignmentMethod is pinned to "linear"
    // for TOPP parity (see makeIrtDetectionParam_). So the risk of cutting per_bin is loss of
    // redundancy and RT-region coverage, NOT under-determining a spline.
    registerIntOption_("calibration_nonlinear_per_bin", "<n>", 10,
                       "Nonlinear iRT anchors sampled per RT bin (OpenMS default 50). Anchor "
                       "QUALITY, not count, limits the fit -- 67,429 candidates yielded 136 pairs. "
                       "0 would silently disable nonlinear calibration upstream, so it is rejected.",
                       false, true);
    setMinInt_("calibration_nonlinear_per_bin", 1);
    registerDoubleOption_("calibration_nonlinear_rt_window", "<s>", 2400.0,
                          "FULL-width RT window (i.e. +/-1200 s) for finding nonlinear iRT anchors, "
                          "applied around the LINEAR-transformed RT (OpenMS default 600). Too narrow "
                          "and anchors are not found; too wide and the feature finder can lock onto "
                          "the wrong peak, which corrupts the transform rather than merely weakening "
                          "it. Size it from the linear phase's own residual spread, not from the raw "
                          "library prediction error.", false, true);
    setMinFloat_("calibration_nonlinear_rt_window", 1.0);   // 0 or negative silently means "whole run"
    setMinInt_("calibration_nonlinear_bins", 1);
    // The calibration's estimated window is fitted to the anchors that SURVIVED both the discovery
    // window and outlier removal, so it is in-sample, censored and survivor-selected all at once.
    // On this data it produced 746.9 s (OSW: 375.3 s) while the true predicted-vs-observed residual
    // p95 on unseen precursors is ~1720 s -- i.e. the estimate understates the tail it is being used
    // to cover, which is how pass 1 comes to extract ~35% of true peaks and score noise for the
    // rest. This floor exists so the estimate can never shrink the window below what the library's
    // own RT error demands. 0 = off (OSW-parity behaviour, estimate used as-is).
    // RT warping between a PREDICTED library and a real gradient is curved, not affine, and the
    // binned-median fit resolves it only as coarsely as its <=60 bins. LOESS fits a local line at
    // each control point over all anchors, then PAVA restores monotonicity (elution order is
    // physics; a non-monotone RT map is wrong however well it fits). Span is the fraction of
    // anchors in each local neighbourhood: smaller follows curvature more closely but is noisier,
    // and with few anchors a large span degenerates towards the global line.
    registerDoubleOption_("rt_fit_loess_span", "<frac>", 0.3,
                          "LOESS neighbourhood as a fraction of anchors for the RT recalibration "
                          "fit. 0 = use the historical binned-median fit instead.", false, true);
    registerDoubleOption_("rt_calib_min_yield_pct", "<pct>", 10.0,
                          "Minimum calibration anchor yield (percent of candidates that produced a "
                          "usable anchor) before the estimated RT window is trusted. Below this the "
                          "configured -rt_extraction_window is kept. The estimate is also never "
                          "allowed to WIDEN the window: a wider answer means the fit failed.", false, true);
    registerDoubleOption_("rt_extraction_window_min", "<s>", 0.0,
                          "Lower bound (full width, seconds) on the pass-1 RT window when "
                          "-use_estimated_rt_window is true. The calibration's estimate is fitted to "
                          "surviving anchors and understates the tail for unseen precursors; set this "
                          "to ~2x the library's predicted-RT residual p95. 0 = no floor.", false, true);
    registerIntOption_("max_concurrent_swaths", "<n>", -1, "Cap SWATH windows extracted concurrently (-1 = auto).", false);
    registerDoubleOption_("rt_extraction_window_recal", "<s>", 240.0, "Pass-2 (narrow) RT window in seconds, after recalibration.", false);
    // FDR-safe by symmetry (every precursor, both classes, uses its OWN pass-1 apex), but
    // default OFF until confirmed with decoy/entrapment diagnostics on real data. This is the
    // key lever to recover the narrow-pass collapse (predicted-RT residual p95 ~342s >> window).
    registerStringOption_("empirical_rt", "true|false", "false", "Pass-2+: replace predicted RT of every pass-1-seen precursor (target AND decoy) with its OWN measured apex RT (DIA-NN-style empirical library; symmetric -> FDR-safe).", false);
    setValidStrings_("empirical_rt", {"true", "false"});
    registerDoubleOption_("mz_extraction_window", "<ppm>", 30.0, "MS2 m/z window (ppm).", false);
    registerDoubleOption_("mz_extraction_window_ms1", "<ppm>", 30.0, "MS1 m/z window (ppm).", false);
    registerDoubleOption_("ion_mobility_window", "<1/K0>", -1.0,
                          "Ion-mobility extraction window. -1 = OFF, which is required for data "
                          "with no ion-mobility dimension; with -pasef auto this is switched on "
                          "automatically when the SWATH maps carry IM.", false);

    registerIntOption_("recal_passes", "<n>", 2, "1 = single wide pass (OpenSWATH-like); 2 = recalibrated two-pass.", false);

    registerIntOption_("threads", "<n>", 1, "OpenMP worker threads for extraction.", false);
    registerStringOption_("tempDirectory", "<dir>", "", "Scratch dir for cached SWATH data.", false);
    registerStringOption_("pasef", "<mode>", "auto", "diaPASEF IM windowing: auto (on if SWATHs carry IM) | true | false.", false);
    registerStringOption_("score_osw", "<file>", "", "Score an existing raw .osw in-process (LDA FDR) and exit — validation/standalone use.", false);
    registerFlag_("fdr_pi0", "Apply the Storey pi0 correction (pyprophet/DIA-NN parity: more IDs, "
                             "but a nominal 1% FDR is ~2% actual). Default off = honest true-1% FDR.");
    // Error control above the precursor. A 1% precursor FDR is not a 1% peptide or protein FDR
    // (a protein's score is the best of many precursor draws), so these levels get their own
    // target-decoy pass rather than inheriting the precursor q-value. See src/odia_fdr.h.
    registerStringOption_("fdr_context", "global|none", "global",
                          "Also control FDR at the peptide and protein level and write "
                          "SCORE_PEPTIDE / SCORE_PROTEIN. 'none' writes only SCORE_MS2.", false);
    registerStringOption_("picked_protein", "true|false", "true",
                          "Protein-level FDR by picked target-decoy competition (Savitski 2015; "
                          "The 2022): each protein competes with its OWN decoy, removing the "
                          "protein-size bias of ranking all targets against all decoys.", false);
    registerStringOption_("entrapment_tag", "<prefix>", "",
                          "Id prefix marking entrapment sequences in the library. When set, report "
                          "the combined entrapment FDP estimate (Wen 2025) next to the nominal "
                          "q-value — the way to find out whether the reported 1% is real.", false);
    // The learner inside the semi-supervised loop. Everything around it -- folds, training-set
    // selection, fold normalisation, q-values -- is identical either way, so this is a like-for-like
    // swap. GBT expresses interactions between sub-scores that a single hyperplane cannot; measured
    // on synthetic interaction data through the full loop: 530 -> 615 IDs at a controlled FDR.
    registerStringOption_("classifier", "lda|gbt", "lda",
                          "Semi-supervised learner: lda (linear, default) or gbt (histogram "
                          "gradient-boosted trees; captures sub-score interactions).", false);
    setValidStrings_("classifier", {"lda", "gbt"});
    registerFlag_("selftest", "Run the recalibration-fit self-check and exit (no data needed).");
  }

  // --- inlined loaders (call libOpenMS directly; not inherited) ----------------
  /// Cache path for a TSV library.
  ///
  /// Cache path for a TSV library.
  ///
  /// PQP (SQLite), not parquet -- measured, not preference. OSWPQ cannot round-trip a library this
  /// size in EITHER direction on this OpenMS version:
  ///   * write: each string column used a single `arrow::StringBuilder` (32-bit offsets, 2 GB cap);
  ///     78.5M transitions with long `traml_id`s overflowed it. Fixed here by chunking the columns.
  ///   * read: `ParquetFile::readTable` calls `CombineChunks()` whenever a column has >1 chunk, and
  ///     `getColumn` then returns `chunk(0)` alone while the row loop runs to `num_rows` -- so a
  ///     chunked column is read past its end and dies in `basic_string::_M_create`. Making that
  ///     safe needs a chunk-aware cursor across all 13 transition columns.
  /// Until the reader is fixed, PQP delivers the same win (no serial re-parse) with a path that is
  /// already exercised at this scale. `-convert_library <file>.oswpq` still writes parquet if asked.
  /// NOTE the cache is parquet, not PQP. An interim version used ".cache.pqp"; reading such a cache
  /// back at proteome scale dies with SIGBUS, and because the cache is picked up silently on the
  /// next run, a stale one from that experiment crashed an unrelated calibration run 5 minutes in.
  /// The extension is part of the fix: a cache written by an older build must not be mistaken for a
  /// current one.
  static std::string libraryCachePath_(const std::string& tr) { return tr + ".oswpq"; }

  /// Is the parquet cache present and NOT older than the library it was built from?
  /// Conservative: any doubt (missing, unreadable timestamp, older) means "no".
  static bool parquetCacheFresh_(const std::string& tr, const std::string& cache)
  {
    std::error_code ec;
    if (!std::filesystem::exists(cache, ec) || ec) { return false; }
    const auto src = std::filesystem::last_write_time(tr, ec);
    if (ec) { return false; }
    const auto dst = std::filesystem::last_write_time(cache, ec);
    if (ec) { return false; }
    return dst >= src;
  }



  /// Continuous RSS probe with phase attribution, plus allocator and component accounting.
  ///
  /// WHY THIS EXISTS. docs/OpenDIAlyzer-memory-review.md section 2 modelled peak RSS as
  /// library + chromatograms + features. Its own section 6 REFUTED that model: removing 62% of the
  /// chromatogram term moved peak RSS by 4.6%, and the document concludes the dominant term is
  /// unidentified and needs "a phase-resolved measurement, not another model". This is it.
  ///
  /// Three things the previous entry/exit sampling could not do:
  ///   1. Catch a peak INSIDE a phase. Extraction's high-water mark is mid-phase; sampling at the
  ///      boundaries reports the quiet ends and misses it entirely.
  ///   2. Separate live data from allocator retention. glibc returns freed blocks to an arena, not
  ///      to the OS, so RSS can stay high with nothing live. mallinfo2's uordblks (in use) against
  ///      fordblks (free but retained) says which.
  ///   3. Attribute anything to a named structure. A delta says WHICH phase grew; a component size
  ///      says WHAT.
  struct MemProbe
  {
    static MemProbe& instance() { static MemProbe m; return m; }

    static double rssGB()
    {
      std::FILE* f = std::fopen("/proc/self/statm", "r");
      if (!f) { return -1.0; }
      long long total = 0, resident = 0;
      const int n = std::fscanf(f, "%lld %lld", &total, &resident);
      std::fclose(f);
      if (n != 2) { return -1.0; }
      return static_cast<double>(resident) * static_cast<double>(::sysconf(_SC_PAGESIZE)) / 1073741824.0;
    }

    void start()
    {
      if (running_.exchange(true)) { return; }
      th_ = std::thread([this] {
        while (running_.load())
        {
          const double r = rssGB();
          if (r >= 0.0)
          {
            std::lock_guard<std::mutex> g(mu_);
            const std::string ph = stack_.empty() ? std::string("(outside any phase)") : stack_.back();
            if (r > peak_) { peak_ = r; peak_phase_ = ph; }
            double& p = phase_peak_[ph];
            if (r > p) { p = r; }
            // 27% of wall (537 s) and 27% of CPU sat outside every phase, and the global RSS
            // peak landed there. Un-instrumented time is not free time -- charge it to the
            // phase it follows so the gap has a name instead of being invisible.
            if (stack_.empty()) { gap_secs_["after " + last_phase_] += 0.2; }
          }
          std::this_thread::sleep_for(std::chrono::milliseconds(200));
        }
      });
    }

    void stop()
    {
      if (!running_.exchange(false)) { return; }
      if (th_.joinable()) { th_.join(); }
    }

    void push(const std::string& n) { std::lock_guard<std::mutex> g(mu_); stack_.push_back(n); }
    void pop()
    {
      std::lock_guard<std::mutex> g(mu_);
      if (!stack_.empty()) { last_phase_ = stack_.back(); stack_.pop_back(); }
    }

    /// Allocator state. `in_use` is live; `retained` is freed-but-not-returned, i.e. RSS the process
    /// holds for nothing. `mmapped` is large blocks, which glibc DOES return on free.
    static void logAllocator(const char* where)
    {
#if defined(__GLIBC__) && defined(__GLIBC_MINOR__) && (__GLIBC__ > 2 || __GLIBC_MINOR__ >= 33)
      const struct mallinfo2 mi = mallinfo2();
      const double g = 1073741824.0;
      OPENMS_LOG_INFO << "OpenDIAlyzer[mem/alloc] " << where << ": in_use "
                      << std::fixed << std::setprecision(2) << mi.uordblks / g << " GB, retained "
                      << mi.fordblks / g << " GB, mmapped " << mi.hblkhd / g << " GB, arena "
                      << mi.arena / g << " GB; RSS " << rssGB() << " GB" << std::endl;
#else
      OPENMS_LOG_INFO << "OpenDIAlyzer[mem/alloc] " << where << ": RSS " << std::fixed
                      << std::setprecision(2) << rssGB() << " GB (mallinfo2 unavailable)" << std::endl;
#endif
    }

    void report() const
    {
      std::lock_guard<std::mutex> g(mu_);
      OPENMS_LOG_INFO << "OpenDIAlyzer[mem] ---- peak RSS by phase (sampled at 5 Hz) ----" << std::endl;
      std::vector<std::pair<std::string, double>> v(phase_peak_.begin(), phase_peak_.end());
      std::sort(v.begin(), v.end(), [](const auto& a, const auto& b) { return a.second > b.second; });
      for (const auto& kv : v)
      {
        OPENMS_LOG_INFO << "OpenDIAlyzer[mem]   " << std::fixed << std::setprecision(2)
                        << std::setw(8) << kv.second << " GB  " << kv.first << std::endl;
      }
      OPENMS_LOG_INFO << "OpenDIAlyzer[mem] GLOBAL PEAK " << std::fixed << std::setprecision(2)
                      << peak_ << " GB, reached during: " << peak_phase_ << std::endl;

      std::vector<std::pair<std::string, double>> gv(gap_secs_.begin(), gap_secs_.end());
      std::sort(gv.begin(), gv.end(), [](const auto& a, const auto& b) { return a.second > b.second; });
      double gtot = 0.0;
      for (const auto& kv : gv) { gtot += kv.second; }
      OPENMS_LOG_INFO << "OpenDIAlyzer[mem] ---- un-phased wall time, total "
                      << std::fixed << std::setprecision(1) << gtot << " s ----" << std::endl;
      for (const auto& kv : gv)
      {
        if (kv.second < 1.0) { continue; }
        OPENMS_LOG_INFO << "OpenDIAlyzer[mem]   " << std::fixed << std::setprecision(1)
                        << std::setw(8) << kv.second << " s  " << kv.first << std::endl;
      }
    }

  private:
    MemProbe() = default;
    ~MemProbe() { stop(); }
    std::thread th_;
    std::atomic<bool> running_{false};
    mutable std::mutex mu_;
    std::vector<std::string> stack_;
    std::map<std::string, double> phase_peak_;
    std::map<std::string, double> gap_secs_;
    std::string last_phase_ = "(startup)";
    double peak_ = 0.0;
    std::string peak_phase_ = "(none)";
  };

  /// Bytes actually held by a LightTargetedExperiment, capacity-based (capacity, not size, is what
  /// the process holds) and including the heap each std::string owns beyond its SSO buffer.
  static std::size_t libraryBytes_(const OpenSwath::LightTargetedExperiment& e)
  {
    auto strb = [](const std::string& s) -> std::size_t {
      return sizeof(std::string) + (s.capacity() > 15 ? s.capacity() + 1 : 0);   // libstdc++ SSO = 15
    };
    std::size_t b = e.compounds.capacity() * sizeof(OpenSwath::LightCompound)
                  + e.transitions.capacity() * sizeof(OpenSwath::LightTransition)
                  + e.proteins.capacity() * sizeof(OpenSwath::LightProtein);
    for (const auto& c : e.compounds)
    {
      b += strb(c.id) + strb(c.sequence) + strb(c.peptide_group_label) + strb(c.gene_name);
      b += c.protein_refs.capacity() * sizeof(std::string);
      for (const auto& r : c.protein_refs) { b += strb(r); }
    }
    for (const auto& t : e.transitions) { b += strb(t.transition_name) + strb(t.peptide_ref); }
    return b;
  }

  /// Wall+CPU stopwatch for the phases OpenMS does not instrument.
  ///
  /// The four OpenMS "Progress of ..." phases accounted for only 49% of a 37-minute run; the other
  /// 51% ran at ~7 of 224 cores and was invisible. Reporting CPU alongside wall makes a serial
  /// phase self-identifying: avg cores = cpu/wall.
  struct PhaseTimer
  {
    const char* name;
    std::chrono::steady_clock::time_point t0;
    double c0;
    double r0;
    static double cpuSeconds()
    {
      struct rusage ru{};
      getrusage(RUSAGE_SELF, &ru);
      return ru.ru_utime.tv_sec + ru.ru_utime.tv_usec * 1e-6
           + ru.ru_stime.tv_sec + ru.ru_stime.tv_usec * 1e-6;
    }

    /// CURRENT resident set in GB. ru_maxrss is a high-water mark and so cannot attribute growth to
    /// a phase; /proc/self/statm's resident field is the live value, which is what a per-phase
    /// delta needs.
    ///
    /// This exists because the memory decomposition in docs/OpenDIAlyzer-memory-review.md was a
    /// MODEL, and its own A/B refuted it: cutting 62% of the chromatogram term moved peak RSS by
    /// 4.6%. The dominant term is unidentified, and the document says finding it needs a
    /// phase-resolved measurement rather than another model. This is that measurement.
    static double rssGB()
    {
      std::FILE* f = std::fopen("/proc/self/statm", "r");
      if (!f) { return -1.0; }                       // not Linux; report nothing rather than a guess
      long long total = 0, resident = 0;
      const int n = std::fscanf(f, "%lld %lld", &total, &resident);
      std::fclose(f);
      if (n != 2) { return -1.0; }
      return static_cast<double>(resident) * static_cast<double>(::sysconf(_SC_PAGESIZE)) / 1073741824.0;
    }
    explicit PhaseTimer(const char* n)
      : name(n), t0(std::chrono::steady_clock::now()), c0(cpuSeconds()), r0(rssGB())
    {
      MemProbe::instance().push(n);
    }
    ~PhaseTimer()
    {
      const double w = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
      const double c = cpuSeconds() - c0;
      const double r1 = rssGB();
      MemProbe::instance().pop();
      OPENMS_LOG_INFO << "OpenDIAlyzer[phase] " << name << ": " << std::fixed << std::setprecision(1)
                      << w << " s wall, " << c << " s cpu, " << (w > 0.01 ? c / w : 0.0)
                      << " avg cores";
      if (r0 >= 0.0 && r1 >= 0.0)
      {
        OPENMS_LOG_INFO << "; RSS " << std::setprecision(2) << r0 << " -> " << r1 << " GB ("
                        << std::showpos << (r1 - r0) << std::noshowpos << ")";
      }
      OPENMS_LOG_INFO << std::endl;
    }
  };

  OpenSwath::LightTargetedExperiment loadLibrary_(const std::string& tr)
  {
    OpenSwath::LightTargetedExperiment exp;
    const FileTypes::Type t = FileHandler::getTypeByFileName(tr);
    if (t == FileTypes::PQP)
    {
      // legacy_traml_id=true is REQUIRED here, not cosmetic. With the default (false) the reader
      // uses PRECURSOR.ID -- a row number -- as the compound id, so every id comes back as "0",
      // "1", "2", ... and the original "DECOY_<target id>" naming is gone. Decoy-to-target pairing
      // in prefilterLibrary_() is by that id (MRMDecoy.cpp:842), so a PQP library read the default
      // way pairs zero decoys and the FDR has no null distribution left. Measured on the 7,149,966-
      // precursor library: 3,546,541 decoys, 0 paired. TRAML_ID keeps the source ids.
      TransitionPQPFile().convertPQPToTargetedExperiment(tr.c_str(), exp, /*legacy_traml_id=*/true);
    }
    else if (t == FileTypes::TSV)
    {
      // Parsing the TSV is the single largest SERIAL cost in the whole tool: measured 7:39 wall,
      // one core, on the 18.5 GB / 78.5M-transition library -- 25% of a 10 ppm run's wall time and
      // a hard Amdahl ceiling no amount of extraction parallelism can get under. The parsed form is
      // identical every time, so parse once and cache it as parquet next to the library.
      const std::string cache = libraryCachePath_(tr);
      const bool use_cache = getStringOption_("library_cache") != "false";
      if (use_cache && parquetCacheFresh_(tr, cache))
      {
        OPENMS_LOG_INFO << "OpenDIAlyzer: loading library from parquet cache " << cache
                        << " (skipping the serial TSV parse)." << std::endl;
        TransitionParquetFile().convertParquetToTargetedExperiment(cache, exp);
        return exp;
      }
      TransitionTSVFile reader;
      reader.setParameters(TransitionTSVFile().getDefaults());
      reader.convertTSVToTargetedExperiment(tr.c_str(), t, exp);
      if (use_cache)
      {
        // Best-effort: a cache that cannot be written must not fail the run. Write to a temporary
        // and rename, so an interrupted write cannot leave a truncated cache that later looks
        // "fresh" by timestamp and gets loaded as if complete.
        const std::string tmp = cache + ".partial";
        try
        {
          OPENMS_LOG_INFO << "OpenDIAlyzer: writing parquet library cache " << cache
                          << " (subsequent runs skip the TSV parse)." << std::endl;
          std::error_code ec;
          std::filesystem::remove_all(tmp, ec);
          TransitionParquetFile().convertLightTargetedExperimentToParquet(tmp, exp);
          std::filesystem::rename(tmp, cache, ec);
          if (ec)
          {
            OPENMS_LOG_WARN << "OpenDIAlyzer: could not finalise the library cache: " << ec.message()
                            << std::endl;
            std::filesystem::remove(tmp, ec);
          }
        }
        catch (const std::exception& e)
        {
          std::error_code ec;
          std::filesystem::remove(tmp, ec);
          OPENMS_LOG_WARN << "OpenDIAlyzer: could not write the library cache (" << e.what()
                          << "); continuing with the parsed TSV." << std::endl;
        }
      }
    }
    else if (t == FileTypes::OSWPQ)
    {
      TransitionParquetFile().convertParquetToTargetedExperiment(tr, exp);
    }
    else
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                       "Library must be OpenSWATH TSV, PQP or OSWPQ.");
    }
    return exp;
  }

#ifdef WITH_MZPEAK
  // Owns the mzPeak metadata index for the whole run: every SwathMap shares it, so it must
  // outlive them.
  odia::MzPeakIndexPtr mzpeak_index_;
#endif

  // Whether this run carries ion mobility. ONE definition, used by calibration, extraction and
  // the IM extraction window, so they cannot disagree -- they did, and requesting IM extraction
  // on non-IM data aborts the run outright.
  bool detectPasef_(const std::vector<OpenSwath::SwathMap>& maps) const
  {
    const std::string pm = getStringOption_("pasef");
    if (pm == "true")  { pasef_ = true;  return true; }
    if (pm == "false") { pasef_ = false; return false; }
    pasef_ = false;
    for (const auto& sm : maps) { if (!sm.ms1 && sm.imLower >= 0 && sm.imUpper >= 0) { pasef_ = true; break; } }
    return pasef_;
  }
  mutable bool pasef_ = false;                 // set by detectPasef_ once maps are loaded
  // Set only by the (currently empty) mass-accuracy hook. <= 0 means "not inferred", in which case
  // makeChromParams_ keeps the configured -mz_extraction_window.
  mutable double inferred_mz_window_ms2_ = -1.0;
  mutable double inferred_mz_window_ms1_ = -1.0;

  bool loadDIARun_(const std::string& in, std::shared_ptr<ExperimentalSettings>& exp_meta,
                   std::vector<OpenSwath::SwathMap>& swath_maps, const std::string& tmp)
  {
#ifdef WITH_MZPEAK
    // mzPeak is the streaming input: the adapter decodes only the (frame, window) groups the
    // extractor asks for, so peak memory is O(one group) rather than O(run) -- measured
    // 154 MB on a 13.7 GB diaPASEF run against 500 GB - 2.0 TB for the SwathFile path.
    // FileHandler has no mzPeak type, so dispatch on the extension.
    if (in.size() > 7 && in.compare(in.size() - 7, 7, ".mzpeak") == 0)
    {
      swath_maps = odia::loadMzPeakSwathMaps(in, mzpeak_index_);
      // exp_meta stays null: mzPeak carries no OpenMS ExperimentalSettings. Downstream only
      // uses it for provenance in the output, and OpenSwathWorkflow tolerates a null here.
      return !swath_maps.empty();
    }
#endif
    const FileTypes::Type t = FileHandler::getTypeByFileName(in);
    // IN-MEMORY reads (the wave scheduler needs them; 'cache' was disk-backed -> the Dl
    // I/O stalls that made pass 1 crawl). Bruker: 'normal' = RegularSwathFileConsumer
    // (in RAM); mzML/mzXML: 'cacheWorkingInMemory'.
    // -readOptions drives this; previously hardcoded, which made the option a no-op AND
    // passed a TOPP-level name into SwathFile, which rejects it.
    std::string readopts; bool dummy_load = false;
    resolveReadOptions_(in, readopts, dummy_load);
    SwathFile sf;
    sf.setLogType(log_type_);
    if (t == FileTypes::MZML)        { swath_maps = TargetedDataFileLoader::loadFile(in, tmp, exp_meta, readopts, nullptr); }
    else if (t == FileTypes::MZXML)  { swath_maps = sf.loadMzXML(in, tmp, exp_meta, readopts); }
    else if (t == FileTypes::SQMASS) { swath_maps = sf.loadSqMass(in, exp_meta); }
#ifdef WITH_OPENTIMS
    else if (t == FileTypes::BRUKER_TDF) { swath_maps = sf.loadBrukerTdf(in, tmp, exp_meta, readopts); }
#endif
    else { return false; }
    return !swath_maps.empty();
  }

  // --- recalibration fit (shared by the engine and the selftest) --------------
  // Pool-Adjacent-Violators: the optimal L2 monotone-nondecreasing fit to y (weighted by w).
  // Replaces the old forward-max clamp, which pinned y to a running max and so created flat
  // plateaus wherever one noisy bin dipped (distorting the map on both sides). PAVA instead
  // *averages* adjacent violators, the isotonic regression DIA-NN/Calib-RT use. O(n).
  static std::vector<double> pava_(const std::vector<double>& y, const std::vector<double>& w)
  {
    assert(y.size() == w.size());                             // parallel; positive finite weights
    const size_t n = y.size();
    std::vector<double> val(n), wt(n);
    std::vector<size_t> len(n);
    size_t m = 0;                                             // active pooled blocks
    for (size_t i = 0; i < n; ++i)
    {
      val[m] = y[i]; wt[m] = w[i]; len[m] = 1;
      while (m > 0 && val[m] < val[m - 1])                    // downward violation -> pool
      {
        const double nw = wt[m - 1] + wt[m];
        val[m - 1] = (val[m - 1] * wt[m - 1] + val[m] * wt[m]) / nw;
        wt[m - 1] = nw; len[m - 1] += len[m];
        --m;
      }
      ++m;
    }
    std::vector<double> out; out.reserve(n);
    for (size_t b = 0; b < m; ++b) { for (size_t k = 0; k < len[b]; ++k) { out.push_back(val[b]); } }
    return out;
  }

  // LOESS: locally weighted linear regression with tricube weights, evaluated at xout.
  //
  // Why this and not the binned median it replaces: binning quantises the anchor set into <=60
  // fixed groups and takes a median per group, so the fit resolution is capped by bin count and a
  // bin straddling a curvature change averages across it. LOESS instead fits a LOCAL line at every
  // evaluation point, weighting neighbours by distance, so it follows gradient curvature (which is
  // exactly what RT warping between a predicted library and a real run looks like) without
  // imposing a global polynomial. It is what DIA-NN and OpenSWATH's own `lowess` alignment use.
  //
  // Output is NOT guaranteed monotone -- local fits can cross -- so callers must still run PAVA.
  // Elution order is physics: a non-monotone RT map is always wrong, however well it fits.
  static std::vector<double> loess_(const std::vector<double>& x, const std::vector<double>& y,
                                    const std::vector<double>& w, const std::vector<double>& xout,
                                    double span, int robust_iters = 2)
  {
    const size_t n = x.size();
    std::vector<double> out(xout.size(), 0.0);
    if (n == 0) { return out; }
    if (n < 3)
    {
      for (size_t k = 0; k < xout.size(); ++k) { out[k] = y[n / 2]; }
      return out;
    }
    // Span rule. A FIXED fraction is wrong across the range of anchor counts this sees (19..2000):
    // at small n, 0.3*n is a handful of points and the local line is mostly variance. Require an
    // absolute floor of ~15 neighbours as well, so the effective span GROWS as n shrinks.
    const double eff_span = std::max(span, 15.0 / static_cast<double>(n));
    size_t q = static_cast<size_t>(std::llround(eff_span * static_cast<double>(n)));
    q = std::max<size_t>(std::min<size_t>(15, n), std::min(q, n));
    q = std::max<size_t>(3, q);

    // Robustness weights (Cleveland's bisquare loop). WITHOUT this it is not LOESS, just one pass
    // of local regression -- and that matters here specifically: the anchors are pass-1 peak groups
    // hunted in windows narrower than the RT error, so gross outliers (residuals of thousands of
    // seconds) are guaranteed in the set. A single tricube pass has breakdown ~1/q, i.e. ~2%,
    // whereas the binned-MEDIAN path it replaces has 50%. Skipping the robust loop would swap a
    // high-breakdown estimator for a low-breakdown one in exactly the dirty-anchor regime.
    std::vector<double> rob(n, 1.0);
    std::vector<double> fit_at_x(n, 0.0);
    for (int it = 0; it <= robust_iters; ++it)
    {
      // Evaluate at xout on the final pass, at the anchor x on earlier passes (to get residuals).
      const std::vector<double>& grid = (it == robust_iters) ? xout : x;
      std::vector<double> res(grid.size(), 0.0);
      for (size_t k = 0; k < grid.size(); ++k)
      {
        const double x0 = grid[k];
        size_t lo = static_cast<size_t>(std::lower_bound(x.begin(), x.end(), x0) - x.begin());
        size_t hi = lo;
        if (lo > 0) { --lo; } else { hi = std::min(n, q); }
        while (hi - lo < q && (lo > 0 || hi < n))
        {
          const bool take_left = (hi >= n) || (lo > 0 && (x0 - x[lo - 1]) <= (x[hi] - x0));
          if (take_left && lo > 0) { --lo; } else if (hi < n) { ++hi; } else { break; }
        }
        const double dmax = std::max(std::abs(x0 - x[lo]), std::abs(x[hi - 1] - x0));
        double sw = 0, sx = 0, sy = 0, sxx = 0, sxy = 0;
        for (size_t i = lo; i < hi; ++i)
        {
          double u = (dmax > 0.0) ? std::abs(x[i] - x0) / dmax : 0.0;
          if (u >= 1.0) { u = 1.0; }
          const double tri = std::pow(1.0 - u * u * u, 3.0);
          const double wi = tri * rob[i] * (w.empty() ? 1.0 : w[i]);
          if (wi <= 0.0) { continue; }
          const double dx = x[i] - x0;                     // centre at x0 -> intercept IS the fit
          sw += wi; sx += wi * dx; sy += wi * y[i];
          sxx += wi * dx * dx; sxy += wi * dx * y[i];
        }
        if (sw <= 0.0) { res[k] = y[std::min(lo, n - 1)]; continue; }
        const double det = sw * sxx - sx * sx;
        // RELATIVE degeneracy test: det has units of weight^2 * x^4, so an absolute epsilon is
        // scale-dependent and would be miscalibrated by ~1e8 if RT were normalised to [0,1].
        if (std::abs(det) <= 1e-12 * std::abs(sw * sxx)) { res[k] = sy / sw; continue; }
        const double b = (sw * sxy - sx * sy) / det;
        res[k] = (sy - b * sx) / sw;
      }
      if (it == robust_iters) { out = res; break; }
      fit_at_x = res;
      // Bisquare reweight from the median absolute residual.
      std::vector<double> ar(n);
      for (size_t i = 0; i < n; ++i) { ar[i] = std::abs(y[i] - fit_at_x[i]); }
      std::vector<double> tmp = ar;
      std::nth_element(tmp.begin(), tmp.begin() + tmp.size() / 2, tmp.end());
      const double mad = tmp[tmp.size() / 2];
      const double s = 6.0 * mad;
      for (size_t i = 0; i < n; ++i)
      {
        if (s <= 0.0) { rob[i] = 1.0; continue; }          // perfect fit -> no downweighting
        const double u = ar[i] / s;
        rob[i] = (u >= 1.0) ? 0.0 : std::pow(1.0 - u * u, 2.0);
      }
    }
    return out;
  }

  // Write the nearest-rank 95th-pct |anchor residual| of a fitted map (diagnostic).
  static void setP95Resid_(const TransformationDescription& td,
                           const std::vector<std::pair<double, double>>& pts, double* p95_resid)
  {
    if (!p95_resid) { return; }
    if (pts.empty()) { *p95_resid = 0.0; return; }
    std::vector<double> r; r.reserve(pts.size());
    for (const auto& pr : pts) { r.push_back(std::abs(pr.second - td.apply(pr.first))); }
    std::sort(r.begin(), r.end());
    const size_t nr = r.size();
    size_t k = (nr * 95 + 99) / 100;                         // ceil(0.95*nr) ...
    k = (k == 0 ? 0 : k - 1);                                // ... nearest-rank index
    *p95_resid = r[std::min(k, nr - 1)];
  }

  static TransformationDescription identityTrafo_()
  {
    TransformationDescription td; td.fitModel("identity"); return td;
  }

  // Fit a robust, monotone library_RT -> observed_RT map. Binned medians (outlier-robust) +
  // PAVA isotonic (monotone) + piecewise-linear interpolation. If p95_resid != nullptr, the
  // 95th-pct |anchor residual| after the fit is written there (diagnostic only; this is the
  // IN-SAMPLE residual of the ANCHORS, not the predictive residual on unseen peptides).
  //
  // Robustness (codex review): OpenMS's interpolated model throws below 3 UNIQUE x-values even
  // for linear interpolation (TransformationModelInterpolated::preprocessDataPoints_). Small or
  // degenerate anchor sets are real (few IDs on a small library), so we guarantee >=3 unique x
  // or fall back to identity -- never let fitModel throw.
  // loess_span <= 0 selects the historical binned-median path; >0 fits LOESS at the bin centres
  // and then isotonises. Both end in PAVA + an interpolated model, so the only thing that changes
  // is how the control-point y values are estimated.
  static TransformationDescription fitTrafo_(std::vector<std::pair<double, double>> pts,
                                             double* p95_resid = nullptr,
                                             double loess_span = 0.0)
  {
    if (pts.empty()) { if (p95_resid) { *p95_resid = 0.0; } return identityTrafo_(); }
    std::sort(pts.begin(), pts.end());
    const size_t n = pts.size();
    const int NB = std::max(1, std::min<int>(60, static_cast<int>(n / 10)));  // >=~10 anchors/bin (C10)
    std::vector<double> bx, by, bw;                           // per-bin median x, median y, count
    for (int b = 0; b < NB; ++b)
    {
      const size_t lo = n * b / NB, hi = n * (b + 1) / NB;
      if (hi - lo < 3) { continue; }                          // need >=3 for a robust median
      std::vector<double> xs, ys;
      xs.reserve(hi - lo); ys.reserve(hi - lo);
      for (size_t i = lo; i < hi; ++i) { xs.push_back(pts[i].first); ys.push_back(pts[i].second); }
      std::nth_element(xs.begin(), xs.begin() + xs.size() / 2, xs.end());
      std::nth_element(ys.begin(), ys.begin() + ys.size() / 2, ys.end());
      bx.push_back(xs[xs.size() / 2]);
      by.push_back(ys[ys.size() / 2]);
      bw.push_back(static_cast<double>(hi - lo));
    }
    // Coalesce tied bin-x (weighted) BEFORE isotonic so no bin is silently dropped, and the
    // control-point x-values are strictly increasing for the interpolated model.
    std::vector<double> ux, uy, uw;
    for (size_t i = 0; i < bx.size(); ++i)
    {
      if (!ux.empty() && bx[i] <= ux.back() + 1e-6)
      {
        const double nw = uw.back() + bw[i];
        uy.back() = (uy.back() * uw.back() + by[i] * bw[i]) / nw;   // weighted-mean y for tied x
        uw.back() = nw;
      }
      else { ux.push_back(bx[i]); uy.push_back(by[i]); uw.push_back(bw[i]); }
    }
    // LOESS path. Two things had to change from the first attempt:
    //
    // 1. GATING. Binning gives NB = min(60, n/10) bins, so at n=19 -- the low end actually observed
    //    on this data -- NB is 1, ux has ONE element, and the old `ux.size() >= 2` guard skipped
    //    LOESS silently. The feature was inert at exactly the anchor count that motivated it. Below
    //    ~50 anchors a local linear fit is mostly variance anyway, so refuse it EXPLICITLY and say
    //    so, rather than appearing to be on while doing nothing.
    // 2. GRID. Evaluating only at the <=60 bin centres keeps the old shape resolution: LOESS then
    //    improves the y at each knot but never the knot DENSITY, and curvature between knots stays
    //    linearly interpolated. Use a uniform grid over the anchor x-range instead (never beyond
    //    it -- extrapolation is where local fits are worst). Cost is trivial.
    const size_t n_anchor = pts.size();
    if (loess_span > 0.0 && n_anchor >= 50)
    {
      std::vector<double> rx, ry, rw;
      rx.reserve(n_anchor); ry.reserve(n_anchor);
      for (const auto& pr : pts) { rx.push_back(pr.first); ry.push_back(pr.second); }
      const double x_lo = rx.front(), x_hi = rx.back();
      if (x_hi > x_lo)
      {
        const size_t NG = std::min<size_t>(128, std::max<size_t>(8, n_anchor / 2));
        std::vector<double> gx(NG);
        for (size_t i = 0; i < NG; ++i)
        {
          gx[i] = x_lo + (x_hi - x_lo) * static_cast<double>(i) / static_cast<double>(NG - 1);
        }
        const std::vector<double> gy = loess_(rx, ry, rw, gx, loess_span);
        // Replace the binned control points wholesale; weights become uniform because they are no
        // longer bin counts (carrying the stale counts into PAVA would weight a LOESS value by how
        // many anchors happened to fall in an unrelated bin).
        ux = gx; uy = gy; uw.assign(NG, 1.0);
      }
    }
    else if (loess_span > 0.0)
    {
      OPENMS_LOG_INFO << "OpenDIAlyzer: only " << n_anchor << " RT anchors -- using the binned-median "
                      << "fit rather than LOESS (local regression needs >=50 to beat a median)."
                      << std::endl;
    }
    const std::vector<double> iso = (uy.size() >= 2) ? pava_(uy, uw) : uy;   // monotone fit
    TransformationDescription::DataPoints cps;
    for (size_t i = 0; i < ux.size(); ++i) { cps.emplace_back(ux[i], iso[i]); }

    // If binning collapsed to <2 unique x, rebuild from RAW anchors grouped by x -- per-group
    // MEDIAN y (outlier-robust) + PAVA, consistent with the main path (NOT first-y + a forward
    // clamp, which would reintroduce the plateau artefact PAVA exists to avoid). pts is sorted
    // by (x,y), so equal-x anchors are contiguous.
    if (cps.size() < 2)
    {
      std::vector<double> gx, gy, gw;
      size_t i = 0;
      while (i < pts.size())
      {
        size_t j = i;
        while (j < pts.size() && pts[j].first <= pts[i].first + 1e-6) { ++j; }
        std::vector<double> ys; ys.reserve(j - i);
        for (size_t k = i; k < j; ++k) { ys.push_back(pts[k].second); }
        std::nth_element(ys.begin(), ys.begin() + ys.size() / 2, ys.end());
        gx.push_back(pts[i].first);
        gy.push_back(ys[ys.size() / 2]);
        gw.push_back(static_cast<double>(j - i));
        i = j;
      }
      const std::vector<double> giso = (gy.size() >= 2) ? pava_(gy, gw) : gy;
      cps.clear();
      for (size_t k = 0; k < gx.size(); ++k) { cps.emplace_back(gx[k], giso[k]); }
    }
    if (cps.size() < 2)                                       // one unique x -> no map possible
    {
      TransformationDescription td = identityTrafo_();
      setP95Resid_(td, pts, p95_resid);
      return td;
    }
    if (cps.size() == 2)                                      // interpolated needs 3: add collinear midpoint
    {
      const double mx = cps[0].first + 0.5 * (cps[1].first - cps[0].first);   // overflow-safe
      const double my = cps[0].second + 0.5 * (cps[1].second - cps[0].second);
      if (!(mx > cps[0].first && mx < cps[1].first))          // FP degenerate -> can't interpolate
      {
        TransformationDescription td = identityTrafo_();
        setP95Resid_(td, pts, p95_resid);
        return td;
      }
      cps.insert(cps.begin() + 1, std::make_pair(mx, my));
    }
    TransformationDescription td;
    td.setDataPoints(cps);
    Param p;
    p.setValue("interpolation_type", "linear");
    p.setValue("extrapolation_type", "four-point-linear");
    td.fitModel("interpolated", p);
    setP95Resid_(td, pts, p95_resid);
    return td;
  }

  // Rewrite each library precursor's RT through `t`. The ChromatogramExtractor centers
  // the RT window on compound.rt DIRECTLY (the trafo argument to performExtraction is not
  // used for the coordinate path), so we pre-scale the library RT into run seconds here.
  static void rescaleLibraryRT_(OpenSwath::LightTargetedExperiment& exp, const TransformationDescription& t)
  {
    for (auto& c : exp.compounds) { c.rt = t.apply(c.rt); }
  }

  // Empirical library (DIA-NN-style): replace the PREDICTED RT of precursors identified in
  // pass 1 with their MEASURED apex RT. Call AFTER rescaleLibraryRT_ so unseen precursors
  // keep the global-map RT while seen ones snap to the measured pass-1 apex (residual -> measurement
  // floor -> a same-width window is now centred exactly, so the peak sits mid-window rather
  // than near a truncating edge). Returns how many compounds were overwritten.
  static size_t applyEmpiricalRT_(OpenSwath::LightTargetedExperiment& exp,
                                  const std::map<std::string, double>& empirical_rt)
  {
    size_t n = 0;
    for (auto& c : exp.compounds)
    {
      auto it = empirical_rt.find(c.id);
      if (it != empirical_rt.end() && std::isfinite(it->second)) { c.rt = it->second; ++n; }
    }
    return n;
  }

  // Peak-group sub-scores read from an .osw for the in-process LDA.
  struct OswRows
  {
    std::vector<std::vector<double>> feats;   // rows x VAR_ columns (NaN->0, constant cols dropped)
    std::vector<int> labels;                  // 1=target, 0=decoy
    std::vector<long long> group;             // precursor id
    std::vector<long long> feature_id;
    std::vector<double> library_rt, exp_rt;
    std::vector<std::string> traml_id;        // PRECURSOR.TRAML_ID (== library compound.id)

    /// Put the rows in an order that depends on the DATA, not on the order extraction produced them.
    ///
    /// Two separate order-dependencies fed off row position. Fold assignment used a group's
    /// first-occurrence index (fixed in odia_lda.h). GBT's histogram reduction partitions rows into
    /// chunks by index -- lo = n_rows*c/nchunk -- so permuting rows regroups the partial sums, and
    /// the resulting floating-point difference flips split points and diverges the trees. That one
    /// is invisible without OpenMP and identical at 1/8/64 threads, so it reads as deterministic:
    /// measured max |delta| 1.965 on permuted rows, against 6.2e-14 once the order is canonical.
    ///
    /// Sort on physical quantities first (precursor, then peak apex RT) so the key means something
    /// even if feature ids are themselves assigned in extraction order; feature_id only breaks ties.
    void canonicalize()
    {
      const std::size_t n = feats.size();
      if (n < 2) { return; }
      std::vector<std::size_t> idx(n);
      for (std::size_t i = 0; i < n; ++i) { idx[i] = i; }
      std::sort(idx.begin(), idx.end(), [&](std::size_t a, std::size_t b) {
        if (group[a] != group[b]) { return group[a] < group[b]; }
        if (exp_rt[a] != exp_rt[b]) { return exp_rt[a] < exp_rt[b]; }
        return feature_id[a] < feature_id[b];
      });
      // Both loaders fill all seven columns once per row, so a size mismatch is a bug, not a case
      // to tolerate. Skipping the odd one out would leave it misaligned against the rest -- scores
      // silently attached to the wrong precursor, which no test would catch and no output would
      // flag. Refuse to reorder anything rather than reorder some of it.
      for (const std::size_t sz : {labels.size(), group.size(), feature_id.size(),
                                   library_rt.size(), exp_rt.size(), traml_id.size()})
      {
        if (sz != n)
        {
          OPENMS_LOG_ERROR << "OpenDIAlyzer: OswRows column length " << sz << " != " << n
                           << "; refusing to canonicalise (rows left in extraction order, so this "
                              "run is not reproducible)" << std::endl;
          return;
        }
      }
      const auto apply = [&](auto& v) {
        std::decay_t<decltype(v)> out;
        out.reserve(n);
        for (std::size_t i = 0; i < n; ++i) { out.push_back(std::move(v[idx[i]])); }
        v = std::move(out);
      };
      apply(feats); apply(labels); apply(group); apply(feature_id);
      apply(library_rt); apply(exp_rt); apply(traml_id);
    }
  };

  // Read all FEATURE_MS2 VAR_* sub-scores + label/RT for every candidate peak group.
  // for_anchors=true excludes RT-derived sub-scores so RT-recalibration anchors are NOT
  // selected by the very RT agreement they are meant to correct (C8, anti-circularity).
  /// Anchor yield (percent) behind the calibration's RT window estimate. Set where the estimate is
  /// produced, read where it is applied -- different functions, so it lives here rather than as a
  /// local that only looked like it was in scope.
  double calib_yield_pct_ = -1.0;

  /// Set when -out names a .oswpq: features stay in memory, nothing is written to sqlite, and the
  /// run's output is a parquet bundle. Decided once from the output extension.
  bool parquet_out_ = false;
  /// Features retained by the parquet path, handed straight to scoring instead of a disk round trip.
  FeatureMap pass_features_;
  UInt64 run_id_ = 0;
  /// What the sqlite path obtained by joining PRECURSOR / PEPTIDE / PROTEIN. Without a database the
  /// same facts have to come from the library -- which is where they originated, so this is the
  /// shorter path, not a workaround.
  struct PrecursorMeta
  {
    long long id = 0;
    int decoy = 0;
    std::string sequence;                    ///< modified sequence: the peptide-level rollup key
    std::vector<std::string> proteins;       ///< accessions: the protein-level rollup keys
    double library_rt = std::numeric_limits<double>::quiet_NaN();  ///< the LIBRARY's predicted RT
  };
  std::map<std::string, PrecursorMeta> precursor_index_;

  /// Build the group-id -> (precursor id, decoy) map from the searched library. Decoy status lives
  /// on the TRANSITION (LightCompound has no decoy field), so it is derived per peptide ref exactly
  /// as prefilterLibrary_ does -- one definition of "decoy" for the whole tool.
  void buildPrecursorIndex_(const OpenSwath::LightTargetedExperiment& exp)
  {
    precursor_index_.clear();
    std::unordered_set<std::string> decoy_refs;
    for (const auto& t : exp.getTransitions())
    {
      if (t.getDecoy()) { decoy_refs.insert(t.getPeptideRef()); }
    }
    long long next = 0;
    for (const auto& c : exp.getCompounds())
    {
      PrecursorMeta m;
      m.id = next++;
      m.decoy = decoy_refs.count(c.id) ? 1 : 0;
      m.sequence = c.sequence;
      m.library_rt = c.rt;                   // PRECURSOR.LIBRARY_RT in the sqlite schema
      m.proteins.assign(c.protein_refs.begin(), c.protein_refs.end());
      precursor_index_[c.id] = std::move(m);
    }
  }

  OswRows loadOswScores_(const std::string& osw, bool for_anchors = false)
  {
    OswRows R;
    sqlite3* db = nullptr;
    if (sqlite3_open_v2(osw.c_str(), &db, SQLITE_OPEN_READONLY, nullptr) != SQLITE_OK) { sqlite3_close(db); return R; }

    std::vector<std::string> vcols;              // discover VAR_* columns
    sqlite3_stmt* ps = nullptr;
    if (sqlite3_prepare_v2(db, "PRAGMA table_info(FEATURE_MS2)", -1, &ps, nullptr) == SQLITE_OK)
    {
      while (sqlite3_step(ps) == SQLITE_ROW)
      {
        const unsigned char* nm = sqlite3_column_text(ps, 1);
        if (!nm) { continue; }
        std::string s(reinterpret_cast<const char*>(nm));
        if (s.rfind("VAR_", 0) != 0) { continue; }
        if (for_anchors && s == "VAR_NORM_RT_SCORE") { continue; }   // circular for RT anchors
        vcols.push_back(s);
      }
      sqlite3_finalize(ps);
    }
    OPENMS_LOG_INFO << "OpenDIAlyzer: sqlite scoring found " << vcols.size()
                    << " VAR_ sub-scores in FEATURE_MS2." << std::endl;
    if (vcols.empty()) { sqlite3_close(db); return R; }

    // Schema-robust join (C11): our own OSWs store FEATURE.PRECURSOR_ID as the library
    // string id (== PRECURSOR.TRAML_ID); a canonical OpenSWATH OSW remaps it to the integer
    // PRECURSOR.ID. Detect which and join accordingly; group on the integer p.ID either way.
    std::string join_col = "p.TRAML_ID";
    {
      sqlite3_stmt* ts = nullptr;
      if (sqlite3_prepare_v2(db, "SELECT typeof(PRECURSOR_ID) FROM FEATURE LIMIT 1", -1, &ts, nullptr) == SQLITE_OK
          && sqlite3_step(ts) == SQLITE_ROW)
      {
        const unsigned char* t = sqlite3_column_text(ts, 0);
        if (t && std::string(reinterpret_cast<const char*>(t)) == "integer") { join_col = "p.ID"; }
      }
      sqlite3_finalize(ts);
    }
    std::string q = "SELECT f.ID, p.ID, f.EXP_RT, p.DECOY, p.LIBRARY_RT, p.TRAML_ID";
    const int vc0 = 6;                                        // first VAR_* column index in the row
    for (const auto& c : vcols) { q += ", m." + c; }
    q += " FROM FEATURE f JOIN PRECURSOR p ON f.PRECURSOR_ID=" + join_col + " JOIN FEATURE_MS2 m ON m.FEATURE_ID=f.ID";
    sqlite3_stmt* st = nullptr;
    if (sqlite3_prepare_v2(db, q.c_str(), -1, &st, nullptr) != SQLITE_OK) { sqlite3_close(db); return R; }
    const int nv = static_cast<int>(vcols.size());
    while (sqlite3_step(st) == SQLITE_ROW)
    {
      if (sqlite3_column_type(st, 3) == SQLITE_NULL) { continue; }    // DECOY NULL -> skip (C9)
      const long long dec = sqlite3_column_int64(st, 3);
      if (dec != 0 && dec != 1) { continue; }                        // require DECOY in {0,1}
      std::vector<double> x(nv);
      for (int j = 0; j < nv; ++j)
      {
        if (sqlite3_column_type(st, vc0 + j) == SQLITE_NULL) { x[j] = std::numeric_limits<double>::quiet_NaN(); continue; }
        const double v = sqlite3_column_double(st, vc0 + j);
        x[j] = std::isinf(v) ? std::numeric_limits<double>::quiet_NaN() : v;   // NULL/inf -> missing
      }
      R.feature_id.push_back(sqlite3_column_int64(st, 0));
      R.group.push_back(sqlite3_column_int64(st, 1));
      R.exp_rt.push_back(sqlite3_column_type(st, 2) == SQLITE_NULL     // NULL EXP_RT -> NaN, not 0.0
                         ? std::numeric_limits<double>::quiet_NaN()
                         : sqlite3_column_double(st, 2));
      R.labels.push_back(dec == 0 ? 1 : 0);
      R.library_rt.push_back(sqlite3_column_type(st, 4) == SQLITE_NULL
                             ? std::numeric_limits<double>::quiet_NaN()
                             : sqlite3_column_double(st, 4));
      const unsigned char* tid = sqlite3_column_text(st, 5);
      R.traml_id.push_back(tid ? std::string(reinterpret_cast<const char*>(tid)) : std::string());
      R.feats.push_back(std::move(x));
    }
    sqlite3_finalize(st);
    sqlite3_close(db);
    // Drop columns with NO usable signal: all-missing (e.g. VAR_IM_XCORR_SHAPE when IM
    // xcorr is absent -> all NULL) or constant among finite values. This is an
    // UNSUPERVISED structural filter (uses no labels -> not CV leakage) and prevents an
    // all-NaN column from propagating NaN through the discriminant (-> constant d-scores).
    // Per-fold z-scoring/imputation still happens in the LDA.
    const std::size_t before = R.feats.empty() ? 0 : R.feats[0].size();
    dropUninformativeColumns_(R);
    OPENMS_LOG_INFO << "OpenDIAlyzer[scoreload/sqlite] rows " << R.feats.size()
                    << "; columns " << before << " -> "
                    << (R.feats.empty() ? 0 : R.feats[0].size())
                    << " after dropping uninformative; targets "
                    << std::count(R.labels.begin(), R.labels.end(), 1) << std::endl;
    R.canonicalize();
    return R;
  }

  // Drop columns with NO usable signal: all-missing, or constant among the finite values. This is
  // an UNSUPERVISED structural filter (uses no labels, so it is not CV leakage) and it prevents an
  // all-NaN column from propagating NaN through the discriminant and flattening every d-score.
  // Shared by both score sources so the parquet and sqlite paths cannot drift apart on it.
  static void dropUninformativeColumns_(OswRows& R)
  {
    if (R.feats.empty()) { return; }
    const int nc = static_cast<int>(R.feats[0].size());
    std::vector<char> keep(nc, 0);
    for (int j = 0; j < nc; ++j)
    {
      double first = 0.0; bool have_first = false, varies = false, any_finite = false;
      for (const auto& row : R.feats)
      {
        const double v = row[j];
        if (!std::isfinite(v)) { continue; }
        any_finite = true;
        if (!have_first) { first = v; have_first = true; }
        else if (std::abs(v - first) > 1e-12) { varies = true; break; }
      }
      keep[j] = (any_finite && varies) ? 1 : 0;
    }
    for (auto& row : R.feats)
    {
      std::vector<double> nr; nr.reserve(nc);
      for (int j = 0; j < nc; ++j) { if (keep[j]) { nr.push_back(row[j]); } }
      row = std::move(nr);
    }
  }

  // Score rows straight out of the in-memory FeatureMap, so a run never has to write its features
  // to disk and read them back just to score them. Same OswRows contract as loadOswScores_ --
  // deliberately, so finalScore_/recalibration do not care which source produced them.
  //
  // The sub-scores that loadOswScores_ finds as VAR_* COLUMNS are meta values on the feature here,
  // so discovery is over meta keys instead of a table schema. Keys are collected from the first
  // feature and then applied to all of them: OpenSWATH writes the same score set for every feature
  // of a run, and taking the union instead would let one odd feature add a column that is missing
  // (NaN) everywhere else -- which dropUninformativeColumns_ would then delete anyway.
  OswRows loadScoresFromFeatureMap_(const FeatureMap& fmap,
                                    const std::map<std::string, PrecursorMeta>& prec_of_id,
                                    bool for_anchors = false)
  {
    OswRows R;
    if (fmap.empty()) { return R; }

    // UNION over a sample, not fmap[0] alone. The first version of this took the key set from the
    // first feature and applied it to all of them, justified by "OpenSWATH writes the same score set
    // for every feature" -- an assumption asserted, not checked. If feature 0 happens to be missing
    // any score, that column vanishes for the whole run and the classifier silently trains on a
    // smaller feature space than the sqlite path, which reads the FEATURE_MS2 table schema and
    // therefore always sees all 29 VAR_ columns.
    //
    // A sample rather than all 2M features: the union saturates almost immediately, and scanning
    // every feature's key list to build a set costs more than it can possibly add.
    std::set<std::string> keyset;
    const std::size_t probe = std::min<std::size_t>(fmap.size(), 4096);
    for (std::size_t i = 0; i < probe; ++i)
    {
      std::vector<std::string> keys;
      fmap[i].getKeys(keys);
      for (const auto& k : keys) { keyset.insert(std::string(k)); }
    }
    // MS1 SUB-SCOPE. loadOswScores_ reads FEATURE_MS2, so it sees MS2-level sub-scores only;
    // OpenSWATH keeps the MS1-level ones (var_ms1_*) in a separate FEATURE_MS1 table. This path has
    // both on the Feature and can use either scope.
    //
    // THE EVIDENCE THAT ORIGINALLY JUSTIFIED DROPPING THEM WAS CONFOUNDED. It was: "identical
    // settings, only the -out extension differing: sqlite (29 sub-scores) gave 6,607 identifications,
    // this path (35-36 sub-scores) gave 4,913 -- 26% fewer for having MORE features", concluding the
    // MS1 columns were sparse noise the fit wasted capacity on.
    //
    // That 4,913 is now explained by a different defect entirely: this path read `library_rt` from
    // the feature's norm_RT (the OBSERVED rt in iRT space) instead of the library's PREDICTION, so
    // recalibrate_ fitted a function of exp_rt against exp_rt and pass 2 extracted on a corrupted RT
    // axis. With that fixed the same path reaches 6,522 -- WITHOUT any MS1 scores. So the MS1
    // columns were never shown to hurt; they were blamed for a deficit another bug caused.
    //
    // MEASURED 2026-08-01, controlled A/B on the Astral benchmark, both arms on one node:
    //     MS2 scope   (24 features)  6,506 IDs   5,699 peptides   583 proteins
    //     MS1+MS2     (36 features)  6,574 IDs   5,709 peptides   628 proteins
    //   +68 precursors against a measured noise floor of +/-83
    // i.e. NO MEASURABLE EFFECT. The MS1 sub-scores neither help nor hurt. Default stays 'false'
    // because 12 more features cost compute for nothing -- which is the right conclusion for a
    // reason the original comment got wrong.
    // -ms1_scores true hands the classifier ~12 more features (isotope correlation/overlap, mass
    // deviation, MS1 xcorr shape/coelution): a real precursor has the right isotope envelope at MS1,
    // which is information the MS2-only scope discards. The GBT treats missing as its own category
    // (odia_gbt_test T5), so sparsity is handled rather than imputed.
    const bool use_ms1_scores = getFlag_("ms1_scores");
    std::vector<std::string> vkeys, dropped;
    for (const auto& s : keyset)
    {
      std::string up = s;
      std::transform(up.begin(), up.end(), up.begin(), ::toupper);
      if (up.rfind("VAR_", 0) != 0) { continue; }
      if (!use_ms1_scores && up.rfind("VAR_MS1_", 0) == 0) { dropped.push_back(s); continue; }
      if (for_anchors && up == "VAR_NORM_RT_SCORE") { continue; }   // circular for RT anchors (C8)
      vkeys.push_back(s);
    }
    OPENMS_LOG_INFO << "OpenDIAlyzer: in-memory scoring found " << vkeys.size() << " VAR_ sub-scores ("
                    << (use_ms1_scores ? "MS1+MS2 scope" : "MS2 scope; " + std::to_string(dropped.size())
                                                           + " var_ms1_* excluded, -ms1_scores to include")
                    << ") across " << probe << " probed features (of " << fmap.size() << ")." << std::endl;
    if (vkeys.empty()) { return R; }

    R.feats.reserve(fmap.size());
    std::size_t skip_nogid = 0, skip_unmapped = 0, skip_baddecoy = 0;
    for (const Feature& f : fmap)
    {
      // The precursor this feature belongs to, and its decoy flag, come from the library rather
      // than from the feature: the feature only knows its transition group id.
      const std::string gid =
        f.metaValueExists("PeptideRef") ? std::string(f.getMetaValue("PeptideRef").toString())
                                        : std::string();
      if (gid.empty()) { ++skip_nogid; continue; }
      const auto it = prec_of_id.find(gid);
      if (it == prec_of_id.end()) { ++skip_unmapped; continue; }   // unmapped -> cannot be labelled
      const int dec = it->second.decoy;
      if (dec != 0 && dec != 1) { ++skip_baddecoy; continue; }

      std::vector<double> x(vkeys.size(), std::numeric_limits<double>::quiet_NaN());
      for (std::size_t j = 0; j < vkeys.size(); ++j)
      {
        if (!f.metaValueExists(vkeys[j])) { continue; }
        const DataValue& dv = f.getMetaValue(vkeys[j]);
        if (dv.isEmpty()) { continue; }
        const double v = static_cast<double>(dv);
        if (std::isfinite(v)) { x[j] = v; }              // NULL/inf -> missing, as in the sqlite path
      }
      R.feature_id.push_back(static_cast<long long>(f.getUniqueId()));
      R.group.push_back(it->second.id);
      R.exp_rt.push_back(f.getRT());
      R.labels.push_back(dec == 0 ? 1 : 0);
      // The LIBRARY's predicted RT, from the library -- NOT the feature's norm_RT.
      //
      // norm_RT is scores.normalized_experimental_rt (MRMFeatureFinderScoring.cpp:1018): the
      // OBSERVED retention time mapped into iRT space. The sqlite path reads PRECURSOR.LIBRARY_RT,
      // which is the library's PREDICTION. recalibrate_ fits library_rt -> exp_rt to learn the
      // library-to-run warp; feeding it norm_RT instead fits a function of exp_rt against exp_rt,
      // which is close to the identity and teaches the calibration nothing. Pass 2 then extracts on
      // an uncorrected RT axis.
      R.library_rt.push_back(it->second.library_rt);
      R.traml_id.push_back(gid);
      R.feats.push_back(std::move(x));
    }
    const std::size_t before_mem = R.feats.empty() ? 0 : R.feats[0].size();
    dropUninformativeColumns_(R);
    OPENMS_LOG_INFO << "OpenDIAlyzer[scoreload/memory] rows " << R.feats.size() << "/" << fmap.size()
                    << " (skipped: no group id " << skip_nogid << ", unmapped " << skip_unmapped
                    << ", bad decoy " << skip_baddecoy << "); columns " << before_mem << " -> "
                    << (R.feats.empty() ? 0 : R.feats[0].size())
                    << " after dropping uninformative; targets "
                    << std::count(R.labels.begin(), R.labels.end(), 1) << std::endl;
    R.canonicalize();
    return R;
  }

  // best-scoring row per precursor group (max d-score); returns group -> row index
  static std::unordered_map<long long, size_t> bestPerGroup_(const OswRows& R, const std::vector<double>& dscore)
  {
    std::unordered_map<long long, size_t> best;
    for (size_t i = 0; i < R.feats.size(); ++i)
    {
      auto it = best.find(R.group[i]);
      if (it == best.end() || dscore[i] > dscore[it->second]) { best[R.group[i]] = i; }
    }
    return best;
  }

  // Recalibrate: in-process LDA over the pass-1 .osw, take target peak groups at
  // q<0.01 as RT anchors, fit library_RT -> observed_RT.
  static bool hasBothClasses_(const OswRows& R)
  {
    bool t = false, d = false;
    for (int l : R.labels) { (l == 1 ? t : d) = true; if (t && d) { return true; } }
    return false;
  }

  // Fit the global library_RT -> observed_RT map from confident pass-1 target anchors.
  // Also (optionally) emit, per confident target, its OBSERVED apex RT keyed by TRAML_ID
  // (empirical_rt) so the caller can replace predicted RT with measured RT for seen
  // precursors (DIA-NN-style empirical library), and the 95th-pct anchor residual (p95_resid)
  // for residual-driven window sizing.
  /// Filled by recalibrate_(): the transition-group ids of the best-scoring TARGET peak groups from
  /// the pass just finished. These are the only precursors in the run with evidence behind them, and
  /// they are what the mass calibration should measure -- see calibrateMassFromPass_.
  std::vector<std::string> pass_confident_ids_;

  TransformationDescription recalibrate_(const std::string& osw, int& n_anchors,
                                         std::map<std::string, double>* empirical_rt = nullptr,
                                         double* p95_resid = nullptr)
  {
    n_anchors = 0;
    pass_confident_ids_.clear();
    OswRows R = parquet_out_
                  ? loadScoresFromFeatureMap_(pass_features_, precursor_index_, /*for_anchors=*/true)
                  : loadOswScores_(osw, /*for_anchors=*/true);   // exclude RT scores (C8)
    if (R.feats.size() < 50 || R.feats[0].empty() || !hasBothClasses_(R))
    {
      OPENMS_LOG_ERROR << "OpenDIAlyzer: too few scored features / missing a class for LDA recalibration." << std::endl;
      return TransformationDescription();
    }
    odia::ScoredGroups s = odia::scoreSemiSupervisedLDA(R.feats, R.labels, R.group);
    auto best = bestPerGroup_(R, s.dscore);
    // Anchors = the TOP target peak groups by d-score (NOT the honest q<0.01 set, which is
    // empty on the WIDE bootstrap pass). Approximate RT anchors are enough to fit the
    // transform; the binned-median fit tolerates the fraction that are wrong. The final
    // FDR (finalScore_) stays strict at q<0.01 -- only the RT bootstrap is lenient.
    std::vector<std::pair<double, size_t>> tgt;                 // (d-score, row)
    for (const auto& kv : best)
    {
      const size_t i = kv.second;
      if (R.labels[i] == 1 && std::isfinite(R.library_rt[i]) && std::isfinite(R.exp_rt[i])) { tgt.emplace_back(s.dscore[i], i); }
    }
    std::sort(tgt.begin(), tgt.end(), [](const auto& a, const auto& b) { return a.first > b.first; });
    const size_t take = std::min<size_t>(tgt.size(), 2000);    // top-2000 highest-confidence targets
    std::vector<std::pair<double, double>> pts;
    pts.reserve(take);
    for (size_t k = 0; k < take; ++k) { const size_t i = tgt[k].second; pts.emplace_back(R.library_rt[i], R.exp_rt[i]); }
    // Same ranked set, kept as ids: these precursors have real evidence in THIS run, which is
    // exactly what a mass calibration needs and what the iRT anchor list does not provide.
    pass_confident_ids_.reserve(take);
    for (size_t k = 0; k < take; ++k)
    {
      const size_t i = tgt[k].second;
      if (i < R.traml_id.size() && !R.traml_id[i].empty()) { pass_confident_ids_.push_back(R.traml_id[i]); }
    }
    n_anchors = static_cast<int>(pts.size());
    if (pts.size() < 20)
    {
      OPENMS_LOG_ERROR << "OpenDIAlyzer: only " << pts.size() << " target anchors for recalibration." << std::endl;
      return TransformationDescription();
    }
    // Empirical library (DIA-NN-style, FDR-SAFE by symmetry): EVERY precursor -- target AND
    // decoy -- that had a pass-1 peak gets ITS OWN observed apex RT (the best-scoring pass-1
    // peak group's EXP_RT). Pass 2 then centres each precursor's window on the peak it actually
    // found, collapsing the RT residual that a smooth calibration cannot remove (measured p95
    // ~342 s -- far wider than the narrow window, which is why pure recalibration collapses).
    //
    // Why this is FDR-safe where the old target-only scheme was not (codex): the procedure is
    // IDENTICAL for both classes and each uses its OWN apex. A false target locked onto its own
    // pass-1 fluke is mirrored by a decoy locked onto ITS own fluke -> the decoy null captures
    // the winner's-curse inflation -> target/decoy calibration is preserved. (The prior scheme
    // gave the decoy the TARGET's selected RT, not its own -> asymmetric.) No pairing / "DECOY_"
    // assumption needed; the final FDR (finalScore_) still validates via target-decoy on pass 2.
    if (empirical_rt)
    {
      empirical_rt->clear();
      size_t n_t = 0, n_d = 0;
      for (const auto& kv : best)
      {
        const size_t i = kv.second;
        if (R.traml_id[i].empty() || !std::isfinite(R.exp_rt[i])) { continue; }
        (*empirical_rt)[R.traml_id[i]] = R.exp_rt[i];          // own apex; both classes, same rule
        (R.labels[i] == 1 ? n_t : n_d)++;
      }
      OPENMS_LOG_INFO << "OpenDIAlyzer: empirical library -> " << n_t << " targets + " << n_d
                      << " decoys get their OWN observed apex RT (symmetric procedure -> FDR-safe)." << std::endl;
    }
    OPENMS_LOG_INFO << "OpenDIAlyzer: recalibration from top-" << pts.size() << " target anchors by d-score." << std::endl;
    return fitTrafo_(std::move(pts), p95_resid, getDoubleOption_("rt_fit_loess_span"));
  }

  // Write pyprophet-compatible SCORE_MS2 (FEATURE_ID, SCORE, RANK, PVALUE, QVALUE, PEP).
  void writeScoreMs2_(const std::string& osw, const OswRows& R, const odia::ScoredGroups& s)
  {
    sqlite3* db = nullptr;
    if (sqlite3_open(osw.c_str(), &db) != SQLITE_OK) { sqlite3_close(db); return; }
    sqlite3_exec(db, "DROP TABLE IF EXISTS SCORE_MS2;", nullptr, nullptr, nullptr);
    sqlite3_exec(db, "CREATE TABLE SCORE_MS2(FEATURE_ID INT, SCORE REAL, RANK INT, PVALUE REAL, QVALUE REAL, PEP REAL);", nullptr, nullptr, nullptr);
    auto best = bestPerGroup_(R, s.dscore);
    sqlite3_exec(db, "BEGIN;", nullptr, nullptr, nullptr);
    sqlite3_stmt* ins = nullptr;
    sqlite3_prepare_v2(db, "INSERT INTO SCORE_MS2(FEATURE_ID,SCORE,RANK,PVALUE,QVALUE,PEP) VALUES(?,?,?,?,?,?);", -1, &ins, nullptr);
    for (size_t i = 0; i < R.feature_id.size(); ++i)
    {
      const int rank = (best[R.group[i]] == i) ? 1 : 2;
      sqlite3_bind_int64(ins, 1, R.feature_id[i]);
      sqlite3_bind_double(ins, 2, s.dscore[i]);
      sqlite3_bind_int(ins, 3, rank);
      // PVALUE / QVALUE / PEP are three DIFFERENT statistics. All three used to be bound to the
      // q-value, so every row in every .osw this tool has ever written had them identical -- and
      // anything downstream reading PEP (IPF, protein-level inference) was consuming a q-value.
      sqlite3_bind_double(ins, 4, s.pvalue.empty() ? s.qvalue[i] : s.pvalue[i]);
      sqlite3_bind_double(ins, 5, s.qvalue[i]);
      sqlite3_bind_double(ins, 6, s.pep.empty() ? s.qvalue[i] : s.pep[i]);
      sqlite3_step(ins);
      sqlite3_reset(ins);
    }
    sqlite3_finalize(ins);
    sqlite3_exec(db, "COMMIT;", nullptr, nullptr, nullptr);
    sqlite3_close(db);
  }

  // Roll the precursor d-scores up to peptides and proteins and control FDR THERE too, writing
  // SCORE_PEPTIDE / SCORE_PROTEIN in the schema pyprophet produces (so anything downstream reads
  // them without knowing which tool wrote the file).
  //
  // This is not a nicety. A 1% precursor FDR is NOT a 1% peptide or protein FDR: a protein takes
  // the best of many independent precursor draws, so its score is an extreme-value statistic and
  // its error rate is strictly worse. Reporting a precursor q-value as though it were a protein
  // q-value is how a DIA result most commonly overstates its confidence, and until now that is
  // exactly what this tool's output invited a reader to do -- it wrote SCORE_MS2 and nothing else.
  //
  // See src/odia_fdr.h for the estimators and docs/OpenDIAlyzer-scoring-fdr-backlog.md (Area 2).
  void contextFdr_(const std::string& osw, const OswRows& R, const odia::ScoredGroups& s)
  {
    if (getStringOption_("fdr_context") == "none") { return; }

    // precursor id -> (modified sequence, peptide id) and peptide id -> protein accession.
    // Decoy status is taken from PRECURSOR (already validated at the precursor level) rather
    // than from PEPTIDE/PROTEIN, so the three levels cannot disagree about what a decoy is.
    std::unordered_map<long long, std::pair<long long, std::string>> pep_of_prec;
    std::unordered_map<long long, std::vector<std::string>> prot_of_pep;
    {
      sqlite3* db = nullptr;
      if (sqlite3_open(osw.c_str(), &db) != SQLITE_OK) { sqlite3_close(db); return; }
      sqlite3_stmt* st = nullptr;
      if (sqlite3_prepare_v2(db,
            "SELECT m.PRECURSOR_ID, p.ID, p.MODIFIED_SEQUENCE FROM PRECURSOR_PEPTIDE_MAPPING m "
            "JOIN PEPTIDE p ON m.PEPTIDE_ID = p.ID;", -1, &st, nullptr) == SQLITE_OK)
      {
        while (sqlite3_step(st) == SQLITE_ROW)
        {
          const char* seq = reinterpret_cast<const char*>(sqlite3_column_text(st, 2));
          pep_of_prec[sqlite3_column_int64(st, 0)] =
            {sqlite3_column_int64(st, 1), seq ? seq : ""};
        }
      }
      sqlite3_finalize(st);
      st = nullptr;
      if (sqlite3_prepare_v2(db,
            "SELECT m.PEPTIDE_ID, pr.PROTEIN_ACCESSION FROM PEPTIDE_PROTEIN_MAPPING m "
            "JOIN PROTEIN pr ON m.PROTEIN_ID = pr.ID;", -1, &st, nullptr) == SQLITE_OK)
      {
        while (sqlite3_step(st) == SQLITE_ROW)
        {
          const char* acc = reinterpret_cast<const char*>(sqlite3_column_text(st, 1));
          if (acc && *acc) { prot_of_pep[sqlite3_column_int64(st, 0)].push_back(acc); }
        }
      }
      sqlite3_finalize(st);
      sqlite3_close(db);
    }
    if (pep_of_prec.empty())
    {
      OPENMS_LOG_WARN << "OpenDIAlyzer: no PRECURSOR_PEPTIDE_MAPPING in the .osw; skipping "
                      << "peptide/protein-level FDR (precursor-level q-values are unaffected)."
                      << std::endl;
      return;
    }

    const bool use_pi0 = getFlag_("fdr_pi0");
    const std::string dtag = getStringOption_("decoy_tag");

    // ---- peptide level ---------------------------------------------------------------
    // Key on the MODIFIED sequence, so a peptide's charge states collapse into one entity.
    // Decoy peptides are keyed by the tag + sequence: their shuffled sequence could otherwise
    // collide with a real target sequence and silently merge a decoy into a target entity.
    std::vector<std::string> pep_key;
    std::vector<double> pep_score;
    std::vector<int> pep_label;
    std::unordered_map<std::string, long long> pep_id_of_key;
    pep_key.reserve(R.group.size());
    for (std::size_t i = 0; i < R.group.size(); ++i)
    {
      const auto it = pep_of_prec.find(R.group[i]);
      if (it == pep_of_prec.end() || it->second.second.empty())
      {
        pep_key.emplace_back();
        pep_score.push_back(0.0);
        pep_label.push_back(R.labels[i]);
        continue;
      }
      std::string key = (R.labels[i] == 1 ? "" : dtag) + it->second.second;
      pep_id_of_key.emplace(key, it->second.first);
      pep_key.push_back(std::move(key));
      pep_score.push_back(s.dscore[i]);
      pep_label.push_back(R.labels[i]);
    }
    auto peptides = odia::fdr::rollUp(pep_key, pep_score, pep_label);
    odia::fdr::assignQValues(peptides, use_pi0);

    // ---- protein level ---------------------------------------------------------------
    // Roll peptides up to proteins. A shared peptide contributes to every protein it maps to;
    // this is the "protein" context, not a parsimony/protein-group inference, and the log says
    // so rather than letting the column name imply more than was computed.
    std::vector<std::string> prot_key;
    std::vector<double> prot_score;
    std::vector<int> prot_label;
    for (const auto& p : peptides)
    {
      const auto id_it = pep_id_of_key.find(p.id);
      if (id_it == pep_id_of_key.end()) { continue; }
      const auto acc_it = prot_of_pep.find(id_it->second);
      if (acc_it == prot_of_pep.end()) { continue; }
      for (const auto& acc : acc_it->second)
      {
        // PROTEIN_ACCESSION already carries the tag for decoys (MRMDecoy.cpp:507); only add it
        // if it is missing, so the picked pairing below sees exactly one tag either way.
        const bool tagged = acc.size() > dtag.size() && acc.compare(0, dtag.size(), dtag) == 0;
        prot_key.push_back((p.label == 1 || tagged) ? acc : dtag + acc);
        prot_score.push_back(p.score);
        prot_label.push_back(p.label);
      }
    }
    auto proteins = odia::fdr::rollUp(prot_key, prot_score, prot_label);
    std::size_t n_paired = 0;
    const bool picked = getStringOption_("picked_protein") != "false";
    if (picked) { proteins = odia::fdr::pickedCompetition(proteins, dtag, &n_paired); }
    odia::fdr::assignQValues(proteins, use_pi0);

    writeScoreLevel_(osw, "SCORE_PEPTIDE", "PEPTIDE_ID", peptides, pep_id_of_key);
    writeScoreProtein_(osw, proteins);

    OPENMS_LOG_INFO << "OpenDIAlyzer: context FDR -> " << odia::fdr::countAtQ(peptides, 0.01)
                    << " peptides and " << odia::fdr::countAtQ(proteins, 0.01)
                    << " proteins at q<0.01"
                    << (picked ? " (picked target-decoy competition, " + std::to_string(n_paired)
                                   + " pairs)"
                               : " (plain target-decoy; -picked_protein false)")
                    << ". Protein level is per-accession, NOT a parsimony protein-group "
                    << "inference." << std::endl;

    // ---- entrapment check (only when the library was salted) --------------------------
    const std::string etag = getStringOption_("entrapment_tag");
    if (!etag.empty())
    {
      std::size_t db_ent = 0, db_tar = 0, rep = 0, rep_ent = 0;
      for (const auto& p : peptides)
      {
        if (p.label != 1) { continue; }
        const bool is_ent = p.id.size() > etag.size() && p.id.compare(0, etag.size(), etag) == 0;
        (is_ent ? db_ent : db_tar)++;
        if (p.qvalue <= 0.01) { ++rep; if (is_ent) { ++rep_ent; } }
      }
      const auto e = odia::fdr::entrapmentFdp(rep, rep_ent, db_tar, db_ent);
      if (!e.valid)
      {
        OPENMS_LOG_WARN << "OpenDIAlyzer: -entrapment_tag '" << etag << "' matched " << db_ent
                        << " of " << (db_ent + db_tar) << " target peptides; cannot estimate an "
                        << "entrapment FDP (need entrapment sequences AND discoveries)."
                        << std::endl;
      }
      else
      {
        OPENMS_LOG_INFO << "OpenDIAlyzer: ENTRAPMENT check at a nominal 1% peptide q-value: "
                        << rep_ent << " of " << rep << " discoveries are entrapments (r="
                        << e.ratio << ") -> estimated true FDP " << (100.0 * e.fdp) << "%. "
                        << (e.fdp > 0.02
                              ? "That is more than double the nominal rate -- the reported "
                                "q-values are OPTIMISTIC."
                              : "Consistent with the nominal rate.")
                        << std::endl;
      }
    }
  }

  // SCORE_PEPTIDE, in pyprophet's schema. `id_of_key` maps the entity key back to PEPTIDE.ID.
  void writeScoreLevel_(const std::string& osw, const std::string& table,
                        const std::string& id_column,
                        const std::vector<odia::fdr::Entity>& entities,
                        const std::unordered_map<std::string, long long>& id_of_key)
  {
    sqlite3* db = nullptr;
    if (sqlite3_open(osw.c_str(), &db) != SQLITE_OK) { sqlite3_close(db); return; }
    sqlite3_exec(db, ("DROP TABLE IF EXISTS " + table + ";").c_str(), nullptr, nullptr, nullptr);
    sqlite3_exec(db, ("CREATE TABLE " + table + "(CONTEXT TEXT, RUN_ID INT, " + id_column +
                      " INT, SCORE REAL, PVALUE REAL, QVALUE REAL, PEP REAL);").c_str(),
                 nullptr, nullptr, nullptr);
    sqlite3_exec(db, "BEGIN;", nullptr, nullptr, nullptr);
    sqlite3_stmt* ins = nullptr;
    sqlite3_prepare_v2(db, ("INSERT INTO " + table + "(CONTEXT,RUN_ID," + id_column +
                            ",SCORE,PVALUE,QVALUE,PEP) VALUES('global',NULL,?,?,?,?,?);").c_str(),
                       -1, &ins, nullptr);
    for (const auto& e : entities)
    {
      const auto it = id_of_key.find(e.id);
      if (it == id_of_key.end()) { continue; }
      sqlite3_bind_int64(ins, 1, it->second);
      sqlite3_bind_double(ins, 2, e.score);
      sqlite3_bind_double(ins, 3, e.pvalue);
      sqlite3_bind_double(ins, 4, e.qvalue);
      sqlite3_bind_double(ins, 5, e.pep);
      sqlite3_step(ins);
      sqlite3_reset(ins);
    }
    sqlite3_finalize(ins);
    sqlite3_exec(db, "COMMIT;", nullptr, nullptr, nullptr);
    sqlite3_close(db);
  }

  // SCORE_PROTEIN. Resolved by accession rather than by a precomputed id map, because picked
  // competition drops one member of every pair and the surviving accession is not known until
  // after it has run.
  void writeScoreProtein_(const std::string& osw,
                          const std::vector<odia::fdr::Entity>& proteins)
  {
    sqlite3* db = nullptr;
    if (sqlite3_open(osw.c_str(), &db) != SQLITE_OK) { sqlite3_close(db); return; }
    std::unordered_map<std::string, long long> id_of_accession;
    sqlite3_stmt* st = nullptr;
    if (sqlite3_prepare_v2(db, "SELECT ID, PROTEIN_ACCESSION FROM PROTEIN;", -1, &st, nullptr)
        == SQLITE_OK)
    {
      while (sqlite3_step(st) == SQLITE_ROW)
      {
        const char* acc = reinterpret_cast<const char*>(sqlite3_column_text(st, 1));
        if (acc && *acc) { id_of_accession[acc] = sqlite3_column_int64(st, 0); }
      }
    }
    sqlite3_finalize(st);
    sqlite3_exec(db, "DROP TABLE IF EXISTS SCORE_PROTEIN;", nullptr, nullptr, nullptr);
    sqlite3_exec(db, "CREATE TABLE SCORE_PROTEIN(CONTEXT TEXT, RUN_ID INT, PROTEIN_ID INT, "
                     "SCORE REAL, PVALUE REAL, QVALUE REAL, PEP REAL);",
                 nullptr, nullptr, nullptr);
    sqlite3_exec(db, "BEGIN;", nullptr, nullptr, nullptr);
    sqlite3_stmt* ins = nullptr;
    sqlite3_prepare_v2(db, "INSERT INTO SCORE_PROTEIN(CONTEXT,RUN_ID,PROTEIN_ID,SCORE,PVALUE,"
                           "QVALUE,PEP) VALUES('global',NULL,?,?,?,?,?);",
                       -1, &ins, nullptr);
    for (const auto& e : proteins)
    {
      const auto it = id_of_accession.find(e.id);
      if (it == id_of_accession.end()) { continue; }
      sqlite3_bind_int64(ins, 1, it->second);
      sqlite3_bind_double(ins, 2, e.score);
      sqlite3_bind_double(ins, 3, e.pvalue);
      sqlite3_bind_double(ins, 4, e.qvalue);
      sqlite3_bind_double(ins, 5, e.pep);
      sqlite3_step(ins);
      sqlite3_reset(ins);
    }
    sqlite3_finalize(ins);
    sqlite3_exec(db, "COMMIT;", nullptr, nullptr, nullptr);
    sqlite3_close(db);
  }

  // SCORE_MS2 as parquet, written into the run's bundle directory. Same columns as the sqlite
  // table so anything that reads one can read the other; the difference is that this is a typed
  // columnar file rather than rows in a database that exists only to be read back.
  void writeScoresParquet_(const std::string& out_dir, const OswRows& R, const odia::ScoredGroups& s)
  {
    auto best = bestPerGroup_(R, s.dscore);
    arrow::Int64Builder feature_id;
    arrow::DoubleBuilder score, pvalue, qvalue, pep;
    arrow::Int32Builder rank;
    for (std::size_t i = 0; i < R.feature_id.size(); ++i)
    {
      ParquetFile::appendOrThrow(feature_id.Append(R.feature_id[i]), "feature_id");
      ParquetFile::appendOrThrow(score.Append(s.dscore[i]), "score");
      ParquetFile::appendOrThrow(rank.Append(best[R.group[i]] == i ? 1 : 2), "rank");
      // pvalue / qvalue / pep are three DIFFERENT statistics; writing the q-value into all three
      // (as this tool once did for sqlite) silently feeds a q-value to anything reading PEP.
      ParquetFile::appendOrThrow(pvalue.Append(s.pvalue.empty() ? s.qvalue[i] : s.pvalue[i]), "pvalue");
      ParquetFile::appendOrThrow(qvalue.Append(s.qvalue[i]), "qvalue");
      ParquetFile::appendOrThrow(pep.Append(s.pep.empty() ? s.qvalue[i] : s.pep[i]), "pep");
    }
    auto schema = arrow::schema({arrow::field("feature_id", arrow::int64()),
                                 arrow::field("score", arrow::float64()),
                                 arrow::field("rank", arrow::int32()),
                                 arrow::field("pvalue", arrow::float64()),
                                 arrow::field("qvalue", arrow::float64()),
                                 arrow::field("pep", arrow::float64())});
    auto table = arrow::Table::Make(schema,
      {ParquetFile::finishArray(feature_id, "feature_id"), ParquetFile::finishArray(score, "score"),
       ParquetFile::finishArray(rank, "rank"), ParquetFile::finishArray(pvalue, "pvalue"),
       ParquetFile::finishArray(qvalue, "qvalue"), ParquetFile::finishArray(pep, "pep")});
    // .oswpq is a ZIP archive (library/*, runs/run_id=*/features.parquet), not a directory -- the
    // first version of this wrote out_dir + "/score_ms2.parquet" and got "not writable for the
    // current user" after a 47-minute run, because that path is INSIDE the archive file. Write to a
    // temp file and add the entry, which is what the TOPP tool does for its own tables.
    addTableToBundle_(out_dir, "runs/run_id=" + std::to_string(run_id_) + "/score_ms2.parquet", table);
  }

  /// Write an arrow table into the .oswpq archive under @p entry_name.
  void addTableToBundle_(const std::string& archive, const std::string& entry_name,
                         const std::shared_ptr<arrow::Table>& table)
  {
    File::TempDir tmp;
    const std::string local = tmp.getPath() + "/" + File::basename(entry_name);
    ParquetFile::writeTable(table, local);
    ZipArchiveFile::addOrReplaceFromFile(archive, entry_name, local);
  }

  // Peptide- and protein-level FDR for the parquet path. Identical estimators to contextFdr_ --
  // odia::fdr::rollUp / assignQValues / pickedCompetition, so the two output formats cannot drift
  // apart statistically. What differs is only where the mappings come from and where the result
  // goes: precursor->peptide and peptide->protein are read from the LIBRARY (LightCompound carries
  // the modified sequence and the protein refs) rather than from PRECURSOR_PEPTIDE_MAPPING /
  // PEPTIDE_PROTEIN_MAPPING, and the rows are written as parquet keyed by sequence/accession
  // instead of by the integer ids a database needed.
  void contextFdrParquet_(const std::string& out_dir, const OswRows& R, const odia::ScoredGroups& s)
  {
    if (getStringOption_("fdr_context") == "none") { return; }
    const bool use_pi0 = getFlag_("fdr_pi0");
    const std::string dtag = getStringOption_("decoy_tag");

    // group id -> its library entry, by the integer id the rows carry
    std::unordered_map<long long, const PrecursorMeta*> by_id;
    by_id.reserve(precursor_index_.size());
    for (const auto& kv : precursor_index_) { by_id[kv.second.id] = &kv.second; }

    // ---- peptide level: key on the MODIFIED sequence so charge states collapse -------------
    // Decoys are prefixed, because a shuffled decoy sequence can collide with a real target
    // sequence and would otherwise be merged into the target's entity.
    std::vector<std::string> pep_key;
    std::vector<double> pep_score;
    std::vector<int> pep_label;
    pep_key.reserve(R.group.size());
    for (std::size_t i = 0; i < R.group.size(); ++i)
    {
      const auto it = by_id.find(R.group[i]);
      if (it == by_id.end() || it->second->sequence.empty())
      {
        pep_key.emplace_back(); pep_score.push_back(0.0); pep_label.push_back(R.labels[i]);
        continue;
      }
      pep_key.push_back((R.labels[i] == 1 ? std::string() : dtag) + it->second->sequence);
      pep_score.push_back(s.dscore[i]);
      pep_label.push_back(R.labels[i]);
    }
    auto peptides = odia::fdr::rollUp(pep_key, pep_score, pep_label);
    odia::fdr::assignQValues(peptides, use_pi0);

    // ---- protein level ---------------------------------------------------------------------
    std::unordered_map<std::string, const PrecursorMeta*> meta_of_pep;
    for (std::size_t i = 0; i < pep_key.size(); ++i)
    {
      if (pep_key[i].empty()) { continue; }
      const auto it = by_id.find(R.group[i]);
      if (it != by_id.end()) { meta_of_pep.emplace(pep_key[i], it->second); }
    }
    std::vector<std::string> prot_key;
    std::vector<double> prot_score;
    std::vector<int> prot_label;
    for (const auto& pep : peptides)
    {
      const auto it = meta_of_pep.find(pep.id);
      if (it == meta_of_pep.end()) { continue; }
      for (const auto& acc : it->second->proteins)
      {
        if (acc.empty()) { continue; }
        const bool tagged = acc.size() > dtag.size() && acc.compare(0, dtag.size(), dtag) == 0;
        prot_key.push_back((pep.label == 1 || tagged) ? acc : dtag + acc);
        prot_score.push_back(pep.score);
        prot_label.push_back(pep.label);
      }
    }
    auto proteins = odia::fdr::rollUp(prot_key, prot_score, prot_label);
    std::size_t n_paired = 0;
    const bool picked = getStringOption_("picked_protein") != "false";
    if (picked) { proteins = odia::fdr::pickedCompetition(proteins, dtag, &n_paired); }
    odia::fdr::assignQValues(proteins, use_pi0);

    writeEntitiesParquet_(out_dir, "runs/run_id=" + std::to_string(run_id_) + "/score_peptide.parquet",
                          "modified_sequence", peptides);
    writeEntitiesParquet_(out_dir, "runs/run_id=" + std::to_string(run_id_) + "/score_protein.parquet",
                          "protein_accession", proteins);

    OPENMS_LOG_INFO << "OpenDIAlyzer: context FDR -> " << odia::fdr::countAtQ(peptides, 0.01)
                    << " peptides and " << odia::fdr::countAtQ(proteins, 0.01)
                    << " proteins at q<0.01"
                    << (picked ? " (picked target-decoy competition, " + std::to_string(n_paired)
                                   + " pairs)"
                               : " (plain target-decoy; -picked_protein false)")
                    << ". Protein level is per-accession, NOT a parsimony protein-group inference."
                    << std::endl;
  }

  // One entity table (peptides or proteins) as parquet. Keyed by the entity's own name rather than
  // by an integer id: without a database there is no id to join back to, and the name is what a
  // reader actually wants.
  void writeEntitiesParquet_(const std::string& archive, const std::string& entry_name,
                             const std::string& key_column,
                             const std::vector<odia::fdr::Entity>& entities)
  {
    arrow::StringBuilder key;
    arrow::BooleanBuilder decoy;
    arrow::DoubleBuilder score, pvalue, qvalue, pep;
    for (const auto& e : entities)
    {
      ParquetFile::appendOrThrow(key.Append(e.id), "key");
      ParquetFile::appendOrThrow(decoy.Append(e.label != 1), "decoy");
      ParquetFile::appendOrThrow(score.Append(e.score), "score");
      ParquetFile::appendOrThrow(pvalue.Append(e.pvalue), "pvalue");
      ParquetFile::appendOrThrow(qvalue.Append(e.qvalue), "qvalue");
      ParquetFile::appendOrThrow(pep.Append(e.pep), "pep");
    }
    auto schema = arrow::schema({arrow::field(key_column, arrow::utf8()),
                                 arrow::field("decoy", arrow::boolean()),
                                 arrow::field("score", arrow::float64()),
                                 arrow::field("pvalue", arrow::float64()),
                                 arrow::field("qvalue", arrow::float64()),
                                 arrow::field("pep", arrow::float64())});
    auto table = arrow::Table::Make(schema,
      {ParquetFile::finishArray(key, key_column), ParquetFile::finishArray(decoy, "decoy"),
       ParquetFile::finishArray(score, "score"), ParquetFile::finishArray(pvalue, "pvalue"),
       ParquetFile::finishArray(qvalue, "qvalue"), ParquetFile::finishArray(pep, "pep")});
    addTableToBundle_(archive, entry_name, table);
  }

  // Final FDR: in-process LDA on the last pass, write SCORE_MS2, report IDs@1%.
  int finalScore_(const std::string& osw)
  {
    if (parquet_out_)
    {
      // pass_features_ is the other structure that OUTLIVES the phase that fills it. Meta values
      // dominate it: ~30 named sub-scores per feature, each a string key plus a DataValue.
      // Count what actually costs, not just the top-level array. sizeof(Feature) is 296 B, but
      // MetaInfoInterface is an 8-BYTE POINTER to a separately heap-allocated MetaInfo holding a
      // flat_map<UInt, DataValue> -- so the meta values are invisible in sizeof() and are the
      // larger term. Subordinates are a vector<Feature> per feature (the per-transition
      // sub-features in MRM), each 296 B with a MetaInfo of its own.
      if (getFlag_("mem_components"))
      {
      PhaseTimer pt_mem("mem_components");
      std::size_t mv = 0, subs = 0, sub_mv = 0, str_vals = 0;
      std::vector<std::string> meta_keys;
      std::vector<std::string> sk;                 // hoisted: was constructed 30.2M times
      for (const Feature& f : pass_features_)
      {
        meta_keys.clear();
        f.getKeys(meta_keys);
        mv += meta_keys.size();
        for (const auto& k : meta_keys)
        {
          if (f.getMetaValue(k).valueType() == DataValue::STRING_VALUE) { ++str_vals; }
        }
        subs += f.getSubordinates().size();
        for (const Feature& sf : f.getSubordinates())
        {
          sk.clear();
          sf.getKeys(sk);
          sub_mv += sk.size();
        }
      }
      const double gb = 1073741824.0;
      const std::size_t pair_b = sizeof(UInt) + sizeof(DataValue) + 4;   // flat_map element, padded
      OPENMS_LOG_INFO << "OpenDIAlyzer[mem/component] feature map detail: "
                      << std::fixed << std::setprecision(2)
                      << (pass_features_.size() * sizeof(Feature)) / gb << " GB top-level Features; "
                      << subs << " subordinates = " << (subs * sizeof(Feature)) / gb << " GB; "
                      << mv << "+" << sub_mv << " meta values ~= "
                      << ((mv + sub_mv) * pair_b) / gb << " GB in flat_maps ("
                      << str_vals << " are strings, each an extra allocation); "
                      << (pass_features_.size() + subs) << " MetaInfo allocations" << std::endl;
      OPENMS_LOG_INFO << "OpenDIAlyzer[mem/component] pass_features_: " << pass_features_.size()
                      << " features, " << mv << " meta values, "
                      << std::fixed << std::setprecision(2)
                      << (pass_features_.size() * sizeof(Feature)) / 1073741824.0
                      << " GB in Feature objects alone (meta values are extra and not counted here)"
                      << std::endl;
      }
      MemProbe::logAllocator("before score_load");
    }
    OswRows R;
    { PhaseTimer pt("score_load");
      R = parquet_out_ ? loadScoresFromFeatureMap_(pass_features_, precursor_index_)
                       : loadOswScores_(osw); }
    if (R.feats.empty() || R.feats[0].empty() || !hasBothClasses_(R))
    {
      OPENMS_LOG_ERROR << "OpenDIAlyzer: no scored features or missing target/decoy class for final FDR (cannot control FDR)." << std::endl;
      return 0;
    }
    odia::LDAParams p;
    p.use_pi0 = getFlag_("fdr_pi0");
    p.n_folds = getIntOption_("lda_folds");
    p.n_iter = getIntOption_("lda_iterations");
    p.train_fdr_initial = getDoubleOption_("lda_train_fdr_initial");
    p.train_fdr = getDoubleOption_("lda_train_fdr");
    p.classifier = (getStringOption_("classifier") == "gbt") ? odia::LDAParams::Classifier::GBT
                                                             : odia::LDAParams::Classifier::LDA;
    odia::ScoredGroups s;
    { PhaseTimer pt(p.classifier == odia::LDAParams::Classifier::GBT ? "classifier_fit_gbt"
                                                                     : "classifier_fit_lda");
      s = odia::scoreSemiSupervisedLDA(R.feats, R.labels, R.group, p); }
    auto best = bestPerGroup_(R, s.dscore);
    if (s.n_iterations_trained == 0)
    {
      OPENMS_LOG_WARN << "OpenDIAlyzer: the semi-supervised LDA never fitted a discriminant ("
                      << s.n_iterations_skipped << " iterations skipped for want of confident "
                      << "positives). Scores come from the single-feature initialisation, NOT a "
                      << "learned model. Raise -lda_train_fdr_initial." << std::endl;
    }
    else
    {
      OPENMS_LOG_INFO << "OpenDIAlyzer: LDA fitted " << s.n_iterations_trained
                      << " iteration(s), skipped " << s.n_iterations_skipped << "." << std::endl;
    }
    int ids = 0;
    for (const auto& kv : best) { const size_t i = kv.second; if (R.labels[i] == 1 && s.qvalue[i] < 0.01) { ++ids; } }
    if (parquet_out_) { PhaseTimer pt("write_scores"); writeScoresParquet_(osw, R, s); }
    else               { PhaseTimer pt("write_scores_sqlite"); writeScoreMs2_(osw, R, s); }
    OPENMS_LOG_INFO << "OpenDIAlyzer: in-process LDA FDR -> " << ids << " target precursors at q<0.01." << std::endl;
    { PhaseTimer pt("context_fdr");
      if (parquet_out_) { contextFdrParquet_(osw, R, s); }
      else               { contextFdr_(osw, R, s); } }
    return ids;
  }

  // Feature-finder configuration COPIED VERBATIM from TOPP OpenSwathWorkflow's "Scoring"
  // defaults (src/topp/OpenSwathWorkflow.cpp:394-450). Previously we set only the two MS1
  // flags and inherited MRMFeatureFinderScoring's raw defaults for everything else -- which
  // silently disabled 6 discriminating sub-scores (VAR_MI_*, VAR_IM_*: all NULL in our .osw
  // vs populated in OpenSWATH's) and left the peak picker unconfigured. Matching this block
  // is a prerequisite for reproducing OpenSWATH. Shared by calibration and extraction.
  static Param makeFeatureFinderParam_()
  {
    Param ff = MRMFeatureFinderScoring().getDefaults();
    ff.remove("rt_extraction_window");
    ff.setValue("stop_report_after_feature", 5);
    ff.setValue("rt_normalization_factor", 100.0);   // iRT peptides live on a ~0..100 scale
    ff.setValue("Scores:use_ms1_correlation", "true");
    ff.setValue("Scores:use_ms1_fullscan", "true");
    ff.setValue("Scores:use_ms1_mi", "true");
    ff.setValue("Scores:use_mi_score", "true");

    ff.setValue("TransitionGroupPicker:min_peak_width", -1.0);
    ff.setValue("TransitionGroupPicker:recalculate_peaks", "true");
    ff.setValue("TransitionGroupPicker:compute_peak_quality", "false");
    ff.setValue("TransitionGroupPicker:minimal_quality", -1.5);
    ff.setValue("TransitionGroupPicker:background_subtraction", "none");
    ff.setValue("TransitionGroupPicker:compute_peak_shape_metrics", "false");
    ff.remove("TransitionGroupPicker:stop_after_intensity_ratio");
    ff.setValue("TransitionGroupPicker:recalculate_peaks_max_z", 0.75);

    ff.setValue("TransitionGroupPicker:PeakPickerChromatogram:use_gauss", "false");
    ff.setValue("TransitionGroupPicker:PeakPickerChromatogram:sgolay_polynomial_order", 3);
    ff.setValue("TransitionGroupPicker:PeakPickerChromatogram:sgolay_frame_length", 11);
    ff.setValue("TransitionGroupPicker:PeakPickerChromatogram:peak_width", -1.0);
    ff.setValue("TransitionGroupPicker:PeakPickerChromatogram:remove_overlapping_peaks", "true");
    ff.setValue("TransitionGroupPicker:PeakPickerChromatogram:write_sn_log_messages", "false");
    ff.setValue("TransitionGroupPicker:PeakPickerChromatogram:method", "corrected");
    ff.setValue("TransitionGroupPicker:PeakPickerChromatogram:signal_to_noise", 0.1);
    ff.setValue("TransitionGroupPicker:PeakPickerChromatogram:gauss_width", 30.0);
    ff.setValue("uis_threshold_sn", -1);
    ff.setValue("uis_threshold_peak_area", 0);
    ff.remove("TransitionGroupPicker:PeakPickerChromatogram:sn_win_len");
    ff.remove("TransitionGroupPicker:PeakPickerChromatogram:sn_bin_count");
    ff.remove("TransitionGroupPicker:PeakPickerChromatogram:stop_after_feature");

    // EMG scoring is very CPU-intensive -> off, as in the TOPP tool.
    ff.remove("Scores:use_elution_model_score");
    ff.setValue("EMGScoring:max_iteration", 10);
    ff.remove("EMGScoring:interpolation_step");
    ff.remove("EMGScoring:tolerance_stdev_bounding_box");
    ff.remove("EMGScoring:deltaAbsError");
    ff.remove("EMGScoring:statistics:mean");
    ff.remove("EMGScoring:statistics:variance");
    return ff;
  }

  // OpenSWATH-equivalent RT calibration: auto-sample iRT anchor peptides from the library
  // (prioritising the built-in CiRT kit), extract and score ONLY those, outlier-filter, and
  // fit the iRT -> run-RT transform (OpenMS CalibrationWorkflow -- the very code standalone
  // OpenSwathWorkflow runs). This replaces our data-free linear bootstrap, which was the
  // last measured difference vs OpenSWATH (same precursors and feature counts, but different
  // peaks selected because the expected RT fed into scoring was wrong).
  // Returns false if calibration could not be performed (caller falls back to bootstrap).
  // ---- robust location/scale for calibration residuals ------------------------------------
  //
  // Half-sample mode (Bickel & Fruhwirth 2006). Recursively keep the SHORTEST interval containing
  // half the remaining points; the limit is the densest region, i.e. the mode. Chosen over a KDE
  // mode-finder for three reasons that matter here:
  //   * no bandwidth to select -- and the literature is explicit that Silverman's critical-bandwidth
  //     calibration is not correct for general modality work, so a bandwidth rule would itself need
  //     defending;
  //   * it is reported more outlier-resistant than KDE-based low-bias mode estimators, and our
  //     anchor sample is contaminated by construction (~99.8% of library targets are absent);
  //   * the dense core it converges to is also where the scale must be measured (see
  //     localScaleAboutMode_), which is the other number the hook has to return.
  // O(n log n) after the sort. `v` must be sorted ascending.
  static double halfSampleMode_(const std::vector<double>& v)
  {
    const std::size_t n = v.size();
    if (n == 0) { return 0.0; }
    if (n == 1) { return v[0]; }
    if (n == 2) { return 0.5 * (v[0] + v[1]); }
    if (n == 3)
    {
      const double d1 = v[1] - v[0], d2 = v[2] - v[1];
      if (d1 < d2) { return 0.5 * (v[0] + v[1]); }
      if (d2 < d1) { return 0.5 * (v[1] + v[2]); }
      return v[1];
    }
    const std::size_t k = (n + 1) / 2;                 // ceil(n/2)
    std::size_t best = 0;
    double best_w = std::numeric_limits<double>::max();
    for (std::size_t i = 0; i + k - 1 < n; ++i)
    {
      const double w = v[i + k - 1] - v[i];
      if (w < best_w) { best_w = w; best = i; }
    }
    return halfSampleMode_(std::vector<double>(v.begin() + best, v.begin() + best + k));
  }

  /// Robust scale of the signal component, estimated LOCALLY about a known mode.
  ///
  /// The obvious choice -- the shorth (width of the shortest interval holding half the points,
  /// /1.349) -- is wrong here, and the selftest caught it: the shorth has 50% breakdown, so once
  /// the signal is a MINORITY the shortest half is forced to include background and the estimate
  /// blows up (measured 11.6 ppm for a 1 ppm signal at 40% purity). The mode itself survives that
  /// because halfSampleMode_ recurses into the dense core; the scale must be estimated the same
  /// way -- locally.
  ///
  /// So: restrict to a shrinking neighbourhood of the mode and take the MAD there. Inside a tight
  /// neighbourhood the signal is the majority even when it is a global minority, which is exactly
  /// the condition MAD needs. Two or three passes converge; more would start chasing the tail.
  static double localScaleAboutMode_(const std::vector<double>& v, double mode, double init_window)
  {
    if (v.size() < 8 || !(init_window > 0.0)) { return 0.0; }
    double w = init_window;
    double sigma = 0.0;
    for (int pass = 0; pass < 3; ++pass)
    {
      std::vector<double> dev;
      dev.reserve(v.size());
      for (double x : v) { if (std::abs(x - mode) <= w) { dev.push_back(std::abs(x - mode)); } }
      if (dev.size() < 8) { break; }                    // neighbourhood too thin to trust
      std::nth_element(dev.begin(), dev.begin() + dev.size() / 2, dev.end());
      const double mad = dev[dev.size() / 2];
      if (!(mad > 0.0)) { break; }
      sigma = 1.4826 * mad;                             // MAD -> Gaussian-equivalent sigma
      const double next = 3.0 * sigma;
      if (next >= w) { break; }                         // converged / not shrinking
      w = next;
    }
    return sigma;
  }

  /// Is the residual sample actually PEAKED, or is it flat noise?
  ///
  /// This is the guard that stops the estimator widening itself to death. In the regime this tool
  /// runs in, most library targets are absent, so a wide search window returns the nearest unrelated
  /// centroid and those errors are ~UNIFORM over the window. A uniform sample still has a perfectly
  /// well-defined mode and scale -- they are just meaningless, and acting on them produces a window
  /// that grows until it swallows everything. Real measurement error is centrally peaked; noise is
  /// not. Compare density near zero against density out at the edge: ~1 for uniform, >>1 for a
  /// genuine error distribution.
  static double peakednessRatio_(const std::vector<double>& abs_err, double window)
  {
    if (window <= 0.0 || abs_err.empty()) { return 0.0; }
    std::size_t c = 0, e = 0;
    for (double a : abs_err)
    {
      if (a <= 0.2 * window) { ++c; }
      else if (a >= 0.6 * window && a <= 0.8 * window) { ++e; }
    }
    // Both bands are 0.2*window wide, so raw counts are already comparable densities.
    if (e == 0) { return c > 0 ? 1e9 : 0.0; }          // no edge mass at all -> maximally peaked
    return static_cast<double>(c) / static_cast<double>(e);
  }

  /// Measure the run's real MS1/MS2 mass accuracy from the precursors that pass 1 actually found,
  /// and use it for the next pass.
  ///
  /// The hook already existed and already worked; what it lacked was anchors. It was being handed
  /// the iRT anchor list, and the peakedness gate correctly refused to fit: inverting the statistic
  /// (centre |err| <= 0.2W against edge 0.6-0.8W over a 50 ppm bootstrap, ratio = 1 + 5N/M) puts the
  /// observed 1.14 at **2.8% real anchors** for MS2 and 1.69 at 12% for MS1. A mode and a width
  /// fitted to 97% noise would be meaningless, so declining was right -- and narrowing the bootstrap
  /// does not help, because it shrinks signal and noise together.
  ///
  /// Pass-1's top target peak groups by d-score are a different population entirely: they have
  /// measured evidence. This is what DIA-NN does, and on this run DIA-NN reports 1.66 ppm (2.0) and
  /// 1.71 ppm (1.7.12) while OpenDIAlyzer extracts at the configured 10 ppm -- six times wider than
  /// the instrument's real error, admitting interference into every correlation sub-score.
  void calibrateMassFromPass_(const std::vector<OpenSwath::SwathMap>& swath_maps,
                              const OpenSwath::LightTargetedExperiment& searched,
                              const TransformationDescription& rt_trafo)
  {
    if (pass_confident_ids_.empty()) { return; }
    const std::unordered_set<std::string> keep(pass_confident_ids_.begin(), pass_confident_ids_.end());

    OpenSwath::LightTargetedExperiment anchors;
    for (const auto& c : searched.getCompounds()) { if (keep.count(c.id)) { anchors.compounds.push_back(c); } }
    for (const auto& t : searched.getTransitions())
    {
      if (keep.count(t.getPeptideRef())) { anchors.transitions.push_back(t); }
    }
    if (anchors.compounds.empty()) { return; }

    double ms2 = inferMassAccuracyPpm_(swath_maps, anchors, rt_trafo, false);
    double ms1 = inferMassAccuracyPpm_(swath_maps, anchors, rt_trafo, true);

    // CALIBRATION MAY ONLY NARROW. The configured window is the user's assertion about the
    // instrument; the point of measuring is to discover the instrument is BETTER than that and
    // tighten accordingly. A fit that says "actually, use a wider window" is telling you the fit
    // failed, not that the instrument is bad -- and acting on it is strictly worse than doing
    // nothing, because it admits interference the user explicitly excluded.
    //
    // This is not hypothetical. On the first run with pass-1 anchors the MS1 fit returned
    // offset 24.6 ppm, sigma 24.2 ppm, peakedness exactly 3.0 -- noise filling the 50 ppm bootstrap
    // window, squeaking past the gate -- and set the MS1 window to 72.6 ppm on an instrument
    // measured at 1.66 ppm. Identifications fell from 6,798 to 4,496.
    const double cfg_ms2 = getDoubleOption_("mz_extraction_window");
    const double cfg_ms1 = getDoubleOption_("mz_extraction_window_ms1") > 0.0
                             ? getDoubleOption_("mz_extraction_window_ms1") : cfg_ms2;
    if (ms2 > cfg_ms2)
    {
      OPENMS_LOG_WARN << "OpenDIAlyzer: MS2 mass calibration returned " << ms2 << " ppm, WIDER than "
                      << "the configured " << cfg_ms2 << " ppm -- rejecting it. Calibration may only "
                      << "narrow; a wider answer means the anchors did not support a fit."
                      << std::endl;
      ms2 = -1.0;
    }
    if (ms1 > cfg_ms1)
    {
      OPENMS_LOG_WARN << "OpenDIAlyzer: MS1 mass calibration returned " << ms1 << " ppm, WIDER than "
                      << "the configured " << cfg_ms1 << " ppm -- rejecting it (same rule)."
                      << std::endl;
      ms1 = -1.0;
    }
    if (ms2 > 0.0) { inferred_mz_window_ms2_ = ms2; }
    if (ms1 > 0.0) { inferred_mz_window_ms1_ = ms1; }
    if (ms2 > 0.0 || ms1 > 0.0)
    {
      OPENMS_LOG_INFO << "OpenDIAlyzer: mass accuracy inferred from " << anchors.compounds.size()
                      << " pass-1 confident precursors -- MS2 " << ms2 << " ppm, MS1 " << ms1
                      << " ppm (overriding the configured window for the next pass)." << std::endl;
    }
    else
    {
      OPENMS_LOG_INFO << "OpenDIAlyzer: mass accuracy still not inferable from "
                      << anchors.compounds.size() << " pass-1 anchors; keeping the configured window."
                      << std::endl;
    }
  }

  // ---- calibration inference hooks -------------------------------------------------------
  //
  // Both are DELIBERATELY EMPTY. They exist so the two data-driven calibrations DIA-NN performs and
  // this tool does not have a defined place to land, with the call sites, inputs and return
  // contract already wired.
  //
  // Contract for both: return a POSITIVE value to override the configured setting, or <= 0 to mean
  // "not inferred -- keep what the user configured". A hook must never silently return a plausible
  // default; that is how an unimplemented step gets mistaken for a working one.
  //
  // WHY these two specifically (measured on astral.mzML, same library, 2026-07-31):
  //   * Mass accuracy. DIA-NN logs `Calibrating with mass accuracies 30 (MS1), 20 (MS2)` and then
  //     `Optimised mass accuracy: 6 ppm`. It starts wide, measures its own error, and narrows.
  //     This tool used a fixed 30 ppm against an instrument delivering 1.63 ppm (0.93 corrected),
  //     admitting ~20x more interfering signal than necessary. Effect of narrowing by hand:
  //     30 ppm -> 0 identifications, 27.3M features, 881 GB;
  //     10 ppm -> 2,957 identifications, 2.07M features, 189 GB.
  //     Automating this is the single highest-value calibration change available.
  //   * RT window. The window currently comes from `estimateWindow()` over the anchors that
  //     survived both the discovery window and outlier removal -- in-sample, censored and
  //     survivor-selected. It produced 746.9-1135.5 s where the true predicted-vs-observed residual
  //     on unseen precursors is p95 ~1720 s. A held-out estimate belongs here instead.
  //
  /// HOOK (empty): infer the m/z extraction window in ppm from observed mass error.
  /// Intended implementation: extract at the configured (wide) window, take confident anchor
  /// fragments, build the observed-minus-theoretical m/z distribution, and return ~5x its robust
  /// sd -- separately for MS1 and MS2, which have different error scales.
  /// @param ms1 true for the MS1 window, false for MS2.
  /// @return ppm window to use, or <= 0 to keep the configured value.
  double inferMassAccuracyPpm_(const std::vector<OpenSwath::SwathMap>& swath_maps,
                               const OpenSwath::LightTargetedExperiment& anchors,
                               const TransformationDescription& rt_trafo,
                               bool ms1) const
  {
    // The bootstrap is a HALF-width, so the configured default of 50 searched +/-50 ppm = 100 ppm
    // full width -- 10x the 10 ppm extraction window whose error it is trying to measure. At this
    // instrument's MS2 peak density (~1.37 centroids/Da, so ~0.07 peaks in a +/-50 ppm window at
    // m/z 500) that window is overwhelmingly interference, and "most intense within it" selects an
    // interferent whenever one outshines the true fragment. The residual distribution then comes
    // out FLAT -- which is exactly what the peakedness gate reported (1.14-2.64 against a
    // threshold of 3), correctly refusing to fit noise.
    //
    // The configured extraction window is the user's assertion that the true error lies inside it.
    // Searching wider than that assertion cannot find signal it excludes; it can only add noise.
    // Search 3x the asserted half-width: wide enough to see the distribution's shape and its
    // shoulders, narrow enough that the sample is mostly real.
    const double cfg_full = ms1 && getDoubleOption_("mz_extraction_window_ms1") > 0.0
                              ? getDoubleOption_("mz_extraction_window_ms1")
                              : getDoubleOption_("mz_extraction_window");
    const double boot = std::min(getDoubleOption_("mz_calib_bootstrap_ppm"),
                                 cfg_full > 0.0 ? 1.5 * cfg_full : 50.0);   // 3x the half-width
    const int    min_n = getIntOption_("mz_calib_min_anchors");
    const double kmul = getDoubleOption_("mz_calib_sigma_multiple");
    const double floor_ppm = getDoubleOption_("mz_calib_floor_ppm");
    const double min_peak = getDoubleOption_("mz_calib_min_peakedness");
    if (boot <= 0.0) { return -1.0; }

    // Index the MS2 maps by isolation window so a fragment is only ever sought in the window that
    // actually isolated its precursor; searching all maps would manufacture interference.
    std::vector<const OpenSwath::SwathMap*> maps;
    for (const auto& sm : swath_maps)
    {
      if (!sm.sptr) { continue; }
      if (ms1 == sm.ms1) { maps.push_back(&sm); }
    }
    if (maps.empty()) { return -1.0; }

    // Signed relative errors, in ppm.
    std::vector<double> err;
    err.reserve(4096);

    const auto& transitions = anchors.getTransitions();
    // Group transitions by precursor so we can cap the work per precursor and keep the sample
    // spread across the library rather than dominated by a few fragment-rich assays.
    std::unordered_map<std::string, std::vector<const OpenSwath::LightTransition*>> by_prec;
    for (const auto& t : transitions) { by_prec[t.getPeptideRef()].push_back(&t); }

    const std::size_t max_prec = static_cast<std::size_t>(std::max(1, getIntOption_("mz_calib_max_anchors")));
    std::size_t used = 0;
    for (const auto& c : anchors.getCompounds())
    {
      if (used >= max_prec) { break; }
      auto it = by_prec.find(c.id);
      if (it == by_prec.end()) { continue; }

      // Where in RT to look. rt_trafo maps library RT onto run RT; without it we would be searching
      // the whole gradient and the sample would be dominated by interference.
      const double rt = rt_trafo.getDataPoints().empty() ? c.rt : rt_trafo.apply(c.rt);

      // Precursor m/z lives on the TRANSITION, not on LightCompound.
      const double prec_mz = it->second.empty() ? 0.0 : it->second.front()->precursor_mz;
      if (prec_mz <= 0.0) { continue; }

      // Pick the map whose isolation window contains this precursor (MS1: the single survey map).
      const OpenSwath::SwathMap* map = nullptr;
      for (const auto* m : maps)
      {
        if (ms1) { map = m; break; }
        if (prec_mz >= m->lower && prec_mz <= m->upper) { map = m; break; }
      }
      if (!map || !map->sptr) { continue; }

      const std::size_t ns = map->sptr->getNrSpectra();
      if (ns == 0) { continue; }
      // Nearest spectrum in time -- binary search would need sorted RTs; a coarse scan is fine at
      // these anchor counts and avoids assuming ordering.
      std::size_t best_i = 0; double best_dt = std::numeric_limits<double>::max();
      for (std::size_t i = 0; i < ns; ++i)
      {
        const double d = std::abs(map->sptr->getSpectrumMetaById(static_cast<int>(i)).RT - rt);
        if (d < best_dt) { best_dt = d; best_i = i; }
      }
      OpenSwath::SpectrumPtr sp = map->sptr->getSpectrumById(static_cast<int>(best_i));
      if (!sp || sp->getMZArray()->data.empty()) { continue; }
      const auto& mz = sp->getMZArray()->data;
      const auto& in = sp->getIntensityArray()->data;

      int taken = 0;
      for (const auto* t : it->second)
      {
        // MS1 has ONE observable per precursor: prec_mz. The loop body below is keyed on the
        // transition, so at MS1 level every iteration computed the same theo, matched the same
        // peak and pushed the same residual again -- 400 real measurements reported as 2,400,
        // clearing the min_anchors gate on duplicates and distorting the peakedness histogram
        // (every unique value entered it six times).
        if (ms1 && taken >= 1) { break; }
        if (taken >= 6) { break; }                       // cap fragments per precursor
        const double theo = ms1 ? prec_mz : t->product_mz;
        if (theo <= 0.0) { continue; }
        const double tol = theo * boot * 1e-6;
        // MOST INTENSE within the window, not nearest. On this instrument class the interferent
        // population is numerous but individually weak (most centroids carry very few ions), so
        // "nearest" preferentially selects single-ion noise sitting close to the query centre,
        // while "most intense" selects the real peak whenever it is above that floor.
        double best_mz = -1.0, best_int = -1.0;
        for (std::size_t j = 0; j < mz.size(); ++j)
        {
          if (mz[j] < theo - tol) { continue; }
          if (mz[j] > theo + tol) { break; }             // m/z arrays are ascending
          if (j < in.size() && in[j] > best_int) { best_int = in[j]; best_mz = mz[j]; }
        }
        if (best_mz > 0.0)
        {
          err.push_back((best_mz - theo) / theo * 1e6);
          ++taken;
        }
      }
      if (taken > 0) { ++used; }
    }

    const char* lvl = ms1 ? "MS1" : "MS2";
    if (static_cast<int>(err.size()) < min_n)
    {
      OPENMS_LOG_INFO << "OpenDIAlyzer[mzcal] " << lvl << ": only " << err.size()
                      << " anchor errors (need " << min_n << ") -- not inferred." << std::endl;
      return -1.0;
    }

    std::sort(err.begin(), err.end());
    const double offset = halfSampleMode_(err);
    // Scale is taken about the MODE, not the mean: a contaminated sample has a long tail that would
    // drag a mean-centred scale outward.
    std::vector<double> dev;
    dev.reserve(err.size());
    for (double e : err) { dev.push_back(std::abs(e - offset)); }
    std::sort(dev.begin(), dev.end());
    const double sigma = localScaleAboutMode_(err, offset, boot);

    const double peak = peakednessRatio_(dev, boot);
    if (peak < min_peak)
    {
      OPENMS_LOG_WARN << "OpenDIAlyzer[mzcal] " << lvl << ": residuals are FLAT (peakedness "
                      << peak << " < " << min_peak << ") over " << err.size() << " anchors. That is "
                      << "what a mostly-noise anchor set looks like -- a mode and a width can still "
                      << "be computed from it and would be meaningless. Not inferred; keeping the "
                      << "configured window." << std::endl;
      return -1.0;
    }
    if (!(sigma > 0.0) || !std::isfinite(sigma) || !std::isfinite(offset))
    {
      OPENMS_LOG_INFO << "OpenDIAlyzer[mzcal] " << lvl << ": degenerate scale -- not inferred." << std::endl;
      return -1.0;
    }

    // Window covers the RANDOM error only; the systematic offset belongs to SwathMapMassCorrection,
    // which already fits MZTrafoModel and rewrites the spectra. Adding |offset| back in here would
    // re-pay for a bias that is corrected elsewhere.
    const double win = std::max(kmul * sigma, floor_ppm);
    OPENMS_LOG_INFO << "OpenDIAlyzer[mzcal] " << lvl << ": " << err.size() << " anchor errors, "
                    << "offset " << offset << " ppm, sigma " << sigma << " ppm, peakedness "
                    << peak << " -> window " << win << " ppm (" << kmul << " sigma, floor "
                    << floor_ppm << ")." << std::endl;
    return win;
  }

  /// HOOK (empty): infer the pass-1 RT extraction window in seconds.
  /// Intended implementation: hold out a fraction of the anchors from the transform fit, measure
  /// |observed - predicted| on those, and return ~2x their p95 -- an out-of-sample estimate, unlike
  /// `estimateWindow()` which is fitted to the survivors it is derived from.
  /// @param n_anchor_pairs anchors that survived outlier removal (fit quality context).
  /// @return full-width window in seconds, or <= 0 to keep the configured/estimated value.
  double inferRtWindowSeconds_(const std::vector<OpenSwath::SwathMap>& /*swath_maps*/,
                               const OpenSwath::LightTargetedExperiment& /*anchors*/,
                               const TransformationDescription& /*rt_trafo*/,
                               std::size_t /*n_anchor_pairs*/) const
  {
    return -1.0;   // not inferred
  }

  bool calibrateCiRT_(std::vector<OpenSwath::SwathMap>& swath_maps,
                      OpenSwath::LightTargetedExperiment& transition_exp,
                      bool pasef, double rt_window, const std::string& in_file,
                      TransformationDescription& rt_trafo, double& estimated_rt_window)
  {
    try
    {
      CalibrationWorkflow cw;
      cw.setLogType(log_type_);
      {
        Param cwp = cw.getParameters();
        cwp.setValue("qc:min_rsq", getDoubleOption_("calibration_min_rsq"));
        cwp.setValue("qc:min_coverage", getDoubleOption_("calibration_min_coverage"));
        // These were previously left at OpenMS defaults, which is how the calibration came to cost
        // 32.7 CPU-hours for 136 anchors. See the option comments for the measurement.
        cwp.setValue("auto_irt:irt_bins_nonlinear", getIntOption_("calibration_nonlinear_bins"));
        cwp.setValue("auto_irt:irt_peptides_per_bin_nonlinear", getIntOption_("calibration_nonlinear_per_bin"));
        cwp.setValue("auto_irt:irt_nonlinear_rt_extraction_window", getDoubleOption_("calibration_nonlinear_rt_window"));
        cw.setParameters(cwp);
      }

      // iRT anchors sampled once from the target library itself (auto_irt). Priority
      // sequences (CiRT/iRT kits) are optional: the sampler falls back to binning the
      // library across the iRT range when none are supplied.
      std::vector<std::string> priority_peptides = loadPriorityIrtSequences_();
      CalibrationWorkflow::IrtExperiments irt =
        cw.prepareIrtExperiments(IrtStrategy::SAMPLE_ONCE, transition_exp, priority_peptides);
      if (!irt.is_prepared || irt.linear_irt.getTransitions().empty())
      {
        OPENMS_LOG_WARN << "OpenDIAlyzer: could not prepare iRT anchors for calibration." << std::endl;
        return false;
      }
      OPENMS_LOG_INFO << "OpenDIAlyzer: CiRT calibration using "
                      << irt.linear_irt.getCompounds().size() << " linear iRT anchor compounds ("
                      << irt.linear_irt.getTransitions().size() << " transitions); nonlinear anchors: "
                      << irt.nonlinear_irt.getCompounds().size() << "." << std::endl;

      ChromExtractParams cp = makeChromParams_(rt_window);
      ChromExtractParams cp_ms1 = makeMs1ChromParams_(cp);
      ChromExtractParams cp_irt = cp;
      cp_irt.rt_extraction_window = -1;                 // iRT anchors: search the whole run
      // The irt_* options exist because the calibration extraction is NOT the analyte
      // extraction: anchors are few and must be found before any m/z correction, so TOPP uses
      // a deliberately wider window here. Registering them without applying them would make
      // the help text a lie and silently reuse the analyte window.
      cp_irt.mz_extraction_window = getDoubleOption_("irt_mz_extraction_window");
      cp_irt.ppm = (getStringOption_("irt_mz_extraction_window_unit") == "ppm");
      cp_irt.im_extraction_window = getDoubleOption_("irt_im_extraction_window");

      Param ff = makeFeatureFinderParam_();
      ff.setValue("Scores:use_ion_mobility_scores", pasef ? "true" : "false");
      const Param irt_detection = makeIrtDetectionParam_();
      // These subsystems need their own defaults; passing empty Params makes
      // performCalibration throw ("the element 'mz_extraction_window' could not be found").
      const Param calibration_param = SwathMapMassCorrection().getDefaults();
      Param mrm_mapping_param = MRMMapping().getDefaults();       // iRT variant, as in the TOPP tool
      mrm_mapping_param.setValue("precursor_tolerance", irt_detection.getValue("irt_precursor_tolerance"));
      mrm_mapping_param.setValue("product_tolerance", irt_detection.getValue("irt_product_tolerance"));
      mrm_mapping_param.setValue("map_multiple_assays", irt_detection.getValue("map_multiple_assays"));
      mrm_mapping_param.setValue("error_on_unmapped", irt_detection.getValue("error_on_unmapped"));

      // load_into_memory was hardcoded true here, which made the whole -readOptions setting
      // cosmetic for memory: calibration forced full residency regardless, so the process peak
      // RSS was set before extraction even began. Resolve it the SAME way every other site does.
      std::string eff_ro; bool cal_load_into_memory = false;
      resolveReadOptions_(in_file, eff_ro, cal_load_into_memory);
      CalibrationWorkflow::CalibrationResult res =
        cw.performCalibration(swath_maps, transition_exp, cp, cp_ms1, irt, ff, cp_irt,
                              irt_detection, calibration_param, mrm_mapping_param,
                              pasef, cal_load_into_memory);
      rt_trafo = res.rt_trafo;
      estimated_rt_window = res.estimated_rt_window;
      if (rt_trafo.getDataPoints().empty())
      {
        OPENMS_LOG_WARN << "OpenDIAlyzer: calibration produced an empty transform." << std::endl;
        return false;
      }

      // ---- data-driven inference hooks (currently empty; see their declarations) -------------
      // Wired here so that implementing either one needs no plumbing changes. Both no-op today and
      // say so once, rather than pretending the calibration is doing something it is not.
      {
        // SAME DIRECTION RULE as calibrateMassFromPass_. res.rt_trafo is the calibration's output in
        // OpenSWATH's convention -- RUN rt -> iRT (see the comment at the `native_trafo` declaration)
        // -- while irt.nonlinear_irt's compound.rt values are library iRT. inferMassAccuracyPpm_ needs
        // to go the other way: from a library rt to the RUN rt at which to look for the peak.
        //
        // Passing rt_trafo unmodified fed a run->iRT map a value already in iRT. The result is not a
        // retention time, so the "nearest spectrum in time" search landed on an arbitrary spectrum and
        // the peak search sampled interference -- the FLAT residual distribution this site has been
        // reporting (peakedness 1.17 over 1066 MS2 anchors).
        TransformationDescription mz_trafo = rt_trafo;
        mz_trafo.invert();                                  // run -> iRT  becomes  iRT -> run
        const double mz_ms2 = inferMassAccuracyPpm_(swath_maps, irt.nonlinear_irt, mz_trafo, false);
        const double mz_ms1 = inferMassAccuracyPpm_(swath_maps, irt.nonlinear_irt, mz_trafo, true);
        if (mz_ms2 > 0.0) { inferred_mz_window_ms2_ = mz_ms2; }
        if (mz_ms1 > 0.0) { inferred_mz_window_ms1_ = mz_ms1; }
        if (mz_ms2 > 0.0 || mz_ms1 > 0.0)
        {
          OPENMS_LOG_INFO << "OpenDIAlyzer: inferred m/z windows -- MS2 " << mz_ms2 << " ppm, MS1 "
                          << mz_ms1 << " ppm (overriding the configured values)." << std::endl;
        }
        else
        {
          OPENMS_LOG_INFO << "OpenDIAlyzer: mass-accuracy inference not implemented -- using the "
                          << "configured -mz_extraction_window (" << getDoubleOption_("mz_extraction_window")
                          << " ppm). NOTE this instrument's measured MS2 accuracy is ~1.6 ppm; a "
                          << "window far wider than the real error admits interference and destroys "
                          << "target/decoy separation." << std::endl;
        }

        const double rt_win = inferRtWindowSeconds_(swath_maps, irt.nonlinear_irt, rt_trafo,
                                                    rt_trafo.getDataPoints().size());
        if (rt_win > 0.0)
        {
          OPENMS_LOG_INFO << "OpenDIAlyzer: inferred RT window " << rt_win << " s (replacing the "
                          << "calibration estimate " << estimated_rt_window << " s)." << std::endl;
          estimated_rt_window = rt_win;
        }
      }
      // Report the YIELD, not just the count. 136 pairs looks fine in isolation; "136 from 67,429
      // candidates" is obviously broken, and the difference between those two log lines is 32
      // wasted CPU-hours nobody noticed. A low yield means the anchors' true peaks are not inside
      // the window they were searched in -- widen -calibration_nonlinear_rt_window, do not sample more.
      const std::size_t n_cand = irt.nonlinear_irt.getCompounds().empty()
                                   ? irt.linear_irt.getCompounds().size()
                                   : irt.nonlinear_irt.getCompounds().size();
      const std::size_t n_pairs = rt_trafo.getDataPoints().size();
      const double yield = n_cand > 0 ? 100.0 * static_cast<double>(n_pairs) / static_cast<double>(n_cand) : 0.0;
      OPENMS_LOG_INFO << "OpenDIAlyzer: CiRT calibration OK -- " << n_pairs << " anchor pairs from "
                      << n_cand << " candidates (" << yield << "% yield), estimated RT window = "
                      << estimated_rt_window << " s." << std::endl;
      // A healthy yield is tens of percent, not single digits, so warn well above 5%.
      if (yield < 15.0)
      {
        // n_pairs is the count AFTER outlier removal, so a low yield has TWO causes that this
        // number cannot separate, and they demand opposite fixes:
        //   (a) anchors never found  -> peaks outside the search window -> WIDEN the window
        //   (b) wrong peaks found, then discarded by outlier removal -> NARROW it
        // Prescribing "widen" unconditionally is actively harmful under (b), because a wrong anchor
        // corrupts the transform whereas a missing one merely weakens it. Outlier removal may
        // discard down to qc:min_coverage of the matched pairs, which bounds how many were matched
        // before removal -- report that bound rather than guessing which mode we are in.
        const double mc = getDoubleOption_("calibration_min_coverage");
        const double matched_upper = mc > 0.0 ? static_cast<double>(n_pairs) / mc : 0.0;
        calib_yield_pct_ = yield;
        OPENMS_LOG_WARN << "OpenDIAlyzer: calibration anchor yield is only " << yield << "% ("
                        << n_pairs << "/" << n_cand << "). At most " << static_cast<std::size_t>(matched_upper)
                        << " candidates matched before outlier removal (bound from qc:min_coverage="
                        << mc << "). If that bound is also small the anchors were NOT FOUND -- widen "
                        << "-calibration_nonlinear_rt_window. If it is large, wrong peaks were found "
                        << "and discarded -- NARROW it instead. The estimated window ("
                        << estimated_rt_window << " s) is fitted to the survivors and is optimistic "
                        << "for unseen precursors either way." << std::endl;
      }
      // Floor-hugging on min_coverage is the signature of a transform fitted to a self-selected
      // clique: the outlier remover wanted to discard more and only the QC floor stopped it.
      if (getDoubleOption_("calibration_min_coverage") <= 0.35 && getDoubleOption_("calibration_min_rsq") <= 0.75)
      {
        OPENMS_LOG_WARN << "OpenDIAlyzer: BOTH calibration QC gates are relaxed (min_coverage="
                        << getDoubleOption_("calibration_min_coverage") << ", min_rsq="
                        << getDoubleOption_("calibration_min_rsq") << "). Relaxing one for a "
                        << "predicted library is defensible; relaxing both lets outlier removal "
                        << "discard most matched anchors to manufacture the R^2 it is checked against."
                        << std::endl;
      }
      return true;
    }
    catch (const Exception::BaseException& e)
    {
      OPENMS_LOG_WARN << "OpenDIAlyzer: CiRT calibration failed (" << e.what()
                      << "); falling back to the linear bootstrap." << std::endl;
      return false;
    }
  }

  // "Calibration:RTNormalization" defaults, replicated VERBATIM from TOPP OpenSwathWorkflow
  // (MRMRTNormalizer is not a DefaultParamHandler, so the tool builds this by hand).
  // Keep these exactly at the TOPP defaults: deviating (lowess + estimateBestPeptides=true)
  // made MRMRTNormalizer's quality filter discard so many anchors that calibration aborted
  // with "insufficient RT coverage after outlier removal", while stock OpenSwathWorkflow
  // calibrated the same library fine.
  static Param makeIrtDetectionParam_()
  {
    Param p;
    p.setValue("alignmentMethod", "linear");
    p.setValidStrings("alignmentMethod", {"linear", "interpolated", "lowess", "b_spline"});
    p.setValue("lowess:auto_span", "true");
    p.setValidStrings("lowess:auto_span", {"true", "false"});
    p.setValue("lowess:span", 0.05);
    p.setValue("lowess:auto_span_min", 0.15);
    p.setValue("lowess:auto_span_max", 0.80);
    p.setValue("lowess:auto_span_grid", "0.005,0.01,0.05,0.15,0.25,0.30,0.50,0.70,0.90");
    p.setValue("b_spline:num_nodes", 5);
    p.setValue("useIterativeChauvenet", "false");
    p.setValidStrings("useIterativeChauvenet", {"true", "false"});
    p.setValue("RANSACMaxIterations", 1000);
    p.setValue("RANSACMaxPercentRTThreshold", 3);
    p.setValue("RANSACSamplingSize", 10);
    p.setValue("estimateBestPeptides", "false");
    p.setValidStrings("estimateBestPeptides", {"true", "false"});
    p.setValue("InitialQualityCutoff", 0.5);
    p.setValue("OverallQualityCutoff", 5.5);
    p.setValue("NrRTBins", 10);
    p.setValue("MinPeptidesPerBin", 1);
    p.setValue("MinBinsFilled", 8);
    p.setValue("precursor_tolerance", 0.9);
    p.setValue("product_tolerance", 1.2);
    p.setValue("irt_precursor_tolerance", 1.5);
    p.setValue("irt_product_tolerance", 1.5);
    p.setValue("map_multiple_assays", "false");
    p.setValidStrings("map_multiple_assays", {"true", "false"});
    p.setValue("error_on_unmapped", "false");
    p.setValidStrings("error_on_unmapped", {"true", "false"});
    return p;
  }

  // Sequences of the built-in iRT/CiRT kits, used to prioritise anchor sampling.
  static std::vector<std::string> loadPriorityIrtSequences_()
  {
    std::vector<std::string> seqs;
    for (const char* fn : {"CHEMISTRY/cirtkit.tsv", "CHEMISTRY/irtkit.tsv"})
    {
      std::string path;
      try { path = File::find(fn); } catch (...) { continue; }
      std::ifstream in(path);
      if (!in) { continue; }
      std::string line;
      bool header = true;
      while (std::getline(in, line))
      {
        if (header) { header = false; continue; }
        std::vector<std::string> f;
        size_t start = 0, tab;
        while ((tab = line.find('\t', start)) != std::string::npos) { f.push_back(line.substr(start, tab - start)); start = tab + 1; }
        f.push_back(line.substr(start));
        if (f.size() > 6 && !f[6].empty()) { seqs.push_back(f[6]); }
      }
    }
    OPENMS_LOG_INFO << "OpenDIAlyzer: " << seqs.size() << " priority iRT/CiRT sequences loaded." << std::endl;
    return seqs;
  }

  // Chromatogram-extraction parameters (shared by calibration and the extraction passes).
  ChromExtractParams makeChromParams_(double rt_window, double mz_window_override = -1.0) const
  {
    ChromExtractParams cp;
    cp.min_upper_edge_dist = getDoubleOption_("min_upper_edge_dist");
    // An inferred window (from the mass-accuracy hook, once implemented) wins over the configured
    // one; it is expressed in ppm by contract, so the unit flag is forced to match.
    // mz_window_override lets the PREFILTER screen at a different width than extraction scores at
    // -- see -prefilter_mz_extraction_window.
    cp.mz_extraction_window = (mz_window_override > 0.0)     ? mz_window_override
                            : (inferred_mz_window_ms2_ > 0.0) ? inferred_mz_window_ms2_
                                                              : getDoubleOption_("mz_extraction_window");
    cp.ppm = (inferred_mz_window_ms2_ > 0.0) || (getStringOption_("mz_extraction_window_unit") == "ppm");
    cp.rt_extraction_window = rt_window;                       // < 0 -> whole RT range
    // Requesting IM extraction on data WITHOUT an ion-mobility array is a hard abort in
    // ChromatogramExtractorAlgorithm.cpp:307 ("Requested ion mobility extraction but no ion
    // mobility array found") -- SIGABRT, zero features, on every TripleTOF/Orbitrap run.
    // The window must therefore follow the DATA, not a default: -pasef auto already detects
    // whether the SWATH maps carry IM, so gate on that. An explicit positive value from the
    // user is still honoured, but only when the data can support it.
    const double imw = getDoubleOption_("ion_mobility_window");
    cp.im_extraction_window = (pasef_ && imw > 0.0) ? imw : -1.0;
    cp.extraction_function = getStringOption_("extraction_function");
    cp.extra_rt_extract = getDoubleOption_("extra_rt_extraction_window");
    return cp;
  }

  // MS1 extraction params. Deliberately NOT `cp_ms1 = cp` with the window overwritten: that
  // inherits the MS2 m/z UNIT, so selecting Th for MS2 would reinterpret a 30 ppm MS1 window
  // as 30 Th. TOPP keeps the units separate; so do we.
  ChromExtractParams makeMs1ChromParams_(const ChromExtractParams& cp) const
  {
    ChromExtractParams cp_ms1 = cp;
    cp_ms1.mz_extraction_window = (inferred_mz_window_ms1_ > 0.0)
                                    ? inferred_mz_window_ms1_
                                    : getDoubleOption_("mz_extraction_window_ms1");
    cp_ms1.ppm = (inferred_mz_window_ms1_ > 0.0)
                 || (getStringOption_("mz_extraction_window_ms1_unit") == "ppm");
    cp_ms1.im_extraction_window = getDoubleOption_("im_extraction_window_ms1");
    return cp_ms1;
  }

  // TOPP semantics (OpenSwathWorkflow.cpp:1020-1030): 'cacheWorkingInMemory' and
  // 'workingInMemory' are TOPP-LEVEL names. SwathFile itself only accepts normal|cache|split
  // and THROWS on anything else -- so passing 'cacheWorkingInMemory' straight through, as this
  // tool did, made every mzML/mzXML run abort with "Unknown or unsupported option". Map first.
  void resolveReadOptions_(const std::string& in_file, std::string& swathfile_opt,
                           bool& load_into_memory) const
  {
    std::string opt = getStringOption_("readOptions");
    if (opt == "auto")
    {
      // Established per-format behaviour: Bruker .d streams natively; everything else was
      // meant to work in memory (which is also what the wave scheduler needs).
      opt = (FileHandler::getTypeByFileName(in_file) == FileTypes::BRUKER_TDF)
                ? "normal" : "cacheWorkingInMemory";
    }
    load_into_memory = false;
    if (opt == "cacheWorkingInMemory") { swathfile_opt = "cache";  load_into_memory = true; }
    else if (opt == "workingInMemory") { swathfile_opt = "normal"; load_into_memory = true; }
    else                               { swathfile_opt = opt; }
  }

  // Does this decoy id name a target that survived the evidence screen? Pure and static so the
  // selftest can exercise it without a run: this is the pairing rule that was wrong, and getting
  // it wrong disables the prefilter silently rather than loudly.
  /// Selftest-only since the decoy pass was inverted (see "pass 2" below). Kept because it is the
  /// written-down statement of the pairing rule: decoy id == decoy_tag + target id, tag stripped
  /// exactly ONCE and only as a prefix. The inverted loop cannot violate those cases by
  /// construction -- it only ever forms `decoy_tag + <a kept target id>` -- so the selftest now
  /// guards the RULE rather than the code path.
  static bool decoyPartnerKept_(const std::string& decoy_id, const std::string& decoy_tag,
                                const std::unordered_set<std::string>& kept_targets)
  {
    return decoy_id.size() > decoy_tag.size()
           && decoy_id.compare(0, decoy_tag.size(), decoy_tag) == 0
           && kept_targets.count(decoy_id.substr(decoy_tag.size())) > 0;
  }

  // Screen the library against evidence in the run, then keep the surviving TARGETS *and their
  // paired DECOYS*.
  //
  // Result::filtered_targets is documented target-only ("for auto-iRT sampling"), so using it
  // directly as the search library would drop every decoy and silently destroy FDR -- the one
  // thing in this tool that has actually been validated against pyprophet. So the evidence is
  // used only to choose which peptides survive; the experiment is rebuilt keeping decoys whose
  // target partner survived, preserving the ~1:1 target:decoy ratio the FDR estimate needs.
  //
  // Why this exists: without it, extraction runs over EVERY library precursor. On the
  // whole-proteome library that is 7,149,966 candidates and, with permissive peak-picking,
  // 27,344,163 features / 1753 GB / 2.9 h -- against OpenSwathWorkflow's 34,165 retained.
  bool prefilterLibrary_(const std::vector<OpenSwath::SwathMap>& swath_maps,
                         OpenSwath::LightTargetedExperiment& transition_exp,
                         const ChromExtractParams& cp, const ChromExtractParams& cp_ms1) const
  {
    if (getStringOption_("prefilter") != "true") { return false; }
    const std::size_t n_before = transition_exp.getCompounds().size();

    TransitionListEvidenceFilter filt;
    // The filter was previously constructed with PURE DEFAULTS, and its default
    // evidence_sources="hybrid" means supported_ms1 OR supported_ms2 -- where supported_ms1 is
    // just `ms1_hit_count > 0`: ONE peak anywhere in the run within the m/z tolerance, with no
    // intensity floor, no RT coherence and no SWATH-window consistency. At proteome scale that
    // saturates. Measured on 12_80.mzML (1083 MS1 spectra, 30 ppm, 3,603,425 targets):
    //     MS1: 3,150,365 (87.4%)   MS2: 514,349 (14.3%)   both: 474,593   hybrid(OR): 3,190,121
    // i.e. the MS1 arm rubber-stamps almost everything and the OR hides the selective MS2 arm.
    // MS2 evidence (>= ms2_min_fragment_hits distinct library fragments in ONE spectrum) is the
    // criterion that actually discriminates, and "both" (474,593) is barely tighter than MS2
    // alone -- so requiring MS2 captures essentially all the available selectivity without
    // needing the AND that OpenMS does not expose.
    Param pf = filt.getDefaults();
    // Only ms2_top_transitions_per_precursor fragments are indexed per precursor, so demanding more
    // distinct hits than that can never be satisfied -- it would silently support nothing and (via
    // the guards below) fall back to the full library. Fail loudly instead of mysteriously.
    const int min_frag = getIntOption_("prefilter_min_fragments");
    const int top_tr   = static_cast<int>(pf.getValue("ms2_top_transitions_per_precursor"));
    if (min_frag > top_tr)
    {
      OPENMS_LOG_ERROR << "OpenDIAlyzer: -prefilter_min_fragments " << min_frag << " exceeds the "
                       << top_tr << " transitions indexed per precursor; no precursor could ever "
                       << "be supported." << std::endl;
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        "prefilter_min_fragments > ms2_top_transitions_per_precursor");
    }
    pf.setValue("evidence_sources", getStringOption_("prefilter_evidence"));
    pf.setValue("ms2_min_fragment_hits", min_frag);
    pf.setValue("ms1_top_peaks_per_spectrum", getIntOption_("prefilter_top_peaks"));
    pf.setValue("ms2_top_peaks_per_spectrum", getIntOption_("prefilter_top_peaks"));
    pf.setValue("ms2_min_qualifying_spectra", getIntOption_("prefilter_min_spectra"));
    // Leave "enabled" false: that flag only arms an internal min_supported_precursors THROW.
    // We prefer our own soft guards below (keep-full-library) over aborting the run.
    filt.setParameters(pf);

    TransitionListEvidenceFilter::Result res;
    { PhaseTimer pt("prefilter/scan_targets");
      res = filt.filter(swath_maps, transition_exp, cp_ms1, cp, pasef_, getIntOption_("threads")); }

    // Decoy status lives on the TRANSITION (getDecoy()), not on LightCompound, so derive it per
    // peptide ref first -- the evidence split below needs it.
    std::unordered_set<std::string> decoy_refs;
    { PhaseTimer pt("prefilter/build_sets");
      for (const auto& t : transition_exp.getTransitions())
      {
        if (t.getDecoy()) { decoy_refs.insert(t.getPeptideRef()); }
      }
    }

    // Peptide refs with evidence. Two things this must NOT do, both previously wrong:
    //
    // 1. It must honour evidence_sources. Hardcoding `supported_ms1 || supported_ms2` re-implements
    //    the filter's "hybrid" rule no matter what was configured, so -prefilter_evidence would
    //    change only the LOGGED count while keeping exactly the same precursors -- a setting that
    //    silently does nothing is worse than no setting.
    // 2. It must key on compound_id only. compound_id is always populated (the filter sets it from
    //    getPeptideRef(), TransitionListEvidenceFilter.cpp:380/390), so also inserting the bare
    //    sequence adds nothing except collisions: evidence for PEPTIDEK_2 would then keep
    //    PEPTIDEK_3 -- a different precursor, different charge, zero evidence -- and pass 2 would
    //    attach its decoy too. Symmetric, so not an FDR bug, but it dilutes exactly the
    //    selectivity this filter exists to provide.
    const std::string ev = getStringOption_("prefilter_evidence");
    const std::string dtag = getStringOption_("decoy_tag");
    std::unordered_set<std::string> keep;        // supported TARGETS, by compound id
    for (const auto& e : res.evidence)
    {
      const bool supported = (ev == "ms1") ? e.supported_ms1
                           : (ev == "ms2") ? e.supported_ms2
                                           : (e.supported_ms1 || e.supported_ms2);
      if (supported && !e.compound_id.empty()) { keep.insert(e.compound_id); }
    }

    // Check the pairing PRECONDITION before spending the two extraction passes below. Pairing is by
    // id ("DECOY_" + target id); a library whose decoy ids do not carry the tag cannot be paired at
    // all, and the ratio guard at the end would then abort after ~6 min of scanning with only a
    // number to go on. This says what is actually wrong, in seconds. The observed cause: reading a
    // .pqp with legacy_traml_id=false, which renames every compound to its row number ("0", "1", ...).
    if (!decoy_refs.empty())
    {
      const std::size_t tagged = std::count_if(decoy_refs.begin(), decoy_refs.end(),
        [&](const std::string& s) { return s.compare(0, dtag.size(), dtag) == 0; });
      if (tagged == 0)
      {
        OPENMS_LOG_ERROR << "OpenDIAlyzer: none of the " << decoy_refs.size() << " decoy precursors "
                         << "has an id starting with -decoy_tag ('" << dtag << "') -- e.g. '"
                         << *decoy_refs.begin() << "'. Decoys are paired to their targets BY ID, so "
                         << "no pair can be formed and the FDR would have no null distribution. "
                         << "Numeric ids mean the library lost its original names in a format "
                         << "round-trip; rebuild it from the TSV/TraML source." << std::endl;
        throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                          "library decoy ids do not carry the decoy tag");
      }
    }

    // ---- score the DECOYS with the identical criterion -------------------------------------
    //
    // Selection must be label-symmetric: a false target retained BECAUSE it had coincidental
    // evidence, paired with a decoy that faced no such test, makes retained decoys score
    // systematically lower than retained false targets -- an anti-conservative FDR. So both classes
    // are scored and a pair is retained when EITHER member passes.
    //
    // The filter refuses decoys twice -- on getDecoy() AND on a hardcoded "DECOY"-prefix id test --
    // so presenting them for scoring needs the flag cleared AND the ids aliased. Measured cost:
    // 43.5 s to build this ~21M-transition view plus 63.5 s for the second spectrum pass. That is
    // the price of not modifying OpenMS, and an invalid FDR costs more.
    std::unordered_map<std::string, std::string> alias_of, id_of_alias;
    OpenSwath::LightTargetedExperiment dview;
    { PhaseTimer pt("prefilter/build_decoy_view");
      for (const auto& c : transition_exp.getCompounds())
      {
        if (!decoy_refs.count(c.id)) { continue; }
        const std::string a = "PFD" + std::to_string(alias_of.size());
        alias_of[c.id] = a;
        id_of_alias[a] = c.id;
        auto cc = c;
        cc.id = a;
        dview.compounds.push_back(cc);
      }
      for (const auto& t : transition_exp.getTransitions())
      {
        const auto it = alias_of.find(t.getPeptideRef());
        if (it == alias_of.end()) { continue; }
        auto tt = t;
        tt.peptide_ref = it->second;
        tt.setDecoy(false);
        dview.transitions.push_back(tt);
      }
    }

    std::unordered_set<std::string> keep_decoy;   // keyed by the decoy's TARGET id
    if (!dview.compounds.empty())
    {
      TransitionListEvidenceFilter::Result res_d;
      { PhaseTimer pt("prefilter/scan_decoys");
        res_d = filt.filter(swath_maps, dview, cp_ms1, cp, pasef_, getIntOption_("threads")); }
      for (const auto& e : res_d.evidence)
      {
        const bool sup = (ev == "ms1") ? e.supported_ms1
                       : (ev == "ms2") ? e.supported_ms2
                                       : (e.supported_ms1 || e.supported_ms2);
        if (!sup) { continue; }
        const auto it = id_of_alias.find(e.compound_id);
        // Stored under the TARGET id (tag stripped once here) so the pair-union test below is a
        // plain lookup instead of building `dtag + c.id` for all 7.1M compounds.
        if (it != id_of_alias.end() && it->second.compare(0, dtag.size(), dtag) == 0)
        {
          keep_decoy.insert(it->second.substr(dtag.size()));
        }
      }
    }
    OPENMS_LOG_INFO << "OpenDIAlyzer[prefilter] evidence: " << keep.size() << " targets, "
                    << keep_decoy.size() << " decoys pass the same criterion (label-symmetric "
                    << "selection; a large asymmetry here means the criterion is not discriminating)."
                    << std::endl;

    // Fail LOUDLY, not open. The old behaviour here was to fall back to the unfiltered library,
    // which is the single most expensive thing this program can do AND -- measured on astral.mzML,
    // same library -- yields 27,344,193 features whose best q-value is 0.508, i.e. zero
    // identifications at any usable FDR. "Silently do the useless 1753 GB thing" is not a safe
    // default for a mis-set window. Searching the full library must be an explicit -prefilter false.
    if (keep.empty() && keep_decoy.empty())
    {
      OPENMS_LOG_ERROR << "OpenDIAlyzer: prefilter found NO supported precursors in either class. "
                       << "Check the m/z / RT / IM windows and -prefilter_evidence. Pass "
                       << "-prefilter false to search the unfiltered library deliberately."
                       << std::endl;
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        "prefilter supported no precursors");
    }

    // A decoy survives iff its target partner does -- and the link is the ID, not the sequence.
    // MRMDecoy derives every decoy id as decoy_tag + target id ("computed deterministically from
    // the target id", MRMDecoy.cpp:642; also :557/:694) while SHUFFLING the sequence. So matching a
    // decoy to its target BY SEQUENCE can only ever fire on an accidental shuffle collision.
    // Measured on the real 7,149,966-precursor library (12_80.mzML, 2026-07-30): it kept
    // 3,546,918 targets and 23 decoys. That is the dangerous failure mode -- 23 is not 0, so the
    // then-current `n_dec == 0` guard did NOT fire, and the run proceeded with a 154,000:1
    // target:decoy ratio, i.e. an FDR estimate with no null distribution left to calibrate
    // against. That guard is now a ratio band (see below), so partial pairing failure is caught too.
    PhaseTimer t_rebuild("prefilter/pair_and_rebuild");
    std::unordered_set<std::string> kept_ids;
    for (const auto& c : transition_exp.getCompounds())        // pass 1: targets, on evidence
    {
      if (decoy_refs.count(c.id)) { continue; }                // decoys handled in pass 2
      // Pair-union: the pair is retained when EITHER member has evidence. This is what makes the
      // selection invariant under swapping the two labels.
      if (keep.count(c.id) || keep_decoy.count(c.id)) { kept_ids.insert(c.id); }
    }
    // Pass 2 walks the KEPT TARGETS and asks whether each has a decoy, rather than walking all
    // 7.1M compounds and asking whether each decoy has a kept target. Same set, because alias_of's
    // key set is exactly {compound ids in decoy_refs} -- the predicate the old loop filtered on.
    // One concatenation per kept target (~4e5) instead of one substr() per decoy compound (~3.5M).
    std::size_t n_dec = 0;
    {
      std::vector<std::string> add;
      add.reserve(kept_ids.size());
      for (const auto& id : kept_ids)              // holds TARGETS only at this point
      {
        std::string d = dtag + id;
        if (alias_of.count(d)) { add.push_back(std::move(d)); }
      }
      n_dec = add.size();
      for (auto& d : add) { kept_ids.insert(std::move(d)); }   // insert AFTER iterating
    }

    OpenSwath::LightTargetedExperiment out;                    // rebuild, preserving library order
    for (const auto& c : transition_exp.getCompounds())
    {
      if (kept_ids.count(c.id)) { out.compounds.push_back(c); }
    }
    for (const auto& t : transition_exp.getTransitions())
    {
      if (kept_ids.count(t.getPeptideRef())) { out.transitions.push_back(t); }
    }
    out.proteins = transition_exp.proteins;
    OPENMS_LOG_INFO << "OpenDIAlyzer[prefilter] " << n_before << " -> " << out.compounds.size()
                    << " precursors (" << (out.compounds.size() - n_dec) << " target / " << n_dec
                    << " decoy), " << out.transitions.size() << " transitions; evidence: "
                    << res.supported_precursors << "/" << res.total_target_precursors
                    << " targets supported." << std::endl;
    // Guard on the RATIO, not on n_dec != 0. The 23-decoy run is precisely why: 23 is not zero, so
    // an `n_dec == 0` test passes it, and the search then runs with a 154,000:1 target:decoy ratio
    // and a meaningless FDR. Correct pairing yields n_dec/n_target ~= the library-wide decoy ratio
    // by construction, so checking that IS a self-test of the pairing logic -- any partial pairing
    // failure (suffix tags, a tag mismatch, charge-suffix drift) shows up here instead of silently
    // producing q-values nobody can trust. Orphan decoys (target absent from the library) legitimately
    // depress the ratio slightly, hence a band rather than equality.
    const std::size_t n_tgt = out.compounds.size() - n_dec;
    const std::size_t lib_dec = decoy_refs.size();
    const std::size_t lib_tgt = n_before > lib_dec ? n_before - lib_dec : 0;
    const double r_lib  = lib_tgt > 0 ? static_cast<double>(lib_dec) / static_cast<double>(lib_tgt) : 0.0;
    const double r_kept = n_tgt   > 0 ? static_cast<double>(n_dec)   / static_cast<double>(n_tgt)   : 0.0;
    if (n_dec == 0 || (r_lib > 0.0 && r_kept < 0.8 * r_lib))
    {
      OPENMS_LOG_ERROR << "OpenDIAlyzer[prefilter] decoy pairing looks broken: kept target:decoy ratio "
                       << r_kept << " vs library ratio " << r_lib << " (" << n_dec << " decoys for "
                       << n_tgt << " targets). FDR would be uncalibrated. Check that -decoy_tag ('"
                       << dtag << "') matches how the library was built." << std::endl;
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        "prefilter decoy pairing produced an implausible target:decoy ratio");
    }
    transition_exp = std::move(out);

    // Hand the freed library back to the OS. The move above destroys ~94% of the library (7,149,966
    // precursors down to ~423,000), but that library is ~140 MILLION small string allocations -- two
    // per transition, both past the 15-char SSO threshold -- and glibc keeps freed blocks that small
    // in its arena free-lists rather than returning the pages. Measured on this library: 38.79 GB
    // loaded, 27.67 GB still resident after the move, 17.80 GB after malloc_trim(0). Nearly 10 GB,
    // for one call, held for the entire extraction that follows.
    //
    // It does NOT recover everything: ~15 GB stays resident against ~2.3 GB of live data, because
    // the surviving 6% is interleaved with the freed 94% in the same pages and no amount of trimming
    // can return a page with one live object on it. That residue needs the library not to be built
    // out of 140M separate allocations in the first place; this is the cheap half of the fix.
#ifdef __GLIBC__
    malloc_trim(0);
#endif
    return true;
  }

  // --- one in-process extraction pass -----------------------------------------
  ExitCodes extractPass_(const std::vector<OpenSwath::SwathMap>& swath_maps,
                         const std::shared_ptr<ExperimentalSettings>& exp_meta,
                         const OpenSwath::LightTargetedExperiment& transition_exp,
                         const TransformationDescription& trafo,
                         double rt_window, const std::string& in_file, const std::string& osw_path)
  {
    // The sqlite path seeds the .osw with the library tables (PRECURSOR / TRANSITION / PEPTIDE) so
    // the feature rows have something to join against. The parquet writer takes the assay library
    // as an argument instead, so there is nothing to seed.
    if (!parquet_out_)
    {
      TransitionPQPFile().convertLightTargetedExperimentToPQP(osw_path.c_str(), transition_exp);
    }

    ChromExtractParams cp = makeChromParams_(rt_window);
    ChromExtractParams cp_ms1 = makeMs1ChromParams_(cp);

    Param ff = makeFeatureFinderParam_();
    // Progress reporting. OpenSwathWorkflow prints nothing between "will analyse N peptides"
    // and completion, so a long run is indistinguishable from a hung one -- during this
    // session that ambiguity repeatedly cost time. Emit the shape of the work and a wall-clock
    // stamp at every phase boundary, so elapsed time can be attributed to a phase.
    {
      std::size_t n_prec = transition_exp.getCompounds().size();
      std::size_t n_tr   = transition_exp.getTransitions().size();
      std::size_t n_ms2  = 0;
      for (const auto& sm : swath_maps) { if (!sm.ms1) { ++n_ms2; } }
      OPENMS_LOG_INFO << "OpenDIAlyzer[progress] extraction start: " << n_prec << " precursors, "
                      << n_tr << " transitions, " << n_ms2 << " MS2 windows, "
                      << getIntOption_("threads") << " threads, RT window "
                      << rt_window << " s" << std::endl;
    }
    const auto t_extract0 = std::chrono::steady_clock::now();

    // PARQUET path: hand the OSW writer an EMPTY filename so it stays inactive, and ask
    // performExtraction to retain the features instead. OpenSwathWorkflow.cpp:911 gates this
    // explicitly (`if (!osw_writer.isActive() && store_features)`) -- features go either to sqlite
    // or to the FeatureMap, never both. Keeping them means the run never writes 2.07M features to
    // disk, rebuilds their PRECURSOR_ID with a SQL hash join, and reads them back, purely to score
    // them in the same process. Measured: that round trip was 121 s of SYSTEM time per scoring
    // call, 78% of the phase's CPU, and it happened once per pass.
    //
    // NOTE: with sqlite the writer had to be DESTROYED before remapFeaturePrecursorIds_ could run
    // (it held an open connection and the remap needs exclusive access). Neither exists here.
    {
    OpenSwathOSWWriter oswwriter(parquet_out_ ? std::string() : osw_path, /*uis*/ false);
    if (!parquet_out_) { oswwriter.writeHeader(); }             // create RUN / FEATURE / FEATURE_MS2
    const UInt64 cur_run = UniqueIdGenerator::getUniqueId();
    Interfaces::IMSDataConsumer* chrom = new NoopMSDataWritingConsumer("");   // no chromatogram output
    if (!parquet_out_) { oswwriter.addRun(cur_run, in_file); }
    oswwriter.setRunId(cur_run);
    run_id_ = cur_run;

    // Detect diaPASEF (SWATHs carry IM bounds). Enabling pasef is required for correct
    // diaPASEF IM windowing (C2). outer_loop_threads = -1 selects the OpenMS 3.6 SWATH
    // WAVE SCHEDULER (measured 2.7x faster than the legacy nested path that >=0 requests);
    // it needs batchSize 0 + in-memory reads (handled in loadDIARun_/load_into_memory).
    bool has_im = false;
    for (const auto& sm : swath_maps)
    {
      if (!sm.ms1 && sm.imLower >= 0 && sm.imUpper >= 0) { has_im = true; break; }
    }
    const std::string pmode = getStringOption_("pasef");      // auto|true|false (isolation control)
    const bool pasef = (pmode == "true") ? true : (pmode == "false") ? false : has_im;
    // Ion-mobility sub-scores (VAR_IM_*): the TOPP tool defaults these to "auto" and resolves
    // auto->true for PASEF. Without this they are silently NULL -- on diaPASEF data that
    // discards four of the most discriminating features.
    ff.setValue("Scores:use_ion_mobility_scores", pasef ? "true" : "false");
    const bool use_ms1 = true, use_ms1_im = pasef, prm = false, mrm = false;
    const int outer_loop_threads = -1;                        // wave scheduler (NOT legacy)
    FeatureMap fmap;
    OpenSwathWorkflow wf(use_ms1, use_ms1_im, prm, pasef, mrm, outer_loop_threads);
    wf.setLogType(log_type_);
    Param mrm_map;
    // load_into_memory was hardcoded true, which materialises the working SWATH in RAM --
    // exactly the whole-run residency the streaming input exists to avoid, and it would have
    // silently defeated an mzPeak-backed run. Now driven by -readOptions.
    // Must AGREE with what loadDIARun_ actually did, or the workflow is told the maps are
    // resident when they are not -- so both sites go through the same resolver.
    std::string eff; bool load_into_memory = false;
    resolveReadOptions_(in_file, eff, load_into_memory);
    if (!load_into_memory && getIntOption_("batchSize") <= 0)
    {
      // INFO, not WARN: this is the default path now. The old WARN asserted "throughput will drop";
      // that is not established -- but neither is the reverse, since the only comparison we have is
      // confounded by thread count (see the -readOptions registration). State the mechanism, claim
      // nothing about speed.
      OPENMS_LOG_INFO << "OpenDIAlyzer: readOptions='" << eff << "' -- streaming working SWATH map "
                         "(bounded memory, legacy inner-batch scheduling, as OpenSwathWorkflow does "
                         "by default). The SWATH wave scheduler needs full residency; pass "
                         "-readOptions cacheWorkingInMemory to use it and compare on your data."
                      << std::endl;
    }
    wf.performExtraction(swath_maps, trafo, cp, cp_ms1, ff, transition_exp,
                         fmap, /*store_features*/ parquet_out_, oswwriter, chrom,
                         getIntOption_("batchSize"), getIntOption_("ms1_isotopes"),
                         load_into_memory,
                         mrm_map, /*mobilogram*/ nullptr, getIntOption_("innerBatchSize"),
                         getIntOption_("max_concurrent_swaths"));
    {
      const double secs = std::chrono::duration<double>(
          std::chrono::steady_clock::now() - t_extract0).count();
      OPENMS_LOG_INFO << "OpenDIAlyzer[progress] extraction done in " << (int) secs << " s ("
                      << (int) (secs / 60.0) << " min)" << std::endl;
    }
    if (parquet_out_)
    {
      OPENMS_LOG_INFO << "OpenDIAlyzer[progress] retained " << fmap.size()
                      << " features in memory (no sqlite round trip)." << std::endl;
      { PhaseTimer pt_fm("setup/retain_features");
        pass_features_ = std::move(fmap); }   // must happen while fmap is still in scope
    }
    { PhaseTimer pt_ch("setup/free_chromatograms");
    delete chrom; }
    }                                         // <- oswwriter closes its sqlite connection here
    // remapFeaturePrecursorIds_ is a SQLITE-ONLY fixup: OpenSwathOSWWriter emits whatever string id
    // the library carried into a column declared INT NOT NULL. The parquet schema declares
    // PRECURSOR_ID as int64 outright (ArrowSchemaRegistry OSWFeatureSchema), so there is nothing to
    // repair and the per-pass hash-join rebuild of a 2.07M-row table simply does not happen.
    if (!parquet_out_) { remapFeaturePrecursorIds_(osw_path); }
    return EXECUTION_OK;
  }

  // Canonical OSW requires FEATURE.PRECURSOR_ID to be the INTEGER PRECURSOR.ID (the column is
  // declared INT NOT NULL). OpenSwathOSWWriter writes whatever string id the library carried:
  // a PQP-loaded library yields the integer id already, but our TSV-loaded library yields the
  // TransitionGroupId string -- so pyprophet/duckdb rejects our .osw with a type mismatch
  // ("column declared as integer, found ... text"). Remap once after extraction. Idempotent:
  // rows that don't join (already integer) are left untouched.
  static void remapFeaturePrecursorIds_(const std::string& osw)
  {
    sqlite3* db = nullptr;
    if (sqlite3_open(osw.c_str(), &db) != SQLITE_OK) { sqlite3_close(db); return; }
    // Nothing to do if the writer already produced integers.
    bool is_text = false;
    sqlite3_stmt* ts = nullptr;
    if (sqlite3_prepare_v2(db, "SELECT typeof(PRECURSOR_ID) FROM FEATURE LIMIT 1", -1, &ts, nullptr) == SQLITE_OK
        && sqlite3_step(ts) == SQLITE_ROW)
    {
      const unsigned char* t = sqlite3_column_text(ts, 0);
      is_text = (t && std::string(reinterpret_cast<const char*>(t)) == "text");
    }
    sqlite3_finalize(ts);
    if (!is_text) { sqlite3_close(db); return; }

    // A correlated `UPDATE ... (SELECT ...)` runs one indexed lookup PER FEATURE ROW and took
    // >80 min on a 48 GB / 4.15M-feature .osw. Rebuilding the table with a single hash join is
    // one pass over each table instead, and lets SQLite pick a join strategy.
    char* err = nullptr;
    auto run = [&](const char* sql) {
      if (sqlite3_exec(db, sql, nullptr, nullptr, &err) != SQLITE_OK)
      {
        OPENMS_LOG_ERROR << "OpenDIAlyzer: PRECURSOR_ID remap failed on [" << sql << "]: "
                         << (err ? err : "?") << std::endl;
        if (err) { sqlite3_free(err); err = nullptr; }
        return false;
      }
      if (err) { sqlite3_free(err); err = nullptr; }
      return true;
    };
    // Build the column list from the live schema: hardcoding it would silently drop columns
    // if the OSW schema gains fields (it has 11 today, incl. EXP_IM_LEFT/RIGHTWIDTH).
    std::string cols;
    std::string decl;                      // explicit "NAME TYPE [NOT NULL]" list for CREATE TABLE
    {
      sqlite3_stmt* ps = nullptr;
      if (sqlite3_prepare_v2(db, "PRAGMA table_info(FEATURE)", -1, &ps, nullptr) == SQLITE_OK)
      {
        while (sqlite3_step(ps) == SQLITE_ROW)
        {
          const unsigned char* nm = sqlite3_column_text(ps, 1);
          if (!nm) { continue; }
          const std::string c(reinterpret_cast<const char*>(nm));
          const unsigned char* ty = sqlite3_column_text(ps, 2);
          std::string t = ty ? reinterpret_cast<const char*>(ty) : "";
          const int notnull = sqlite3_column_int(ps, 3);
          // PRECURSOR_ID currently holds TEXT (the TRAML_ID); after the remap it is an integer key.
          if (c == "PRECURSOR_ID" || t.empty()) { t = "INT"; }
          if (!cols.empty()) { cols += ", "; decl += ", "; }
          cols += (c == "PRECURSOR_ID") ? "COALESCE(p.ID, f.PRECURSOR_ID) AS PRECURSOR_ID"
                                        : ("f." + c + " AS " + c);
          decl += c + " " + t + (notnull ? " NOT NULL" : "");
        }
        sqlite3_finalize(ps);
      }
    }
    if (cols.empty()) { OPENMS_LOG_ERROR << "OpenDIAlyzer: cannot read FEATURE schema for remap." << std::endl; sqlite3_close(db); return; }
    // CAST(... AS TEXT) is REQUIRED: FEATURE.PRECURSOR_ID is declared INT, so SQLite applies
    // numeric affinity when comparing it to the TEXT column PRECURSOR.TRAML_ID and the join
    // matches nothing (observed: every row fell through COALESCE and stayed text).
    // CREATE TABLE ... AS SELECT would be shorter, and was what this did -- but SQLite gives a CTAS
    // column NO DECLARED TYPE unless the select-item is a bare column reference. `COALESCE(...) AS
    // PRECURSOR_ID` is not, so the rebuilt table ended up with `PRECURSOR_ID` (no type), i.e. NONE
    // affinity. SQLite does not care; pyprophet reads the .osw through duckdb, which does, and dies
    // with `Unimplemented type for cast (BIGINT -> BLOB) when casting from source column ID`.
    // Declare the schema explicitly and INSERT into it so the column types survive.
    const std::string create_sql = "CREATE TABLE FEATURE_REMAP(" + decl + ");";
    const std::string insert_sql =
      "INSERT INTO FEATURE_REMAP SELECT " + cols +
      " FROM FEATURE f LEFT JOIN PRECURSOR p ON p.TRAML_ID = CAST(f.PRECURSOR_ID AS TEXT);";
    const bool ok =
      run("PRAGMA journal_mode=OFF;") &&
      run("PRAGMA synchronous=OFF;") &&
      run("CREATE INDEX IF NOT EXISTS ix_precursor_traml ON PRECURSOR(TRAML_ID);") &&
      run(create_sql.c_str()) &&
      run(insert_sql.c_str()) &&
      run("DROP TABLE FEATURE;") &&
      run("ALTER TABLE FEATURE_REMAP RENAME TO FEATURE;") &&
      run("CREATE INDEX IF NOT EXISTS ix_feature_precursor ON FEATURE(PRECURSOR_ID);");
    if (ok)
    {
      // Verify BOTH the stored value type AND the DECLARED column type. Checking only the former
      // is how this shipped broken: `typeof(PRECURSOR_ID)` returns "integer" because the values are
      // integers, so the old check passed and logged "pyprophet-compatible" on a table whose column
      // had no declared type at all -- which is precisely what pyprophet's duckdb reader rejects.
      std::string now = "?";
      sqlite3_stmt* vs = nullptr;
      if (sqlite3_prepare_v2(db, "SELECT typeof(PRECURSOR_ID) FROM FEATURE LIMIT 1", -1, &vs, nullptr) == SQLITE_OK
          && sqlite3_step(vs) == SQLITE_ROW)
      {
        const unsigned char* t = sqlite3_column_text(vs, 0);
        if (t) { now = reinterpret_cast<const char*>(t); }
      }
      sqlite3_finalize(vs);
      std::string declared = "";
      {
        sqlite3_stmt* ds = nullptr;
        if (sqlite3_prepare_v2(db, "PRAGMA table_info(FEATURE)", -1, &ds, nullptr) == SQLITE_OK)
        {
          while (sqlite3_step(ds) == SQLITE_ROW)
          {
            const unsigned char* nm = sqlite3_column_text(ds, 1);
            if (nm && std::string(reinterpret_cast<const char*>(nm)) == "PRECURSOR_ID")
            {
              const unsigned char* ty = sqlite3_column_text(ds, 2);
              declared = ty ? reinterpret_cast<const char*>(ty) : "";
              break;
            }
          }
          sqlite3_finalize(ds);
        }
      }
      if (declared.empty())
      {
        OPENMS_LOG_ERROR << "OpenDIAlyzer: FEATURE.PRECURSOR_ID has NO DECLARED TYPE after the "
                         << "remap. SQLite tolerates this; pyprophet (via duckdb) does not and will "
                         << "fail with a BIGINT->BLOB cast error. The .osw is NOT scorable."
                         << std::endl;
        now = "untyped";
      }
      if (now == "integer")
      {
        OPENMS_LOG_INFO << "OpenDIAlyzer: remapped FEATURE.PRECURSOR_ID to integer PRECURSOR.ID "
                        << "(hash-join rebuild) -- .osw is now pyprophet-compatible." << std::endl;
      }
      else
      {
        OPENMS_LOG_ERROR << "OpenDIAlyzer: PRECURSOR_ID remap did NOT take effect (typeof=" << now
                         << "); the .osw will not be readable by pyprophet." << std::endl;
      }
    }
    sqlite3_close(db);
  }



  /// Global-row-range access to an arrow::ChunkedArray without per-row shared_ptr traffic.
  ///
  /// Two facts make the obvious approaches wrong. Columns do NOT share a chunk layout -- string
  /// columns split at arrow's 2 GB array limit while an int64 column of the same table does not --
  /// so a chunk index from one column cannot address another. And ChunkedColumn::resolve() returns
  /// shared_ptr BY VALUE, so per-row access costs atomic refcount updates on control blocks shared
  /// by every thread; measured, that turned a 20.1 s serial loop into 34.5 s on 64 threads.
  ///
  /// This resolves once per (column, segment) and hands back a raw pointer plus the local offset.
  struct ChunkCursor
  {
    const arrow::ChunkedArray* col = nullptr;
    std::vector<int64_t> start;                 ///< global row where each chunk begins
    explicit ChunkCursor(const std::shared_ptr<arrow::ChunkedArray>& c) : col(c.get())
    {
      if (!col) { return; }
      start.resize(col->num_chunks() + 1, 0);
      for (int i = 0; i < col->num_chunks(); ++i) { start[i + 1] = start[i] + col->chunk(i)->length(); }
    }
    /// Chunk index containing global row r (binary search: once per segment, not per row).
    int chunkOf(int64_t r) const
    {
      return int(std::upper_bound(start.begin(), start.end(), r) - start.begin()) - 1;
    }
    const arrow::Array* chunk(int i) const { return col->chunk(i).get(); }
    int64_t local(int64_t r, int i) const { return r - start[i]; }
    int64_t chunkEnd(int i) const { return start[i + 1]; }
    bool valid() const { return col != nullptr && col->num_chunks() > 0; }
  };

  /// Text of row `r`, handling both utf8 and large_utf8 -- the writer picks either depending on
  /// column size, and a dynamic_cast to only one of them silently yields nothing. That is exactly
  /// how a probe reported "0 distinct fragment annotations" where there were 102.
  static std::string_view arrowText(const arrow::Array* a, int64_t i)
  {
    if (!a || a->IsNull(i)) { return {}; }
    if (const auto* s = dynamic_cast<const arrow::StringArray*>(a)) { return std::string_view(s->GetView(i)); }
    if (const auto* l = dynamic_cast<const arrow::LargeStringArray*>(a)) { return std::string_view(l->GetView(i)); }
    return {};
  }



  /// Materialise a LightTargetedExperiment from the compact library, using SHORT synthetic ids.
  ///
  /// This is the point of the whole exercise. LightTransition owns two std::strings and
  /// LightCompound four; the library's real ids ("DECOY_PEPTIDEK_2_y7_1") exceed libstdc++'s
  /// 15-character SSO buffer, so 78.6M transitions cost ~471M heap allocations and leave the arena
  /// fragmented -- measured at 83-88% of a 186 GB peak. Synthetic ids
  /// ("p1z141z3" / "DECOY_p1z141z3", base-36, <=14 chars) stay INSIDE the SSO buffer, so the
  /// materialised experiment allocates nothing per row.
  ///
  /// The "DECOY_" prefix is preserved deliberately: prefilterLibrary_ pairs decoys to targets by
  /// exactly that convention, so the pairing logic needs no change. The real ids live on in the
  /// compact library and are restored for output.
  void materializeFromCompact_(const odia::CompactLibrary& clib,
                               OpenSwath::LightTargetedExperiment& exp) const
  {
    using CL = odia::CompactLibrary;
    exp.compounds.reserve(clib.peptideCount());
    exp.transitions.reserve(clib.transitionCount());

    std::unordered_map<std::string, std::size_t> prot_seen;
    for (std::size_t i = 0; i < clib.peptideCount(); ++i)
    {
      const auto p = CL::Peptide(std::uint32_t(i));
      OpenSwath::LightCompound c;
      c.id = CL::syntheticId(std::uint32_t(i), clib.isDecoy(p));      // SSO: no allocation
      c.sequence = std::string(clib.sequence(p));                     // real sequence: scoring needs it
      c.charge = clib.charge(p);
      c.rt = clib.rt(p);
      c.drift_time = clib.driftTime(p);
      const auto parent = clib.parent(p);
      if (parent != CL::no_protein)
      {
        std::string acc(clib.accession(parent));
        if (!acc.empty()) { c.protein_refs.push_back(std::move(acc)); }
      }
      exp.compounds.push_back(std::move(c));
    }
    for (std::size_t i = 0; i < clib.transitionCount(); ++i)
    {
      const auto t = CL::Transition(std::uint32_t(i));
      const auto pep = clib.peptideOf(t);
      const std::uint32_t pi = static_cast<std::uint32_t>(pep);
      const std::uint8_t fl = clib.transitionFlags(t);
      OpenSwath::LightTransition tr;
      // "t" + base-36 index: <= 8 chars, always SSO.
      tr.transition_name = "t" + CL::syntheticId(std::uint32_t(i), false).substr(1);
      tr.peptide_ref = CL::syntheticId(pi, clib.isDecoy(pep));
      tr.precursor_mz = clib.precursorMz(pep);
      tr.product_mz = clib.productMz(t);
      tr.library_intensity = clib.intensity(t);
      tr.precursor_im = clib.driftTime(pep);
      tr.fragment_charge = clib.fragmentCharge(t);
      tr.setDecoy((fl & CL::Decoy) != 0);
      tr.setDetectingTransition((fl & CL::Detecting) != 0);
      tr.setIdentifyingTransition((fl & CL::Identifying) != 0);
      tr.setQuantifyingTransition((fl & CL::Quantifying) != 0);
      const auto ann = clib.annotation(t);
      if (!ann.empty()) { tr.setFragmentType(std::string(ann)); }
      exp.transitions.push_back(std::move(tr));
    }
    // Proteins: one entry per distinct accession actually referenced.
    for (std::size_t i = 0; i < clib.proteinCount(); ++i)
    {
      OpenSwath::LightProtein pr;
      pr.id = std::string(clib.accession(CL::Protein(std::uint32_t(i))));
      exp.proteins.push_back(std::move(pr));
    }
  }

  /// Load the library into CompactLibrary and report what it costs. A PROBE, not the production
  /// path: it exists so the compact representation can be compared against the measured 38.79 GB
  /// of the LightTargetedExperiment path on the same file and the same machine, before anything is
  /// restructured around it.
  ///
  /// The .oswpq carries protein ACCESSIONS but not protein SEQUENCES, so the substring win needs a
  /// FASTA. With one supplied, a peptide costs (protein, offset, length) and no characters;
  /// without, sequences are interned -- still no per-string malloc and no duplication, but every
  /// distinct peptide keeps its own characters. Reporting both is the point.
  ExitCodes compactProbe_(const std::string& lib, const std::string& fasta)
  {
    odia::CompactLibrary clib;
    std::unordered_map<std::string, odia::CompactLibrary::Protein> prot_by_acc;

    if (!fasta.empty())
    {
      PhaseTimer pt("compact_probe/fasta");
      std::vector<FASTAFile::FASTAEntry> entries;
      FASTAFile().load(fasta, entries);
      std::size_t chars = 0;
      for (const auto& e : entries) { chars += e.sequence.size(); }
      clib.reserve(entries.size(), 0, 0, chars);
      for (const auto& e : entries)
      {
        // OpenSWATH accessions are the FASTA id; index by both the raw id and its middle field
        // ("sp|P02768|ALBU_HUMAN" -> "P02768") because libraries differ in which they carry.
        const auto h = clib.addProtein(e.identifier, e.sequence);
        prot_by_acc[e.identifier] = h;
        const auto bar1 = e.identifier.find('|');
        if (bar1 != std::string::npos)
        {
          const auto bar2 = e.identifier.find('|', bar1 + 1);
          if (bar2 != std::string::npos) { prot_by_acc[e.identifier.substr(bar1 + 1, bar2 - bar1 - 1)] = h; }
        }
      }
      // INDEX THE PROTEOME ONCE, before any peptide is added. Every peptide is then looked up in
      // the whole FASTA rather than in whichever protein its accession names -- accessions go
      // missing, get renamed or versioned, and peptides are shared between proteins, and in each
      // of those cases an accession-keyed lookup stores characters for a sequence that is right
      // there in the file.
      clib.indexProteins();
      OPENMS_LOG_INFO << "OpenDIAlyzer[compact] FASTA: " << entries.size() << " proteins, "
                      << std::fixed << std::setprecision(2) << chars / 1048576.0
                      << " MB of sequence, k-mer indexed" << std::endl;
    }

    std::unique_ptr<File::TempDir> temp_dir;
    // ZipRandomAccessFile is internal to OpenMS and not in the installed prefix, so this takes the
    // documented fallback: extract the entry to a temp file and read it. Slower and it touches
    // disk, which is acceptable for a probe -- the measurement here is MEMORY, not load time.
    auto open_entry = [&](const std::string& entry) -> std::shared_ptr<arrow::Table> {
      return ParquetFile::readTable(ZipArchiveFile::extractEntryToTempFile(lib, entry, temp_dir));
    };

    std::unordered_map<long long, odia::CompactLibrary::Peptide> pep_by_id;
    {
      PhaseTimer pt("compact_probe/precursors");
      auto tbl = open_entry("library/precursors.parquet");
      auto id_c  = ParquetFile::getColumn(tbl, OSWPrecursorSchema::PRECURSOR_ID);
      auto ch_c  = ParquetFile::getColumn(tbl, OSWPrecursorSchema::CHARGE);
      auto mod_c = ParquetFile::getOptionalColumn(tbl, OSWPrecursorSchema::MODIFIED_SEQUENCE);
      auto unm_c = ParquetFile::getOptionalColumn(tbl, OSWPrecursorSchema::UNMODIFIED_SEQUENCE);
      auto acc_c = ParquetFile::getOptionalColumn(tbl, OSWPrecursorSchema::PROTEIN_ACCESSIONS);
      const int64_t n = tbl->num_rows();
      pep_by_id.reserve(n);
      clib.reserve(prot_by_acc.size(), n, 0, 0);
      clib.resizePeptides(std::size_t(n));

      // PASS 1, PARALLEL: resolve each peptide's span. locate() only READS the store, so it is
      // safe from many threads; interning a miss MUTATES it, so misses are only recorded here and
      // stored serially afterwards. The plain (non-chunked) Arrow arrays are read-only and shared.
      std::vector<odia::SequenceStore::Span> spans(n);
      std::vector<odia::CompactLibrary::Protein> parents(n, odia::CompactLibrary::no_protein);
      std::vector<int> charges(n, 0);
      std::vector<std::string> misses(n);          // sequence, only where locate() failed
      std::vector<long long> ids(n, 0);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
      for (long long r = 0; r < n; ++r)
      {
        const std::string acc = acc_c ? ParquetFile::getString(acc_c, r) : std::string();
        const std::string seq_m = mod_c ? ParquetFile::getString(mod_c, r) : std::string();
        const std::string seq_u = unm_c ? ParquetFile::getString(unm_c, r) : std::string();
        // Substring matching must use the UNMODIFIED sequence: a modified one carries "(UniMod:4)"
        // and occurs in no protein. The peptide identity for scoring is still the modified form,
        // which is why a library needs both.
        const std::string& seq = !seq_u.empty() ? seq_u : seq_m;
        auto parent = odia::CompactLibrary::no_protein;
        if (!acc.empty())
        {
          const auto first = acc.substr(0, acc.find('/'));
          const auto it = prot_by_acc.find(first);      // read-only after the FASTA pass
          if (it != prot_by_acc.end()) { parent = it->second; }
        }
        parents[r] = parent;
        charges[r] = int(ParquetFile::getInt64(ch_c, r, 0, true));
        ids[r] = ParquetFile::getInt64(id_c, r, 0, false);
        const auto sp = clib.locateOrNull(seq);        // proteome-wide, not accession-keyed
        if (sp.valid()) { spans[r] = sp; } else { misses[r] = seq; }
      }
      // PASS 2, SERIAL and only for the misses -- on this library that is the decoys, which are
      // shuffled and occur in no protein.
      // Serial pass for the misses -- interning MUTATES the store. Decoys are folded INTO the
      // index in batches as they are stored, so a decoy that shares a subsequence with one already
      // seen resolves to a span instead of paying for its own characters. Batched because each
      // indexInterned() sorts the new tail and merges: per-peptide would be quadratic.
      constexpr long long kIndexEvery = 200000;
      long long since_index = 0, late_hits = 0;
      for (long long r = 0; r < n; ++r)
      {
        bool derived = misses[r].empty();
        if (!derived)
        {
          // Retry against everything indexed since this peptide was first looked up.
          const auto retry = clib.locateOrNull(misses[r]);
          if (retry.valid()) { spans[r] = retry; derived = true; ++late_hits; }
          else
          {
            spans[r] = clib.internPeptideSequence(misses[r]);
            if (++since_index >= kIndexEvery) { clib.indexInterned(); since_index = 0; }
          }
        }
        clib.setPeptide(std::size_t(r), parents[r], spans[r], charges[r], derived);
        pep_by_id[ids[r]] = odia::CompactLibrary::Peptide(std::uint32_t(r));
      }
      clib.indexInterned();
      OPENMS_LOG_INFO << "OpenDIAlyzer[compact] " << late_hits
                      << " peptides resolved against previously-interned sequences (decoys included "
                      << "in the index); index " << std::fixed << std::setprecision(2)
                      << clib.indexBytes() / 1073741824.0 << " GB" << std::endl;
    }
    {
      PhaseTimer pt("compact_probe/transitions");
      auto tbl = open_entry("library/transitions.parquet");
      using CC = ParquetFile::ChunkedColumn;
      const int64_t n = tbl->num_rows();
      clib.resizeTransitions(std::size_t(n));

      // PARALLEL OVER ROW RANGES, with each column resolved once per segment.
      std::unordered_map<std::string, odia::CompactLibrary::AnnotationId> ann_id;
      std::mutex ann_mu;
      std::atomic<std::size_t> unmapped{0};
      {
        ChunkCursor pidc(tbl->GetColumnByName(OSWTransitionSchema::PRECURSOR_ID));
        ChunkCursor mzc (tbl->GetColumnByName(OSWTransitionSchema::PRODUCT_MZ));
        ChunkCursor inc (tbl->GetColumnByName(OSWTransitionSchema::LIBRARY_INTENSITY));
        ChunkCursor annc(tbl->GetColumnByName(OSWTransitionSchema::ANNOTATION));
        if (!pidc.valid() || !mzc.valid() || !inc.valid())
        {
          OPENMS_LOG_ERROR << "OpenDIAlyzer[compact] transitions table is missing a required column."
                           << std::endl;
          return INPUT_FILE_CORRUPT;
        }
        // Fixed block count, so the partition depends on the DATA rather than the thread count.
        const int64_t nblk = std::max<int64_t>(1, std::min<int64_t>(4096, n / 20000));
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1)
#endif
        for (int64_t b = 0; b < nblk; ++b)
        {
          const int64_t lo = n * b / nblk, hi = n * (b + 1) / nblk;
          std::unordered_map<std::string_view, odia::CompactLibrary::AnnotationId> local_ann;
          std::size_t local_unmapped = 0;
          for (int64_t r = lo; r < hi; )
          {
            // Resolve every column ONCE for this segment, then walk it.
            const int ci = pidc.chunkOf(r), cm = mzc.chunkOf(r), cn = inc.chunkOf(r);
            const int ca = annc.valid() ? annc.chunkOf(r) : -1;
            int64_t seg = std::min({pidc.chunkEnd(ci), mzc.chunkEnd(cm), inc.chunkEnd(cn), hi});
            if (ca >= 0) { seg = std::min(seg, annc.chunkEnd(ca)); }
            const auto* pa = static_cast<const arrow::Int64Array*>(pidc.chunk(ci));
            const auto* ma = static_cast<const arrow::DoubleArray*>(mzc.chunk(cm));
            const auto* ia = static_cast<const arrow::DoubleArray*>(inc.chunk(cn));
            const arrow::Array* aa = ca >= 0 ? annc.chunk(ca) : nullptr;
            for (int64_t g = r; g < seg; ++g)
            {
              const auto it = pep_by_id.find(pa->Value(pidc.local(g, ci)));
              if (it == pep_by_id.end()) { ++local_unmapped; continue; }
              odia::CompactLibrary::AnnotationId aid = odia::SequenceStore::npos;
              if (aa)
              {
                const std::string_view av = arrowText(aa, annc.local(g, ca));
                if (!av.empty())
                {
                  const auto f = local_ann.find(av);
                  if (f != local_ann.end()) { aid = f->second; }
                  else
                  {
                    std::lock_guard<std::mutex> lk(ann_mu);
                    const std::string key(av);
                    auto g2 = ann_id.find(key);
                    if (g2 == ann_id.end()) { g2 = ann_id.emplace(key, clib.internAnnotation(av)).first; }
                    aid = g2->second;
                    local_ann.emplace(av, aid);
                  }
                }
              }
              clib.setTransition(std::size_t(g), it->second,
                                 ma->Value(mzc.local(g, cm)), float(ia->Value(inc.local(g, cn))), aid);
            }
            r = seg;
          }
          unmapped += local_unmapped;
        }
      }
      OPENMS_LOG_INFO << "OpenDIAlyzer[compact] " << ann_id.size() << " distinct fragment annotations"
                      << (unmapped ? ", " + std::to_string(unmapped.load()) + " transitions with no precursor" : "")
                      << std::endl;
    }
    clib.finalize();

    const auto st = clib.stats();
    const double gb = 1073741824.0;
    OPENMS_LOG_INFO << "OpenDIAlyzer[compact] " << st.proteins << " proteins, " << st.peptides
                    << " peptides (" << st.peptides_from_fasta << " as FASTA substrings, "
                    << st.peptides_standalone << " standalone), " << st.transitions << " transitions"
                    << std::endl;
    OPENMS_LOG_INFO << "OpenDIAlyzer[compact] CompactLibrary " << std::fixed << std::setprecision(2)
                    << st.bytes / gb << " GB  vs  string-bearing objects " << st.bytes_as_objects / gb
                    << " GB  (" << (st.bytes ? double(st.bytes_as_objects) / double(st.bytes) : 0.0)
                    << "x)" << std::endl;
    MemProbe::logAllocator("compact probe done");
    return EXECUTION_OK;
  }

  ExitCodes main_(int, const char**) override
  {
    if (getFlag_("selftest")) { return selftest_(); }
    const std::string score_only = getStringOption_("score_osw");
    if (!score_only.empty()) { finalScore_(score_only); return EXECUTION_OK; }

    const std::string in = getStringOption_("in");
    const std::string tr = getStringOption_("tr");
    const std::string out = getStringOption_("out");
    const int passes = getIntOption_("recal_passes");

    // Convert-and-exit: needs -tr only, so it runs without -in/-out and before their check.
    const std::string convert_to = getStringOption_("convert_library");
    if (!convert_to.empty())
    {
      if (tr.empty())
      {
        OPENMS_LOG_ERROR << "OpenDIAlyzer: -convert_library needs -tr." << std::endl;
        return MISSING_PARAMETERS;
      }
      const auto t0 = std::chrono::steady_clock::now();
      OpenSwath::LightTargetedExperiment exp = loadLibrary_(tr);
      const double t_read = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
      const auto t1 = std::chrono::steady_clock::now();
      // Dispatch on the requested extension: .pqp -> SQLite, anything else -> parquet bundle.
      // NOTE parquet cannot currently be READ BACK at whole-proteome scale (chunked columns vs a
      // reader that combines/truncates them -- see libraryCachePath_), so .pqp is the working path.
      if (FileHandler::getTypeByFileName(convert_to) == FileTypes::PQP)
      {
        TransitionPQPFile().convertLightTargetedExperimentToPQP(convert_to.c_str(), exp);
      }
      else
      {
        TransitionParquetFile().convertLightTargetedExperimentToParquet(convert_to, exp);
      }
      const double t_write = std::chrono::duration<double>(std::chrono::steady_clock::now() - t1).count();
      OPENMS_LOG_INFO << "OpenDIAlyzer: converted " << exp.getCompounds().size() << " precursors / "
                      << exp.getTransitions().size() << " transitions -> " << convert_to
                      << " (read " << t_read << " s, write " << t_write << " s)." << std::endl;
      return EXECUTION_OK;
    }

    // Standalone probe: needs a library and nothing else, so it must be handled BEFORE the
    // -in/-tr/-out requirement, next to the other modes that do not run a search.
    if (!getStringOption_("compact_probe").empty())
    {
      MemProbe::instance().start();
      MemProbe::logAllocator("startup");
      const ExitCodes rc = compactProbe_(getStringOption_("compact_probe"),
                                         getStringOption_("compact_probe_fasta"));
      MemProbe::instance().stop();
      MemProbe::instance().report();
      return rc;
    }

    if (in.empty() || tr.empty() || out.empty())
    {
      OPENMS_LOG_ERROR << "OpenDIAlyzer: -in, -tr and -out are required." << std::endl;
      return MISSING_PARAMETERS;
    }
    if (passes < 1 || passes > 2)                              // C17
    {
      OPENMS_LOG_ERROR << "OpenDIAlyzer: recal_passes must be 1 (single) or 2 (recalibrated)." << std::endl;
      return ILLEGAL_PARAMETERS;
    }
#ifdef _OPENMP
    omp_set_num_threads(getIntOption_("threads"));
    // Extraction parallelism is capped by the number of SWATH windows (24 on this data), so a
    // 224-core machine sits ~90% idle in that phase. OpenSwathWorkflow now also splits the
    // coordinate list inside each SWATH, which needs a SECOND active OpenMP level.
    omp_set_max_active_levels(2);
#endif

    // load library + DIA run ONCE (reused across passes)
    // .oswpq output selects the parquet path: features stay in memory and nothing is written to
    // sqlite. Chosen by extension, exactly as TOPP OpenSwathWorkflow does it.
    MemProbe::instance().start();
    MemProbe::logAllocator("startup");

    // RAII, not a call before `return`: main_ has six return sites and an early one would have been
    // missed. This reports on every exit path, including the error ones -- which are exactly the
    // runs where knowing the memory profile matters most.
    struct MemReport
    {
      ~MemReport()
      {
        MemProbe::logAllocator("final");
        MemProbe::instance().stop();
        MemProbe::instance().report();
      }
    } mem_report_guard;
    (void)mem_report_guard;
    parquet_out_ = (FileHandler::getTypeByFileName(out) == FileTypes::OSWPQ);
    if (parquet_out_)
    {
      OPENMS_LOG_INFO << "OpenDIAlyzer: parquet output -- features are scored in memory; no .osw is "
                      << "written and no PRECURSOR_ID remap is needed." << std::endl;
    }

    OpenSwath::LightTargetedExperiment transition_exp;
    { PhaseTimer pt("library_load"); transition_exp = loadLibrary_(tr); }
    {
      const double gb = libraryBytes_(transition_exp) / 1073741824.0;
      OPENMS_LOG_INFO << "OpenDIAlyzer[mem/component] library AS LOADED: "
                      << transition_exp.getCompounds().size() << " compounds, "
                      << transition_exp.getTransitions().size() << " transitions, "
                      << std::fixed << std::setprecision(2) << gb << " GB accounted" << std::endl;
      MemProbe::logAllocator("after library_load");
    }
    std::shared_ptr<ExperimentalSettings> exp_meta;
    std::vector<OpenSwath::SwathMap> swath_maps;
    const std::string tmp = getStringOption_("tempDirectory");
    // Detect IM immediately after loading, NOT inside calibration: with -rt_calibration none
    // calibration never runs, and pasef_ would stay false -- silently disabling ion-mobility
    // extraction on genuine diaPASEF data.
    bool run_ok = false;
    { PhaseTimer pt("dia_run_load"); run_ok = loadDIARun_(in, exp_meta, swath_maps, tmp); }
    if (!run_ok)
    {
      OPENMS_LOG_ERROR << "OpenDIAlyzer: failed to load DIA run " << in << std::endl;
      return INPUT_FILE_CORRUPT;
    }
    {
      detectPasef_(swath_maps);                // sets pasef_ before ANY extraction params are built
      OPENMS_LOG_INFO << "OpenDIAlyzer: ion mobility " << (pasef_ ? "PRESENT" : "absent")
                      << "; IM extraction window "
                      << ((pasef_ && getDoubleOption_("ion_mobility_window") > 0.0) ? "on" : "OFF")
                      << "." << std::endl;

      // Screen the library against evidence in THIS run before extracting anything. Must come
      // after detectPasef_ (it needs the IM flag) and before any extraction pass, since every
      // pass and the calibration all consume transition_exp.
      // The prefilter and the extraction want OPPOSITE things from the m/z window, and until now
      // they shared it. Screening asks "could this precursor be here at all" -- a narrow window
      // makes that test miss real precursors, and they are then gone before any scorer sees them.
      // Extraction asks "what is the signal" -- a wide window admits interference that corrupts the
      // correlation sub-scores. Measured on this benchmark: at 5 ppm the prefilter kept 95,676
      // precursors, at 10 ppm it kept 423,079, and a comparison against DIA-NN's identifications
      // showed 924 of its IDs (11.9%) fall OUTSIDE the survivors -- an unreachable ceiling that no
      // classifier can recover. Default: screen at 2x the extraction window unless told otherwise.
      const double pf_mz = getDoubleOption_("prefilter_mz_extraction_window");
      const ChromExtractParams cp_pf     = makeChromParams_(getDoubleOption_("rt_extraction_window"),
                                                            pf_mz);
      const ChromExtractParams cp_pf_ms1 = makeMs1ChromParams_(cp_pf);
      { PhaseTimer pt("prefilter"); prefilterLibrary_(swath_maps, transition_exp, cp_pf, cp_pf_ms1); }
      // Build AFTER prefiltering: the index must describe the library actually searched, or the
      // integer precursor ids in the output will not line up with the features.
      {
        const double gb = libraryBytes_(transition_exp) / 1073741824.0;
        OPENMS_LOG_INFO << "OpenDIAlyzer[mem/component] library AFTER prefilter: "
                        << transition_exp.getCompounds().size() << " compounds, "
                        << transition_exp.getTransitions().size() << " transitions, "
                        << std::fixed << std::setprecision(2) << gb << " GB accounted" << std::endl;
        MemProbe::logAllocator("after prefilter");
      }
      if (parquet_out_) { PhaseTimer pt("precursor_index"); buildPrecursorIndex_(transition_exp); }

      // Dump-and-exit. Tuning the prefilter otherwise costs a full extraction (hours, ~1 TB) per
      // threshold tried, which makes it untunable in practice. With this, the survivor list is
      // obtainable in the time it takes to read the library and scan the spectra, so no-loss can
      // be checked against a reference ID list directly.
      const std::string pf_out = getStringOption_("prefilter_out");
      if (!pf_out.empty())
      {
        std::set<std::string> dec;
        for (const auto& t : transition_exp.getTransitions())
        {
          if (t.getDecoy()) { dec.insert(t.getPeptideRef()); }
        }
        std::ofstream os(pf_out);
        if (!os) { OPENMS_LOG_ERROR << "OpenDIAlyzer: cannot write " << pf_out << std::endl; return CANNOT_WRITE_OUTPUT_FILE; }
        os << "id\tsequence\tdecoy\n";
        for (const auto& c : transition_exp.getCompounds())
        {
          os << c.id << '\t' << c.sequence << '\t' << (dec.count(c.id) ? 1 : 0) << '\n';
        }
        OPENMS_LOG_INFO << "OpenDIAlyzer: wrote " << transition_exp.getCompounds().size()
                        << " surviving precursors to " << pf_out << " (-prefilter_out set, exiting "
                        << "before extraction)." << std::endl;
        return EXECUTION_OK;
      }
    }

    // Linear bootstrap calibration. The library RT is normalized to ~[0,1]; map it to
    // the run's actual RT range so pass 1 can use a BOUNDED window. (rt_win=-1 / whole
    // range exploded the candidate-peak count ~100x and choked the single-threaded OSW
    // writer.) recalibrate_ then refines this map nonlinearly from confident IDs.
    double rt_min = std::numeric_limits<double>::max(), rt_max = std::numeric_limits<double>::lowest();
    // 271.8 s sat un-phased between precursor_index and pass 1 -- 14% of the run in code nothing
    // was watching. Two 300-call metadata probes look free on a resident map and are not
    // necessarily free on a streaming one; time them rather than assume.
    { PhaseTimer pt_rtr("setup/run_rt_range");
    for (const auto& sm : swath_maps)                          // any map (MS1 or MS2), not only MS1 (C18)
    {
      if (!sm.sptr) { continue; }
      const size_t ns = sm.sptr->getNrSpectra();
      if (ns == 0) { continue; }
      const double a = sm.sptr->getSpectrumMetaById(0).RT;     // spectra are RT-ordered
      const double b = sm.sptr->getSpectrumMetaById(static_cast<int>(ns) - 1).RT;
      rt_min = std::min(rt_min, std::min(a, b));
      rt_max = std::max(rt_max, std::max(a, b));
    }
    }
    if (!(rt_max > rt_min) || !std::isfinite(rt_min) || !std::isfinite(rt_max))
    {
      OPENMS_LOG_ERROR << "OpenDIAlyzer: could not determine a finite run RT range from the input." << std::endl;
      return INCOMPATIBLE_INPUT_DATA;
    }
    // interpolated+linear (proven in fitTrafo_); fitModel("linear") degenerated to identity.
    const std::string calib_mode = getStringOption_("rt_calibration");
    const bool calib_none = (calib_mode == "none");
    bool calib_done = false;
    TransformationDescription pass_map;                       // pass 1: library RT -> run seconds
    // OpenSWATH convention: performExtraction takes a trafo mapping RUN RT -> normalized iRT.
    // It inverts that itself for the extraction windows, and the scorer uses it directly
    // (normalized_experimental_rt = trafo.apply(exp_rt)) to compare against the library's iRT.
    // Feeding it identity while pre-scaling compound.rt into seconds gives correct WINDOWS but
    // a wrong NORM_RT scale, which changes which candidate peaks are reported. So when we have
    // a real calibration we hand OpenMS the transform in its own convention and leave the
    // library in iRT units -- exactly what standalone OpenSwathWorkflow does.
    TransformationDescription native_trafo;                   // run RT -> iRT (OpenMS convention)
    bool use_native_trafo = false;
    double calib_rt_window = -1.0;                            // window estimated by the calibration
    const bool use_estimated_window = getStringOption_("use_estimated_rt_window") != "false";
    if (calib_mode == "cirt")
    {
      const bool pasef_c = detectPasef_(swath_maps);
      double est_win = -1.0;
      if (calibrateCiRT_(swath_maps, transition_exp, pasef_c,
                         getDoubleOption_("rt_extraction_window"), in, pass_map, est_win))
      {
        calib_rt_window = est_win;
        // Direction: OpenSWATH's calibration transform maps RUN RT -> normalized iRT (the
        // extraction code inverts it before use), but we need LIBRARY RT -> run seconds to
        // pre-scale compound.rt. Pick the orientation whose image of the library RT range
        // actually covers the run's RT range -- a midpoint-only test is not enough (an
        // iRT->iRT map lands near 0, which can still sit inside a lenient tolerance).
        double lo = std::numeric_limits<double>::max(), hi = std::numeric_limits<double>::lowest();
        for (const auto& c : transition_exp.compounds) { lo = std::min(lo, c.rt); hi = std::max(hi, c.rt); }
        const double run_span = rt_max - rt_min;
        auto score_dir = [&](const TransformationDescription& t) {   // lower is better; <0 = invalid
          const double a = t.apply(lo), b = t.apply(hi);
          if (!std::isfinite(a) || !std::isfinite(b)) { return -1.0; }
          const double out_span = std::abs(b - a);
          const double out_mid = 0.5 * (a + b);
          if (out_span < 0.25 * run_span) { return -1.0; }           // collapsed -> wrong direction
          if (out_mid < rt_min - 0.5 * run_span || out_mid > rt_max + 0.5 * run_span) { return -1.0; }
          return std::abs(out_span - run_span) / run_span + std::abs(out_mid - 0.5 * (rt_min + rt_max)) / run_span;
        };
        const TransformationDescription original = pass_map;
        TransformationDescription inv = pass_map;
        inv.invert();
        const double s_fwd = score_dir(pass_map), s_inv = score_dir(inv);
        if (s_inv >= 0.0 && (s_fwd < 0.0 || s_inv < s_fwd))
        {
          OPENMS_LOG_INFO << "OpenDIAlyzer: inverting calibration transform (run RT -> iRT becomes library RT -> run RT)." << std::endl;
          pass_map = inv;
          native_trafo = original;         // original was run RT -> iRT: exactly what OpenMS wants
          calib_done = true;
        }
        else if (s_fwd >= 0.0)
        {
          native_trafo = inv;              // pass_map is library -> run, so the inverse is run -> iRT
          calib_done = true;
        }
        else
        {
          OPENMS_LOG_WARN << "OpenDIAlyzer: calibration transform maps the library RT range ["
                          << lo << ", " << hi << "] outside the run range [" << rt_min << ", " << rt_max
                          << "] in both orientations; falling back to the linear bootstrap." << std::endl;
          pass_map = TransformationDescription();
        }
        if (calib_done)
        {
          use_native_trafo = !native_trafo.getDataPoints().empty();
          OPENMS_LOG_INFO << "OpenDIAlyzer: CiRT calibration applied (library RT " << lo << ".." << hi
                          << " -> " << pass_map.apply(lo) << ".." << pass_map.apply(hi)
                          << " s; run range " << rt_min << ".." << rt_max << "); "
                          << (use_native_trafo ? "passing the run-RT -> iRT transform to OpenMS (library kept in iRT units)."
                                               : "pre-scaling the library RT.") << std::endl;
        }
      }
    }
    if (calib_none)
    {
      pass_map.fitModel("identity");                          // use library RT verbatim (== OSW w/o -tr_irt)
      OPENMS_LOG_INFO << "OpenDIAlyzer: rt_calibration=none -> library RT used verbatim; run RT range ["
                      << rt_min << ", " << rt_max << "] s." << std::endl;
    }
    else if (!calib_done)
    {
      pass_map.setDataPoints(std::vector<std::pair<double, double>>{
          {0.0, rt_min}, {0.5, 0.5 * (rt_min + rt_max)}, {1.0, rt_max}});
      Param bp;
      bp.setValue("interpolation_type", "linear");
      bp.setValue("extrapolation_type", "four-point-linear");
      pass_map.fitModel("interpolated", bp);
      OPENMS_LOG_INFO << "OpenDIAlyzer: run RT range [" << rt_min << ", " << rt_max
                      << "] s; linear bootstrap apply(0.5)=" << pass_map.apply(0.5) << "s (expect ~"
                      << 0.5 * (rt_min + rt_max) << ")." << std::endl;
    }
    TransformationDescription identity;
    identity.fitModel("identity");
    const bool use_empirical = getStringOption_("empirical_rt") != "false";
    std::map<std::string, double> empirical_rt;               // TRAML_ID -> observed apex RT (pass-1 IDs)
    for (int p = 1; p <= passes; ++p)
    {
      // Two ways to place the RT window:
      //  (a) native (a real calibration exists): leave the library in iRT units and hand
      //      OpenMS the run-RT -> iRT transform, as standalone OpenSwathWorkflow does. Both
      //      the extraction window AND the NORM_RT score are then computed the same way.
      //  (b) fallback (bootstrap / none): no trustworthy transform, so pre-scale compound.rt
      //      into run seconds and extract with identity.
      if (!use_native_trafo) { rescaleLibraryRT_(transition_exp, pass_map); }
      const bool narrow = (p > 1);
      // Empirical library on recalibrated passes: overwrite RT of pass-1-identified
      // precursors with their measured apex RT (centres their window on the measured pass-1 apex).
      if (narrow && use_empirical && !empirical_rt.empty())
      {
        const size_t n_emp = applyEmpiricalRT_(transition_exp, empirical_rt);
        OPENMS_LOG_INFO << "OpenDIAlyzer: empirical RT applied to " << n_emp << " library precursors." << std::endl;
      }

      // OpenSWATH narrows the extraction window to the one its calibration estimates from the
      // anchor residuals ("[Estimated] RT window applied: ..."), instead of keeping the wide
      // search window. Doing the same is what finally makes the candidate-peak sets comparable
      // (a 14x wider window yields far more candidates, so the reported top-5 diverge).
      double rt_win = narrow ? getDoubleOption_("rt_extraction_window_recal")
                             : getDoubleOption_("rt_extraction_window");
      // A CALIBRATION MAY ONLY NARROW -- the same rule the mass calibration now follows, and for the
      // same reason. The configured window is the user's assertion; measuring exists to discover the
      // run is TIGHTER than that. An estimate that comes back WIDER is reporting that its anchors did
      // not support a fit, and widening admits interference into every correlation sub-score.
      //
      // Measured, on this benchmark: the calibration widened 600 s -> 1435.1 s (61% of a 2333 s
      // gradient) from an anchor yield of 1.8% -- 70 anchors out of 3,897. Identical inputs and
      // classifier at 600 s gave 6,798 identifications; at 1435.1 s they gave 4,897. That one
      // parameter was worth 39% of the result, and the widening was never gated.
      const double min_yield = getDoubleOption_("rt_calib_min_yield_pct");
      const bool yield_ok = (calib_yield_pct_ < 0.0) || (calib_yield_pct_ >= min_yield);
      if (!narrow && use_estimated_window && calib_rt_window > 0.0 && yield_ok
          && calib_rt_window <= getDoubleOption_("rt_extraction_window"))
      {
        const double floor_win = getDoubleOption_("rt_extraction_window_min");
        const double est = calib_rt_window;
        rt_win = std::max(est, floor_win);
        OPENMS_LOG_INFO << "OpenDIAlyzer: using calibration-estimated RT window " << est
                        << " s (was " << getDoubleOption_("rt_extraction_window") << " s)." << std::endl;
        if (rt_win > est)
        {
          OPENMS_LOG_INFO << "OpenDIAlyzer: RAISED to " << rt_win
                          << " s by -rt_extraction_window_min." << std::endl;
        }
        if (floor_win <= 0.0)
        {
          OPENMS_LOG_INFO << "OpenDIAlyzer: no -rt_extraction_window_min floor set. The estimate is "
                          << "fitted to surviving anchors and understates the tail for unseen "
                          << "precursors; if pass 1 finds few confident peaks, that is the first "
                          << "thing to raise." << std::endl;
        }
      }
      if (!narrow && use_estimated_window && calib_rt_window > 0.0 && rt_win != calib_rt_window)
      {
        OPENMS_LOG_WARN << "OpenDIAlyzer: REJECTED the calibration's RT window estimate ("
                        << calib_rt_window << " s) -- "
                        << (calib_yield_pct_ >= 0.0 && calib_yield_pct_ < min_yield
                              ? "anchor yield " + std::to_string(calib_yield_pct_) + "% is below -rt_calib_min_yield_pct"
                              : std::string("it is WIDER than the configured window"))
                        << ". Keeping " << rt_win << " s. Calibration may only narrow." << std::endl;
      }
      const std::string osw = (p == passes) ? out
                                            : (out + ".pass" + std::to_string(p) +
                                               (parquet_out_ ? ".oswpq" : ".osw"));
      OPENMS_LOG_INFO << "[pass " << p << "/" << passes << (narrow ? " NARROW rt_win=" : " WIDE rt_win=")
                      << rt_win << "s] -> " << osw << std::endl;

      // The GLOBAL peak RSS lands here -- 186.28 GB on the 2026-08-01 profile -- and it was the one
      // phase with no timer, so the sampler could only attribute it to "(outside any phase)".
      // Named per pass, because pass 1 (wide window) and pass 2 (narrow) are different workloads.
      ExitCodes rc = EXECUTION_OK;
      {
        PhaseTimer pt(p == 1 ? "extract_pass1_wide" : "extract_pass2_narrow");
        rc = extractPass_(swath_maps, exp_meta, transition_exp,
                          use_native_trafo ? native_trafo : identity, rt_win, in, osw);
      }
      MemProbe::logAllocator(p == 1 ? "after extract pass1" : "after extract pass2");
      if (rc != EXECUTION_OK) { return rc; }

      if (p < passes)
      {
        int n_anchors = 0;
        double p95_resid = -1.0;
        // current run RT -> observed RT; also harvest empirical RTs + fit-residual diagnostic
        pass_map = recalibrate_(osw, n_anchors, use_empirical ? &empirical_rt : nullptr, &p95_resid);
        if (n_anchors < 20)
        {
          OPENMS_LOG_ERROR << "OpenDIAlyzer: recalibration produced too few anchors; aborting." << std::endl;
          return INCOMPATIBLE_INPUT_DATA;
        }
        // recalibrate_ fits library RT -> observed RT. On the native path the library is still
        // in iRT units, so that IS the iRT -> run map; invert it for OpenMS's convention and
        // keep the library untouched. On the fallback path pass_map is applied to compound.rt.
        if (use_native_trafo)
        {
          TransformationDescription t = pass_map;
          t.invert();
          native_trafo = t;
        }
        // Re-measure the run's mass accuracy on the precursors this pass actually found, so the
        // NEXT pass extracts at the instrument's real error instead of the configured guess.
        //
        // Placed HERE deliberately, and both parts of that matter: after the n_anchors >= 20 check,
        // because a pass that could not even fit an RT transform has no trustworthy confident set to
        // calibrate from; and after native_trafo is updated, because inferMassAccuracyPpm_ uses the
        // transform to locate each anchor's spectra and the pre-update value belongs to the previous
        // pass. An earlier version of this call had both wrong.
        // DIRECTION MATTERS, and it is not the transform handed to OpenMS. native_trafo is
        // run RT -> iRT (the convention OpenMS wants, line ~3262) and on this path the library is
        // deliberately kept in iRT units. inferMassAccuracyPpm_ starts from a LIBRARY rt and needs
        // the RUN rt at which to look for the peak, i.e. iRT -> run: the INVERSE.
        //
        // Applying native_trafo directly fed a run->iRT map an iRT value. The result is not a
        // retention time at all, so "nearest spectrum in time" picked an arbitrary spectrum and
        // "most intense peak within the window" then sampled interference -- which is exactly the
        // FLAT residual distribution the peakedness gate kept reporting (1.17 over 1066 anchors).
        // On the fallback path the library was already rescaled into run seconds, so identity is
        // still correct there.
        TransformationDescription mass_cal_trafo = identity;
        if (use_native_trafo)
        {
          mass_cal_trafo = native_trafo;
          mass_cal_trafo.invert();                 // run -> iRT  becomes  iRT -> run
        }
        calibrateMassFromPass_(swath_maps, transition_exp, mass_cal_trafo);
        // Diagnostic only: this is the IN-SAMPLE anchor residual (small); it is NOT the
        // predictive residual on unseen peptides, so it must NOT auto-size the window
        // (would collapse it and lose unseen precursors). Window sizing stays fixed until a
        // held-out predictive-residual estimate exists (backlog).
        OPENMS_LOG_INFO << "OpenDIAlyzer: pass-" << p << " calibration p95 anchor residual = "
                        << p95_resid << " s (diagnostic; window stays fixed)." << std::endl;
      }
    }
    if (parquet_out_)
    {
      // The assay library + scored features, written ONCE, from memory. This is the only time the
      // features touch disk in a parquet run.
      OpenSwathOSWParquetWriter pw;
      pw.write(out, transition_exp, pass_features_, run_id_, in, /*uis*/ false);
      OPENMS_LOG_INFO << "OpenDIAlyzer: wrote parquet bundle " << out << " ("
                      << pass_features_.size() << " features)." << std::endl;
    }
    const int ids = finalScore_(out);            // in-process LDA FDR on the final pass
    OPENMS_LOG_INFO << "OpenDIAlyzer: done -> " << out << " (" << ids << " target precursors at q<0.01)" << std::endl;
    return EXECUTION_OK;
  }

  // ponytail: one runnable check on the non-trivial fit (binning + apply, outlier reject,
  // PAVA monotonicity, residual reporting). Catches a broken isotonic fit on the cluster build.
  ExitCodes selftest_()
  {
    // (a) outlier rejection + accuracy: y=2x+10 with two gross outliers.
    std::vector<std::pair<double, double>> pts;
    for (int i = 0; i < 500; ++i) { pts.emplace_back(i, 2.0 * i + 10.0); }
    pts.emplace_back(250, 9999.0); pts.emplace_back(120, -9999.0);
    double p95 = -1.0;
    TransformationDescription td = fitTrafo_(pts, &p95);
    const double got = td.apply(300.0), want = 610.0;
    if (std::abs(got - want) > 5.0)
    {
      OPENMS_LOG_ERROR << "selftest FAILED: apply(300)=" << got << " expected ~" << want << std::endl;
      return UNEXPECTED_RESULT;
    }
    if (!(p95 >= 0.0) || p95 > 20.0)                          // clean line -> tiny residual
    {
      OPENMS_LOG_ERROR << "selftest FAILED: p95 residual=" << p95 << " (expected small, finite)" << std::endl;
      return UNEXPECTED_RESULT;
    }
    // (b) PAVA monotonicity: an underlying-monotone map with heavy noise that creates local
    // dips in the raw anchors. The fitted map MUST be non-decreasing across a sweep.
    std::vector<std::pair<double, double>> noisy;
    for (int i = 0; i < 600; ++i)
    {
      const double x = i;
      const double dip = ((i / 30) % 2 == 0) ? 120.0 : -120.0;   // alternating local bias
      noisy.emplace_back(x, 1.5 * x + 5.0 + dip);
    }
    TransformationDescription tn = fitTrafo_(std::move(noisy));
    double prev = tn.apply(0.0);
    for (int x = 5; x <= 595; x += 5)
    {
      const double cur = tn.apply(static_cast<double>(x));
      if (cur < prev - 1e-6)
      {
        OPENMS_LOG_ERROR << "selftest FAILED: fitted map not monotone at x=" << x
                         << " (" << cur << " < " << prev << ")" << std::endl;
        return UNEXPECTED_RESULT;
      }
      prev = cur;
    }
    // (c) degenerate/small-anchor sets must NOT throw (OpenMS interpolated needs >=3 unique x)
    // and must stay monotone. These are the cases codex flagged as the runtime blocker.
    {
      // 25 anchors (NB=2), a clean line: must fit and recover the slope, not throw.
      std::vector<std::pair<double, double>> few;
      for (int i = 0; i < 25; ++i) { few.emplace_back(i, 3.0 * i + 7.0); }
      double p = -1.0;
      TransformationDescription t = fitTrafo_(few, &p);
      if (std::abs(t.apply(10.0) - 37.0) > 3.0)
      {
        OPENMS_LOG_ERROR << "selftest FAILED: 25-anchor fit apply(10)=" << t.apply(10.0) << " expected ~37" << std::endl;
        return UNEXPECTED_RESULT;
      }
    }
    {
      // exactly two unique x -> raw-x fallback (NB=1 collapses the single bin). Per-group
      // MEDIAN y must win: cluster medians are 205 and 605, so the collinear midpoint gives
      // apply(200)=405 and apply(100)=205. A min-y fallback would give 200/400 -- assert
      // tightly enough (within 2 s) to catch that regression.
      std::vector<std::pair<double, double>> two{{100, 200}, {100, 210}, {100, 205}, {300, 600}, {300, 610}, {300, 605}};
      TransformationDescription t = fitTrafo_(two);
      if (std::abs(t.apply(100.0) - 205.0) > 2.0 || std::abs(t.apply(200.0) - 405.0) > 2.0
          || std::abs(t.apply(300.0) - 605.0) > 2.0)
      {
        OPENMS_LOG_ERROR << "selftest FAILED: two-unique-x fit apply(100/200/300)=" << t.apply(100.0)
                         << "/" << t.apply(200.0) << "/" << t.apply(300.0) << " expected ~205/405/605" << std::endl;
        return UNEXPECTED_RESULT;
      }
    }
    {
      // single unique x -> identity (no map possible); must not throw and apply(x)=x.
      std::vector<std::pair<double, double>> one{{50, 111}, {50, 112}, {50, 110}};
      TransformationDescription to = fitTrafo_(one);
      if (std::abs(to.apply(42.0) - 42.0) > 1e-6)
      {
        OPENMS_LOG_ERROR << "selftest FAILED: single-x fallback not identity (apply(42)=" << to.apply(42.0) << ")" << std::endl;
        return UNEXPECTED_RESULT;
      }
      // empty input -> identity + p95 == 0.
      double pe = -1.0;
      TransformationDescription te = fitTrafo_({}, &pe);
      if (std::abs(te.apply(42.0) - 42.0) > 1e-6 || pe != 0.0)
      {
        OPENMS_LOG_ERROR << "selftest FAILED: empty-input fallback not identity (apply(42)=" << te.apply(42.0) << ")" << std::endl;
        return UNEXPECTED_RESULT;
      }
    }
    {
      // LOESS recalibration. Measured against the TRUE warp at HELD-OUT points, not against the
      // anchors: anchor residual rewards overfitting noise, so "beats the other fit on the anchors"
      // is the wrong criterion and would pass a map that has simply memorised the scatter.
      //
      // Setup is the real one: FEW anchors (150 -- this run produced 19 to 136), NOISY, and a warp
      // that is monotone but curved, which is what a predicted library looks like against a real
      // gradient. Deterministic pseudo-noise so the check is reproducible.
      auto warp = [](double x) { return x + 300.0 * std::sin(3.14159265358979 * x / 3000.0); };
      std::vector<std::pair<double, double>> noisy;
      for (int i = 0; i < 150; ++i)
      {
        const double x = 20.0 * i;                               // 0 .. 2980
        const double jitter = 60.0 * std::sin(12.9898 * i) * std::cos(4.1414 * i);   // ~+/-60 s
        noisy.emplace_back(x, warp(x) + jitter);
      }
      TransformationDescription t_bin = fitTrafo_(noisy, nullptr, 0.0);
      TransformationDescription t_loess = fitTrafo_(noisy, nullptr, 0.3);
      double err_bin = 0.0, err_loess = 0.0;
      int n_eval = 0;
      for (double x = 10.0; x <= 2970.0; x += 20.0)              // midpoints: never an anchor x
      {
        err_bin += std::abs(t_bin.apply(x) - warp(x));
        err_loess += std::abs(t_loess.apply(x) - warp(x));
        ++n_eval;
      }
      err_bin /= n_eval; err_loess /= n_eval;
      if (!(err_loess < err_bin))
      {
        OPENMS_LOG_ERROR << "selftest FAILED: LOESS mean |error vs TRUE warp| " << err_loess
                         << " s did not beat the binned fit " << err_bin << " s on held-out points"
                         << std::endl;
        return UNEXPECTED_RESULT;
      }
      // NOTE: no monotonicity assert here. PAVA + piecewise-linear interpolation is monotone BY
      // CONSTRUCTION, so such a check passes even if LOESS returned a sine wave of crossings -- it
      // would verify that PAVA exists, which we already know. The checks below test things that
      // can actually break.

      // ROBUSTNESS. The anchors are pass-1 peak groups hunted in windows narrower than the RT
      // error, so gross outliers are guaranteed. A single tricube pass has breakdown ~1/q; the
      // binned MEDIAN it replaces has 50%. This asserts the bisquare loop actually recovers that:
      // with 5% of anchors thrown +2000 s off, LOESS must still not be worse than the binned fit.
      std::vector<std::pair<double, double>> dirty = noisy;
      for (size_t i = 0; i < dirty.size(); i += 20) { dirty[i].second += 2000.0; }   // 5% gross
      TransformationDescription d_bin = fitTrafo_(dirty, nullptr, 0.0);
      TransformationDescription d_loess = fitTrafo_(dirty, nullptr, 0.3);
      double d_eb = 0.0, d_el = 0.0; int d_n = 0;
      for (double x = 10.0; x <= 2970.0; x += 20.0)
      {
        d_eb += std::abs(d_bin.apply(x) - warp(x));
        d_el += std::abs(d_loess.apply(x) - warp(x));
        ++d_n;
      }
      d_eb /= d_n; d_el /= d_n;
      if (d_el > d_eb * 1.10)
      {
        OPENMS_LOG_ERROR << "selftest FAILED: with 5% gross outliers LOESS (" << d_el
                         << " s) is materially worse than the binned median (" << d_eb
                         << " s) -- the bisquare robustness loop is not working" << std::endl;
        return UNEXPECTED_RESULT;
      }

      // GATING. Below 50 anchors LOESS must fall back to the binned path rather than silently
      // no-op: at n=19 the binning yields ONE bin, which is how the first version came to be
      // inert at exactly the anchor count this data produces.
      std::vector<std::pair<double, double>> few;
      for (int i = 0; i < 19; ++i)
      {
        const double x = 150.0 * i;
        few.emplace_back(x, warp(x));
      }
      TransformationDescription f_on = fitTrafo_(few, nullptr, 0.3);
      TransformationDescription f_off = fitTrafo_(few, nullptr, 0.0);
      for (double x = 100.0; x <= 2600.0; x += 100.0)
      {
        if (std::abs(f_on.apply(x) - f_off.apply(x)) > 1e-6)
        {
          OPENMS_LOG_ERROR << "selftest FAILED: at 19 anchors the LOESS path should fall back to "
                           << "the binned fit, but the maps differ at x=" << x << std::endl;
          return UNEXPECTED_RESULT;
        }
      }
    }
    {
      // m/z calibration statistics. Three properties, each of which can actually break.
      auto lcg = [](uint64_t& s) { s = s * 6364136223846793005ULL + 1442695040888963407ULL;
                                   return static_cast<double>((s >> 11) & ((1ULL << 53) - 1)) / 9007199254740992.0; };
      auto gauss = [&](uint64_t& s, double mu, double sd) {
        const double u1 = std::max(1e-12, lcg(s)), u2 = lcg(s);
        return mu + sd * std::sqrt(-2.0 * std::log(u1)) * std::cos(6.283185307179586 * u2);
      };

      // (1) RECOVERY: signal at a known offset, buried in 60% uniform contamination -- the regime
      // this estimator exists for. HSM must find the offset; a mean would be dragged to ~0.
      uint64_t seed = 12345;
      std::vector<double> mix;
      for (int i = 0; i < 400; ++i) { mix.push_back(gauss(seed, 4.0, 1.0)); }          // 40% signal
      for (int i = 0; i < 600; ++i) { mix.push_back(-50.0 + 100.0 * lcg(seed)); }      // 60% uniform
      std::sort(mix.begin(), mix.end());
      const double mode = halfSampleMode_(mix);
      if (std::abs(mode - 4.0) > 1.0)
      {
        OPENMS_LOG_ERROR << "selftest FAILED: half-sample mode " << mode
                         << " did not recover the +4 ppm offset under 60% contamination" << std::endl;
        return UNEXPECTED_RESULT;
      }
      const double sc = localScaleAboutMode_(mix, mode, 50.0);
      if (!(sc > 0.3 && sc < 4.0))
      {
        OPENMS_LOG_ERROR << "selftest FAILED: shorth scale " << sc << " implausible for sd=1 signal"
                         << std::endl;
        return UNEXPECTED_RESULT;
      }

      // (2) THE SAFETY PROPERTY: pure uniform noise must be REJECTED. A uniform sample still has a
      // perfectly computable mode and width -- acting on them is how the window widens to death.
      std::vector<double> noise;
      for (int i = 0; i < 1000; ++i) { noise.push_back(std::abs(-50.0 + 100.0 * lcg(seed))); }
      const double pk_noise = peakednessRatio_(noise, 50.0);
      if (pk_noise >= 3.0)
      {
        OPENMS_LOG_ERROR << "selftest FAILED: uniform noise scored peakedness " << pk_noise
                         << " (>=3) -- the guard would accept a meaningless calibration" << std::endl;
        return UNEXPECTED_RESULT;
      }

      // (3) ...and a real error distribution must be ACCEPTED, or the guard blocks everything.
      std::vector<double> real;
      for (int i = 0; i < 1000; ++i) { real.push_back(std::abs(gauss(seed, 0.0, 2.0))); }
      const double pk_real = peakednessRatio_(real, 50.0);
      if (pk_real < 3.0)
      {
        OPENMS_LOG_ERROR << "selftest FAILED: a genuine sd=2 ppm error distribution scored "
                         << "peakedness " << pk_real << " (<3) -- the guard is too strict"
                         << std::endl;
        return UNEXPECTED_RESULT;
      }
    }
    {
      // Prefilter decoy pairing. The regression this guards: matching a SHUFFLED decoy to its
      // target by sequence keeps only accidental collisions (23 of ~3.5M on the real library),
      // which slips past the n_dec == 0 guard and silently destroys the target:decoy ratio the
      // FDR depends on. Pairing is by id prefix only -- sequences must never be consulted.
      const std::unordered_set<std::string> kept{"PEPTIDEA_2", "PEPTIDEK_3"};
      struct { const char* id; bool want; } cases[] = {
        {"DECOY_PEPTIDEA_2", true },   // partner survived
        {"DECOY_PEPTIDEK_3", true },
        {"DECOY_PEPTIDEZ_2", false},   // partner did not survive
        {"PEPTIDEA_2",       false},   // untagged id is not a decoy of anything
        {"DECOY_",           false},   // tag alone -> empty target id
        {"DEC_PEPTIDEA_2",   false},   // wrong tag
        {"DECOY_DECOY_PEPTIDEA_2", false}, // nested tag: strip ONCE, decoy-of-decoy is not a partner
        {"XDECOY_PEPTIDEA_2", false},  // tag must be a prefix, not merely contained
      };
      for (const auto& c : cases)
      {
        if (decoyPartnerKept_(c.id, "DECOY_", kept) != c.want)
        {
          OPENMS_LOG_ERROR << "selftest FAILED: decoyPartnerKept_(" << c.id << ") != " << c.want << std::endl;
          return UNEXPECTED_RESULT;
        }
      }
      // An empty tag must not make every decoy look paired.
      if (decoyPartnerKept_("DECOY_PEPTIDEA_2", "", kept))
      {
        OPENMS_LOG_ERROR << "selftest FAILED: empty decoy_tag paired a decoy" << std::endl;
        return UNEXPECTED_RESULT;
      }
    }
    OPENMS_LOG_INFO << "OpenDIAlyzer selftest OK: apply(300)=" << got << " (~610), p95_resid="
                    << p95 << ", isotonic monotone, degenerate/small-anchor fits robust." << std::endl;
    return EXECUTION_OK;
  }
};

int main(int argc, const char** argv)
{
  TOPPOpenDIAlyzer tool;
  return tool.main(argc, argv);
}
