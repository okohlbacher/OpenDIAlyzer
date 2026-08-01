# OpenDIALibGen → OpenDIAlyzer merger: design and deep self-review

Status: design plan only. This document is based on the current working tree as read on
2026-07-27; several of the cited files are currently modified or untracked, so the
citations refer to the working-tree versions, not necessarily `HEAD`. No source change is
part of this plan.

## Executive decision

Proceed with the merger, but do not treat these three capabilities as one milestone:

1. **On-demand, in-process library construction:** feasible and worth doing. Extract the
   generator core into a reusable `LibraryBuilder`, return
   `OpenSwath::LightTargetedExperiment` in memory, and retain the PeptDeep inference
   objects for the duration of the run.
2. **Per-run recalibration and re-prediction:** feasible with the current inference-only
   build, provided “recalibration” means fitting a run-specific calibration layer around
   immutable PeptDeep outputs. Implement RT first; add MS2 and IM only after the engine
   exposes trustworthy labels for them.
3. **Neural-weight fine-tuning:** **not feasible in the current build**. The present
   OpenMS wrappers expose inference sessions and prediction methods, not a training
   session, optimizer, gradients, checkpoints, or mutable weights. This phase is a
   separate research and build project and is a no-go until a small RT-only proof of
   concept passes held-out and FDR-safety gates.

The retirement of `OpenDIALibGen` should follow library parity and integrated-search
parity. It must not be coupled to the success of neural fine-tuning.

## 1. Current-state map

### 1.1 What `OpenDIALibGen` actually does

The comments at the top of `src/opendialibgen.cpp` are stale: they say FASTA digestion,
modifications, CCS, assay refinement, and decoys are out of scope
(`src/opendialibgen.cpp:1-14`), but the current implementation contains all of them.
The working pipeline is:

| Stage | Current implementation | Current representation |
|---|---|---|
| Peptide-list input | Parses `PeptideSequence`, `PrecursorCharge`, and optional `ModifiedPeptideSequence` | `std::vector<Precursor>` (`src/opendialibgen.cpp:100-106`, `src/opendialibgen.cpp:319-385`) |
| FASTA input | `FASTAFile` → tryptic `ProteaseDigestion`; enumerates charge states and fixed/variable peptidoforms | `std::vector<Precursor>` (`src/opendialibgen.cpp:395-436`) |
| RT prediction | `PeptDeepRTInference::predictRT` | `std::vector<float>` (`src/opendialibgen.cpp:189-192`) |
| MS2 prediction | `PeptDeepMS2Inference::predictMS2` with charge, NCE, and instrument | `std::vector<std::vector<float>>` (`src/opendialibgen.cpp:189-194`) |
| CCS/IM prediction | `PeptDeepCCSInference::predictCCS`, followed by CCS → 1/K0 conversion | `std::vector<float>`, then one `double` per precursor (`src/opendialibgen.cpp:191-197`, `src/opendialibgen.cpp:223-227`) |
| Assay assembly | Theoretical b/y fragment m/z values are paired to predicted intensity channels | A single in-memory TSV `std::string` (`src/opendialibgen.cpp:169-303`) |
| Assay refinement | TSV string is written to a temporary file, parsed into `TargetedExperiment`, then reannotated, m/z-filtered, and reduced to 6–12 detecting transitions | `TargetedExperiment` (`src/opendialibgen.cpp:588-624`) |
| Decoys | `MRMDecoy::generateDecoys` creates shuffled decoys and appends them to the heavyweight experiment | `TargetedExperiment` (`src/opendialibgen.cpp:626-635`) |
| Output | Converts the heavyweight experiment back to OpenSWATH TSV | File on disk (`src/opendialibgen.cpp:637`) |

There is no `LightTargetedExperiment` output from the generator today. `generate()`
returns a TSV string; `refine_and_decoy()` creates a temporary TSV and a heavyweight
`TargetedExperiment`; `main()` ultimately writes TSV (`src/opendialibgen.cpp:171-303`,
`src/opendialibgen.cpp:598-646`, `src/opendialibgen.cpp:704-724`).

The two input modes also have different current semantics. A FASTA build is refined and
decoyed by default, while a peptide-list build returns raw targets only because
`main()` exits early whenever `fasta_path.empty()` (`src/opendialibgen.cpp:708-718`).
The merged search path must remove that accidental difference: a search library needs
the same refinement and target/decoy policy regardless of whether its precursor list
came from FASTA or TSV.

### 1.2 Predictor classes, APIs, and lifetime

The three predictors are OpenMS classes backed by `ONNXPredictorBase`:

| Predictor | Public prediction API | Owned state |
|---|---|---|
| `OpenMS::PeptDeepRTInference` | `predictRT(vector<string>)` and `predictRT(vector<AASequence>)` → `vector<float>` | Private `ONNXPredictorBase model_` and batch size (`ext/OpenMS/src/openms/include/OpenMS/ML/PEPTDEEP/PeptDeepRTInference.h:20-51`) |
| `OpenMS::PeptDeepMS2Inference` | `predictMS2(peptides, charges, nces, instruments)` → `vector<vector<float>>` | Private `ONNXPredictorBase model_` and batch size (`ext/OpenMS/src/openms/include/OpenMS/ML/PEPTDEEP/PeptDeepMS2Inference.h:25-60`) |
| `OpenMS::PeptDeepCCSInference` | `predictCCS(peptides, charges)` → `vector<float>` | Private `ONNXPredictorBase model_` and batch size (`ext/OpenMS/src/openms/include/OpenMS/ML/PEPTDEEP/PeptDeepCCSInference.h:21-50`) |

`ONNXPredictorBase` owns `Ort::SessionOptions`, `Ort::Session`, and CPU
`Ort::MemoryInfo`; it exposes session inspection and execution support, but no optimizer,
checkpoint, gradient, training graph, or weight-update API
(`ext/OpenMS/src/openms/include/OpenMS/ML/ONNX/ONNXPredictorBase.h:27-69`).
The constructor creates an optimized inference `Ort::Session`
(`ext/OpenMS/src/openms/source/ML/ONNX/ONNXPredictorBase.cpp:21-48`).

`OpenDIALibGen::generate()` constructs all three predictor objects as local variables.
They are destroyed when that one call returns (`src/opendialibgen.cpp:171-197`). This is
the exact lifetime problem the merger must fix.

Prediction conditions are also embedded in source rather than represented as runtime
configuration: NCE is `35.0` and instrument index is `2` (timsTOF)
(`src/opendialibgen.cpp:52-76`). The merged module must put model paths, instrument,
collision energy, ONNX threads, and prediction batch size into a recorded
`PredictionOptions` object. Otherwise the engine cannot reproduce or deliberately
change a run’s initial library.

### 1.3 What `OpenDIAlyzer` does now

`TOPPOpenDIAlyzer` derives directly from `TOPPBase`, loads the DIA run and library using
libOpenMS, and calls `OpenSwathWorkflow::performExtraction` in process
(`src/opendialyzer.cpp:1-11`, `src/opendialyzer.cpp:64-80`,
`src/opendialyzer.cpp:404-460`).

Its library boundary is already the desired search-engine type:

```cpp
OpenSwath::LightTargetedExperiment loadLibrary_(const std::string& tr);
ExitCodes extractPass_(...,
  const OpenSwath::LightTargetedExperiment& transition_exp,
  ...);
```

`loadLibrary_()` converts TSV, PQP, or OSWPQ into a
`LightTargetedExperiment` (`src/opendialyzer.cpp:99-123`). `extractPass_()` consumes that
light experiment directly (`src/opendialyzer.cpp:405-409`) and writes the library tables
into the output OSW through
`TransitionPQPFile::convertLightTargetedExperimentToPQP`
(`src/opendialyzer.cpp:411-412`). Therefore an intermediate transition TSV is not
required by any downstream extraction API.

The current two-pass flow is:

1. Require `-in`, `-tr`, and `-out` (`src/opendialyzer.cpp:469-477`).
2. Load the library and DIA run once (`src/opendialyzer.cpp:487-496`).
3. Bootstrap predicted RT to the run’s RT range (`src/opendialyzer.cpp:498-530`).
4. Mutate every light compound’s RT in place, extract pass 1, fit a new RT transform,
   mutate the same library again, and extract pass 2
   (`src/opendialyzer.cpp:182-188`, `src/opendialyzer.cpp:533-559`).
5. Score the final OSW using the in-process LDA (`src/opendialyzer.cpp:561-562`).

`odia_lda.h` is the relevant anchor/FDR component. It provides deterministic,
group-cross-validated semi-supervised LDA and q-values, with a fixed default seed of 42
(`src/odia_lda.h:50-72`, `src/odia_lda.h:249-294`). `odia_score.h` is a standalone,
header-only scoring kernel and is not included by `opendialyzer.cpp`; today it is built
only into `odia-score-test` (`src/odia_score.h:1-16`, `CMakeLists.txt:38-41`). It may
later help derive MS2 training labels, but it is not part of the library-merger seam.

One important current-state mismatch must be corrected before model adaptation:
`recalibrate_()` is documented as selecting q < 0.01 anchors, but it actually selects up
to the top 2,000 target groups by d-score even when the honest q < 0.01 set is empty
(`src/opendialyzer.cpp:313-355`). That lenient set can be acceptable for a robust RT
mapping, but it is not safe training data for neural fine-tuning.

### 1.4 Exact seam

The join point is the initialization of `transition_exp`:

```cpp
// current
OpenSwath::LightTargetedExperiment transition_exp = loadLibrary_(tr);

// target
LibrarySession library_session =
  tr.empty() ? library_builder.build(source_options)
             : LibrarySession::fromExternal(loadLibrary_(tr), adaptation_options);

OpenSwath::LightTargetedExperiment& transition_exp =
  library_session.state().searchLibrary();
```

This replaces only `src/opendialyzer.cpp:487-488`; `extractPass_()` already accepts the
resulting type. The `LibrarySession`—not a temporary `generate()` call—must stay alive
through both extraction passes so that its predictor bundle, canonical precursor
catalog, base predictions, target/decoy mapping, and mutable search library also stay
alive.

## 2. Target architecture

### 2.1 Architectural decisions

1. **One reusable module, two temporary front ends.** Move generation out of
   `opendialibgen.cpp` into `odia_library_builder.*`. During migration, both executables
   call it. After parity, remove the `OpenDIALibGen` front end.
2. **The search-facing result is always `LightTargetedExperiment`.** No TSV is required
   between prediction and extraction.
3. **Keep an immutable canonical prediction state.** Never repeatedly transform the
   already transformed `LightCompound::rt`. Every pass is materialized from base model
   outputs plus the current run adaptation. This prevents compounding transforms and
   makes rollback possible.
4. **Keep predictor ownership separate from training capability.** A live inference
   session is useful for re-prediction but does not imply that weights are trainable.
5. **Use stable string IDs across passes.** OSW integer precursor IDs are local to an
   output database. Adaptation and library replacement must key on `TRAML_ID` /
   `LightCompound::id`.
6. **Adapt targets and decoys symmetrically.** Any RT, MS2, or IM transformation applied
   to targets must also be applied to decoys. A target-only adaptive model would change
   the null distribution and invalidate FDR.
7. **Do not retain two complete library representations longer than necessary.**
   Initially use `TargetedExperiment` for behavior parity, convert once to the light
   representation, then release it before loading the full DIA run.

### 2.2 New and changed files

| File | Change |
|---|---|
| `src/odia_library_builder.h` | New public data types and `LibraryBuilder`,
  `PeptDeepPredictors`, `LibraryState`, and `LibrarySession` declarations |
| `src/odia_library_builder.cpp` | New digestion, peptidoform enumeration, chunked
  prediction, transition assembly, assay refinement, decoy generation, heavy→light
  conversion, indexing, and optional TSV/PQP export |
| `src/odia_run_adaptation.h` | New `AnchorSet`, `AnchorSelector`,
  `RunLibraryAdapter`, calibration reports, and future trainable-backend interface |
| `src/odia_run_adaptation.cpp` | New RT calibration/rebuild logic first; MS2/IM
  calibration and training backends only in later phases |
| `src/opendialyzer.cpp` | Optional library source CLI, dispatch, persistent
  `LibrarySession`, stable-ID anchor extraction, pass-specific library materialization,
  and adaptation between passes |
| `src/opendialibgen.cpp` | During migration, reduce to a compatibility CLI that calls
  `LibraryBuilder` and exports TSV; delete after retirement |
| `CMakeLists.txt` | Add shared module target, enforce one ONNX-enabled OpenMS, link both
  executables during migration, add model/runtime checks and tests, then remove the
  legacy target |
| `src/odia_library_builder_test.cpp` | New deterministic assembly and semantic-parity
  tests |
| `src/odia_run_adaptation_test.cpp` | New stable-ID, no-compounding, target/decoy
  symmetry, calibration, rollback, and unsupported-training tests |

Do not put the new reusable code in another `main()` translation unit or expose anonymous
namespace types from `opendialibgen.cpp`; the current `Precursor`, `DigestOpts`,
`ModPolicy`, and free functions are private to that file
(`src/opendialibgen.cpp:50-167`, `src/opendialibgen.cpp:305-437`).

### 2.3 Proposed types and function signatures

The names may be adjusted to project style, but the ownership and boundaries should stay
as follows:

```cpp
// src/odia_library_builder.h
namespace odia
{
struct PrecursorSpec
{
  OpenMS::AASequence sequence;       // exact peptidoform
  std::string unmodified_sequence;
  std::string protein_id;
  std::string compound_id;           // stable across passes
  int charge;
  bool decoy;
  std::string paired_compound_id;    // target↔decoy
};

struct DigestionOptions
{
  std::string enzyme{"Trypsin"};
  int missed_cleavages{1};
  int min_length{7}, max_length{30};
  int min_charge{2}, max_charge{3};
};

struct ModificationOptions
{
  bool carbamidomethyl_c{true};
  bool oxidation_m{true};
  bool protein_nterm_acetyl{true};
  int max_variable_mods{2};
};

struct PredictionOptions
{
  std::filesystem::path model_dir;
  float nce{35.0f};
  std::int64_t instrument_index{2};
  int intra_op_threads{1};
  std::size_t batch_size{500};
};

struct AssayOptions
{
  double product_mz_min{150.0}, product_mz_max{2000.0};
  int min_transitions{6}, max_transitions{12};
  bool generate_decoys{true};
  std::uint64_t decoy_seed{42};
};

struct PredictionBatch
{
  std::vector<float> rt;
  std::vector<std::vector<float>> ms2;
  std::vector<float> ccs;
};

class PeptDeepPredictors final
{
public:
  PeptDeepPredictors(const PredictionOptions&);
  PredictionBatch predict(std::span<const PrecursorSpec>,
                          PredictionMask fields = PredictionMask::All);

  OpenMS::PeptDeepRTInference& rt();
  OpenMS::PeptDeepMS2Inference& ms2();
  OpenMS::PeptDeepCCSInference& ccs();

  bool supportsWeightTraining(PredictionMask) const noexcept { return false; }

private:
  // unique_ptr is deliberate: the OpenMS wrappers own non-copyable sessions.
  std::unique_ptr<OpenMS::PeptDeepRTInference> rt_;
  std::unique_ptr<OpenMS::PeptDeepMS2Inference> ms2_;
  std::unique_ptr<OpenMS::PeptDeepCCSInference> ccs_;
};

struct BasePrediction
{
  double model_rt;
  double one_over_k0;
  // Retain enough unpruned fragment predictions to re-rank an assay after MS2
  // adaptation; otherwise a previously pruned transition cannot be restored.
  std::vector<PredictedFragment> fragments;
};

class LibraryState final
{
public:
  OpenSwath::LightTargetedExperiment& searchLibrary();
  const OpenSwath::LightTargetedExperiment& searchLibrary() const;

  const PrecursorSpec& precursor(std::string_view compound_id) const;
  const BasePrediction& basePrediction(std::string_view compound_id) const;
  void restoreBasePredictions(PredictionMask);
  void replaceAssays(std::span<const RebuiltAssay>);
  void rebuildIndex(); // mandatory after transition/compound vector replacement

private:
  OpenSwath::LightTargetedExperiment library_;
  std::vector<PrecursorSpec> precursor_catalog_;
  std::vector<BasePrediction> base_predictions_;
  std::unordered_map<std::string, std::size_t> compound_index_;
  std::unordered_map<std::string, std::vector<std::size_t>> transition_index_;
};

struct LibraryBuildResult
{
  LibraryState state;
  std::shared_ptr<PeptDeepPredictors> predictors;
  LibraryBuildReport report;
};

class LibraryBuilder final
{
public:
  LibraryBuilder(std::shared_ptr<PeptDeepPredictors>,
                 DigestionOptions,
                 ModificationOptions,
                 AssayOptions);

  LibraryBuildResult buildFromFasta(const std::filesystem::path&);
  LibraryBuildResult buildFromPeptideList(const std::filesystem::path&);

  static void writeTSV(const LibraryState&, const std::filesystem::path&);
  static void writePQP(const LibraryState&, const std::filesystem::path&);
};

class LibrarySession final
{
public:
  static LibrarySession generated(LibraryBuildResult&&);
  static LibrarySession external(OpenSwath::LightTargetedExperiment&&,
                                 std::shared_ptr<PeptDeepPredictors> = {});
  LibraryState& state();
  PeptDeepPredictors* predictors() noexcept;
};
} // namespace odia
```

For the first implementation, `LibraryBuilder` should assemble a
`TargetedExperiment` directly, run the same heavyweight `MRMAssay` and `MRMDecoy`
operations as today, call
`OpenSwathDataAccessHelper::convertTargetedExp(const TargetedExperiment&,
LightTargetedExperiment&)`, and immediately release the heavyweight object. That
conversion API already exists
(`ext/OpenMS/src/openms/include/OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/DataAccessHelper.h:47-57`)
and copies proteins, compounds, transitions, intensity, charge, decoy state, and
precursor IM (`ext/OpenMS/src/openms/source/ANALYSIS/OPENSWATH/DATAACCESS/DataAccessHelper.cpp:103-180`).

This is the smallest parity-preserving removal of the TSV round-trip. A later memory
optimization can assemble light structures directly and use the already available
`MRMAssay::reannotateTransitionsLight`, `restrictTransitionsLight`, and
`detectingTransitionsLight`
(`ext/OpenMS/src/openms/include/OpenMS/ANALYSIS/OPENSWATH/MRMAssay.h:189-245`) plus
`MRMDecoy::generateDecoysLight`
(`ext/OpenMS/src/openms/include/OpenMS/ANALYSIS/OPENSWATH/MRMDecoy.h:106-144`).
Do not make that representation switch part of the first parity milestone.

The adaptation API should explicitly distinguish calibration from training:

```cpp
// src/odia_run_adaptation.h
namespace odia
{
struct Anchor
{
  std::string compound_id;
  double qvalue;
  double model_rt;
  double observed_rt;
  std::optional<double> model_im;
  std::optional<double> observed_im;
  std::vector<float> observed_fragment_intensities;
};

struct AnchorSelectionOptions
{
  double max_qvalue{0.01};
  std::size_t min_rt_anchors{100};
  std::size_t max_rt_anchors{5000};
  unsigned split_seed{42};
};

class AnchorSelector final
{
public:
  AnchorSet fromPass1OSW(const std::filesystem::path&,
                        const LibraryState&,
                        const AnchorSelectionOptions&) const;
};

class RunLibraryAdapter final
{
public:
  RunLibraryAdapter(LibraryState&, PeptDeepPredictors*);

  AdaptationReport calibrateRT(const AnchorSet&);
  AdaptationReport calibrateMS2(const AnchorSet&); // later
  AdaptationReport calibrateIM(const AnchorSet&);  // later

  AdaptationReport repredict(std::span<const std::string> compound_ids,
                             PredictionMask);

  FineTuneReport fineTuneRT(const AnchorSet&, const FineTuneOptions&);
  FineTuneReport fineTuneMS2(const AnchorSet&, const FineTuneOptions&);
  FineTuneReport fineTuneIM(const AnchorSet&, const FineTuneOptions&);
};
} // namespace odia
```

Until a trainable backend is linked, all `fineTune*()` methods must return a typed
`UnsupportedBackend` result. They must not silently relabel spline fitting as neural
fine-tuning.

### 2.4 Library-state invariants

The following invariants are required for correctness:

- `LightCompound::id` and `LightTransition::transition_name` remain stable across passes.
- Every target and decoy compound can be mapped back to an exact `AASequence`, charge,
  model inputs, base prediction, and paired target/decoy.
- Base RT, MS2, and CCS predictions are immutable. A pass-specific value is always
  derived from base prediction + current run adaptation.
- RT rebuild updates `LightCompound::rt`.
- IM rebuild updates both `LightCompound::drift_time` and each associated
  `LightTransition::precursor_im`; the existing heavy→light conversion does the same
  propagation (`DataAccessHelper.cpp:139-144`).
- MS2 rebuild can change both intensities and the selected 6–12 transitions. Therefore
  the state must retain the unpruned predicted fragments or be able to regenerate them;
  keeping only the already-pruned light transitions is insufficient.
- Transition and compound indexes are rebuilt after vector replacement. No pointer or
  index into those vectors survives a rebuild.
- Adaptation is applied to targets and decoys using the same learned mapping/backend.
- Shipped model files are immutable. Any future fine-tuned checkpoint/model is
  run-scoped and written only under the run’s temporary directory.
- Every output records source type, FASTA/peptide-list hash or input-library hash, model
  hashes, digestion/modification/prediction options, adaptation mode, seed, anchor
  counts, and fallback/rollback reason.

### 2.5 New end-to-end flow

```text
validate CLI
  |
  +-- -tr ---------------------------------> load LightTargetedExperiment
  |                                           reconstruct adaptation catalog where possible
  |
  +-- -fasta / -peptides
        -> construct one live PeptDeepPredictors bundle
        -> digest/read precursor list
        -> predict RT/MS2/CCS in bounded chunks
        -> assemble/refine/decoy in memory
        -> convert once to LightTargetedExperiment
        -> retain predictor bundle + canonical precursor/base-prediction state
  |
load DIA run once
  |
materialize pass-1 RT from canonical prediction + bootstrap mapping
  -> extract pass 1
  -> in-process group-CV LDA
  -> select stable-ID confident target anchors
  -> fit run RT calibration
  -> [optional future] fine-tune trainable RT backend on training anchors
  -> validate on held-out anchors; rollback on failure
  -> re-predict all affected target and decoy compounds
  -> rebuild affected RT/MS2/IM library fields from canonical state
  -> extract pass 2
  -> final in-process LDA/FDR
```

The generated-library branch should build before the full DIA run is loaded so the
temporary heavyweight library and prediction tensors do not coexist with all in-memory
SWATH maps. Instrument/NCE auto-detection should use a lightweight metadata read where
possible; until that exists, use explicit CLI values with the current timsTOF/NCE 35
defaults. Do not load the full run merely to infer those two fields and then hold it
during library construction.

### 2.6 Pass-1 anchors and stable IDs

`loadOswScores_()` currently selects `p.ID` and uses that database-local integer as the
group key (`src/opendialyzer.cpp:226-265`). Add `p.TRAML_ID` to the query and retain it as
`OswRows::compound_id`. Continue using integer `p.ID` for LDA grouping inside one OSW,
but use the string ID to join an accepted anchor back to `LibraryState`.

Use two anchor policies:

- **RT calibration:** q < 0.01 when enough anchors exist; an explicitly named
  `robust_top_n_fallback` may use high d-score targets for the existing median-spline
  calibration only. Record when the fallback is used.
- **Any weight update:** no top-N fallback. Require q-controlled anchors, a minimum
  count, train/validation separation by precursor (preferably by unmodified sequence or
  protein family), and a held-out improvement gate.

Anchor selection must continue to exclude RT-derived scores when producing RT anchors;
the current code already excludes `VAR_NORM_RT_SCORE`
(`src/opendialyzer.cpp:200-220`). MS2 or IM adaptation needs analogous feature-exclusion
rules so that the labels are not selected mainly by the same agreement score being
optimized.

### 2.7 What “re-predict / fine-tune” means by phase

**Phase 2, feasible now:**

- Keep the three inference sessions alive.
- Fit a monotone run-specific RT calibration from base PeptDeep RT to observed seconds.
- Materialize all target and decoy RT values from base prediction + calibration.
- Optionally fit small, explicit calibration heads for CCS→1/K0 and MS2 intensity only
  after suitable labels exist.
- Re-run inference when model inputs or a future backend state changes; do not re-run it
  merely as ceremony when the immutable model and inputs are unchanged.

**Phase 3, true neural fine-tuning:**

- Introduce a separate `TrainablePredictorBackend`, initially RT-only.
- Bundle release-generated training model, evaluation model, optimizer model, and
  checkpoint artifacts for the exact PeptDeep model hash.
- Train only a constrained subset such as the RT output head at first.
- Export or otherwise materialize an updated inference model, reconstruct the
  `PeptDeepRTInference` session, predict all target and decoy precursors, and rebuild RT
  from canonical inputs.
- Keep the base inference sessions/models available for rollback.

MS2 and CCS/IM fine-tuning are later independent subphases. Pass-1 OSW already contains
transition-level intensity tables in OpenSWATH, but the current `OswRows` reader only
loads `FEATURE_MS2` aggregate sub-scores and RT (`src/opendialyzer.cpp:200-265`).
Before MS2 tuning, define and validate how `FEATURE_TRANSITION` apex/area intensities are
normalized into fragment labels. Before IM tuning, validate that pass-1 exposes an
observed precursor IM with enough coverage and excludes IM-derived selection scores.

## 3. CLI changes

### 3.1 Supported invocations

Keep `-in` reserved for the DIA run. Do not overload it with the old
`OpenDIALibGen -in peptides.tsv` meaning.

```bash
# Existing-library mode
OpenDIAlyzer -in run.d -tr library.tsv -out result.osw

# On-demand FASTA mode
OpenDIAlyzer -in run.d -fasta proteins.fasta -out result.osw \
  -model_dir /path/to/share/OpenMS/models

# On-demand peptide-list mode
OpenDIAlyzer -in run.d -peptides precursors.tsv -out result.osw \
  -model_dir /path/to/share/OpenMS/models
```

`-in` and `-out` remain required. Exactly one of `-tr`, `-fasta`, or `-peptides` is
required:

| Selection | Dispatch |
|---|---|
| `-tr` only | Call existing `loadLibrary_()`; no library build. Construct a catalog from
  light compounds only if adaptation needs it. |
| `-fasta` only | Instantiate live predictors, digest/enumerate, predict, refine, decoy,
  and return the light experiment in memory. |
| `-peptides` only | Instantiate live predictors, parse exact precursors, then use the
  same refine/decoy pipeline as FASTA mode. |
| None or more than one | `ILLEGAL_PARAMETERS` with the mutually exclusive inputs named. |

For a supplied library, RT calibration can work without PeptDeep. Re-prediction or
future weight tuning requires `-model_dir` and successful conversion of compounds back
to supported `AASequence` inputs. Unsupported external-library entries must be reported
as adaptation coverage; they retain original values rather than being silently dropped.

### 3.2 Options migrated from `OpenDIALibGen`

| Old option | Merged behavior |
|---|---|
| `-model_dir` | Migrate. Optional when `-tr` is used with calibration only; required for
  generated libraries and model re-prediction. Default to the installed OpenMS
  `share/OpenMS/models`, not a fragile path relative to `OpenMS_DIR`. |
| `-missed_cleavages` | Migrate for `-fasta`. |
| `-min_len`, `-max_len` | Migrate for `-fasta`. Rename internally to unambiguous
  `min_peptide_length` / `max_peptide_length`. |
| `-min_charge`, `-max_charge` | Migrate for `-fasta`. |
| `-max_var_mods` | Migrate for `-fasta`. |
| `-no_mods` | Migrate initially for compatibility; later prefer explicit fixed/variable
  modification lists. |
| `-threads` | Do not reuse one number for two runtimes. Keep current `-threads` for
  extraction and add `-onnx_threads` (default 1) for predictor intra-op threads. |
| `-no_decoys` | Do **not** allow in normal integrated search mode. Final LDA refuses to
  control FDR without both target and decoy classes
  (`src/opendialyzer.cpp:386-400`). Retain only in the temporary compatibility exporter
  or an explicitly non-scoring diagnostic mode. |
| `-raw` | Do not migrate to normal search. It means pre-refinement/pre-decoy TSV output,
  which is not a valid integrated search library. |
| `--selftest` | Move generator assertions into CTest. Extend OpenDIAlyzer’s existing
  self-test only with cheap model discovery/build checks. |

Add:

- `-peptides <file>`
- `-onnx_threads <n>`
- `-prediction_batch_size <n>`
- `-instrument <QE|Lumos|timsTOF|SciexTOF|ThermoTOF|index>`
- `-nce <float>`
- `-adaptation <off|calibrate|finetune>`; `finetune` must fail clearly when no
  trainable backend is compiled, never fall back silently.
- `-adapt_models <rt|rt,ms2|rt,ms2,im>` with only RT enabled initially.
- `-adaptation_seed <n>`
- optional `-library_out <tsv|pqp>` for debugging/provenance and migration parity. It
  must serialize the already-built in-memory library, not feed extraction through that
  file.

## 4. Build and dependency changes

### 4.1 Current fault line

The top-level project calls `find_package(OpenMS REQUIRED)` once and links both
`OpenDIALibGen` and `OpenDIAlyzer` to that same imported `OpenMS` target
(`CMakeLists.txt:21`, `CMakeLists.txt:49-65`). A single CMake configure therefore cannot
actually give one executable `openms3` and the other `openms-onnx`; the operational split
comes from separate build directories and runtime `LD_LIBRARY_PATH` choices. The
benchmark scripts show `openms3` for the older search child and `openms-onnx` for the
parent/predictor (`experiments/run_odialyzer.sh:16-32`), while the current fully
in-process script already runs OpenDIAlyzer against the ONNX tree
(`experiments/run_odialyzer_inproc.sh:13-22`).

OpenMS’s `WITH_ONNX` option is additive: it appends ONNX Runtime and PeptDeep sources to
the normal `OpenMS` library (`ext/OpenMS/src/openms/CMakeLists.txt:110-139`,
`ext/OpenMS/src/openms/CMakeLists.txt:358-367`). Therefore the correct target state is
not two OpenMS DSOs in one process. It is **one superset OpenMS build** from one source
commit, with:

- every search/input feature required from `openms3` (especially `WITH_OPENTIMS`);
- all local OpenDIAlyzer/OpenMS patches;
- `WITH_ONNX=ON`;
- one ABI, compiler, standard library, dependency set, and install prefix.

Never try to link or `dlopen` both `openms3/libOpenMS` and
`openms-onnx/libOpenMS` into OpenDIAlyzer. They export the same namespaces and library
identity; symbol interposition and ABI drift would make behavior loader-order-dependent.
The single-superset build is a hard gate before merging runtime code.

### 4.2 Proposed CMake shape

During migration:

```cmake
find_package(OpenMS REQUIRED)

add_library(odia_library
  src/odia_library_builder.cpp
  src/odia_run_adaptation.cpp)
target_include_directories(odia_library PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/src")
target_link_libraries(odia_library PUBLIC OpenMS)

add_executable(OpenDIAlyzer src/opendialyzer.cpp)
target_link_libraries(OpenDIAlyzer PRIVATE odia_library OpenMS)

option(ODIA_BUILD_LEGACY_LIBGEN "Build compatibility OpenDIALibGen" ON)
if(ODIA_BUILD_LEGACY_LIBGEN)
  add_executable(OpenDIALibGen src/opendialibgen.cpp)
  target_link_libraries(OpenDIALibGen PRIVATE odia_library)
endif()
```

`OpenMSConfig.cmake` currently exports several feature flags but not `WITH_ONNX`
(`ext/OpenMS/cmake/OpenMSConfig.cmake.in:70-82`). Add an upstream/exported
`OPENMS_WITH_ONNX` flag if this vendored OpenMS can be changed. Independently add a CMake
compile-and-link probe that includes `PeptDeepRTInference.h` and references its
constructor/prediction symbol. Header presence alone is insufficient because the
headers can exist while PeptDeep sources were omitted from `libOpenMS`.

OpenMS’s current finder checks only the core inference header
`onnxruntime_cxx_api.h` and `libonnxruntime`
(`ext/OpenMS/cmake/FindONNXRuntime.cmake:21-45`). That is enough for Phases 1–2 and is
positive evidence that the current build is inference-only, not evidence of training
support.

### 4.3 Models and runtime packaging

Rename `OPENDIALIBGEN_MODEL_DIR` to `ODIA_PEPTDEEP_MODEL_DIR` and derive the default from
the exported `OPENMS_DATA_DIR` / installed share directory. The current top-level
relative default (`CMakeLists.txt:93-94`) is tied to a build-tree layout. The vendored
OpenMS build already lists all three model files and hashes and installs them
(`ext/OpenMS/src/openms/CMakeLists.txt:128-146`,
`ext/OpenMS/src/openms/CMakeLists.txt:204-206`).

At configure/package/test time:

- verify all three model files exist and match expected SHA-256 values;
- package the matching ONNX Runtime shared library and set a correct install RPATH;
- record the OpenMS commit/build ID, ONNX Runtime version, model hashes, compiler, and
  relevant feature flags;
- use `ldd`/`readelf` on Linux or `otool -L` on macOS in CI to assert exactly one
  `libOpenMS` is resolved;
- run both a PeptDeep prediction smoke test and a Bruker/mzML extraction smoke test
  against that same install tree.

Phase 3 needs a distinct dependency contract. Official ONNX Runtime documentation says
the training C/C++ API requires a training-enabled build plus a training model,
checkpoint, optimizer model, and optional eval model; these artifacts are prepared
offline, currently with the Python training-artifact tools
([ORT training C/C++ API](https://onnxruntime.ai/docs/api/c/training_c_cpp_api.html),
[training artifact preparation](https://onnxruntime.ai/docs/api/python/on_device_training/training_artifacts.html),
[training build](https://onnxruntime.ai/docs/build/training.html)). Do not enable
`-adaptation finetune` merely because core ONNX Runtime was found.

## 5. Migration, tests, CI, and retirement

### 5.1 Migration sequence

1. **Build unification prerequisite:** produce and qualify one ONNX-enabled superset
   OpenMS install.
2. **Extract without behavior change:** move parsing, digestion, peptidoforms,
   prediction, transition assembly, refinement, decoys, and export into
   `odia_library`. Keep `OpenDIALibGen` output behavior.
3. **Remove TSV handoff internally:** assemble `TargetedExperiment` directly, convert to
   light, and let the compatibility CLI serialize only at its outer boundary.
4. **Integrate generated-library dispatch:** add `-fasta`/`-peptides` to OpenDIAlyzer and
   pass the in-memory light experiment directly to `extractPass_()`.
5. **Add runtime RT calibration state:** preserve canonical predictions, retain live
   predictor sessions, add stable-ID anchors, and rebuild RT for pass 2.
6. **Deprecate:** emit an `OpenDIALibGen` warning pointing to OpenDIAlyzer on-demand mode;
   keep the target for one compatibility window.
7. **Retire:** default `ODIA_BUILD_LEGACY_LIBGEN=OFF`, then remove
   `src/opendialibgen.cpp`, its target, old self-test, old model-dir variable, and stale
   documentation after the gates below remain green.

### 5.2 Test gates

**Unit tests**

- FASTA digestion and peptidoform enumeration preserve current behavior, including
  protein-N-terminal acetyl gating (`src/opendialibgen.cpp:125-166`,
  `src/opendialibgen.cpp:401-436`).
- Peptide-list and FASTA sources converge on the same downstream refinement/decoy path.
- Chunked prediction equals one-batch prediction within a declared float tolerance.
- Transition channel→ion annotation, m/z, intensity, charge, RT, and 1/K0 match the
  current generator (`src/opendialibgen.cpp:199-300`).
- Direct in-memory assembly has no invalid references and converts to light without
  losing modifications, charge, RT, IM, decoy state, or transition flags.
- Repeated builds with the same seed produce semantically identical libraries.
- Pass materialization never compounds RT transforms.
- Rebuilds preserve compound/transition IDs and target↔decoy pairs.
- `finetune` returns `UnsupportedBackend` in inference-only builds.

**Semantic parity tests**

- Run old `OpenDIALibGen` and the new compatibility front end on a fixed peptide list and
  FASTA. Compare normalized `TargetedExperiment`/light contents rather than raw TSV
  formatting alone.
- Compare target/decoy compound counts, transition counts, m/z, intensity, RT, IM,
  modifications, flags, and stable IDs.
- Because current shuffled decoys are time-seeded
  (`src/opendialibgen.cpp:595-597` and
  `ext/OpenMS/src/openms/source/ANALYSIS/OPENSWATH/MRMDecoy.cpp:591-594`), either add an
  explicit seed to the vendored decoy API or use a deterministic decoy method before
  declaring byte/semantic determinism. A seed field in `AssayOptions` is not enough if
  `MRMDecoy` ignores it.

**Integration tests**

- Tiny FASTA + tiny mzML: invoke OpenDIAlyzer without `-tr`, produce a valid OSW with both
  classes, run both passes, and write final scores.
- Build-and-search in one invocation versus prebuild-then-`-tr`: the same semantic
  library must give equivalent extraction/scoring output within declared tolerances.
- Existing TSV, PQP, and OSWPQ `-tr` paths remain unchanged.
- Unsupported external-library modifications report partial adaptation coverage without
  dropping entries.
- Synthetic RT distortion is recovered by pass-2 calibration; held-out RT error
  improves and fallback behavior is deterministic.
- A small real run confirms IDs/FDR do not regress before making adaptation the default.

**Build/CI tests**

The root project currently has CTest registrations but no root CI workflow. Add a Linux
CPU workflow that builds the vendored OpenMS superset once, caches ONNX Runtime/models by
version and hash, then builds and tests this project. Add a packaging job that checks
runtime resolution and a small installed-tree smoke test. Keep the pure
`odia-score-test` and LDA tests fast and independent, but the main OpenDIAlyzer job must
always use the ONNX-enabled OpenMS.

### 5.3 Retirement gate

Remove the standalone tool only when all are true:

- on-demand FASTA and peptide-list searches pass integration tests;
- generated in-memory library semantics match the compatibility exporter;
- `-tr` behavior is not regressed;
- model discovery and runtime packaging work from an installed tree;
- deterministic target and decoy generation is solved;
- one real-run benchmark meets declared memory/time limits;
- the integrated command can optionally export the built library for audit/debug;
- user documentation and experiment scripts no longer invoke `OpenDIALibGen`.

Neural fine-tuning is deliberately not a retirement prerequisite. If it were, the
inference/training uncertainty could keep a now-redundant CLI alive indefinitely.

## 6. Phased implementation plan

### Phase 0 — one OpenMS runtime

Deliver one OpenMS build with search features plus `WITH_ONNX=ON`; prove one-DSO runtime
resolution, model discovery, predictor smoke tests, and `.d`/mzML extraction.

**Gate:** no merger code lands until the same OpenMS install passes both prediction and
search smoke tests.

### Phase 1 — smallest change: in-memory library handoff

Extract `odia_library`, keep the legacy CLI as a wrapper, replace TSV-string/temp-TSV
handoff with direct in-memory assembly and heavy→light conversion, and add
`-fasta`/`-peptides` dispatch to OpenDIAlyzer. Keep predictors alive in
`LibrarySession`, but do not alter pass-2 predictions yet.

**Success:** parity tests pass; OpenDIAlyzer searches a tiny FASTA without `-tr`; peak
memory and elapsed time are recorded; no external program is executed.

### Phase 2 — runtime re-prediction and calibration

Add stable string IDs to pass-1 rows, canonical base predictions, RT anchor policies,
run-specific RT calibration, target/decoy-symmetric rebuild, provenance, and rollback.
Keep the predictor sessions live and expose re-prediction, but retain immutable model
weights. Add MS2/IM calibration only after label-quality tests.

**Success:** held-out RT error improves; a real-run A/B does not worsen q < 1% IDs or
FDR diagnostics; repeated runs with the same inputs/seed agree.

### Phase 3a — RT neural fine-tuning proof of concept

Build and package a training-enabled ONNX Runtime separately, generate version-locked
training artifacts offline, implement an RT-only trainable backend, freeze most weights,
train on q-controlled anchors, validate on held-out anchors, export/reload an inference
model, and roll back automatically unless it beats calibration alone.

**Success:** the training artifact and optimizer pipeline works on every supported
platform; held-out RT and end-to-end ID metrics improve over Phase 2 without FDR
inflation; memory/time are acceptable.

### Phase 3b — MS2, then IM/CCS fine-tuning

First define trustworthy empirical labels and score-exclusion rules. Implement MS2 and
IM as separate experiments with separate rollback and success criteria. Do not enable
them merely because RT training worked.

### Phase 4 — retirement

After Phase 1 parity and Phase 2 operational stability, deprecate for one compatibility
window, default the legacy target off, then delete it and update CMake/tests/docs/CI.

## 7. Deep self-review

### 7.1 First-principles attack

The plan’s original premise—“keep inference predictors alive, therefore fine-tune
them”—is false. Object lifetime solves only model reload cost and access to prediction
methods. The OpenMS wrappers hide an inference `Ort::Session` behind private members and
provide no training surface. True fine-tuning changes the dependency, model artifacts,
failure modes, memory use, reproducibility, and validation burden. It must not be
smuggled into the merger as a routine method call.

The real minimum needed for goal (a) is much smaller: a reusable builder returning a
light library. The minimum useful version of goal (b) is also smaller: preserve canonical
PeptDeep predictions and fit a per-run calibration layer. That is likely to capture much
of the RT benefit with far less risk. The plan should earn the complexity of weight
updates only after calibration is measured.

The current engine also does not yet provide clean neural-training labels. Its RT
fallback intentionally uses unvalidated top-N targets; its MS2 reader sees aggregate
scores rather than fragment targets; and any same-run adaptation can reinforce false
identifications. A training API alone would not make the science sound.

### 7.2 Failure-mode and pre-mortem risk matrix

Scale: likelihood and impact are 1 (low) to 5 (high). Overall rating is the higher of the
two, elevated when the failure invalidates FDR or the whole process.

| Risk / imagined failure | Likelihood | Impact | Rating | Mitigation and gate |
|---|---:|---:|---|---|
| Two incompatible OpenMS builds resolve in one process; prediction or extraction crashes or silently calls the wrong ABI | 4 | 5 | **Critical 5** | One superset OpenMS build only; inspect dynamic dependencies; prediction + `.d` + mzML tests from the installed tree. Hard no-go before Phase 1. |
| The ONNX-enabled OpenMS tree lacks a patch or search feature present in `openms3` | 4 | 5 | **Critical 5** | Rebuild from one source commit/options manifest; diff patches and feature flags; qualify search parity before merging. |
| “Fine-tuning” is implemented against inference-only ORT and cannot update weights | 5 | 5 | **Critical 5** | Typed `supportsWeightTraining=false`; Phase 3 requires training-enabled ORT and versioned training artifacts. Phase 3 is currently no-go. |
| The shipped PeptDeep ONNX graph cannot generate usable training artifacts or expose a safe subset of trainable parameters | 4 | 5 | **Critical 5** | Offline artifact-generation spike on the exact model hashes; inspect parameter names/dynamic axes; RT-head-only POC before engine integration. |
| Same-run false IDs train the model, which makes those IDs score better in pass 2 and defeats target-decoy assumptions | 4 | 5 | **Critical 5** | q-controlled anchors only for training, score-feature exclusion, group-separated holdout, target/decoy-symmetric application, rollback, independent runs for final validation. |
| Decoys retain target-like old predictions while targets are re-predicted, or decoys are trained/adapted differently | 3 | 5 | **Critical 5** | Persist target↔decoy mapping and exact decoy sequences; rebuild both classes with the same backend/calibrator; assert symmetry and monitor null-score shifts. |
| Current top-2,000 RT fallback is mistaken for “confident IDs” and used as training truth | 4 | 5 | **Critical 5** | Allow fallback only for robust calibration; prohibit it in all weight-update paths. |
| In-place RT transforms compound across passes or adaptation attempts | 3 | 4 | **High 4** | Immutable base predictions; restore/materialize from canonical values; explicit no-compounding test. |
| OSW integer IDs do not match rebuilt library IDs, so wrong precursors are updated | 3 | 5 | **Critical 5** | Carry `TRAML_ID` through `OswRows`; use integer IDs only inside one OSW; stable-ID join tests. |
| MS2 re-ranking cannot restore transitions pruned during initial assay refinement | 4 | 4 | **High 4** | Retain unpruned fragment predictions or regenerate complete assays; replace and re-index the affected assay, not just its surviving intensities. |
| Whole-proteome FASTA expansion plus modifications, three prediction outputs, heavyweight and light libraries, and in-memory SWATH maps exhaust RAM | 4 | 5 | **Critical 5** | Chunk prediction; avoid the TSV string; release prediction batches/heavy experiment before DIA load; benchmark peak RSS; later use light refinement/decoys; consider an optional content-addressed PQP cache. |
| ONNX and OpenMP thread pools oversubscribe the host or make timing nondeterministic | 4 | 3 | **High 4** | Separate `-onnx_threads` from extraction threads; default ONNX to 1; do not run prediction concurrently with extraction initially; benchmark explicitly. |
| Shuffle decoys vary between runs despite a nominal seed option | 5 | 3 | **High 4** | Patch `MRMDecoy` to accept/propagate a seed or choose a deterministic method; require semantic repeatability before retirement. |
| External libraries contain unsupported modifications/RT scales and cannot be reconstructed for re-prediction | 3 | 4 | **High 4** | Adaptation coverage report; preserve unsupported entries; explicit RT-scale metadata/override; require full coverage only for `finetune`. |
| PeptDeep’s fixed NCE/instrument assumptions are wrong for a run, so on-demand prediction is reproducibly wrong | 3 | 4 | **High 4** | Runtime prediction options, provenance, CLI overrides, later metadata-derived defaults; validate instrument mapping/model hash. |
| Fine-tuned model overfits a small run and harms portability/reproducibility | 4 | 4 | **High 4** | Run-scoped models only; base model immutable; constrained head/layers; held-out gate and automatic rollback; never overwrite packaged models. |
| Retirement removes useful library-only workflows before the integrated replacement is auditable | 2 | 3 | **Medium 3** | One compatibility window and optional `-library_out`; keep module-level export APIs even after deleting the executable. |

### 7.3 Determinism attack

The current target library is intended to be deterministic, but the full generated
library is not: the source itself notes time-seeded decoys
(`src/opendialibgen.cpp:595-597`), and `generateDecoys()` calls
`shufflePeptide(..., -1, ...)` (`MRMDecoy.cpp:591-594`). Fine-tuning adds more sources of
variation: anchor folds, optimizer order, floating-point reductions, ORT kernels, thread
counts, and exported checkpoint state.

Required determinism levels should be explicit:

- **Phase 1:** semantic identity for compounds/transitions and fixed seeded decoys;
  byte-identical serialization where practical.
- **Phase 2:** identical anchor selection, calibration knots, and rebuilt library for the
  same inputs/seed/build.
- **Phase 3:** reproducible within declared numeric tolerances on the same platform and
  build. Cross-platform byte identity is not a credible requirement for neural training.

Every run must record seeds, thread counts, library/model/build hashes, anchor IDs, and
adaptation result. A reproducibility manifest is more useful than pretending all
floating-point execution is byte-stable.

### 7.4 Memory and performance attack

Keeping three model sessions alive is not itself the largest memory problem. The current
generator holds all precursor inputs, all RT/MS2/CCS outputs, a complete TSV string, a
temporary parsed heavyweight experiment, and then an output file. The current search
also deliberately loads SWATH data in memory (`src/opendialyzer.cpp:125-142`) and asks
the workflow to operate in memory (`src/opendialyzer.cpp:454-457`).

The merger wins only if it removes the TSV string and temp parse, predicts in bounded
chunks, and releases heavyweight construction state before loading spectra. A naive
shared class that simply retains every current intermediate would make peak RSS worse.
Set an explicit Phase 1 memory budget and measure:

- precursor count after digestion/mod expansion;
- transitions before and after refinement;
- predictor-session RSS;
- peak builder RSS;
- peak full-run RSS;
- time in digestion, prediction by model, refinement/decoy, load, and each pass.

If the heavy→light peak is unacceptable, use the existing light `MRMAssay`/`MRMDecoy`
APIs. If repeated on-demand generation dominates runtime, add a content-addressed cache
keyed by FASTA/peptide-list content, digestion/modification options, prediction options,
model hashes, and builder version. The cache is an optimization, not the primary
handoff.

### 7.5 Is the merge worth it versus a thin in-process call?

Three alternatives were considered:

| Alternative | On-demand/no subprocess | No TSV round-trip | Predictors stay live | Supports adaptation | Verdict |
|---|---:|---:|---:|---:|---|
| Execute `OpenDIALibGen` | No | No | No | No | Rejected by directive and architecture |
| Call current `generate()` from OpenDIAlyzer | Yes | No; returns TSV string | No; predictors are locals | No | Useful throwaway spike only |
| Shared `LibraryBuilder` + `LibrarySession` | Yes | Yes | Yes | Yes, with honest backend limits | Recommended |

A “thin in-process call” is attractive for Phase 1 schedule, but it preserves exactly the
wrong boundaries: serialization and predictor destruction. The shared module is worth
the modest extra design because both requested capabilities depend on its ownership
model. It should still remain narrow—builder, state, adaptation—not become a general
plugin framework before Phase 3 proves a second backend is real.

Even if neural fine-tuning never ships, the merge remains worthwhile for on-demand
generation, consistent FASTA/peptide-list behavior, removal of TSV I/O, provenance, and
run-specific calibration. What would not be justified is claiming that those benefits
alone close the DIA-NN identification gap; that remains an empirical result, not an
architectural consequence.

### 7.6 Phase go/no-go

| Phase | Decision now | Conditions |
|---|---|---|
| 0 — unified ONNX-enabled OpenMS | **GO, mandatory spike** | Must pass prediction and search parity with one runtime. If it fails, all later phases are no-go. |
| 1 — shared builder + in-memory light handoff | **Conditional GO** | Proceed after Phase 0; require semantic parity, deterministic decoys, and memory budget. |
| 2 — live inference + RT recalibration/rebuild | **GO** | Call it calibration, not weight training; require stable IDs, canonical predictions, target/decoy symmetry, held-out RT improvement, and rollback. |
| 2b — MS2/IM calibration | **Conditional GO** | Only after fragment/IM labels and anti-circular selection rules are validated. |
| 3a — RT neural fine-tuning | **NO-GO today** | Reconsider only after exact-model ORT-training artifact POC, supported-platform build, held-out improvement, FDR-safety, and resource gates. |
| 3b — MS2/IM neural fine-tuning | **NO-GO today** | Requires Phase 3a success plus trustworthy modality-specific labels and independent validation. |
| 4 — retire `OpenDIALibGen` | **Conditional GO** | After Phase 1 parity, Phase 2 operational stability, installed-tree packaging, one compatibility window, and updated docs/scripts. Neural fine-tuning is not required. |

## Final recommendation

Start with Phase 0 and Phase 1. The exact seam is already favorable:
`extractPass_()` consumes `LightTargetedExperiment`, so the first valuable merger removes
serialization without disturbing extraction. Then implement Phase 2 as an explicit
run-calibration layer with live inference sessions and immutable base predictions.

Do not schedule true PeptDeep fine-tuning as routine follow-on implementation. Schedule
an RT-only feasibility spike with training-enabled ONNX Runtime and versioned training
artifacts. Until that spike passes, the product language and CLI must distinguish
“calibrate/re-predict” from “fine-tune neural weights.”

---

# Review verdict (Claude, 2026-07-27) — APPROVED, phased, with reconciliations

Deep-reviewed against the actual tree. **Every load-bearing claim verified in code:**
- ONNXPredictorBase exposes NO optimizer/gradient/checkpoint/training symbols
  (`ext/OpenMS/.../ONNX/ONNXPredictorBase.h`) → neural fine-tuning genuinely infeasible in
  this build. **Confirmed.**
- MRMDecoy shuffle calls `shufflePeptide(..., -1, ...)` (`MRMDecoy.cpp:593`) = time-seeded →
  non-deterministic decoys; already documented in-tree as "vendored patch P3"
  (`opendialibgen.cpp:595`). **Confirmed** — real benchmark-reproducibility hazard.
- `opendialibgen.cpp` really does FASTA digestion/mods/CCS/refine/decoy (top comment stale).
  **Confirmed.**
- `extractPass_` already consumes `LightTargetedExperiment` → Phase-1 seam is clean.
  **Confirmed.**

The plan's central correction — *"predictors alive" ≠ "predictors trainable"; the three
capabilities are three milestones, not one* — is correct and is the most important thing in
the document. Approved as the design of record.

## Reconciliations / deltas to apply when implementing

1. **The merger is NOT the DIA-NN ID-gap lever.** Closing the gap is calibration + the
   empirical (observed-value) library — done in Task B WITHOUT the merger, against the
   existing `-tr` library. The merger buys on-demand generation, consistent FASTA/peptide
   behaviour, no TSV round-trip, provenance, and the *ownership model* fine-tuning would one
   day need. Do not gate the ID-gap work on it. (The plan says this in §7.5; stating it up
   front so it drives scheduling.)

2. **"Fine tuning of the model" in the current build = per-run CALIBRATION, not weight
   updates.** Phase 3a (neural fine-tune) is NO-GO until a training-enabled ONNX Runtime +
   version-locked training artifacts exist — a separate spike, not follow-on work. CLI and
   docs must keep `calibrate`/`re-predict` distinct from `fine-tune` (plan §3.2, §7.1).

3. **Task B already down-paid part of Phase 2.** `OswRows` now carries `TRAML_ID` (plan
   §2.6 asked for exactly this) and the engine does an empirical RT overwrite of pass-1 IDs
   (isotonic PAVA transform + `applyEmpiricalRT_`). So Phase 2's RT-calibration slice is
   partly built; the merger's job there is to formalise it behind `RunLibraryAdapter`.

4. **ADD to Phase 2 explicitly: empirical (observed-value) RT/MS2/IM replacement**, not only
   re-prediction of model values. The DIA-NN calibration research ranks the observed-value
   second pass as the #1 lever — bigger than re-predicting the same model. The plan mentions
   "materialize from canonical + calibration" but should name the observed-value path as a
   first-class Phase-2 deliverable.

5. **Ponytail on Phase 1 scope.** Build Phase 1 with the *minimal* types: a thin
   `LibraryBuilder` returning `LightTargetedExperiment` + a `shared_ptr<PeptDeepPredictors>`
   kept alive. Defer `LibraryState`/`BasePrediction`/`RunLibraryAdapter`/`AnchorSelector`
   and the immutable-canonical-state machinery to Phase 2, when they're actually exercised.
   The rich §2.3 hierarchy is the Phase-2/3 end state, not the Phase-1 diff.

6. **Self-critique the plan surfaces about Task B code:** the engine mutates
   `LightCompound::rt` in place across passes (`rescaleLibraryRT_` then `applyEmpiricalRT_`).
   Correct for `passes ≤ 2` (verified; the C17 guard enforces it), but the plan's
   immutable-canonical-state design (§2.1.3) is more robust and should replace the in-place
   mutation in Phase 2.

## Go / no-go tonight

**Do NOT start merger implementation now.** Blocked by:
- **Phase 0** (one ONNX-enabled *superset* OpenMS build: openms3 search features incl.
  WITH_OPENTIMS + WITH_ONNX, one DSO) — a heavy prerequisite that needs the cluster
  (currently under maintenance). This is the correct hard gate.
- Retirement is further gated on Phase-1 parity + Phase-2 stability + the decoy-determinism
  fix + installed-tree packaging.

**Near-term targets once the cluster returns:** Phase 0 (superset build) → Phase 1
(in-memory handoff + `-fasta`/`-peptides`). Phase 2 formalises the Task-B calibration work.
Phase 3 stays NO-GO. `OpenDIALibGen` keeps building until Phase-1 parity + Phase-2 stability.
