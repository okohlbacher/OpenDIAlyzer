# Vendored OpenMS patches

**Base commit:** `d77542de65` of <https://github.com/OpenMS/OpenMS>
(never push to that remote; the local clone's push URL is deliberately disabled).
`scripts/openms/setup_node.sh` parses this exact line, so keep the format.

OpenDIALibGen builds against a **vendored, patched** OpenMS (the user's decision;
see `vault/20-OpenSWATH/Vendored OpenMS patches.md`). `ext/OpenMS` is gitignored,
so the patches themselves are tracked here and copied into the vendored tree.

| Patch | Applied to | Purpose | Validated by |
|---|---|---|---|
| `OpenMS/PeptDeepModX.h` | `src/openms/include/OpenMS/ML/PEPTDEEP/PeptDeepModX.h` (new file) | **P2 featurization**: fills `mod_x` from an `AASequence`'s mods, matching AlphaPeptDeep exactly. Header-only. | `src/odia_modx_test.cpp` (bit-parity vs `src/testdata/peptdeep_modx_reference.json`) |
| `OpenMS/peptdeep-mod-support.patch` | `PeptDeepInput.{h,cpp}`, `PeptDeepRTInference.{h,cpp}`, `PeptDeepMS2Inference.{h,cpp}` (`git apply` in the OpenMS clone) | **P2 wiring**: `buildModified*` batch builders + `predictRT`/`predictMS2` `AASequence` overloads that populate `mod_x`. The string overloads delegate to these, so they are transparent supersets. Requires a libOpenMS rebuild. | `PeptDeepInference_test` (unmodified 1e-7 preserved) + `src/odia_modpredict_test.cpp` (modified prediction differs, correctly) |
| `OpenMS/opendialyzer-openswath.patch` | 11 files under `ANALYSIS/OPENSWATH`, `FORMAT`, `openswathalgo` (`git apply` in the OpenMS clone) | **OpenDIAlyzer's OpenSWATH delta** — see the breakdown below. | `OpenDIAlyzer -selftest`, `MRMScoring_test`, `TransitionListEvidenceFilter_test`, and the Astral benchmark |

To apply on a fresh OpenMS clone: copy `PeptDeepModX.h` into place, then
`git apply vendored-patches/OpenMS/peptdeep-mod-support.patch` from the OpenMS
root, then rebuild the `OpenMS` target.

Build-system patches (P1 FindONNXRuntime module-path, curl brotli/zstd off) are
described in the vault note. They are one-line changes applied by hand at
configure time and are NOT tracked here. (An earlier revision of this file
attributed them to a `build_openms_onnx.sh` that does not exist in the
repository -- if that script is recreated, list it here.)

Still to come: wiring the validated featurization into `PeptDeepInputBuilder`
(`buildModified*`) and the inference classes, which does require an OpenMS
rebuild — deferred until the plumbing is written.

## `opendialyzer-openswath.patch` — contents

**This patch IS the complete OpenSWATH-side delta.** Regenerate it with:

```bash
git -C ext/OpenMS diff -- . ':(exclude)*PEPTDEEP*' > vendored-patches/OpenMS/opendialyzer-openswath.patch
```

`ext/OpenMS` is gitignored, so an edit to the OpenMS core otherwise leaves no trace in this
repository. That is exactly how a set of ad-hoc modifications once accumulated undetected while the
tracked patch documented only a subset. `scripts/openms/regen-patch.sh --check` fails if the tree and this patch disagree — run it before
committing, so forgetting becomes a failed check rather than a silent loss.

After a build, both the installed prefix and `ext/OpenMS/src` are `chmod -R a-w`
(`scripts/openms/setup_clean_tree.sh` checks this), so modifying the core requires a deliberate
`chmod -R u+w`. **All OpenDIAlyzer code belongs in `src/`, never in the OpenMS tree.**

| Change | Why |
|---|---|
| `OpenSwathWorkflow.cpp` — chunked parallel extraction | Upstream issues one `extractChromatograms` call, leaving most cores idle at full scale. Splits coordinates into per-thread chunks, each with its own `lightClone()`. Safe: the merge-join cursor is per-call and each chunk writes only its own output slice. |
| `OpenSwathWorkflow.cpp` — free `chrom_list` after `return_chromatogram` | It otherwise keeps a full duplicate of every chromatogram alive through scoring. |
| `TransitionListEvidenceFilter.{h,cpp}` — `ms2_min_qualifying_spectra` | Upstream decides MS2 support from the best SINGLE spectrum in the run, so on a 3889-cycle gradient every precursor gets 3889 independent chances at an m/z coincidence — measured, shuffled decoys pass at 1.03x the target rate, i.e. no discrimination. Counting recurrence restores the time dimension. |
| `ParquetFile.{h,cpp}`, `TransitionParquetFile.cpp`, `ArrowSchemaRegistry.cpp` — `ChunkedColumn` | **Load-bearing.** A single arrow string array cannot exceed 2 GB and `traml_id` on the 7.1M-precursor library is far past that, so `getColumn()` returning only chunk 0 would silently truncate the library. It now throws, and every column resolves (chunk, offset). |
| `ChromatogramExtractorAlgorithm.cpp` — exact `reserve()` | Metadata pass + binary search instead of growth-by-doubling. Measured 42% capacity slack and ~2.7M reallocations per batch on the Astral benchmark. |
| `ParquetFile.{h,cpp}` — `ChunkedColumn::resolveRaw()` + non-owning value accessors | The owning `resolve()` returns `shared_ptr<arrow::Array>` BY VALUE, and the accessors call it once per (row, column): ~1e9 shared_ptr copies on a 78.6M-row table, each an atomic refcount pair on control blocks shared by every thread. Cost therefore GROWS with thread count — measured, a row loop went 20.1 s serial to 34.5 s on 64 threads. The owning forms now delegate to the raw ones, so there is one implementation of the null/type handling. **Measured: library_load 144.8 s -> 133.9 s.** |
| `TransitionParquetFile.cpp` — exact `reserve()` on the destination vectors | `num_rows()` is already read on the same lines. **Measured effect: ~0** — vector growth was not a real cost, contrary to the estimate that motivated it. Kept because it is free and correct, but it is not the win. |
| `OpenSwathWorkflow.cpp` — region accounting inside `performExtraction` | That phase is **65% of wall and 90% of CPU** (48,120 CPU-s) and had never been measured from inside. Both external profilers are unavailable on the benchmark nodes: `perf_event_paranoid=4` blocks sampling and `ptrace_scope` refuses a gdb attach. Five `ScopedRegion` accumulators (extract / convert / score / write / MS1) report thread-seconds summed across the parallel loop, directly comparable with the phase timers' CPU-seconds. Instrumentation only — no behaviour change. |
| `Scoring.{h,cpp}` + `MRMScoring.cpp` — `normalizedCrossCorrelationMaxPost()` | Allocation-free maximum fusing the cross-correlation loop with the max scan. **Implemented but NOT wired in**: `MRMScoring_test` asserts on the full matrix via `getXCorrMatrix()`, so switching is an API decision. |
