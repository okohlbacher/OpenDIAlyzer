# Vendored OpenMS patches

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
described in the vault note; they are one-line changes applied at build time by
`build_openms_onnx.sh` rather than tracked files.

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
tracked patch documented only a subset. Verify the two agree with
`git -C ext/OpenMS diff --stat`; any file listed there and absent here is undocumented drift.

After a build, both the installed prefix and `ext/OpenMS/src` are `chmod -R a-w`
(`experiments/setup_clean_tree.sh` checks this), so modifying the core requires a deliberate
`chmod -R u+w`. **All OpenDIAlyzer code belongs in `src/`, never in the OpenMS tree.**

| Change | Why |
|---|---|
| `OpenSwathWorkflow.cpp` — chunked parallel extraction | Upstream issues one `extractChromatograms` call, leaving most cores idle at full scale. Splits coordinates into per-thread chunks, each with its own `lightClone()`. Safe: the merge-join cursor is per-call and each chunk writes only its own output slice. |
| `OpenSwathWorkflow.cpp` — free `chrom_list` after `return_chromatogram` | It otherwise keeps a full duplicate of every chromatogram alive through scoring. |
| `TransitionListEvidenceFilter.{h,cpp}` — `ms2_min_qualifying_spectra` | Upstream decides MS2 support from the best SINGLE spectrum in the run, so on a 3889-cycle gradient every precursor gets 3889 independent chances at an m/z coincidence — measured, shuffled decoys pass at 1.03x the target rate, i.e. no discrimination. Counting recurrence restores the time dimension. |
| `ParquetFile.{h,cpp}`, `TransitionParquetFile.cpp`, `ArrowSchemaRegistry.cpp` — `ChunkedColumn` | **Load-bearing.** A single arrow string array cannot exceed 2 GB and `traml_id` on the 7.1M-precursor library is far past that, so `getColumn()` returning only chunk 0 would silently truncate the library. It now throws, and every column resolves (chunk, offset). |
| `ChromatogramExtractorAlgorithm.cpp` — exact `reserve()` | Metadata pass + binary search instead of growth-by-doubling. Measured 42% capacity slack and ~2.7M reallocations per batch on the Astral benchmark. |
| `Scoring.{h,cpp}` + `MRMScoring.cpp` — `normalizedCrossCorrelationMaxPost()` | Allocation-free maximum fusing the cross-correlation loop with the max scan. **Implemented but NOT wired in**: `MRMScoring_test` asserts on the full matrix via `getXCorrMatrix()`, so switching is an API decision. |
