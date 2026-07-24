# Vendored OpenMS patches

OpenDIALibGen builds against a **vendored, patched** OpenMS (the user's decision;
see `vault/20-OpenSWATH/Vendored OpenMS patches.md`). `ext/OpenMS` is gitignored,
so the patches themselves are tracked here and copied into the vendored tree.

| Patch | Applied to | Purpose | Validated by |
|---|---|---|---|
| `OpenMS/PeptDeepModX.h` | `src/openms/include/OpenMS/ML/PEPTDEEP/PeptDeepModX.h` (new file) | **P2**: modification featurization for PeptDeep ONNX — fills `mod_x` from an `AASequence`'s mods, matching AlphaPeptDeep exactly. Header-only (no libOpenMS rebuild to use). | `src/odia_modx_test.cpp` (bit-parity vs `src/testdata/peptdeep_modx_reference.json`) |

Build-system patches (P1 FindONNXRuntime module-path, curl brotli/zstd off) are
described in the vault note; they are one-line changes applied at build time by
`build_openms_onnx.sh` rather than tracked files.

Still to come: wiring the validated featurization into `PeptDeepInputBuilder`
(`buildModified*`) and the inference classes, which does require an OpenMS
rebuild — deferred until the plumbing is written.
