# OpenDIALibGen — Stage 1

Peptide-list TSV in, OpenSWATH transition TSV out, with RT and MS2 fragment
intensities predicted by the vendored OpenMS PeptDeep ONNX bindings.

```
peptide-list TSV ──▶ PeptDeepRTInference (iRT)
                     PeptDeepMS2Inference (fragment intensities)
                     TheoreticalSpectrumGenerator (b/y m/z, paired by annotation)
                   ──▶ OpenSWATH transition TSV (targets only)
```

## Scope

**In:** unmodified uppercase peptide sequences with a precursor charge;
predicted iRT; predicted b/y z1/z2 fragment intensities paired to theoretical
m/z; one OpenSWATH TSV with the exact column set that
`experiments/diann_lib_to_openswath.py` emits (minus `GeneName` and
`PrecursorIonMobility`), so both library arms feed OpenSWATH identically.

**Out (later stages, each with its own gate):** FASTA digest, modifications,
CCS/ion mobility, decoys (use the existing `OpenSwathDecoyGenerator`),
assay refinement, pqp/TraML, non-specific, fine-tuning, CUDA.

**Parity is inherited, not re-tested.** Prediction correctness vs AlphaPeptDeep
Python ONNX-Runtime is already proven by OpenMS's own `PeptDeepInference_test`
(RT max err 5.96e-08, CCS 1.2e-04, MS2 intensity 7.15e-07 against saved Python
reference values). Stage 1 does not re-validate the predictor; `--selftest` is
a shape/sanity check of the assembly logic this repo owns (column set, finite
non-zero intensities, base-peak normalization, byte-identical repeat).

## Fixed prediction conditions (provenance)

The MS2 model is conditioned on (charge, NCE, instrument); wrong values give
plausible-but-wrong spectra. Stage 1 pins:

- **NCE = 30.0** — AlphaPeptDeep's default for the pretrained ms2 model
  (`nce=30, instrument='QE'`), and the central value of OpenMS's parity
  reference data (`proteomicsml_test_data_ms2_spectra.csv`, nce 25.03–32.5).
- **instrument_index = 0 (QE)** — the same reference data pairs
  `instrument_index 0` with `instrument "QE"` and passes
  `PeptDeepInference_test` against AlphaPeptDeep Python ONNX output, so
  `0 → QE` is the empirically verified encoding. The docstring in
  `PeptDeepMS2Inference.h` claims `0=Lumos`; the parity data contradicts it
  and is the only source backed by a passing test. **TODO (reviewer):**
  cross-check against AlphaPeptDeep's `instrument_dict` for the pinned
  checkpoint before Stage 2.

## Intensity ↔ m/z pairing

`predictMS2` returns a flat `(nAA-1)*8` array per peptide, fragment-major:
`index = fragment_position*8 + ion_index` with channels
`0=b_z1, 1=b_z2, 2=y_z1, 3=y_z2, 4..7=b/y modloss` (channel order read from
`proteomicsml_test_data_ms2_predicted_intensities.csv`). Modloss channels are
zero for unmodified peptides and are dropped. Pairing to
`TheoreticalSpectrumGenerator` annotations (`b3+`, `y5++` via
`add_metainfo=true`): b-series number *s* ↔ fragment position *s−1*; y-series
number *s* ↔ fragment position *nAA−1−s*. `add_first_prefix_ion=true` so every
predicted channel has a theoretical m/z partner (b1_z1 is generated but never
emitted — `predictMS2` forces its intensity to 0). `predictMS2` already
base-peak normalizes to 1.0 and zeroes negatives and values < 1e-4; rows with
zero intensity are not emitted. Fragment charge is capped at
`min(precursor charge, 2)` (AlphaPeptDeep's library-builder convention).

## Output conventions

- `NormalizedRetentionTime` carries the predicted **iRT-scale** value as-is —
  not seconds. RT/IM alignment ownership (spike-in/iRT normalization) is an
  open question deferred to Stage 4 (see the plan, kimi S7).
- `ProteinId` is `UNKNOWN` (a bare peptide list has no proteins; Stage 2's
  digest supplies them).
- `TransitionGroupId` = `PeptideSequence_PrecursorCharge`, `TransitionId` =
  group + 0-based running index, matching the converter's fallback synthesis.
- `Decoy=0`, `DetectingTransition=1`, `IdentifyingTransition=0`,
  `QuantifyingTransition=1`, as in the converter.
- Deterministic: fixed `printf` precision, transitions sorted by ProductMz,
  peptides in input order, no unseeded randomness. Same input → byte-identical
  output (asserted by `--selftest`).

## Usage

```
OpenDIALibGen -in peptides.tsv -out library.tsv -model_dir <dir> [-threads N]
OpenDIALibGen --selftest -model_dir <dir> [-threads N]
```

- `-in`: TSV with columns `PeptideSequence`, `PrecursorCharge` (unmodified
  uppercase sequences, length ≥ 2, charge ≥ 1).
- `-model_dir`: directory containing `peptdeep_rt_dynamic.onnx` and
  `peptdeep_ms2_dynamic.onnx` — the OpenMS build tree's
  `share/OpenMS/models/` (downloaded by OpenMS's CMake; not installed).
- `-threads`: ONNX intra-op thread count, explicit; default 1. Never derived
  from `omp_get_max_threads()`.

## Build & test

```
cmake -B build -DOpenMS_DIR=/scratch/kohlbach/opendialyzer/build/openms-onnx/lib/cmake/OpenMS
cmake --build build --target OpenDIALibGen
ctest --test-dir build -R OpenDIALibGen
```

The selftest's default model dir is derived from `OpenMS_DIR`
(`<build>/share/OpenMS/models`); override with
`-DOPENDIALIBGEN_MODEL_DIR=...` if the layout differs.
