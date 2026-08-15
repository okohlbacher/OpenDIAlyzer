# OpenDIAlyzer — handoff for a coding agent

Everything needed to pick this up: what the best version is, how each phase works, where the time
goes, and which traps will cost you an afternoon.

| | |
|---|---|
| **Repository** | `https://github.com/okohlbacher/OpenDIAlyzer.git` |
| **Branch** | `odia-engine-and-openms-boundary` |
| **Tag** | **`odia-v0.4.0`** |
| **Commit** | `32c9dcd73d4f7c9c7a35f698b9284a62aec39f91` |

The benchmark binary was built from `17eb853`. Everything after it is documentation plus
`790bce9`, whose every change is gated on `-prefilter_out` and therefore inert in a normal run —
so **search behaviour at the tag is identical to the measured configuration.**

## What it achieves

Astral benchmark (`astral.mzML`, 6.4 GB, 2333 s gradient, 150 SWATH windows), 3,603,425-precursor
library, `-classifier gbt`, 224 threads, **idle** node, 2026-08-07. DIA-NN 2.0 measured on the same
node the same day. *Processed using DIA-NN.*

| | ODIA | DIA-NN 2.0 | ratio |
|---|---:|---:|---:|
| precursors @ q<0.01 | **7,008** | 9,261 | **75.7%** |
| peptides | 6,129 | 7,831 | — |
| proteins | 601 | 943 | — |
| wall | **1,353 s** (22:33) | 331 s (5:31) | 4.1× slower |
| CPU | 89,878 s | 31,174 s | 2.9× more |
| avg cores (of 224) | 66.4 | 94.2 | |
| peak RSS | 82.4 GB | 19.5 GB | 4.2× more |
| IDs per CPU-hour | 281 | 1,069 | 3.8× less |

**Do not read the wall-clock history as progress.** A recorded 1,847 s became 1,353 s purely by
moving to an idle node. DIA-NN is load-insensitive (rerun 331 s vs archived 349 s). Identification
counts *are* load-independent; timings are not. Record `uptime` with every timed run.

Peptides are not strictly comparable — ODIA's are modified-sequence level, DIA-NN's are stripped
sequences. **Precursors are the sound row.**

## Architecture in one paragraph

ODIA is **not** a reimplementation of OpenSWATH. It calls OpenMS's
`OpenSwathWorkflow::performExtraction` and `MRMFeatureFinderScoring` for extraction and
sub-scoring, and owns everything around them: the library prefilter, the input layer, a two-pass
extraction with recalibration between, an in-process classifier (no pyProphet round-trip), the FDR,
and parquet output. Optimisation effort belongs in the **harness**, not the extraction kernel.

OpenMS lives at `ext/OpenMS` — a genuine upstream checkout (`d77542d`) with **exactly 20 modified
files, all documented in `vendored-patches/`**. Do not edit `ext/` without regenerating the patch
(`scripts/openms/regen-patch.sh`); `ext/` is gitignored, so an undocumented edit leaves no trace.

## The phases, in execution order

Timings are from the measured run above. `avg cores` is CPU÷wall over the phase.

| # | phase | wall | cores | what it does |
|---|---|---:|---:|---|
| 1 | `library_load` | 108.6 s | **1.0** | TSV/PQP/OSWPQ → `LightTargetedExperiment`. ~471M allocations (78.6M transitions × ~6). |
| 2 | `dia_run_load` | 105.4 s | 45.6 | mzML/mzPeak/`.d` → `vector<SwathMap>`. On mzPeak: mmap archive + `SpectrumStore`. |
| 3 | `prefilter/*` | 206.0 s | 2.3 | see below |
| 4 | `precursor_index` | 4.0 s | 1.0 | id → precursor metadata map |
| 5 | `setup/cirt_calibration` | 20.2 s | 75.9 | CiRT anchors extracted, RT transform fitted |
| 6 | **`extract_pass1_wide`** | **305.9 s** | **138.2** | wide RT window; features scored and retained |
| 7 | `setup/compact_pass_features` | 27.7 s | 1.0 | `FeatureMap` → `ScoreRows`, OpenMS objects dropped |
| 8 | `setup/mass_calibration` | 1.0 s | 1.0 | ppm windows from pass-1 confident IDs |
| 9 | **`extract_pass2_narrow`** | **273.8 s** | **146.4** | narrowed RT and m/z windows |
| 10 | `write_parquet_bundle` | 140.0 s | **1.0** | `OpenSwathOSWParquetWriter` → `.oswpq` |
| 11 | `score_load` | 7.6 s | 1.0 | features → `OswRows` |
| 12 | `classifier_fit_gbt` | 75.2 s | 4.1 | semi-supervised, 3 iterations, 3-fold **group** CV |
| 13 | `write_scores` + `context_fdr` | 10.4 s | 1.0 | q-values; peptide/protein rollup, picked competition |

### Phase 3 in detail — the prefilter

`TransitionListEvidenceFilter` (a vendored OpenMS patch) keeps a target if **4 of its top-6
predicted fragments co-occur in the top-1000 peaks of one spectrum at ±5 ppm**. Sub-phases:
`scan_targets` 61.4 s → `build_sets` 7.9 s → `build_decoy_view` 30.1 s → `scan_decoys` 42.4 s →
`pair_and_rebuild` 29.5 s.

Two structural facts you need before touching it:

- **`build_decoy_view` (30.1 s) exists only as a workaround.** OpenMS refuses a flagged decoy
  library, so ~21M transitions are copied into an aliased, un-flagged view. A fused scan that
  handles both classes in one pass would delete this phase *and* the second spectrum sweep.
- **The pair-union is load-bearing for FDR.** It keeps a pair if *either* member passes, which is
  what makes selection label-blind. **You may reallocate pair slots, never target slots.** Ranking
  targets and decoys in separate pools and stitching the result breaks exchangeability, and neither
  the ratio guard nor a nominal q will catch it.

### Phase 6/9 in detail — two-pass extraction

Pass 1 extracts on a wide RT window and scores; `recalibrate_` fits library-RT → observed-RT from
the confident anchors; `calibrateMassFromPass_` narrows the ppm window (calibration may only
*narrow*, never widen); pass 2 re-extracts. Extraction itself is OpenMS's.

## Where the time actually goes

**Extraction is the best-behaved part of the tool** — 43% of wall at 138–146 of 224 cores. An
occupancy figure of "44 of 224" in `docs/OpenDIAlyzer-parallel-efficiency.md` was measured under
contention and **understates the tool**; ignore it.

**The bottleneck is the serial tail.** Summing every phase at ~1.0 average cores:

```
library_load            108.6
write_parquet_bundle    140.0
build_decoy_view         30.1
pair_and_rebuild         29.5
compact_pass_features    27.7
build_sets                7.9
score_load/write/context 18.0
precursor_index + misc    6.0
                        ------
                        367.8 s  = 27% of wall, single-threaded
```

> **Those 368 seconds alone exceed DIA-NN's entire 331-second runtime.** Even with infinite cores
> for everything else, ODIA cannot finish faster than DIA-NN does today. That is an Amdahl ceiling
> and it is the single most actionable performance fact in the project.

Ranked by value, none of it scientific work:

1. **`write_parquet_bundle`, 140 s, single-threaded.** Phase B of
   `docs/feature-table-memory-plan.md` replaces this writer anyway. Biggest single win.
2. **`library_load`, 109 s, single-threaded**, dominated by ~471M small allocations rather than
   bytes. `-compact_library` addresses it but has open defects (below).
3. **`build_decoy_view`, 30 s**, pure workaround — deletable by fusing the two prefilter scans.
4. **`compact_pass_features`, 28 s**, my own addition; parallelisable per feature.

Memory: peak **82.4 GB**. The allocator is the story — live data never exceeds ~33 GB while the
glibc arena reaches 186 GB and never shrinks (88% debris). `malloc_trim(0)` recovers ~9.9 GB;
~15 GB survives it as fragmentation. tcmalloc via `LD_PRELOAD` gives −44% peak for +21% wall, with
identical IDs.

## Build and run

OpenMS must be built from this branch's vendored patch — a stock OpenMS lacks
`ms2_qualifying_spectra` and `ChunkedColumn` and will not compile against `opendialyzer.cpp`.
A correctly patched install is at `/scratch/kohlbach/openms` on IBMI node `data`.

```bash
cmake -G Ninja <src> -DOpenMS_DIR=<openms>/lib/cmake/OpenMS -DCMAKE_BUILD_TYPE=Release
ninja OpenDIAlyzer
```

```bash
OpenDIAlyzer -in astral.mzML -tr library.oswpq -out out.oswpq \
  -threads 224 -classifier gbt \
  -mz_extraction_window 10 -mz_extraction_window_ms1 10 \
  -prefilter_mz_extraction_window 10 -tempDirectory $TMPDIR
```

`-in` also accepts `.mzpeak` (streaming, mmap, bounded memory), `.d`, mzXML, sqMass.

## Open defects — read before changing anything nearby

1. **`restoreRealIds_` and `setOriginalId` must be fixed as a pair.** `setOriginalId` is never
   called, so `original_id` is always `npos` and `restoreRealIds_` restores *zero* ids; the two
   defects cancel. Fixing only the obvious half activates a positional mapping that no longer holds
   after the prefilter rebuild, turning a silent no-op into a hard throw **after** a full
   extraction run.
2. **`min_child_rows` is not enforced per child** (`odia_gbt.h:383-390`): `nL` is declared, never
   incremented, discarded with `(void)nL`. Only the parent is checked, so leaves of ~4 rows are
   possible — an overfitting surface in exactly the score tail that sets the FDR threshold.
3. **The composite main score never reaches the classifier.** OpenSWATH's
   `main_var_xx_swath_prelim_score` starts with `MAIN_VAR_` and fails the `VAR_` prefix filter in
   `scoreColumnsOf_`.
4. **`-ms1_scores` needs re-measuring.** Its recorded verdict (+68 against a ±83 "noise floor")
   is stale — that floor was a fold-assignment determinism bug, since fixed.
5. **`-compact_library`**: the source still says "STILL NOT USABLE / Do not enable" while a later
   commit claims it works. Unreconciled; default is `false`.
6. **EMG peak-shape scoring is probably ON, not off** — `ff.remove(...)` does not disable it
   because `setDefaults` restores removed keys to their default, which is `true`.
7. **The pass-1 `ScoreRows` wiring has no end-to-end acceptance test.** Its criterion is
   bit-identical IDs; the closest run gave 7,008/6,129 against a recorded 6,980/6,089 — a
   no-regression check, not an isolation of that change.

## Traps that will cost you time

- **`pgrep -c OpenSwathWorkflow` silently fails** — Linux truncates `comm` to 15 chars. Match the
  full command line. This made me conclude "dead" twice and launch duplicate jobs onto one output.
- **`setsid`/`disown` did not survive SSH disconnect here; plain `nohup … &` did.**
- **OpenSWATH writes a binary progress meter into its log**; strip with
  `tr -d '\000-\010\013\014\016-\037' | tr '\r' '\n'`.
- **A failed `cp` of the binary is silent** — hard-fail if the published mtime did not advance.
- **`ldd | grep "not found"` cannot see a *wrong-version* library.** Use `ldd -r`.
- **Never compare wall clocks across node loads.** Two runs of one configuration differed 576 s vs
  824 s purely from other tenants.

## Where the knowledge is

- `docs/OpenDIAlyzer-scoring-gap-anatomy.md` — why IDs are lost at scoring
- `docs/OpenDIAlyzer-prefilter-loss-anatomy.md` — the prefilter, with two withdrawn sections
- `docs/feature-table-memory-plan.md` — the memory phasing
- `docs/odia-memory-and-io-handoff.md` — mzPeak I/O and the spectrum store
- `vault/` — Obsidian research vault (local only by decision; a copy plus sources is at
  `/ceph/ibmi/abi/dont-backup/kohlbach/odia/reference/vault`, entry point
  `00-MOC/DIA analysis MOC`)
- Reference runs and their provenance:
  `/ceph/ibmi/abi/dont-backup/kohlbach/odia/reference/HANDOFF.md`

## The caveat on every number here

Wen et al., *Nat Methods* 22:1454–1463 (2025): **no DIA search tool consistently controls FDR at
the peptide level**, with DIA-NN's true precursor FDP measured above **2.3%** at a nominal 1%.
Every count above — ours and DIA-NN's — is a *nominal* 1% number. **Gate any change on entrapment
FDP, not on identification counts.** The two most seductive fake gains are enabling `-fdr_pi0`
(the largest one-step ID jump available in this codebase, resting on a uniformity assumption this
data violates) and any decoy-suppression rule tuned on, or selected by, its effect on decoys.
