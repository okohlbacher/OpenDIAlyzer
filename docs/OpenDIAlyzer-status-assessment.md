# Where the reassembly actually stands — adversarial self-assessment (2026-07-29)

Standard applied: a thing is **PROVEN** only if it was measured on real data with the number
recorded. "Builds", "runs without crashing", and "looks right" are NOT proof.

## PROVEN — measured, with evidence

| Claim | Evidence |
|---|---|
| In-process LDA FDR reproduces pyprophet | Task A validation; π0 correction opt-in, honest true-1% default |
| mzPeak reader decodes real diaPASEF | 1,957,381 peaks from S08; m/z 96.5–1702 inside declared acq range 95.0–1705.0; per-peak 1/K0 0.60–1.57 |
| Streaming input is bounded-memory | **154–158 MB peak RSS** decoding from a 13.7 GB run, vs 500 GB–2.0 TB for the SwathFile path |
| All 24 isolation windows recovered | matches vendor `analysis.tdf` exactly (24 windows / 12 WindowGroups); was 12 before the fix |
| Mobility banding is sound | `frame 232979 = windows 232835 + outside-any-window 144`, no double counting, 3/3 frames |
| True window bounds tile exactly | `[0.5991, 0.8999]` and `[0.8999, 1.4007]` share the boundary, matching vendor scans `[34,602)`/`[602,944)` |
| Stack rebuilds on a bare node | `scripts/openms/setup_node.sh` + `vendored-patches/OpenMS/opendialyzer-openswath.patch`; selftest passes |
| Tool runs end-to-end on Bruker `.d` diaPASEF | rc=0, 16,177 features from 3,574 precursors. **Does NOT generalise — see below** |

## BROKEN or FAILING — with the evidence

0. **OpenDIAlyzer produces ZERO features on TripleTOF SWATH data.** First head-to-head on a
   complete, independent dataset (PASS00779, 5.0 GB run, complete 477,336-transition Mtb
   library, 224 threads):

   | | wall | features | precursors | %CPU | peak RSS |
   |---|---|---|---|---|---|
   | OpenSwathWorkflow | **1,070 s** | **491,493** | 95,418 | 8367% (~84 cores) | 15.3 GB |
   | OpenDIAlyzer | 1,223 s | **0** | 0 | — | SIGABRT |

   **This is the regression baseline and we do not currently reach it.** The `.osw` contains only
   the seeded library tables. Calibration WORKED (450 iRT anchors, 2350/2350 non-empty
   chromatograms, 100%), so extraction reaches the data — it is peak-picking/scoring that yields
   nothing. Two leads:
   - `Resampling spacing (0.05) is smaller than the smallest distance between data points (3.425)`,
     repeated. This run's SWATH cycle is ~3.4 s; our feature-finder parameters were ported from a
     diaPASEF setup with sub-second spacing. **Leading hypothesis**, testable in one run.
   - `Read chromatogram while reading SWATH files, did not expect that!` — the mzML carries
     chromatograms alongside spectra and our load path does not expect them.

   **Method lesson:** every prior "it works" came from Bruker diaPASEF subsets of ONE library.
   Those reported rc=0 and plausible counts. The first complete, independent dataset shows the
   tool does not generalise beyond the data it was tuned on. No amount of further subset timing
   would have revealed this.

1. **Runtime is not viable.** Three-point fit on Bruker `.d` (3,574 / 14,299 / 57,199 precursors
   → 586 / 960 / 2,697 s; features per precursor 4.53 / 4.53 / 4.53, i.e. dead linear):
   **fixed cost F ≈ 445 s, per-precursor c ≈ 35–40 ms** ⇒ **≈ 3.3 days per pass** on the
   whole-proteome library (~6.6 days two-pass). NOT the 13.6 days first reported from a single
   datapoint — fixed cost was 79% of the smallest run. Still ~1000× off DIA-NN.
   At full scale the fixed cost is irrelevant (445 s against ~282,000 s, 0.2%), so **I/O and
   loading are NOT the bottleneck** — per-precursor extraction + scoring is. Whether that is
   parallel inefficiency or algorithmic cost is UNMEASURED; do not guess, profile it.
2. **The architectural conflict is unresolved.** OpenSWATH's wave scheduler requires
   `load_into_memory` (`OpenSwathWorkflow.cpp:486`), while the streaming input exists to avoid
   whole-run residency. Bounded memory currently *costs* the parallel scheduler. No flag fixes
   this; it needs the scheduler to work against a streaming `SpectrumAccess`.
3. **Three distinct silent-abort traps** hit in one session, each exiting fast with 0 features:
   `qc:min_rsq` threshold, missing CiRT anchors in sampled subsets, and overlapping extraction
   windows (needs `-force`). Any timing harness MUST treat `features == 0` as a hard failure
   rather than a datapoint.
4. **PXD017703 download stalls** — repeated `curl (18) transfer closed` on the 69 GB zip.
5. **Box `.wiff` conversion unproven.** All prerequisites verified present; no converted file has
   come back yet. Do not count it until one does.

## NOT DONE

- **mzPeak adapter: WIRED (2026-07-29), but unproven in the tool.** `-DMZPEAK_ROOT=...` enables
  it (`WITH_MZPEAK`), `-in foo.mzpeak` dispatches to `odia::loadMzPeakSwathMaps`, and it builds
  clean on `data`. A run got as far as CiRT calibration (675 anchor sequences loaded) but has NOT
  been seen to completion. **The 154 MB streaming figure still comes from a standalone test
  binary and has never been reproduced by OpenDIAlyzer itself.** Do not quote it as a product
  number until it is.
- **OpenSwathWorkflow parity incomplete.** Missing: `out_chrom`, `out_features(_type)`, `out_qc`,
  `out_mobilogram`, `swath_windows_file`, `sort_swath_maps`, `matching_window_only`,
  `use_elution_model_score`, `enable_ipf`, `split_file_input`, `append_oswpq`,
  `keep_cached_files`, `tr_type`, `outer_loop_threads`, `use_ms1_ion_mobility`, and the whole
  `Calibration:*` / `Debugging:*` subsections.
- **CiRT/LOESS recalibration — not started.** This is what was actually asked for. `Calibration:*`
  is unregistered; `lowess` previously failed with "insufficient RT coverage" on a predicted
  library and that failure was never diagnosed, only worked around by reverting to `linear`.
- **Library online finetuning — not started.** Depends on the above.
- **No benchmark result exists.** Not one complete OpenDIAlyzer-vs-OpenSWATH run on the
  whole-proteome library.
- **DIA-NN comparison is currently impossible.** The old 43,336-precursor figure came from
  `agxt_variants.fasta`, which exists only on the unreachable spock. Nothing we have is
  like-for-like.
- **Two-pass recalibration unvalidated at scale** — only the 1-pass path has been exercised on
  real `.d`.

## Claims I made this session that were WRONG

Recorded because the pattern matters more than the individual errors.

| Claim | Reality |
|---|---|
| "9× speedup from window width" | The run had ABORTED (rsq below limit, 0 features). Retracted. |
| "Copy elision makes `Spectra` safe to return" | Only for direct-return. A named local caused a **segfault**. |
| "Partition test proves the split is correct" | Tautological — a split into the WRONG window conserves too. It passed while the split was 3.5% wrong. |
| "Registered the `irt_*` options" | They were never applied; `cp_irt` kept the analyte window. |
| "readOptions default → streaming is safer" | Would have silently disabled the wave scheduler for every mzML run. |
| "Four jobs running" | Three. `pgrep -f` matched my own command line — twice in one session. |

Failure mode in common: **plausible-but-unverified**. Every one produced output that looked
fine. The reviews (codex/kimi) caught 5 of these; running the code caught 2 more.

## Honest bottom line

The **components** are in better shape than they have ever been, and several are now proven on
real data rather than asserted. The **product** is not: the streaming input that justifies the
whole design is not wired into the tool, the runtime is ~1000× off DIA-NN, and the calibration
work that was actually requested has not been started.

Highest-value next steps, in order:
1. **Make OpenDIAlyzer produce features on PASS00779.** It is the regression baseline and we
   score 0 against OpenSWATH's 491,493. Start with the resampling-spacing hypothesis (~3.4 s
   cycle vs our sub-second diaPASEF defaults) — one run tests it.
2. Profile the per-precursor phase to separate parallel inefficiency from algorithmic cost.
   Note OpenSWATH itself used ~84 cores on mzML but only ~1.5 on Bruker `.d`, so the #22
   parallelism problem may be Bruker-specific rather than general.
3. Resolve scheduler-vs-streaming deliberately; it gates any viable runtime.
4. Only then: `Calibration:*`, CiRT/LOESS, and library finetuning.
