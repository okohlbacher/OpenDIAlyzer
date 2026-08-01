# Reference DIA benchmark datasets

Every accession, file count and byte total below was read from a live PRIDE/ProteomeXchange
API call or FTPS listing, not from memory. Two items are explicitly flagged as unverified.

Persistent copies live in `/ceph/ibmi/abi/data/dia_reference/` — our group directory. Never
`/ceph/ibmi/it` (IT's tree, `itstaff`-only) and never `/home` (at quota). Hot working copies
get staged to node-local `/scratch` by the jobs that need them.

## DIA-NN — Demichev et al., Nat Methods 2020
Data availability names PXD014690 (new) plus PXD005573, PXD002952, PXD010529, PXD006722.

| Accession | Role | Instrument | Files | GB | Format |
|---|---|---|---|---|---|
| PXD014690 | this paper | TripleTOF 6600 SWATH | 57 | 28.98 | `.wiff` |
| PXD005573 | HeLa depth (Bruderer 2017) | Q Exactive HF | 239 | 1266.58 | `.raw` |
| PXD002952 | LFQBench | TripleTOF 5600/6600 | 273 | 469.72 | `.wiff` |
| PXD010529 | yeast kinase KO | TripleTOF 5600 | 795 | 396.29 | `.wiff` |
| PXD006722 | Specter comparison | Q Exactive | 225 | 445.23 | `.raw` |

## DIA-NN 1.8 / diaPASEF — Demichev et al., Nat Commun 2022
PXD029836 (new) plus PXD017703, PXD022216, PXD013658.

| Accession | Instrument | Files | GB | Format |
|---|---|---|---|---|
| PXD029836 | timsTOF Pro 2 diaPASEF, HeLa dilution | 5 | 92.12 | `.d.zip` |
| PXD017703 | timsTOF Pro diaPASEF | 16 | 825.89 | `.d.zip` |
| PXD022216 | timsTOF Pro, 50 CLL patients | 55 | 260.94 | `.d.zip` |
| PXD013658 | Q Exactive HF-X | 167 | 1251.85 | `.raw` |

> UNVERIFIED: a secondary source calls PXD013658 a human/*Arabidopsis* two-species FDR
> benchmark; PRIDE's own metadata says human only. Do not rely on the two-species framing.

## OpenSWATH — Röst et al., Nat Biotechnol 2014
**No PXD.** The paper states: data are at the PeptideAtlas repository, **PASS00289**
(207 files, 683.94 GB; FTPS, creds published on the PASS_View page). The M-score/decoy
validation used the SGS data itself — the manual ground truth ships inside as
`OpenSWATH_SM3_GoldStandardManualResults.csv`. There is no separate decoy accession.

SGS design: 422 SIS peptides, 10 dilutions x 3 backgrounds x 3 replicates = 90 runs,
TripleTOF 5600. `/SGS/mzxml/` is 466 GB — but **water background alone is 25.6 GB** and is
already an open format.

**PASS00779** — the official OpenSWATH tutorial data (M. tuberculosis, TripleTOF 5600):
3 x `.mzML.gz` totalling 10.46 GB, ALREADY mzML, with assay libraries and reference results.

## LFQBench — Navarro et al., Nat Biotechnol 2016
PXD002952. 273 files / 469.72 GB, of which raw is 114 files / 230.15 GB; the other 231 GB is
Spectronaut `.htrms` intermediates (skip). Cheapest complete A-vs-B arm:
`HYE110_TTOF6600_64fix`, 6 runs, 14.98 GB.

## diaPASEF — Meier et al., Nat Methods 2020
PXD017703, confirmed via the bidirectional ProteomeXchange<->PMID 33257825 link (paper is
paywalled, so the data-availability sentence itself was not read at the primary source).

The 16 zips are NOT equal value — read from their ZIP central directories:
- `TwoProteome_diaPASEF_raw.zip` 68.92 GB — **6 `.d` runs, KNOWN 3x yeast ratio ground truth**
- `HeLa_Evosep_diaPASEF_RAW.zip` 11.08 GB — 9 short-gradient runs, ~1.2 GB each
- `*_pqp_library.zip` 0.04-0.13 GB — **prebuilt PQP libraries, so the ~650 GB of
  `*_Library_RAW.zip` ddaPASEF runs are unnecessary**
- `*_pyprophet_export.zip` 0.10-0.34 GB — published OpenSWATH/pyprophet results to score against

## PXD028735 — Van Puyvelde et al., Sci Data 2022
LFQBench three-species design re-acquired on five platforms INCLUDING timsTOF diaPASEF —
i.e. HYE124 ground truth in diaPASEF, the single best validation target. 844 files, 3.43 TB;
the diaPASEF arm is 744 GB (~21.8 GB/run), minimal A-vs-B ~131 GB.
> Bug in the deposition: its SDRF `comment[file uri]` points at PXD010000, an unrelated
> de-novo set. Ignore those URIs.

## Fetch plan

**Tier 1 (~106 GB, in progress):**
1. PASS00779 (10.6 GB) — already mzML, zero-friction smoke test with reference results.
2. PXD017703 `TwoProteome_*` (69.3 GB) — diaPASEF, known ratio ground truth, matches our data
   type, ships library + published result.
3. PXD017703 `HeLa_Evosep_*` (11.2 GB) — short gradients, fast iteration loop.
4. PXD002952 `HYE110_TTOF6600_64fix` (15.0 GB) — cheapest complete LFQBench arm.

**Deliberately NOT fetching:** PXD005573 / PXD006722 / PXD013658 / PXD010529 (0.4-1.3 TB each,
all third-party comparisons rather than ground truth); the four `*_Library_RAW.zip` in
PXD017703 (650 GB; the PQP libraries carry the same information); PXD002952's `.htrms` (231 GB).

## Download gotchas (both cost a retry cycle)
- **PeptideAtlas** drops the login INTO the dataset directory, so the accession must NOT be
  appended to the URL; `-k` is required (cert does not match the host); and
  `--ftp-method nocwd` is required because curl otherwise sends a `CWD` the chrooted account
  refuses with "Server denied you to change to the given directory" — even though `LIST` works.
- **PRIDE** v2 API returns HTTP 200 with an EMPTY body; use v3
  (`/ws/archive/v3/projects/{acc}/files`). The archive path is date-based and easy to guess
  wrong — PXD017703 is under `2020/12`, not `2020/11`.

## PXD014690 as the primary benchmark (chosen 2026-07-29)

"Microflow SWATH data analysed with DIA-NN" — Demichev et al., Nat Methods 2020. The DIA-NN
paper's OWN data, so their published numbers are the reference to beat.
Human K562 + human plasma + yeast, TripleTOF 6600, SWATH + gas-phase fractionation.
**57 files / 28.98 GB: 25 runs as SCIEX `.wiff` + `.wiff.scan` pairs, plus 7 TSVs.**

> NOTE: `PXD014960` (one digit transposed) is a real but unrelated deposition — an apple
> fruit-russeting study on an Agilent 6520A Q-TOF. Verified via the PRIDE v3 API before
> fetching anything.

### Conversion path — the hard part
`.wiff` is vendor-closed; the only documented reader is the managed `Clearcore2.*` .NET
assembly family that ProteoWizard ships. mzPeakConverter has a native SciEX lane
(`src/sciex.rs`) that reaches Clearcore2 through a C# shim by runtime reflection, so it can
produce BOTH mzML (`--to mzml`) and mzPeak from one tool. It needs:

1. `SciexGlue.dll` — **already built** at `glue/sciex/bin/Release/net8.0/`.
2. A .NET 8 runtime — installing via conda-forge `dotnet-runtime` (not present on `data`).
3. The Clearcore2 DLLs — from a ProteoWizard install, pointed at by `MZPC_PWIZ_DIR`; the
   reader expects them under `<MZPC_PWIZ_DIR>/vendor_api/ABI`. **Not yet obtained.**
   ProteoWizard's Linux distribution is behind a click-through licence on
   proteowizard.sourceforge.io, so this step needs a human to accept the terms — it is not
   something to automate around.
4. Build the converter with `--features sciex`.

There is no macOS build of the Clearcore2 stack, so this conversion must happen on `data`.
