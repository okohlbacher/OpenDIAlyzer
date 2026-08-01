# OpenDIAlyzer input layer — mzPeak streaming (PREFERRED mode) — backlog

**Directive (user 2026-07-26): streaming from mzPeak is the *preferred* mode of
reading DIA data.** `.d` (Bruker TDF) and mzML via OpenMS `SwathFile` become the
*fallback*, not the default.

## Status
- The mzPeak C++ reader is **built and available**: `ext/mzpeak/build/libmzpeak.{a,so}`
  + `ext/mzpeak/include/mzpeak.h` (Rust core + C++ bindings; `MzPeak::open() → Index
  → Spectra`, `spectrum.h` for mz/intensity/IM, lazy `fetch(i)` decode). Proven for
  indexing (`src/odia_index.cpp`).
- The MVP engine currently loads via OpenMS `SwathFile::loadBrukerTdf` /
  `TargetedDataFileLoader` (full materialisation; the `.d` load alone is ~5 min wall,
  ~11 GB RSS on the 1/8). That was the fastest path to a working extraction — NOT the
  intended input. mzPeak streaming is the target.
- Deferred in the synthesis as "a spike until the OpenSWATH gates pass". The extraction
  MVP is now running end-to-end, so this becomes the **top input-layer item**.

## Integration path (concrete — the interface is clean)
`OpenSwathWorkflow::performExtraction` consumes a `std::vector<OpenSwath::SwathMap>`.
Each `SwathMap` = `{ sptr (OpenSwath::SpectrumAccessPtr = shared_ptr<ISpectrumAccess>),
lower, upper, center, imLower, imUpper, ms1 }`. So streaming from mzPeak needs exactly
one adapter — no changes downstream of the loader:

1. **`class MzPeakSpectrumAccess : public OpenSwath::ISpectrumAccess`** over the mzPeak
   reader. Implement:
   - `getNrSpectra()` — from the mzPeak `Index`.
   - `getSpectrumById(int id)` and the IM-aware `getSpectrumById(id, drift_start,
     drift_end)` — lazy `fetch(i)`; slice by ion mobility. **IM is at `selected_ion`
     level, not `scan`** (requirements §3) — read `precursors()[].selected_ions[]
     .ion_mobility_value` (a one-line reader fix is pending for the null `scan`-level IM).
   - `getSpectraByRT(RT, deltaRT)` and `getSpectrumMetaById(id)` — from a **once**-read
     metadata table (per-spectrum `fetch()` is ~40 ms → batch the metadata read, else
     the source is the new bottleneck; requirements §3).
2. **Build `SwathMap`s** by grouping mzPeak spectra by isolation window (m/z lower/upper)
   and IM band, one `SwathMap` per window with its `sptr` = a `MzPeakSpectrumAccess`
   scoped to that window (+ one `ms1=true` map). Reuse the window/overlap sanity checks
   already inlined in `loadDIARun_`.
3. **`loadDIARun_` prefers mzPeak**: if `in` is `.mzpeak` (or a `.d`/mzML with a cached
   mzpeak sidecar), build SwathMaps via the adapter (bounded memory, no full load);
   else fall back to `SwathFile`. Add `libmzpeak` to the `OpenDIAlyzer` CMake target
   (`ext/mzpeak/include` on the include path, link `libmzpeak` — the CMakeLists already
   notes this pending step).

## Why it matters (benefits over SwathFile)
- **Bounded memory**: hold only a per-SWATH ring buffer instead of materialising the run
  (the N1 "<128 GB vs 1.7 TB" goal; the `.d` full load is the current 11 GB / 5-min tax).
- **No 5-min Bruker load** per run — stream in RT order, demultiplexed by SWATH window.
- Realises F1/F4 + §3 of docs/OpenDIAlyzer-requirements.md (streaming input; local,
  discarded extraction) that the MVP shortcut skipped.

## Caveats / open questions
- Per-spectrum `fetch()` latency → must batch the metadata table read once (measured ~11
  min for 16k spectra otherwise).
- `getSpectraByRT` semantics: `performExtraction` calls it for RT windowing — the adapter
  must answer from indexed metadata without decoding peaks.
- Determinism + thread-safety: multiple extraction threads call one `sptr` per window;
  `MzPeakSpectrumAccess` must be safe for concurrent `getSpectrumById` (or one adapter
  instance per worker).
- Validate parity: mzPeak-streamed extraction must reproduce the `SwathFile` IDs on the
  same run before it becomes the default.

Related: docs/OpenDIAlyzer-requirements.md (F1, §3) · docs/OpenDIAlyzer-plan-synthesized.md.
