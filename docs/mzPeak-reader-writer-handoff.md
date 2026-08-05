# mzPeak reader/writer — what OpenDIAlyzer needs, and what is missing

**Audience.** An engineer or agent implementing mzPeak reader/writer support to a standard-compliant
level, without access to this conversation. You should be able to work from this document plus the
upstream repository and the HUPO-PSI specification.

**Scope.** This is written from the perspective of one demanding consumer — a DIA search engine that
streams whole runs and must not materialise them — so it is biased toward *read* completeness and
toward *diaPASEF* (ion mobility). It says so wherever that bias matters.

---

## 0. Orientation: three trees, and they are not the same

| what | where | state |
|---|---|---|
| Upstream C++ library | `github.com/okohlbacher/mzpeak-openms` | **`trunk` @ `7f54871`** (2026-08-03) — what we build against |
| Reference implementation | HUPO-PSI Rust (`hupo-mzpeak`) | the conformance oracle |
| Local patch | `patches/mzpeak-per-peak-ion-mobility.patch` | **SUBSUMED** — kept for provenance only |

> **RESOLVED 2026-08-03.** The situation this section originally described — a deployed reader that
> was no upstream commit, carrying 871 lines of uncommitted work — is over. That work landed on
> trunk as `458d067` (per-peak ion mobility + isolation-window mobility limits) and `b599ba4`
> (Bruker TDF ims-compact), and trunk additionally closed the chunked/numpress gap (`36226bb`),
> three silent-corruption defects (`46050a8`), a writer seconds-into-a-minutes-column bug
> (`d9dfc03`), and added the selection/batch-read streaming API (`4f3167a`). The local patch is no
> longer applied; it is retained only so the provenance of those features stays traceable.
>
> Verified on this build: **26/26 reader tests pass** (the tag failed 3 on chunked decoding), and
> OpenDIAlyzer reads a 136 MB `.mzpeak` end to end at **0.60 GB peak RSS**.

**The lesson worth keeping**, now that the specific problem is fixed: for a period this project ran
a reader that was no upstream commit, carrying unbacked-up work whose loss would have silently
removed diaPASEF support. Sections 2.1 and 2.2 below are preserved in their original "absent
upstream" form because they are the specification of what was missing and why — useful to anyone
implementing the same capability elsewhere, and the record of how the gap was characterised before
it was closed. Read them as requirements, not as current status.

Format references:
- Specification: <https://github.com/HUPO-PSI/mzPeak>
- Paper: <https://pubs.acs.org/doi/10.1021/acs.jproteome.5c00435>
- Upstream plan with numbered items (RDR-1…25, writer phases): `docs/roadmap.md` in the C++ repo.
  **Use those IDs.** This document maps onto them rather than inventing a parallel taxonomy.

---

## 1. The consumer contract — the entire API surface actually used

A reader implementation is sufficient for this consumer if it provides exactly the following.
Everything else in the library is unused by us and can be judged on other grounds.

```cpp
MzPeak::open(path)                     // -> file handle
MzPeak::Index                          // per-spectrum metadata, WITHOUT decoding peaks
MzPeak::Spectra                        // spectrum collection
MzPeak::Spectrum::ms_level()
MzPeak::Spectrum::retention_time()     // SECONDS at our boundary; see §2.1
MzPeak::Spectrum::precursors()         // -> isolation windows + selected ions
MzPeak::Spectrum::mz()                 // const std::vector<double>&
MzPeak::Spectrum::intensity()          // const std::vector<float>&
MzPeak::Spectrum::ion_mobility()       // scalar, per spectrum   (older slice layout)
MzPeak::Spectrum::ion_mobility_array() // PER PEAK (frame layout)  -- trunk 458d067
SelectedIonInfo::ion_mobility_value
SelectedIonInfo::ion_mobility_lower_limit / _upper_limit             -- trunk 458d067
```

Trunk additionally offers a selection/streaming layer we do NOT yet use, and probably should:
`Spectra::indices_in_time_range()` (seconds), `Spectra::get_spectra_batch()` (ascending-index reads,
caller-order results) and `Spectra::extract_ion_chromatogram()` (the only selection helper that
decodes peaks, and only for spectra surviving the RT/ms-level filters). Note `Spectra` is now
non-copyable AND non-movable by design — the base class binds a fetch callback to `this`, so a copy
would dispatch through the original and dangle. One `Spectra` per thread is therefore enforced
rather than advisory.

Two hard non-functional requirements, both load-bearing:

1. **Metadata must be readable without decoding peaks.** We read every spectrum's RT, MS level,
   precursors and mobility band up front to plan extraction, then decode peak arrays only for the
   groups actually requested. Measured: 156 MB peak RSS on a 13.7 GB run, against 500 GB – 2.0 TB
   for the OpenMS `SwathFile` path that materialises everything. If metadata access forces a peak
   decode, the entire memory argument for using mzPeak collapses.
2. **Lazy decode must be thread-safe.** We call it from 224 OpenMP threads. Upstream commit
   `0de5a7e` ("make the lazy peak decode thread-safe and share it across copies") fixes this; any
   implementation must assume concurrent access to distinct spectra of one file, and to copies of
   the same spectrum handle.

---

## 2. Gap list — what is missing, ordered by what blocks us first

### 2.1 PER-PEAK ION MOBILITY — CLOSED on trunk (`458d067`); kept as the requirement spec

**Status when written:** did not exist, on any branch, at any commit. `Spectrum` exposes only a scalar
`ion_mobility()` plus `ion_mobility_type()`.

**Why a scalar is not sufficient.** mzPeak stores a diaPASEF **frame** as ONE spectrum carrying N
isolation windows over disjoint mobility ranges, plus a per-peak mobility array. (This is deliberate
in the converter: mzML has nowhere to put the mobility dimension so it splits a frame into N
spectra; mzPeak does, so the frame is kept whole with N precursors attached.) Consumers built on the
OpenSWATH model expect ONE spectrum per (frame, isolation window) with parallel m/z / intensity /
drift-time arrays. Collapsing to a scalar gives every peak in the frame the same drift time, which
(a) destroys the mobility dimension and (b) makes the N windows bleed into one another — and if only
the first precursor is read, N−1 of them vanish entirely.

**What to implement:**

```cpp
const std::vector<double>& Spectrum::ion_mobility_array() const;   // parallel to mz()/intensity()
```

**Decoding rule, and the trap in it.** The converter writes
`MeanInverseReducedIonMobilityArray` (**MS:1002816**). The library's `ArrayType` enum does not model
that term, so the column arrives typed as `NonStandard`. Matching only on the modelled
`ArrayType::IonMobility` therefore silently yields an empty array. Match on **both** the modelled
term and the column name:

```cpp
if (dim.array_type == Schema::PSI::ArrayType::IonMobility ||
    dim.name.find("mobility") != std::string::npos) { decoder.decimal(dim, mobility_); }
```

A standard-compliant implementation should extend the enum to cover MS:1002816 (and the related
mobility terms) rather than rely on the name fallback — but keep the fallback, because writers in
the wild emit the non-modelled form.

### 2.2 ION-MOBILITY WINDOW LIMITS on selected ions — CLOSED on trunk (`458d067`)

**Status when written:** `SelectedIonInfo` has `selected_ion_mz`, `charge_state`, `intensity`,
`ion_mobility_value`, `ion_mobility_type`, `parameters`. It has **no** mobility bounds.

**What to add:**

```cpp
std::optional<double> ion_mobility_lower_limit;   // per selected ion / isolation window
std::optional<double> ion_mobility_upper_limit;
```

**Why.** For diaPASEF, `ion_mobility_value` is the mobility **midpoint** of the window's scan range.
To assign each peak of a shared frame to the correct isolation window you need the *band*, not the
midpoint: the adapter gives each window only the peaks whose per-peak mobility falls inside its own
range. Without the limits the windows cannot be separated at all.

Note the asymmetry to preserve: when the limits are absent but the value is present (older files),
fall back to the value; when both are absent, fall back to the spectrum scalar. Our adapter already
encodes that precedence and any replacement must keep it.

### 2.3 SPLIT / FLAT METADATA LAYOUT — CLOSED on trunk (`ea1266e` read, `2c83369` write)

**Status upstream:** commit `b5e7c49` changed behaviour from "read as empty" to "reject loudly", and
`6dd9d0e` recognises split-metadata facets and parses `column_mapping`. Loud rejection is a genuine
improvement over silent emptiness, but it is still not support.

**What the layout is.** Instead of one nested `spectra_metadata.parquet`, the writer emits a flat
parent table plus side tables:

```
spectra_metadata.parquet                    (flat parent: index, ms_level, time, …)
spectra_metadata_scans.parquet              (scan_start_time, ion mobility)
spectra_metadata_precursors.parquet         (isolation_window struct)
spectra_metadata_selected_ions.parquet      (selected ion m/z, charge, mobility value/limits)
```

**Implementation shape** (this is what the local patch does, in four passes):

1. **Pass 1 — flat spectrum columns.** Build the per-spectrum record keyed by `index`.
2. **Pass 2 — scans side table.** Overrides retention time; carries ion mobility.
3. **Pass 3 — precursors side table.** `isolation_window` is a
   `struct<isolation_window_target, _lower_offset, _upper_offset, parameters>`; walk it in **file-row
   order** to stay aligned with `source_index`.
4. **Pass 4 — selected_ions side table.** Mobility value and (per §2.2) the limits.

**Two rules that are easy to get wrong and fail silently:**

- **Join on the index VALUE, not the row position.** Side tables are not row-aligned with the
  parent. A positional join produces plausible, wrong metadata.
- **Chunk boundaries differ per column.** Any accessor must walk chunks per column rather than
  assume a shared layout. (This is the same class of defect as Arrow's 2 GB `StringArray` offset cap
  biting elsewhere in this project.) The local patch adds a chunk-walking accessor that calls
  `fn(row_index_in_file, typed_array, row_in_chunk)` so callers never see chunk boundaries.

**Units.** `time` in the parent and `scan_start_time` in the scans table are **minutes**
(UO_0000031). Upstream commit `c61f752` cites the specification's normative MUST for this. Our
consumer boundary is **seconds** — convert once, at the edge, and say so in the type or the name.
Getting this wrong is a factor-of-60 error that will look like a calibration failure, not a units
failure.

### 2.4 CHUNKED ARRAYS AND NUMPRESS — CLOSED on trunk (`36226bb`, `3b0a04c`, `8fbd040`)

Was: `encoding.h:157` threw `"chunked array decoding is not implemented (MS:1000515)"`, failing 3
of 24 tests. Trunk decodes the chunked layout including delta (MS:1003089) and MS-Numpress
(MS:1002312/1002314), and validates each chunk against its declared `chunk_end`. **26/26 tests now
pass.** The blocker analysis below is retained because it is a good record of how the gap was
characterised before it was closed.

**Upstream's own analysis of the blockers** (from `docs/roadmap.md`, worth reproducing because it is
precise):

- array-index paths like `chunk.mz_chunk_values` do not match the physical leaf `…list.item`
  (`parquet.cpp:243`);
- `decode_array` assumes a single column (`encoding.h:74`);
- `buffer_format_from_string` never parses `"chunk_transform"`, so it silently falls back to `Point`
  (`buffer_format.cpp:37`).

**Split the work** as upstream suggests: (path mapping + multi-column decode) → (delta,
MS:1003089) → (numpress, MS:1002312/1002314 via vendored ms-numpress).

**Relevance to us:** unknown and worth determining early. Our production files are converter output
whose layout we have not audited against this. If they are chunked, the current reader throws rather
than mis-decodes — which is the safe failure, but a failure.

### 2.5 BRUKER TDF "ims-compact" — CLOSED on trunk (`b599ba4`)

This is the least standard and most easily missed item.

**The layout stores NO m/z array.** A nonstandard `tof` **Int32** column stands in for it, and m/z is
reconstructed from the index calibration:

```
mz(tof) = (a + b * tof)^2
```

with `a`, `b` from the index calibration (`ims_a_`, `ims_b_`, plus a validity flag). It also carries
a per-peak `MeanInverseReducedIonMobilityArray`.

Without this, the m/z + intensity filter matches nothing and **every spectrum reads as empty** — a
silent, total data loss that looks like an empty file rather than an unsupported layout.

**Second trap in the same layout:** Bruker TDF stores **intensities as Int32**, and `decimal()`
throws on integer types. Dispatch on the column's declared type:

```cpp
const Util::Type it = dim.type_or_throw();
if (it == Util::Type::Int32 || it == Util::Type::Int64) { decoder.integer(dim, intensity_); }
else                                                    { decoder.decimal(dim, intensity_); }
```

Do **not** "fix" this by widening `decimal()` to accept integers: that instantiates the delta
estimator for 8-bit types, which Boost's `median()` does not compile for. The dispatch is the fix.

### 2.6 `row_count()` vs `record_count()` — upstream RDR-11

`record_count()` is the **entity** count (spectra), taken from a KV shared by the profile and
centroid files. A centroid-only run has a profile file with **0 rows but a non-zero spectrum count**,
so sizing buffers from `record_count()` is wrong. The local patch adds:

```cpp
std::size_t row_count() const;   // actual PEAK ROWS, from parquet footer metadata; reads no data
```

Upstream RDR-11 covers the related "missing count KV → real count, not 0/max" fallback.

### 2.7 Type system — upstream RDR-4, and upstream says it is not really deferrable

Unsigned indices currently work only for small fixture values. Widening is needed for Phase 5
`mz_delta_model` lists, Phase 6 nested/byte arrays, and Phase 7 aux arrays. Pull it earlier than the
"deferred" label suggests — upstream's own review correction says exactly this.

---

## 3. Writer side

### 3.0 REQUIRED: chromatogram writing — the writer cannot do it, and we need it [CORRECTED — see note below]

**This is now a blocking requirement, not an assessment.** OpenDIAlyzer extracts ~4.46M
chromatograms per run and currently emits **none of them** — its bundles carry features and scores
but no traces, so nothing downstream can plot a peak, run QC, or verify an identification by eye.
That capability was dropped only because OpenMS's `ChromatogramPeak` costs 16 B/point (85 GB on the
reference run). At 1 B/point it is ~5.4 GB and affordable, so the storage problem is solved and the
remaining obstacle is that **mzPeak has nowhere to put them.**

Current writer surface (`include/mzpeak/writer.h`):

```cpp
void write_spectra_directory(const std::filesystem::path&, const std::vector<SpectrumData>&);
void write_spectra_archive  (const std::filesystem::path&, const std::vector<SpectrumData>&);
```

`SpectrumData` is `{mz, intensity, centroid, ms_level, retention_time, polarity, id}`. The emitted
tables are `spectra_data`, `spectra_peaks`, `spectra_metadata`, `spectra_metadata_{precursors,
scans,selected_ions}`. **There is no chromatogram table on any path**, and `grep -i chromatogram`
over `writer.h` + `writer.cpp` returns nothing. The reader models chromatograms
(`include/mzpeak/chromatogram.h`, `chromatograms.h`) and the bundled fixtures contain
`chromatograms_data.parquet` / `chromatograms_metadata.parquet`, so **the format supports them and
only the writer is missing.**

**What is needed**, mirroring the spectra path:

```cpp
struct ChromatogramData {
  std::vector<double> time;            // SECONDS at the API boundary; the file stores minutes
  std::vector<float>  intensity;
  std::optional<std::string> id;       // e.g. the transition's native id
  std::optional<double> precursor_mz;  // for a SRM/DIA transition
  std::optional<double> product_mz;
  std::optional<std::string> type;     // TIC / BPC / SRM / SIC
};
void write_chromatograms_directory(const std::filesystem::path&,
                                   const std::vector<ChromatogramData>&);
void write_chromatograms_archive  (const std::filesystem::path&,
                                   const std::vector<ChromatogramData>&);
```

Plus the ability to write chromatograms **alongside** spectra into one bundle, since a search result
wants both.

**Two things the implementation should exploit**, both established in
`~/Downloads/compact-chromatogram-storage-handoff.md`:

1. **Chromatograms from one SWATH window share an identical time axis.** Storing time per
   chromatogram is the single largest waste in the naive form (42.7 GB of 85 GB on our run). A
   columnar format should factor the axis out, exactly as the reader's point layout factors out
   coordinates.
2. **Intensity does not need `double`, or even linear spacing.** Measured on real dynamic range
   (1e6 apex over 1e2 baseline), one byte log-spaced preserves the baseline 680x better than two
   bytes linear, because the scores that consume chromatograms are scale-invariant and read shape.

**Units warning, repeating §2.3 because it bites hardest here:** the file stores retention time in
**minutes** and the reader converts to seconds. Writer commit `d9dfc03` fixed the writer storing
SECONDS into that minutes column — a silent 60x error. Any chromatogram writer inherits that trap.

### 3.1 The rest of the writer

We do not currently write mzPeak spectra, so this part is assessment rather than requirement.

**Where upstream is:** writer P0/P1a/P1b done on `writer_test` — point directory, zip-STORE archive,
`spectra_metadata` — and **T2 cross-implementation PASS** (the Rust reference reads C++ output with
matching values). Remaining: P1c centroid/peaks split (representation flag; centroids to
`spectra_peaks.parquet`, profile to `spectra_data.parquet`, with `number_of_peaks`/representation set
accordingly), P2 null-marking + m/z delta model, P3 chunked + numpress.

**The null-marking model (writer P2 / reader RDR-5) is the one to get right**, because reader and
writer must agree exactly:

- Writer null-marks flanking zero-intensity points, fits a WLS model
  `δmz ~ β0 + β1·mz + β2·mz²`, and stores `mz_delta_model` plus transforms **MS:1003901 / MS:1003902**.
- Reader reconstructs null m/z from the model (segment median + regression), and **null intensity
  reads as 0** — not delta-interpolated. Upstream `ee011cc` fixes precisely that, and
  `c3aee8f` fixes reconstructing each null-marked m/z from *its own* adjacent run rather than a
  neighbouring one. Both are silent-corruption bugs: the file reads, the numbers are wrong.

Until RDR-5 lands, profile m/z **interiors are 0** where null-marked (39,968 nulls in the bundled
`small` fixture). Upstream's Phase-2 tests deliberately assert only centroid values and profile
endpoints, never profile interiors, so nothing passes falsely — a good example to imitate rather
than paper over.

---

## 4. How to validate — do not invent a new harness

Upstream has one (`docs/e2e-testing.md`), and its guiding principle is right: **validate
semantically, not by byte-diffing Parquet.** parquet-cpp and parquet-rs differ in `created_by`, ZSTD
streams and page layout, so byte identity is neither achievable nor the goal.

| # | pipeline | oracle | validates | status |
|---|---|---|---|---|
| T1 | C++ write → C++ read | input data | writer+reader agree | done |
| T2 | C++ write → **Rust** read | Rust reference | writer conformance | **PASS** |
| T3 | C++ read → C++ write → C++ read | first read | reader↔writer idempotence | done |
| T4 | Rust write → C++ read | pyarrow ground truth | **reader conformance** | partial — blocked on the gaps above |
| T5 | mzML → Rust → C++ read → C++ write → Rust → mzML | original mzML | whole stack | future |

**T4 is the one that matters for everything in §2.** The Rust oracle:

```bash
cd hupo-mzpeak
cargo build --release --example read_spectrum
target/release/examples/read_spectrum <file.mzpeak> <index>
```

Note the `convert` example is mzML→mzpeak **only**; `read_spectrum`/`read` are the reader oracles.

Invariants for every change: keep `meson test` green, add a regression test per fix, re-run
`scripts/e2e_cross_impl.sh` after any writer change, compare semantically.

**Add one test class that does not yet exist:** a *concurrency* test. Lazy decode is now shared
across copies; nothing in the suite exercises it from many threads, and our consumer does exactly
that at 224.

---

## 5. Build notes (environment, not design)

Both of these cost time and neither is a real requirement:

- `meson.build` requires `libzip >= 1.11.4`. The bound dates from the initial import; the previous
  reader built and ran against **1.11.2** (the old meson log shows `found: YES 1.11.2` accepted).
  Newer meson enforces it. Either relax the bound or pin a newer libzip — but know it is drift, not
  an API need.
- `include/mzpeak/util/enumerable_proxy.h:58` declares
  `Iterator(const Iterator&&) = default;` and `Iterator& operator=(const Iterator&&) = default;`.
  A defaulted move constructor/assignment **cannot take a const rvalue reference**; the implicit
  signatures are `Iterator(Iterator&&)` / `operator=(Iterator&&)`, so `= default` on the const&&
  overloads is ill-formed. Clang (their macOS CI) accepts it; GCC rejects it, and it only bites via
  the newly-added `WavelengthSpectrum` instantiation. **Report upstream** — the fix is dropping the
  two `const`s.

---

## 6. Suggested order of work — REVISED after trunk @ 7f54871

Items 1–4 of the original list are done upstream. What remains:

1. **Adopt the streaming/selection API.** Our adapter still walks spectra itself; trunk offers
   `indices_in_time_range()` + `get_spectra_batch()` + `extract_ion_chromatogram()` with lazy decode
   and ascending-index reads. This is the largest remaining win and it is ours to take, not
   upstream's to build.
2. **Report the two build issues upstream** (§5): the `const Iterator&&` defaulted move constructor
   (two `const`s to delete, gcc rejects it), and the `libzip >= 1.11.4` bound that 1.11.2 satisfies
   in practice.
3. **RDR-4 minimal types** and **RDR-11 `row_count()`** (§2.6, §2.7) — the remaining reader items.
4. **RDR-5 null reconstruction** paired with writer P2, so the two agree by construction.
5. **T4 reverse-cross conformance**, and add the concurrency test that still does not exist: lazy
   decode is shared across copies and nothing exercises it multi-threaded, while we call it from
   224 threads.

Deliberately *not* recommended: reimplementing the format from the paper. The Rust reference is the
conformance oracle and the C++ roadmap is accurate about its own gaps; both are better inputs than a
fresh reading of the specification.


## RESOLVED (was: BLOCKING) -- the measurements below were taken against a STALE SNAPSHOT
##
## The reader on current trunk (c230228) keeps a 2-deep decoded row-group cache and a typed
## EqualityScan. Re-measured on the same file: 3,000 sequential spectra in 7.86 s (382/s)
## against >900 s before, ~115x. Scattered access is 46/s, 8.3x slower than sequential, which
## is the reader's rule 3 showing up directly. Per-peak throughput matches the maintainer's
## own fixture, so there is no residual to report.
##
## What follows is kept for the record of what the pre-cache cost model looked like, and
## because the ACCESS RULES it motivated are still the rules: ascending order, never random,
## contiguous ranges per worker, one Spectra per worker.
##
## SUPERSEDED: bulk spectrum reads. A point query per spectrum makes a full scan quadratic

Measured on astral.mzpeak (3.09 GB, 307,590 spectra) with the reader at trunk:

    3,000 SEQUENTIAL spectra did not decode in 900 s   -> under 3.3 spectra/s
    the same analysis reading mzML parses the whole file in 104 s

At that rate one pass over the run's spectra is ~15 h single-threaded. Raising concurrency
trades it straight back for memory: 224 decoders peaked at 105.8 GB, 16 decoders at
30.6 GB but did not finish in 40 minutes. Time and memory are in direct opposition and
neither end is usable.

### Why

`spectra_peaks.parquet` holds 512,278,842 rows in 489 row groups, averaging 1,047,605 rows
per group. A spectrum is ~1,666 of those rows.

`Spectra::fetch(index)` issues a POINT QUERY for one spectrum:

    signals_->select(dims_, signals_->index().eq(index_))      // spectrum.cpp:56

and `Signals::select` builds a fresh Planner and Executor for every call
(data/signals.cpp:165-187), consulting parquet statistics and the row-group page index.
So each spectrum costs a plan across 489 row groups plus a page read out of a 1,047,605-row
group -- roughly 630x more decoded rows than the caller asked for -- and nothing is cached
between calls: 3,000 sequential spectra span ~5 row groups but triggered ~3,000 row-group
decodes, about 37 GB of redundant work.

Iteration does not help. `Spectra` is an `EnumerableProxy<Spectrum>` whose Iterator holds
`fetch_t = std::function<V(std::size_t)>` and calls it per element, so iterating is the same
point query. `get_spectra_batch` sorts indices ascending "so file access is sequential" --
an acknowledgement of the cost -- but still calls `fetch()` per index, so it issues N plans.

### The batch accessor was tried, as intended, and is ~2x -- not the ~420x needed

Measured directly, reading in the file's own ascending order in batches of 2,048:

    Spectra::get_spectra_batch, file order      7 spectra/s
    Spectra::operator[], sequential           ~3.3 spectra/s
    ->  307,590 spectra at 7/s = 12.2 HOURS for one pass
    ->  mzML parses the same run in 104 s

Sorting helps by about 2x, which is what avoiding some seek cost buys, and no more. It
cannot help further: get_spectra_batch calls fetch() per index internally
(spectra.cpp:186), so the per-spectrum plan is still issued 307,590 times. The ordering was
never the dominant term; the per-call plan is.

### What is needed, and why it cannot be done by the caller

A bulk sequential read that plans ONCE and streams row groups in order, yielding spectra as
it goes: 489 planned reads instead of 307,590, with peaks arriving already grouped. The
query API can express it -- `index().ge(a).and_then(index().le(b))` -- but `Spectra::data_`
and `Spectra::peaks_` are PRIVATE, so no caller can issue a range select. Either

  1. a public bulk/streaming entry point on Spectra (preferred: preload metadata once, then
     walk row groups and hand back spectra), or
  2. a decoded-row-group cache inside Spectra so sorted access reuses a group across the
     ~630 spectra that share it.

Until one exists, mzPeak input cannot feed a full analysis at this scale, and ODIA
benchmarks stay on mzML. Everything on the ODIA side is already done: one shared mmap
(zero descriptors), one shared Index, a compact 8 B/peak spectrum store, and a bounded
decoder count.


## CORRECTION (maintainer handoff item 5): the writer CAN write chromatograms

The section above claiming the writer cannot emit chromatograms is false as of current
trunk: `write_run_directory(RunContents{...chromatograms...})` exists. What the writer
refuses is PRODUCT-BEARING (SRM/MRM) chromatograms -- a product-bearing type throws at write
time, because the Q3 product facet is unwritable in the format today (`has_unreadable_product`
on read).

That matters for us specifically, because the ~4.46M transition chromatograms this tool would
write ARE product-bearing. They cannot be stored as SRM chromatograms as-is. Three options,
and this is a format-direction decision rather than an implementation one:

  1. Plain chromatograms, transition identity in the `id` ("500.2->184.1") and Q1 in the
     precursor/selected-ion facet. Works with the writer today; Q3 lives only in the id string.
  2. Push the product facet upstream so SRM chromatograms become first-class. Larger: spec,
     reference implementation and C++ writer all have to move together.
  3. Keep transitions in our own sidecar and put only TIC/BPC in mzPeak.

Undecided. It should be settled before the writer is wired, because it is the one place the
stated goal and the format genuinely collide.

## READ PERFORMANCE, settled (mzPeak c2ffd34, tag perf-no-page-index-2026-08-03)

Root cause of the slow reads was that astral carries NO Parquet page index -- verified here
independently: all three columns of spectra_peaks.parquet report column_index_offset=None and
offset_index_offset=None, while the row group DOES declare SortingColumn(column_index=0,
descending=False). Without a page index the reader could not skip pages, so it full-scanned a
1,047,605-row group per spectrum. Trunk now binary-searches the declared-sorted index.

Measured here on astral (307,590 spectra, 512,278,842 peaks, 489 row groups):

    1 thread     361 -> 4561 spectra/s     12.6x
    4 workers          7195 spectra/s      best
    8 workers          6065 spectra/s
   16 workers          3915 spectra/s

One qualification for planning: the maintainer handoff predicts ~7,600 spectra/s
single-thread; measured here is 4,561, about 60% of that. The practical conclusion is
unaffected -- a single thread now materialises the whole run in 67 s against 104 s for the
mzML parse ON 224 THREADS, and 43 s at the 4-worker optimum.

Note the knee MOVED from 8 workers to 4 when the reader changed: with the per-spectrum scan
gone the job becomes memory-bandwidth-bound sooner. A worker count carried over from the old
reader is quietly wrong.
