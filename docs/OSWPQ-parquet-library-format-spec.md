# The `.oswpq` Parquet assay-library format — specification and design record

**Audience.** An engineer or agent who must read, write, or extend `.oswpq` without access to this
conversation. Everything needed to implement a reader from scratch is here, together with the
measurements behind each design decision and the defects that are still live.

**Status of this document.** Descriptive, not aspirational. Every number is measured on the Astral
benchmark library (7,149,966 precursors / 78,569,077 transitions) unless stated otherwise. Where a
thing is broken, it says so.

---

## 0. Provenance — read this first

`.oswpq` is **an OpenMS format, not an OpenDIAlyzer invention.** It was introduced upstream in
OpenMS PR #9684 (`5692b27ba8`, *"full peak-map extraction (w/ ion mobility) and XIPM parquet
support"*). The schema lives in `OpenMS/FORMAT/ArrowSchemaRegistry.h`, the container logic in
`OpenMS/ANALYSIS/OPENSWATH/OpenSwathOSWParquetWriter.cpp` and
`OpenMS/ANALYSIS/OPENSWATH/TransitionParquetFile.h`.

This matters for two reasons:

1. **Compatibility is a shared contract.** Anything written here is read by OpenSwathWorkflow and
   by downstream OpenSWATH tooling. A unilateral schema change breaks other people's files.
2. **OpenDIAlyzer does not vendor OpenMS.** Modifications live as a patch
   (`patches/openms-opendialyzer.patch`), applied at build time, never committed into an OpenMS
   checkout. Section 6 lists exactly what that patch changes about this format — currently **one
   breaking on-disk change and one reader addition**.

---

## 1. Container

`.oswpq` is a **ZIP archive**, not a directory. This is the single most common implementation
mistake, and it fails late: writing to `<bundle>/library/precursors.parquet` as if it were a
filesystem path produces *"not writable for the current user"* — observed here **after a 47-minute
run had already completed**, because the path is inside an archive file. Write to a temp file, then
add the entry.

### 1.1 Entries in a library bundle

Measured layout of `bench/library_ids.oswpq` (1.35 GB total):

| entry | bytes | compression |
|---|---:|---|
| `library/precursors.parquet` | 241,326,844 | stored (0) |
| `library/transitions.parquet` | 1,109,547,075 | stored (0) |
| `library/metadata.json` | 1,058 | stored (0) |
| `<basename>.tmp.idx.json` | 124 | stored (0) |

**Entries are STORED, not DEFLATEd.** Parquet already applies per-column compression and encoding;
a second general-purpose pass over compressed bytes buys almost nothing and costs a full decompress
on every read. Storing also keeps entries individually seekable.

### 1.2 Entries added by a search run

A results bundle additionally carries, under a run-partitioned prefix:

```
runs/run_id=<run_id>/features.parquet
runs/run_id=<run_id>/score_ms2.parquet
```

The `run_id=` prefix is Hive-style partitioning, so a directory of bundles can be read as one
partitioned dataset by any Arrow-aware engine.

### 1.3 Sidecar index

`<basename>.tmp.idx.json` maps entry name → uncompressed size:

```json
{
  "library/precursors.parquet": 241326844,
  "library/transitions.parquet": 1109547075,
  "library/metadata.json": 1058
}
```

Written by `ZipArchiveFile::writeSidecarIndex`. It lets a consumer size buffers before reading the
central directory. **It is an optimisation, not a source of truth** — a reader must tolerate its
absence and must never trust it over the ZIP central directory.

---

## 2. `library/precursors.parquet`

One row per precursor (charge state), targets and decoys in the same table.

| column | Arrow type | notes |
|---|---|---|
| `precursor_id` | `int64` | join key to transitions |
| `precursor_mz` | `float64` | |
| `charge` | `int32` | |
| `library_rt` | `float64` | **the library's PREDICTION**, in iRT space — see §7.1 |
| `library_drift_time` | `float64` | ion mobility; `-1` / null when absent |
| `decoy` | `bool` | |
| `traml_id` | `utf8` | the TraML compound id; the human-facing identifier |
| `modified_sequence` | `utf8` | |
| `unmodified_sequence` | `utf8` | |
| `protein_accessions` | `utf8` | delimited list in one cell, not a `list<utf8>` |

Measured on the benchmark: 7,149,966 rows (3,603,425 target / 3,546,541 decoy).

---

## 3. `library/transitions.parquet`

One row per fragment.

| column | Arrow type | notes |
|---|---|---|
| `transition_id` | `int64` | |
| `precursor_id` | `int64` | FK to precursors |
| `traml_id` | `utf8` | **denormalised** — see §7.2 |
| `product_mz` | `float64` | |
| `charge` | `int32` | fragment charge |
| `type` | `utf8` | `"b"` / `"y"` / … |
| `annotation` | `utf8` | e.g. `"y7"` |
| `ordinal` | `int32` | |
| `detecting` | `bool` | |
| `identifying` | `bool` | |
| `quantifying` | `bool` | |
| `library_intensity` | **`float32`** | **changed from upstream `float64`** — §6.1 |
| `decoy` | `bool` | |

Measured: 78,569,077 rows (39,589,429 target / 38,979,648 decoy).

---

## 4. `library/metadata.json`

Not free-form. A structured provenance and census block under a single `openms` key:

```json
{"openms": {
  "schema_version": 1,
  "generator": "TransitionParquetFile",
  "openms_version": "3.6.0-pre-HEAD-2026-07-22",
  "build_time": "Jul 29 2026, 17:08:07",
  "tool": {"name": "OpenSwathWorkflow", "version": "..."},
  "counts": {"proteins": {...}, "peptides": {...}, "precursors": {...},
             "compounds": {...}, "transitions": {...}},
  "fragment_type_counts": {"target": {"b":…,"y":…,"other":…}, "decoy": {…}},
  "charge_counts": {"precursor": {...}, "transition": {...}}
}}
```

`schema_version` is the compatibility gate — **check it before parsing anything else.**

The `counts` block is worth more than it looks. Every count is split target/decoy, so a reader can
verify its own row counts against the writer's *before* spending an hour on extraction. A
target/decoy imbalance visible here (e.g. more decoys than targets) predicts an unreachable 1% FDR
by construction, and catching that at load time costs nothing.

---

## 5. Reading it correctly

### 5.1 The 2 GB chunk cap — the defect that shapes the whole reader

A single Arrow `StringArray` uses **32-bit offsets** and therefore cannot hold more than 2 GB of
character data. A 78.5M-row transition table with long `traml_id`s exceeds that, so **any large text
column is necessarily written as several chunks.** This is not an edge case; it is the normal state
of a proteome-scale library.

Two independent failures follow, and both were live:

* **Write.** `TransitionParquetFile` built each string column with a single `arrow::StringBuilder`
  and overflowed it. Fixed by chunking the columns.
* **Read.** `ParquetFile::readTable` calls `CombineChunks()` when a column has more than one chunk;
  `getColumn()` then returns `chunk(0)` alone while the row loop runs to `num_rows`. The column is
  read past its end and the process dies inside `basic_string::_M_create`.

**Implementation rule: never assume `chunk(0)` is the whole column.** Resolve `(global row) →
(chunk, local index)` for every access. `getColumn()` in the patched build now *throws* on a
multi-chunk column rather than silently returning a short one — failing loudly beats a truncated
read that looks like a small library.

### 5.2 `ChunkedColumn` and the refcount trap

The patch adds `ParquetFile::ChunkedColumn`, a cursor that caches the current chunk so a sequential
scan costs O(1) amortised per row instead of O(chunks).

It exposes **two** resolvers, and the difference is not stylistic:

* `resolve()` returns `std::shared_ptr<arrow::Array>` — one **atomic refcount pair per call**.
* `resolveRaw()` returns a raw `const arrow::Array*`, valid for the lifetime of the owning `Table`,
  which spans any read loop.

The accessors are called once per (row, column). On a 78.6M-row table with 13 columns that is
**~1×10⁹ `shared_ptr` copies on control blocks shared by every thread**, so the cost *grows* with
thread count instead of shrinking. Measured: a row loop went **20.1 s serial → 34.5 s on 64
threads, burning 1,135 s of CPU**. `resolveRaw()` removes it.

> Generalisable lesson: `shared_ptr` in an inner loop over a shared control block is a
> negative-scaling construct. It does not merely fail to speed up; it slows down as threads are
> added.

### 5.3 Nulls

Columns are nullable. Accessors take an explicit `(default_value, allow_null)` pair rather than
silently substituting zero — for `library_drift_time`, "absent" and "0.0" are different claims about
the data.

---

## 6. What OpenDIAlyzer's patch changes about this format

### 6.1 `library_intensity`: `float64` → `float32` — **a breaking on-disk change**

Library intensities are *relative* fragment intensities (predicted or measured, then normalised).
`float32` gives ~7 significant digits, which is already far more than the quantity carries.

Measured cost of the extra precision: **654 MB on a 78.5M-transition library — 12% of the whole
bundle** — for precision that does not exist in the input.

**This is a format change, and it is stated as one in the code.** A reader expecting `float64` here
must be updated in step. It is currently carried in the ODIA patch and **has not been contributed
upstream**; anyone consuming files written by a patched build outside this project needs to know.

### 6.2 `ChunkedColumn` + non-owning accessors

Additive: `ParquetFile.h` gains `ChunkedColumn`, `resolveRaw()`, and `arrow::Array*` overloads of
`getDouble` / `getInt64` / `getBool` / `getString`. No on-disk effect. Rationale in §5.1–5.2.

---

## 7. Design considerations and their consequences

### 7.1 `library_rt` is a prediction, not an observation

`library_rt` holds the **library's predicted** retention time. A run's *observed* RT mapped into iRT
space is a different quantity (`norm_rt`, `scores.normalized_experimental_rt`), and it lives in the
feature tables, not here.

Conflating them is not a cosmetic error. RT recalibration fits `library_rt → exp_rt` to learn the
library-to-run warp; feeding it `norm_rt` instead fits a function of `exp_rt` against `exp_rt`,
which is close to the identity and teaches the calibration nothing. Pass 2 then extracts on an
uncorrected RT axis. **Measured cost of this exact bug: 6,522 → 4,913 identifications**, and it was
misattributed for some time to the MS1 sub-scores being noisy — they were not.

### 7.2 `traml_id` is denormalised into the transition table

`traml_id` appears in both tables. It is redundant: `precursor_id` already joins them.

The redundancy is expensive. `traml_id` is the column that overflows the 2 GB string cap (§5.1),
and 78.6M copies of ~7.1M distinct strings is the single largest contributor to bundle size.

It is retained because it is the OpenSWATH-compatible layout and downstream tools join on it. **A
future version should carry `precursor_id` only and resolve the label at read time** — this is
listed as an upstream report in the project's task list, alongside the related observation that
OpenSWATH uses native *string* ids as join keys where indices would do.

### 7.3 Why the TSV cache is PQP (SQLite), not `.oswpq`

`libraryCachePath_()` returns `<library>.oswpq` as a *path*, but the cache written there has, at
times, been PQP. The reason is measured, not preference: until §5.1 was fixed, `.oswpq` could not
round-trip a library of this size **in either direction**. PQP delivered the same win — no serial
re-parse of a multi-GB TSV — over a path already exercised at this scale.

**A trap worth inheriting:** an interim version used `.cache.pqp`. Reading such a cache back at
proteome scale dies with `SIGBUS`, and because the cache is picked up *silently* on the next run, a
stale one from that experiment crashed an unrelated calibration run five minutes in. The extension
is part of the fix — a cache from an older build must not be mistakable for a current one. **Any
cache keyed on a path rather than on a format version has this failure mode.**

Cache validity is deliberately conservative: missing, unreadable timestamp, or older than the source
all mean "no".

### 7.4 Columnar-on-disk does not imply cheap in memory

The bundle is 1.35 GB. Materialising the same content as ordinary C++ objects
(`LightTargetedExperiment`) costs **38.79 GB RSS**, of which **11.29 GB is memory glibc had freed
but could not return** — ~471 million small string allocations had fragmented the arena around the
live ones.

The in-memory counterpart is `src/odia_library.h` (`CompactLibrary`), which exploits three
properties of proteomics data:

1. A tryptic peptide is a **substring** of its protein. Store the ~20k FASTA sequences once (~11 MB
   for the human proteome); a peptide costs `(protein, offset, length)` and no characters of its
   own. Non-substring peptides (variants, semi-tryptic, synthetic) are stored explicitly and the
   caller cannot tell the difference.
2. Fragment annotations repeat massively (`"y7"`, `"b3"` across millions of rows). Interned, 4 bytes
   per reference.
3. `peptide_ref` was 78.6M references to ~7.1M distinct values. As a handle: 4 bytes, no allocation.

Measured live requirement: **7.26 GB of the 32.65 GB the library path held**.

> **Live defect — do not enable blind.** `-compact_library true` makes the prefilter support **zero
> precursors** (`"prefilter supported no precursors"`, exit 6) on the Astral benchmark. The library
> loads correctly (7,149,966 peptides / 78,569,077 transitions) and then no target matches. **Not
> diagnosed.** The default is `false` and the failure is recorded in the option's own help text so
> it cannot be enabled without reading it.

### 7.5 Synthetic ids and SSO

libstdc++'s small-string optimisation inlines up to 15 characters. Ids longer than that allocate —
471M times. `CompactLibrary` assigns synthetic **base-36** ids that stay inside the SSO buffer and
allocate nothing; real ids are restored before the library is written out. Any format that uses
long human-readable strings as per-row join keys inherits this cost.

---

## 8. Reproducing every claim here

```bash
python3 -c "
import zipfile
z = zipfile.ZipFile('bench/library_ids.oswpq')
for i in z.infolist(): print(f'{i.file_size:>14,}  compress_type={i.compress_type}  {i.filename}')
print(z.read('library/metadata.json').decode())
"
```

Schema source of truth (upstream OpenMS, exact types):

```bash
sed -n '/struct OPENMS_DLLAPI OSWPrecursorSchema/,/^  };/p' src/openms/include/OpenMS/FORMAT/ArrowSchemaRegistry.h
```

The project's deviations from upstream, in full:

```bash
git -C src/OpenMS diff src/openms/source/FORMAT/ArrowSchemaRegistry.cpp src/openms/include/OpenMS/FORMAT/ParquetFile.h
```

---

## 9. Open items

| item | state |
|---|---|
| `library_intensity` `float32` not contributed upstream | patch-local; breaks naive `float64` readers |
| `traml_id` denormalised into transitions | to report upstream; `precursor_id` would suffice |
| Native string ids used as join keys | to report upstream; indices would do |
| `-compact_library true` → zero prefilter support | **undiagnosed**; default off |
| TSV cache format vs `.oswpq` reader fix | reader now chunk-aware; cache choice not re-measured since |
