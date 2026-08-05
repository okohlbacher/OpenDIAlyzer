# ODIA: reading mzPeak, materialising spectra, compacting the feature tables

Handoff, 2026-08-05. Three connected pieces of work on the input and memory path. Every number
here is measured on the Astral benchmark (`astral.mzpeak`, 307,590 physical spectra,
512,278,842 peaks, 128-core node) unless it says *projected*.

The through-line: **ODIA's memory and runtime problems were mostly not where they looked.** In
each of the three areas the first plausible culprit was wrong, and the measurement that settled
it is recorded alongside the fix so nobody re-litigates it from the same wrong prior.

---

## 1. Reading mzPeak

### The starting position was a scandal, and it was ours

mzPeak is columnar and lazily decoded. It should be the cheap input. It was not:

| input | peak RSS | wall |
|---|---|---|
| mzML (`SwathFile`, readOptions normal) | 74.5 GB | 20:16 |
| mzPeak, as first written | **127.6 GB** | still in the prefilter at 83 min |

71% *more* expensive than XML. Every bit of that was consumer-side.

### Four rules, from the mzPeak maintainer, that the reader actually enforces

1. **Metadata first** — read all spectrum metadata before touching peaks.
2. **Ascending order** — never walk indices backwards.
3. **Never random-access** — a point query with no page index becomes a full row-group scan.
4. **Contiguous ranges, one `Index`/`Spectra` per worker.**

Rule 4 is the one that is easy to get subtly wrong, and we did: **the decoded row-group cache
belongs to the `Index`, i.e. to each `open()` — not to the `Spectra`.** Deriving every worker's
`Spectra` from one shared `Index` gives them all *one* 2-deep cache, and workers reading distant
ranges then evict each other on every group. Measured that way: ~1400% CPU, RSS creeping
0.1 GB/min, `populate()` unfinished at 41 minutes where 14 workers at the measured rate should
have taken about a minute.

### Descriptors: 1,060 → 4, and the fix that did *not* work

At 224 threads the run held **1,060 open descriptors on one file**, blew the 1024 limit, and
aborted with `failed to open zip archive Can't open file: Too many open files`.

Centralising the `Index` — the obvious fix — moved it to **911**. Not a fix. **908 of those 911
descriptors were `astral.mzpeak` itself**, opened by the per-thread `Spectra`, not by the
`Index`. What has to be shared is the *bytes*.

`src/odia_mmap_archive.h` maps the file **once**, parses the ZIP central directory **once**, and
hands out `MzPeak::IO::File`s that are plain spans of that mapping. The descriptor is closed
right after `mmap` — a mapping keeps the file alive on its own, so **steady state is zero
descriptors** however many workers run. Threads share one read-only mapping with no locking; the
page cache does the concurrency.

Two format facts it relies on, both **checked rather than assumed**:
- Members must be **STORED**. A DEFLATE member cannot be a span of a mapping, and quietly
  returning its compressed bytes would surface far away as a parquet parse error. Rejected by name.
- The data offset comes from the **local** header, whose name/extra lengths are permitted to
  differ from the central directory's.

### The contract that hangs the reader if you get it backwards

`File::read` must return an **empty optional at EOF, not `optional{0}`**. `ZipBuffer::read`
returns `{}` whenever `zip_fread` yields ≤ 0, and `index.cpp:94` loops *while the optional is
engaged*. Returning `optional{0}` is always truthy and hangs the reader forever — the access
test went from 1.03 s to a timeout with no output, because the first metadata read never
terminated. This is written on `SpanFile::read` in the file.

### What each step was actually worth

Four things were confounded. Separated:

| | factor | whose |
|---|---|---|
| stale Aug-3 snapshot with no row-group cache | 115× | **not a real problem** |
| no page index → 1M-row scan per point query | 12.6× | maintainer fix `c2ffd34` |
| scattered → ascending access (rule 3) | 8.3× | ours |
| one shared cache → per-worker `Index` | — | ours; required for any scaling |
| worker count 224 → 4 (knee re-measured) | 1.8× over 16 | ours |

**Only two were ODIA defects. Two were me measuring the wrong thing** — a snapshot that was not
a git checkout and had zero cache code in it, and a worker count carried over from a superseded
reader. A "BLOCKING" finding built on the stale snapshot (3.3 spectra/s, "12.2 hours") was
**retracted in full**.

---

## 2. Materialising the spectra

### The five-pass discovery

The run walks the spectra **five times**, not once:

1. `prefilter/scan_targets` — every MS1 + MS2 spectrum
2. `prefilter/scan_decoys` — the same spectra again
3. `setup/cirt_calibration` — anchor extraction over the SWATH windows
4. `extract_pass1_wide`
5. `extract_pass2_narrow`

**mzML concealed this.** It parses once and leaves everything resident, so passes 2–5 were walks
over RAM. mzPeak streams honestly, so each pass genuinely re-decoded parquet. The format did not
get slower — **it stopped hiding four redundant passes.**

### `SpectrumStore` (`src/odia_spectrumstore.h`)

Decode once, into one arena. A peak is float32 m/z + float32 intensity = **8 B**, against
`OpenSwath::Spectrum`'s 16 B in two `vector<double>`, and in **one** arena for the whole run —
two allocations rather than two per spectrum, keeping ~1.2M allocations out of an allocator
whose retained pool is the largest single item at peak.

```
materialised 307,590 physical spectra, 512,278,842 peaks in 30.5 s
(10,086 spectra/s, 4 workers) into 3.83 GB (7.63 GB as double pairs, 1.99x)
```

Against mzML's 104 s parse **on 224 threads** leaving ~9 GB resident. 10,086 spectra/s is above
the maintainer's predicted ~7,600.

What it buys, measured on the prefilter (the first consumer):

| | store | mzML |
|---|---|---|
| `prefilter/scan_decoys` | 29.0 s | 30.7 s |
| prefilter total | 125.2 s | 115.3 s |

i.e. **the same speed as having every spectrum resident, from 3.83 GB instead of ~9 GB.**

### Precision, measured not assumed

| m/z | float32 quantisation | vs the 10 ppm window |
|---|---|---|
| 200 | 0.0381 ppm | 262× finer |
| 500–2000 | 0.0305 ppm | 328× finer |
| worst over 4,000 random peaks | 0.0590 ppm | 169× finer |

Two orders of margin. This is *not* the explanation for the open 22-peptide delta (§4).

### The memory panic that inverted on measurement

`populate()` was blamed for a 105.8 GB peak. Wrong:

| configuration | peak |
|---|---|
| streaming, per-thread `Index` | 127.6 GB |
| streaming, shared `Index` | 123.7 GB |
| **with** the spectrum store | 105.8 GB |
| same, 8 workers | 21.8 GB |

Every *uncached* configuration was already 120+ GB. The store did not cause the peak — it
**reduced** it. What scales is the **thread count**: each decoding worker holds a `thread_local
MzPeak::Spectra` with parquet row-group buffers, and 224 of those is the footprint.
`-mzpeak_decode_threads` (default 4) bounds it for the *entire run*, because `populate()` is the
only phase that decodes — afterwards `decode()` serves from memory and no `Spectra` is ever
constructed again. 105.8 GB → **26.5 GB**, later 30.6 GB with the final worker count.

### Two silent correctness bugs found by adversarial review (Kimi), both now fixed

Both would have bitten on real data and neither would have announced itself.

**The back-filled drift array.** The arena's per-peak mobility `dt_` is back-filled with `-1`
as soon as *any* spectrum supplies drift. `drift(i)` keyed on "the arena is non-empty", so
afterwards **every** spectrum reported a non-null drift array full of `-1`. In `decode()` a
mobility band test (`pd[k] >= lo && pd[k] < hi`) then fails for **every peak of every drift-less
spectrum** against any real 1/K0 band (0.6–1.4) — the whole spectrum contributes zero peaks to
every window, silently, and the run reports success. `decodeRaw()` computes this *per spectrum*
and keeps such spectra unbanded, so the two paths disagreed on the same file. Fixed with a
per-spectrum `Meta::per_peak_drift` bit; regression test in `odia_spectrumstore_test.cpp` §8.

**The swallowing catch.** `populate()` caught decode exceptions with the comment *"reported once
by decodeRaw's path"* — but once `populated_` is set, `decode()` **never reaches** `decodeRaw`,
so that message could never fire. A throwing spectrum was stored **empty** with zero diagnostic,
ever. The known trigger is upstream and real: `Spectrum::intensity()` throws `TypeError` on
int32-stored intensities, which mzML permits and real files use. A total failure would be
obvious (0 peptides); a *partial* one — one row group, one int32 chunk — would not be. Now
counted and always reported, including the zero case, because "0 undecodable" is the evidence
that the run saw every spectrum.

Also fixed: the populated path pushed a hard `-1.0` where `decodeRaw` pushes the spectrum's
scalar `ion_mobility()`, so IM scoring input differed between paths on per-slice-written files;
and window edges are now rounded **outward** (`floorFloat`/`ceilFloat`) rather than to nearest,
which can otherwise move an edge *inward* by half a float ulp and drop a peak that is genuinely
inside the window.

---

## 3. Compacting the feature tables

### The number, and where it comes from

Releasing pass-1's `FeatureMap` drops `mallinfo2` `in_use` from 29.64 GB to 9.09 GB. The
features are **20.55 GB of live data**. (RSS falls 38.14 GB across the same release, but the
difference is fragmented pages `malloc_trim` returns at that moment — an earlier note in this
repo quoted the 35.5 GB RSS figure as the feature size and **that was wrong**.)

The same content written out is **0.95 GB of parquet**. A 21× representation penalty, and not
mysterious: per precursor, 4.89 `Feature` objects at 296 B, 71.3 subordinates at 296 B, 269 meta
values in string-keyed flat_maps, and **76 separate `MetaInfo` heap allocations**.
`MetaInfoInterface` is an 8-byte pointer to a heap `MetaInfo` holding a `flat_map<UInt,
DataValue>`, so the meta values are invisible to `sizeof()` and are the larger term.

### Three measured facts that make the compact form nearly free

Taken from the emitted parquet of a full run, not from the schema:

1. **35 of 44 `feature_transition` columns are entirely NULL.** The declared schema is aspiration.
2. **`run_id` is constant** — one int64, stored 2.07M times.
3. **The join keys are implied by position.** Subordinates are contiguous per feature, so CSR
   offsets recover the parent exactly; the key column is redundant.

`src/odia_scored.h` acts on all three: `ColumnSet` gives a *declared but absent* column no
storage slot (`kAbsent`), constants are hoisted, and `TransitionRows`/`PrecursorRows` are CSR
with no parent key stored.

### Measured, at the real table shapes (`odia_scored_test.cpp`)

```
features            412.0 B/row  (naive 8 B x 65 = 520)
feature_transition   44.7 B/row  (naive 8 B x 44 = 352)   21.9M rows
feature_precursor    15.7 B/row  (naive 8 B x  5 =  40)
-> full run projection: 1.72 GB against 20.55 GB live as OpenMS Features (11.9x)
```

**Losslessly** — every value still `double`, no precision change. The test deliberately covers
the two failures that would be *silent*: an off-by-one in the CSR offsets (every subordinate
attributed to its neighbour — all rows present, all totals matching), and a feature with **no**
subordinates swallowing the next one's.

### Why precision narrowing is deferred

A float32 store with two orders of margin *still* came back 22 peptides short (6,002 vs 6,024)
on a full run, and **that is still unexplained**. Narrowing 41 more columns on the same
reasoning, in the same change as a representation swap, would make an ID regression impossible
to attribute. So: Phase A flat-but-`double` (~12.3 GB saved, acceptance = **bit-identical IDs**,
which is a real test because the path is deterministic — two runs shared 6,089 peptides, 0
churn); Phase B the ODIA-owned parquet writer; Phase C float32, measured alone.

---

## 4. Open, and what I would not trust

- **The 22-peptide delta on the mzPeak path is unexplained.** float32 m/z has no mechanism to
  *lose* peaks — it can only flip inclusion within ~0.03 ppm of a window edge. The two
  mechanisms that *did* delete peaks wholesale were the drift back-fill and the swallowing
  catch, both fixed above, both silent, both populated-path-only. **Re-measure the delta before
  attributing it to precision.** The run was Astral (no IM), so the drift bug cannot be it — the
  swallowing catch is the prime suspect, and it now reports itself.
- **`odia_scored.h` is built and tested but NOT wired in.** So is the parallel
  `materializeFromCompact_` top-k/keep capability — one call site, still the 2-arg form.
- **SRM chromatogram storage is a live conflict**: the mzPeak writer refuses *product-bearing*
  SRM chromatograms, and our 4.46M transitions are exactly that. Needs a decision, not a fix.
- **libzip is pinned 1.11.4 → 1.11.2 for local measurement builds only.** Not to ship.
- `restoreRealIds_` uses a positional mapping and `setOriginalId` is never called.

Two things to report upstream: a needless deep copy of every `Feature` inside a critical
section, and native IDs used as join keys where indices would do.

## 5. Process notes worth keeping

- **`pgrep -f "<cli> <args>"` matches the waiting shell itself.** Five occurrences; one cost
  ~100 minutes and produced a false "stuck process" diagnosis. Wait on a captured **PID**.
- **A stale checkout invalidates everything measured against it.** `/scratch/.../mzpeak-trunk`
  was an Aug-3 *snapshot*, not a git checkout, with zero cache code. Verify the tree is what you
  think before quoting a factor.
- **A guard can have a hole shaped like the bug.** An OOB read returns zeros — finite and
  non-negative — so a `md > 0.0` precondition skipped the ratio test entirely. Test guards by
  *reintroducing* the defect.
- Adversarial review across models is worth it, but **Codex refuses briefs phrased as
  exploit-hunting** ("what input makes it read out of bounds") with a cybersecurity flag after
  burning ~100k tokens. Rephrase as correctness and robustness; same findings, no refusal.
