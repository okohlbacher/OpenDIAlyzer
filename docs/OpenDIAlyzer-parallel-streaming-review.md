# OpenDIAlyzer: bottleneck review and a fully-parallel mzPeak-streaming design

Status: 2026-07-28. Everything in §1 is **measured** on spock (224 cores / 2.2 TB RAM), Bruker
diaPASEF `.d`, 1/8 library (9.8M transitions, 24 SWATH windows, ~1650 s gradient) unless marked
otherwise. Companion to `OpenDIAlyzer-mzpeak-backlog.md` (task #21) and #22.

## 1. What actually costs time and memory (measured, not assumed)

| Phase | Parallelism | Measured | Note |
|---|---|---|---|
| Library parse (TSV) + `.d` → RAM | **1 thread** | ~3 min, RSS 52 → 83 GB | `loadDIARun_` uses `SwathFile::loadBrukerTdf(readOptions="normal")` = full materialisation |
| Extraction | **≤ 24** (`omp parallel for` over `wave.swath_indices`) | est. 33.66 GiB per swath | capped by the *number of SWATH windows*, not cores |
| Scoring | up to **224** (work-stealing `inner_batch_queue`) | R = 90–142 threads observed | but only *after* the whole wave has been extracted |
| OSW write | **1 writer thread** (`OSWBufferedWriter`, 2 GiB queue) | `write_bytes` only 49 MB during a stall | not the bottleneck we assumed |
| `PRECURSOR_ID` remap (our post-step) | 1 thread | **> 80 min** on a 48 GB / 4.15M-feature `.osw` | correlated `UPDATE`; see §5 |

**Whole-run observation:** 225 threads, load average 225, but only **~2800–3600 %CPU ≈ 28–36
cores of useful work**, peak RSS **500 GB – 2.0 TB**.

### 1.1 The real #1 bottleneck was *not* the ≤24 loop — it was the allocator
`/proc` sampling during a stalled run: **208 of 225 threads in state `D`** (uninterruptible),
all with `wchan` = `vm_mmap_pgoff` / `__vm_munmap`. That is contention on the kernel's
**per-process `mmap_lock`**: glibc `malloc` routes large allocations through `mmap`/`munmap`,
and 224 scoring threads allocating chromatogram buffers serialise on that one lock.
`write_bytes` was 49 MB, so it was definitively not I/O.

**Fix (no code change):** `LD_PRELOAD=libtcmalloc_minimal.so` + `TCMALLOC_RELEASE_RATE=0`.
Result: D-state → **0**, scoring ran R = 90–142 threads, CPU peaked **~7200 % (~72 cores)**, and
a run that had never finished completed in **2:27:32**.

### 1.2 Window width is expected to matter — but my headline number was WRONG (retracted)
Total extraction work ∝ (transitions × spectra inside each RT window), so narrowing the window
must reduce work. **However**, the "full-library OpenSWATH arm finished in 10 minutes" figure I
first cited here is **invalid**: that run *aborted*. Its log ends with

    Error: Unexpected internal error (WARNING: rsq: 0.92644790211059 is below limit of 0.95.)

i.e. OpenSWATH's own calibration QC gate (`qc:min_rsq` default 0.95) rejected its fit on this
**predicted** library, and the run produced **0 features**. It was fast because it did nothing.
**Window width therefore remains an untested hypothesis on full-scale data** — it is measured
only on `lib_small2` (where 3600 s → 407 s changed peak agreement, not a clean timing A/B).
Until the controlled experiment in §2 is run, treat §1.3's ranking of window width as
*unverified*.

**Separate finding of its own importance:** stock OpenSwathWorkflow, at default settings,
**cannot process this predicted library at all** — the rsq gate fails (0.926 < 0.95) and it
errors out rather than degrading. Our engine survives because `-rt_calib_min_rsq` is relaxed to
0.7 and it falls back to the bootstrap on failure. Any "OpenSWATH baseline" number must
therefore state the `qc:min_rsq` used, or it is not reproducible.

### 1.3 Revised ranking (wall time)
1. **Allocator contention** — fixed, free, ~10× on the affected runs.
2. **RT window width** — *hypothesis*, not yet measured at scale (see §1.2 retraction).
3. **Serial library parse** — ~3 min fixed cost; matters only for short runs (see §5.3).
4. **≤24-way extraction cap** — *now* the limiter, but it is third-order until 1–2 are in place.
5. Single-threaded OSW writer — only bites at very high feature rates.

**Falsification tests.** If (4) were dominant, useful CPU during extraction would sit at ~24
cores *with D-state ≈ 0* — that is now testable in one run. If (1) were incomplete, D-state
would return under tcmalloc (it did not).

### 1.4 Memory
RSS is dominated by the **whole run resident in RAM** (`readOptions "normal"`), not by the
library: 24 swaths × ~33.7 GiB estimated ≈ 800 GB, observed 500 GB – 2.0 TB. This is the single
strongest argument for streaming: it is a *memory* problem first and a speed problem second.

## 2. Honest caveat: is the mmap story complete?

Alternative explanations I have **not** excluded, and how to settle them cheaply:
- **THP / page-zeroing / first-touch NUMA.** 2 TB of first-touch on a multi-socket box can look
  like kernel-time stalls. Distinguish with `perf top -g` during extraction, `/proc/vmstat`
  (`thp_fault_alloc`, `compact_stall`, `pgmajfault`), and `numactl --hardware` + per-node
  `numastat`. If THP compaction dominates, `transparent_hugepage=madvise` may matter more than
  the allocator.
- **tcmalloc win ≠ proof of `mmap_lock`.** tcmalloc also avoids `munmap` entirely
  (`TCMALLOC_RELEASE_RATE=0`), so the win is consistent with either lock contention *or* pure
  syscall/TLB-shootdown cost. A cheap discriminator: run stock glibc with
  `MALLOC_MMAP_THRESHOLD_=1073741824` and `M_TRIM_THRESHOLD=-1` — if that alone recovers most of
  the win, it was the mmap path; if not, it is arena/lock behaviour.
- **Confound in the 10-minute OpenSWATH run:** that arm also used a different RT window *and* a
  narrower candidate set. Attributing it purely to window width needs one controlled run
  (same tool, same everything, window 3600 vs 407).

## 3. Target architecture: fully parallel + streaming on mzPeak

### 3.1 Unit of parallelism: the (SWATH × RT-tile) tile
The ≤24 cap exists because the loop is over SWATH windows. Break it by adding a **second
dimension**: split the gradient into RT tiles and make the work item a **tile = (swath_index,
rt_lo, rt_hi)**. With 24 swaths × 32 tiles ≈ **768 independent work items** — enough to saturate
224 cores with dynamic scheduling and good load balance.

```
work item  = { swath, rt_lo, rt_hi, assay_slice }
per item   : stream spectra in [rt_lo-Δ, rt_hi+Δ] from mzPeak
             → extract chromatograms for assays whose window intersects the tile
             → pick peaks → score → emit rows
```

Because a tile owns a bounded RT range, its **working set is bounded**: one tile's spectra plus
its assays' chromatograms, not the whole run. That is the memory fix.

### 3.2 The mzPeak adapter, and why the current interface is the hard part
`OpenSwathWorkflow` consumes `std::vector<OpenSwath::SwathMap>`, each holding a
`shared_ptr<ISpectrumAccess>`. Implementing `MzPeakSpectrumAccess : ISpectrumAccess` is
mechanically easy, but `ISpectrumAccess` is a **random-access pull** interface
(`getSpectrumById(int)`), which is a poor match for a columnar streaming reader:

- `lightClone()` is called **per scoring job** (`tmp.back().sptr = context.current_swath_map->lightClone()`)
  — the adapter must be cheap to clone and **thread-safe for concurrent reads**.
- Random `getSpectrumById` over a compressed columnar store means either (a) decode-on-demand
  with a cache, or (b) sequential decode of a whole tile up front. For tiles, **(b) is right**:
  decode the tile once into a compact in-memory buffer, then serve `getSpectrumById` from it.

**Conclusion:** implement the adapter as a **tile-backed materialiser**, not a naive
lazy-per-spectrum shim. `MzPeakTileAccess` = "the spectra of swath *s* in RT [lo,hi], decoded
once, served random-access, discarded when the tile retires." This keeps the OpenMS extraction
code unchanged while making memory O(tile) instead of O(run).

### 3.3 Pipelining
Extraction and scoring are currently **sequential per wave**. With tiles they become a two-stage
pipeline: a bounded queue of decoded tiles (producers = readers, N≈4–8) feeding scoring
consumers (all remaining cores). Back-pressure = queue depth in **bytes**, not items, so memory
stays bounded regardless of tile density.

### 3.4 Output — in priority order (revised after review)
1. **Fix the integer id at row construction.** The writer is handed the *string* assay id as
   `PRECURSOR_ID` (`OpenSwathOSWWriter.cpp:1112`); carry the seeded PQP integer through the
   library index instead. This deletes our 80-minute post-hoc `UPDATE` outright.
2. **Measure the writer at W = 407 s before replacing it.** A single writer thread is irrelevant
   if it is 5 % of wall time — and the earlier "writer stall" turned out to be the allocator.
3. **Then** prefer incremental **Arrow/Parquet shards**. OpenMS already defines an OSW-Parquet
   schema (`OpenSwathOSWParquetWriter.h:23`), though its writer currently expects a resident
   `FeatureMap`, so it needs an incremental RecordBatch interface.
4. If SQLite output is mandatory, **bulk-load the sorted Parquet into a fresh `.osw` once** —
   never `UPDATE` a 48 GB database in place.

Per-thread SQLite is a *second* choice: one writer per database, duplicated static library
tables, id/foreign-key reconciliation, and a serialised final merge. If used at all, shard by a
**stable domain** (SWATH / RT partition), never by incidental thread number.

### 3.5 What must change *above* the adapter (codex)
A tile-backed `ISpectrumAccess` is necessary but **not sufficient** — the workflow above it also
has to change, otherwise `performExtraction`'s own OpenMP team still owns every thread:
- route the library by **SWATH × transformed expected RT** and hand a tile only its assay slice;
- stop wrapping the source in `SpectrumAccessOpenMSInMemory` (that is the full materialisation);
- remove the whole-wave extraction barrier;
- expose extraction / peak-picking / scoring as **task-level** operations instead of one global
  `performExtraction` call.

### 3.6 Memory arithmetic (sanity check)
With `W = 407 s` and a 1.385 s cycle a chromatogram is ~294 points. A batch of 1024 assays ×
(10 fragments + 4 MS1 traces) × 294 points × 4 B ≈ **16.9 MB** of intensity. So even a few
hundred MB per active task is fine — *provided* decoded slabs and task counts are governed by
**byte credits**. It is emphatically not fine if 224 workers each independently materialise a
haloed tile.

**Caveat that mzPeak does not solve:** the *library* is its own memory floor (a full
`LightTargetedExperiment` expansion has been measured in the tens of GB). Streaming fixes raw
input only; the library needs a compact indexed source (PQP/Parquet or a read-only SoA).

## 4. Correctness risks of tiling (the part that will actually bite)

> **Corrected after adversarial review (codex).** My first proposal — halo + "emit only peaks
> whose apex ∈ core" — is **not sufficient**. Suppressing out-of-core apexes changes the *top-N
> competition* for an assay, and smoothing/S-N would still be computed over a different span.
> The invariant design is **assay-ownership**, not peak-ownership:
>
> 1. assign each assay **exactly once**, to the half-open slab containing its transformed
>    expected RT;
> 2. decode enough halo to cover that assay's **complete original extraction interval**
>    (for `W = 407 s`, an assay at a boundary needs ~**203.5 s on each side** — i.e. Δ = W/2,
>    *not* "max peak width");
> 3. build its chromatograms over **exactly the same RT samples** the non-tiled run would use;
> 4. pick + score **once** on that complete chromatogram;
> 5. emit the normal top-N peak groups **even if an apex falls outside the core**;
> 6. clip only at the physical start/end of the run, exactly as the baseline does.
>
> Under this scheme tiling is a pure execution detail and results are bit-identical by
> construction, because every assay sees the identical chromatogram it would have seen.

| Risk | Why it breaks | Mitigation |
|---|---|---|
| Peak straddling a tile edge | chromatogram truncated → wrong apex/area, or a peak found twice | **Assay-ownership + Δ = W/2 halo** (see box). Never suppress by apex position. |
| **S/N window is 1000 s by default** (`PeakPickerChromatogram.cpp:33`) | that is *wider than a plausible tile* — a tile fragment silently changes S/N and therefore peak selection | The halo must cover the assay's full extraction interval; S/N is then initialised from the same chromatogram as the baseline. **This alone rules out tiles smaller than the extraction window.** |
| Smoothing / peak extension | picker copies and smooths the **whole** chromatogram before picking (`PeakPickerChromatogram.cpp:92`); SG frame is 11 points but extension can look farther | Same fix: reproduce the identical chromatogram, do not feed fragments. |
| Full-scan (MS1/DIA) scores | scoring fetches raw spectra around the chosen apex via `getMultipleSpectra` | The tile must retain the **apex-neighbouring cycles**, or copy them into the scoring batch before the slab retires. |
| Calibration drift per tile | a per-tile calibration would make results tiling-dependent | Calibration stays **global** (§CiRT, one fit for the run). Tiles are an execution detail only. |
| Determinism across thread counts | floating-point reduction order, dynamic scheduling | Tiles are **disjoint and independently reduced**; final ordering by (precursor_id, rt) before writing. No cross-tile accumulation ⇒ results independent of thread count. Add a CI gate: same input, `-threads 1` vs `-threads 224`, byte-identical `.osw` after canonical sort. |
| Assay assigned to the wrong tile | assays whose RT window spans tiles | Assign an assay to the tile containing its **expected RT**; the halo guarantees its full window is available. |

**Invariance gate (must-have):** a test that runs the same file with 1, 4 and 32 tiles and
asserts identical features. Without this, tiling silently changes results and every future
benchmark becomes untrustworthy.

## 5. Cheaper wins that do not need the rewrite

1. **Keep tcmalloc + `/dev/shm` + `OMP_WAIT_POLICY=PASSIVE`** — already validated; bake into
   the harness and document as the supported way to run.
2. **Calibration-estimated RT window** — already implemented (`-use_estimated_rt_window`);
   the single largest measured lever.
3. **Kill the `PRECURSOR_ID` post-`UPDATE` (> 80 min).** Assign the integer id at write time
   from the library index instead of a correlated SQL update. Pure win, no risk.
4. **Parse the library once, cache as PQP.** The 2.3 GB TSV parse is single-threaded; a PQP
   (SQLite) library loads far faster and is content-addressable for reuse across runs.
5. **Sharded processes as a stopgap.** Running *k* processes each owning a subset of SWATHs,
   then merging `.osw`, sidesteps both the ≤24 cap and the shared `mmap_lock` — at the cost of
   *k*× the run in RAM, which is exactly what streaming fixes. Only attractive on a big-RAM box
   and only until tiling lands.

## 6. Recommended sequence (with gates)

| # | Step | Effort | Expected | Gate |
|---|---|---|---|---|
| 1 | Bake allocator/tmpfs/window settings into the harness | trivial | already-measured 10× / 9× | full run completes < 30 min |
| 2 | Drop the post-hoc `UPDATE`: set the integer id at row construction (`OpenSwathOSWWriter.cpp:1112`) | S | −80 min on large runs | `.osw` identical, pyprophet-scorable |
| 3 | Profile to settle §2 (perf + vmstat + one 3600-vs-407 controlled run) | S | knowledge | attribute the 10-min run correctly |
| 4 | `MzPeakTileAccess` adapter + **1 tile per swath** (no tiling yet) | M | memory O(run)→O(swath); unblocks #21 | identical features vs current loader |
| 5 | RT tiling, **assay-ownership + Δ=W/2 halo** (never emit-by-apex) | M/L | ≤24 → ~768 work items; memory O(tile) | **invariance gate** (1 vs 32 tiles identical) |
| 6 | Two-stage pipeline (decode ∥ score) with byte-bounded back-pressure | M | hides decode latency | no RSS growth with tile count |
| 7 | Sharded Arrow/Parquet output + single bulk-load | M | removes writer serialisation | byte-identical after canonical sort; **only if step 3 shows the writer matters** |
| 8 | GPU offload of the extraction inner loop | L | speculative | only after 1–7; profile first |

**The one cheapest decisive experiment:** step 3 — a single controlled pair of runs (same tool,
same library, window 3600 vs 407) with `perf` and `/proc/vmstat` captured. It settles whether
the remaining wall time is window-driven work, allocator/page behaviour, or the ≤24 cap, and
therefore whether steps 5–7 are worth their risk at all.

## 7. Bottom line

The ≤24-way extraction loop was the *assumed* bottleneck all along; measurement says it was
third in line behind **allocator contention** and **window width**, both of which are now fixed
for free. The strongest remaining argument for the mzPeak streaming rewrite is **memory**
(500 GB – 2 TB resident → bounded per-tile), with parallelism a welcome second. Do steps 1–4
before committing to 5–7, and never land tiling without the invariance gate.

---

## 8. Adversarial review (codex) — what it changed

Reviewed independently against the tree. It **agreed** on the ranking (allocator + window width
ahead of the ≤24 loop; memory as the real case for streaming) and on the tile-backed adapter,
but corrected three things that would have caused silent wrongness:

1. **My halo design was unsound.** "Emit only peaks whose apex is in the core" changes the
   top-N competition and still computes smoothing/S-N over a different span. Replaced with
   **assay-ownership + Δ = W/2 halo + emit all top-N** (§4 box). This is the difference between
   "probably fine" and "identical by construction".
2. **The S/N window default is 1000 s** (`PeakPickerChromatogram.cpp:33`) — wider than any
   sensible tile. That single fact rules out tiles smaller than the extraction window and would
   have silently changed peak selection had I not been told.
3. **Output ordering.** Fix the integer id at row construction *first* (deletes the 80-min
   `UPDATE` for free), **measure** the writer before replacing it, and use OpenMS's existing
   OSW-Parquet schema rather than per-thread SQLite.

It also supplied the memory arithmetic (~16.9 MB of intensity per 1024-assay batch — so tasks
are cheap if governed by byte credits) and flagged that **mzPeak does not fix the library
memory floor**; that needs a compact indexed library source separately.

**Net effect on the plan:** the sequence in §6 is unchanged, but step 5 is now materially
different (assay-ownership, not apex-filtering) and step 7 is explicitly gated on measurement.
