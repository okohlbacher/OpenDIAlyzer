# Where OpenDIAlyzer's memory goes — an adversarial review

Measured on the Astral plasma benchmark (`astral.mzML`, 6.4 GB; 7,149,966-precursor library) on the
IBMI `data` node (224 cores, 2.2 TB). Every number below is measured unless labelled otherwise.

## 1. The measurements

Object sizes, from a probe that included the real OpenMS headers (`/scratch/kohlbach/tmp/memprobe`):

| type | bytes |
|---|---:|
| `LightTransition` | 128 |
| `LightCompound` | 376 |
| `LightProtein` | 96 |
| `ChromatogramPeak` | 16 |
| **`MSChromatogram`** | **960** (header alone, before any point) |

Full library resident, loaded standalone from the parquet bundle:

```
after load   RSS = 38.84 GB   (load peak 45.89 GB)
compounds = 7,149,966   transitions = 78,569,077
  vector-only:      compounds 2.50 GB + transitions 9.37 GB = 11.87 GB
  string CHARACTERS: transition_name 2.18 GB, peptide_ref 2.00 GB, id 0.18 GB, sequence 0.15 GB
  heap-allocated (>15 chars, SSO exceeded): names 72,177,501/78,569,077  refs 67,486,683/78,569,077
```

**The library holds 16.4 GB of information in 38.8 GB of RSS — 2.4x.** The excess is ~140 million
individual small string allocations: two strings per transition, both usually past the 15-character
small-string threshold, each therefore a separate heap block with its own allocator header and
16-byte rounding.

Whole-run peaks:

| run | precursors | transitions | RT window | wall | peak RSS |
|---|---:|---:|---:|---:|---:|
| mz5 | 95,676 | 988,052 | 600 s | — | 74.2 GB |
| mz10 | 423,079 | 4,463,919 | 600 s | 30:14 | 189.5 GB |
| e2e (all fixes) | 423,079 | 4,463,919 | 1435 s → 240 s | 40:37 | **366.8 GB** |

## 2. A fitted model, and where it fails

mz5 and mz10 are a controlled pair — same binary, same threads, same RT window, differing only in
scale. Fitting peak against transition count:

```
peak ≈ 41.5 GB  +  34.8 KB × n_transitions          (at a 600 s window)
```

Decomposing that by what each term should scale with:

```
peak ≈ 41.5 GB [library]  +  126 KB × n_precursors  +  n_chrom × (960 B + 16 B × window/cycle)
```

The 41.5 GB intercept is independently corroborated: the standalone library load measures 38.8 GB.
That is not a fitted coincidence — it is the same object.

**Honest limitation.** Two points, two parameters, so the fit is exact by construction and proves
nothing on its own. Its out-of-sample prediction for the e2e run was **317 GB against 366.8 GB
observed — 16% low**. The model captures the shape but not the peak, and the missing ~50 GB is
unattributed. Candidates not separated: the parquet load transient (measured at +7 GB standalone),
feature accumulation, the OSW write buffer (up to 2 GB), and the pass-1 `.osw` being written while
pass 2 runs. Do not quote the model as validated.

## 3. Findings

### 3.1 `estimateFeatureMemoryPerCompound()` is 63x low, and budgets the wrong thing

`OpenSwathWorkflow.cpp:46`:

```cpp
Size estimateFeatureMemoryPerCompound()
{
  // Conservative estimate in bytes per compound (features + temporaries)
  return static_cast<Size>(2 * 1024); // 2 KB
}
```

Measured: **126 KB per compound** for the non-chromatogram term. And that term is the *small* one —
chromatograms are 225 GB of the 367 GB at a 1435 s window. So nothing in the batching logic is sized
against the quantity that actually sets peak RSS.

### 3.2 The adaptive batch sizing is inert on a large machine

`calculateInnerBatchSize()` takes 5% of *free system memory* and divides by that 2 KB. On a 2.2 TB
node that is 100 GB / 2 KB = 50 million, so it always saturates the hard clamp:

```cpp
inner = std::max<Size>(2000, std::min<Size>(10000, inner));
```

The batch is therefore always 10,000 compounds, chosen by a constant rather than by the memory
situation. The "5% of available memory" logic never binds — on a small machine it would, but there
the 2 KB unit cost makes it wrong in the dangerous direction.

### 3.3 `batchSize` defaults to 0 — no outer batching at all

`-batchSize 0` means the whole prefiltered library is carried through extraction. This is the lever
that directly bounds the chromatogram term and it is off by default.

### 3.4 The RT window is the single largest knob

Chromatogram memory is linear in the window. mz10 and e2e searched **identical** precursor and
transition sets; the only material difference was 600 s vs a calibration-estimated 1435 s, and peak
went 189.5 → 366.8 GB. The calibration's own log says the estimate "understates the tail", i.e. it
is deliberately generous. 1435 s on a 2333 s gradient is 61% of the entire run.

### 3.5 MS1 isotopes are 27% of all chromatograms

`ms1_isotopes` defaults to 3, i.e. **4 MS1 chromatograms per precursor** — 1,692,316 of the 6,156,235
total, 61.8 GB at a 1435 s window.

### 3.6 The library's 38.8 GB is dead weight for the whole run

After `prefilterLibrary_` the working set is 423,079 of 7,149,966 precursors (5.9%). The full library
*is* correctly freed (`transition_exp = std::move(out)`), but it is ~140M small allocations, and
glibc keeps blocks that small in its arena free-lists rather than returning the pages.

Measured (`/scratch/kohlbach/tmp/trimprobe`, real library, prefilter emulated at the same 5.9%):

```
library loaded (7,149,966 cmp)  RSS = 38.79 GB
after building filtered set     RSS = 39.54 GB
after move-assign (freed 94%)   RSS = 27.67 GB   <-- what ODIA did
after malloc_trim(0)            RSS = 17.80 GB
after full scope exit + trim    RSS =  7.73 GB
```

Two results, one of which corrects an earlier guess in this document's own drafting:

- **`malloc_trim(0)` recovers 9.9 GB**, not the ~39 GB first speculated. glibc had already returned
  ~11 GB unprompted; trim gets the next 9.9. Applied at `opendialyzer.cpp:2275`.
- **~15 GB survives the trim** against ~2.3 GB of live data. That is fragmentation, not free-list
  retention: the surviving 6% of objects sit interleaved with the freed 94% in the same pages, and a
  page holding one live object cannot be returned. No amount of trimming reaches it. Fixing that
  needs the library not to be 140M separate allocations — string interning, an arena, or offsets
  into one character pool. Not attempted.

Note also that 7.73 GB remains after *everything* is freed and trimmed, which bounds how much of
this is recoverable by allocator calls at all.

## 4. The fix built: `src/odia_chromatogram.h`

Replaces `vector<MSChromatogram>` with a store exploiting two facts:

**The retention-time axis is shared, not per-chromatogram.** `ChromatogramExtractorAlgorithm.cpp:316`
loops over *spectra* on the outside and transitions on the inside, pushing that spectrum's
`s_meta.RT` into every chromatogram it touches and skipping only those whose own `[rt_start, rt_end]`
excludes it. So every chromatogram from one SWATH map lies on that map's spectrum-RT grid and
occupies a **contiguous slice** of it. Store the grid once, plus a `(first, count)` pair per
chromatogram, and the RT half of every point disappears. This is an invariant of the extractor, not
an approximation.

**Intensities as float32.** The mzML arrays these come from are routinely 32-bit to begin with, and
every downstream score is relative.

| | header | per point | per chromatogram @2392 pts |
|---|---:|---:|---:|
| `MSChromatogram` | 960 B | 16 B | 39,232 B |
| `ChromatogramStore` | 24 B | 4 B | 9,592 B |

Measured by the self-check (`odia-chromatogram`, ctest):

```
20000 chromatograms x 2392 points:
  MSChromatogram 0.73 GB   compact 0.18 GB   ratio 4.09x
  full run (6,156,235 chrom @ 1435 s window): 225 GB -> 55 GB, saves 170 GB
marginal bytes per point: 4.000
```

The float32 narrowing is verified not to move the scores that consume it, rather than assumed to be
harmless:

```
worst relative intensity error: 5.865e-08   (float32 eps = 1.192e-07)
XIC correlation  0.997954835304 (double) vs 0.997954835027 (float32)  |delta| = 2.775e-10
XIC dot product  relative delta = 4.230e-09
```

Out-of-range spans throw rather than clamp — a clamped span pairs intensities with the wrong
retention times, which is a wrong answer rather than a smaller one.

**Scope, stated plainly.** This is the container and its adapter, not a rewrite of the extraction
path. `MRMFeatureFinderScoring` still consumes `MSChromatogram`. The 170 GB is realised only when the
store owns the run's data and `MSChromatogram` views are materialised per scoring call or per small
batch; materialising all of them at once costs what it costs today *plus* the store. Integrating it
means changing `ChromatogramExtractorAlgorithm` to write into the store and giving the scorer a view
type — that work is not done.

## 5. Ranked levers

| # | lever | est. saving at the e2e shape | status |
|---|---|---:|---|
| 1 | `ChromatogramStore` instead of `vector<MSChromatogram>` | **~170 GB** | built + tested; **not integrated** |
| 2 | RT window 1435 s → 600 s | ~128 GB | one parameter; costs sensitivity in the tail |
| 3 | `-batchSize` > 0 to bound the chromatogram term | bounded by choice | available today, untested here |
| 4 | `malloc_trim(0)` after the prefilter | **9.9 GB (measured)** | **done** — `opendialyzer.cpp:2275` |
| 5 | `ms1_isotopes` 3 → 2 | ~15 GB | costs an isotope for MS1 scoring |
| 6 | Fix `estimateFeatureMemoryPerCompound()` 2 KB → measured | none directly | makes 3 sizeable rather than arbitrary |
| 7 | De-fragment the library (intern strings / arena) | ~15 GB | not attempted; see 3.6 |

Levers 1 and 2 are multiplicative, not additive: the store shrinks the per-point cost that the window
multiplies. Together the chromatogram term at 600 s would be ~23 GB instead of 97 GB.

---

## 6. The A/B result — and what it refutes

Controlled A/B of the two chromatogram fixes (`chrom_list.clear()` after conversion + exact
`reserve()` in the extractor). Identical config, identical library, only libOpenMS changed:

| | wall | CPU | peak RSS |
|---|---:|---:|---:|
| before | 40:37 | 2893% | **366.8 GB** |
| after | 44:43 | 2805% | **349.9 GB** |

**16.9 GB, 4.6%.** The model in §2 predicted the chromatogram term was 225 GB of the 367 GB, and
these two changes should have removed a large fraction of it. They did not. **The A/B is evidence
against the model's decomposition**, and the model should not be quoted for attribution.

### Why the fixes overlap — the 156 GB estimate was wrong

`OpenSwathDataAccessHelper::convertToOpenMSChromatogram` **already reserves exactly**
(`DataAccessHelper.cpp:71`: `chromatogram.reserve(cptr->getTimeArray()->data.size())`). So the
`MSChromatogram` copy — the one that survives into scoring — never had doubling slack at all. The
42% slack existed **only** on the `OpenSwath::Chromatogram` side (the two `vector<double>` grown by
`push_back`), and that side is exactly what `chrom_list.clear()` now frees immediately after
conversion. The two fixes therefore address overlapping memory, and the ~156 GB figure derived from
"42% of the whole chromatogram term" was wrong because it assumed both representations doubled.

What each is worth, per chromatogram at 2392 points:

| | before | after |
|---|---:|---:|
| `chrom_list` (2 × `vector<double>`, doubling → 4096) | 65.5 KB | 38.3 KB (reserved) |
| `MSChromatogram` (always exactly reserved) | 39.2 KB | 39.2 KB |
| **held through scoring** | **104.7 KB** | **39.2 KB** |

A 62% cut of the chromatogram term during the scoring phase produced 4.6% of total peak. **That
means chromatograms are a far smaller share of peak RSS than §2 attributed to them**, and the true
dominant term is unidentified. Finding it needs a phase-resolved measurement (RSS sampled against
the progress log, or heaptrack on a reduced run) — not another model.

### The wall-time number is confounded, do not use it

40:37 → 44:43 looks like a 10% regression, but the second run shared the node with a 10-worker
pyprophet job (~220 GB) and an unrelated 360 GB `dorado` process, while the baseline ran on a
quieter machine. CPU% fell (2893 → 2805) which is consistent with contention rather than with the
code. There is a real candidate for a genuine slowdown — the exact-reserve change adds a
metadata pass over every spectrum (`getSpectrumMetaById`) before extraction, and under streaming
access that may not be free — but this measurement cannot separate the two. Re-run on an idle node
before concluding anything.
