# Memory and identifications are the same problem

**2026-08-02.** Where ODIA's memory goes, and why the answer changes the ID plan rather than sitting
beside it.

Measured on `recal600` (tuned windows, 6,930 IDs), Astral plasma benchmark.

---

## 1. The profile

| phase | delta RSS | cumulative |
|---|---:|---:|
| `library_load` | **+32.65 GB** | 32.71 |
| `dia_run_load` | +6.68 | 39.39 |
| `prefilter` (transient peak 63.74) | -4.31 net | 35.08 |
| `setup/cirt_calibration` | +3.90 | 38.98 |
| **`extract_pass1_wide`** | **+44.99** | 83.97 |
| `extract_pass2_narrow` | +5.99 | 90.26 |
| `write_parquet_bundle` | +1.34 | **peak 92.93** |
| `write_scores` | -8.78 | 82.90 |

Two consumers: **the library at 32.65 GB** and **extraction features at ~51 GB**. Everything else is
noise by comparison.

### 1.1 The library number is the interesting one

We load all **7,149,966** precursors (32.65 GB), prefilter to **423,079** (5.9%), and RSS falls only
from 39.39 to 35.08 GB. The freed 94% does not return to the OS -- consistent with the arena
fragmentation measured earlier (`in_use` 7.26 GB against `retained` 5.77 GB right after load).

A `CompactLibrary` representation already exists and was measured at **2.13 GB against 16.70 GB**
for the `LightTargetedExperiment` equivalent, but it is wired only as `-compact_probe`, not as the
load path.

### 1.2 What the extraction memory is made of

From the component walk (2.07M features):

| | count | bytes |
|---|---:|---:|
| top-level `Feature` | 2,070,081 | 0.57 GB |
| subordinates | 30,173,963 | 8.32 GB |
| meta values in `flat_map` | 113,855,060 | ~2.7 GB |
| **separate `MetaInfo` heap allocations** | **32.2 M** | -- |

Accounted: ~11.6 GB against ~51 GB of extraction growth. The remainder is chromatograms plus
allocator overhead and fragmentation -- which is why `LD_PRELOAD` tcmalloc previously cut peak RSS
from 189 to 106 GB (-44%) without touching a line of code.

---

## 2. The connection to identifications

**The prefilter is a memory/speed optimisation, and it costs 923 IDs.**

Its registered purpose is "screen the library against MS1/MS2 evidence in the run before
extraction". It exists so that extraction handles 423,079 precursors instead of 7,149,966. Measured
recall against DIA-NN's ID list: **88.1%** -- it deletes 923 of DIA-NN's 7,787 identifications before
extraction can see them, which is **38.8% of the total ID gap** and unrecoverable by any downstream
work.

So the two problems are one problem:

```
  more IDs  ->  keep more of the library  ->  more memory
  less memory  ->  prefilter harder  ->  fewer IDs
```

Every ID avenue on the backlog pushes memory the wrong way: relaxing the prefilter, widening RT
windows, exporting and re-searching an empirical library. At 93 GB on a 2.2 TB node there is room,
but the trade is real and currently unmanaged.

### 2.1 DIA-NN does not make this trade

DIA-NN has **no prefilter**. It bounds memory by **batching**: the library is split into random
batches of 2,000 precursors with a fixed seed (`diann.cpp:137-138`, `:8999-9014`), and batches are
processed until enough identifications accumulate (`:10371-10392`).

Batching bounds memory to one batch's worth of features *and* preserves recall, because every
precursor is eventually searched. Our prefilter bounds memory by **deleting** precursors, which
bounds recall at the same time.

That is the whole difference in one sentence: **DIA-NN pays for memory with wall-clock; we pay for
it with identifications.**

It also explains the reported footprints -- ~20 GB for DIA-NN against our 93 GB -- without any
appeal to better data structures.

---

## 3. Avenues, ordered by what they buy

### A. Library batching -- fixes memory AND recall together
Replace "prefilter then extract once" with "extract the full library in batches". Bounds the
resident feature set to one batch, removes the 32.65 GB full-library residency, and makes the
prefilter unnecessary -- recovering the 923 precursors it deletes.

- **Memory:** bounded by batch size rather than library size.
- **IDs:** +923 ceiling raised, plus whatever the relaxed candidate set yields.
- **Cost:** wall-clock, and a scoring question -- q-values must be computed over the union, not per
  batch, or the FDR is wrong. DIA-NN accumulates and rescores globally.
- **Risk:** this is the largest structural change on the list.

### B. Use `CompactLibrary` on the load path
Already written and tested; measured 2.13 GB against 16.70 GB. Currently reachable only through
`-compact_probe`.
- **Memory:** plausibly -25 GB of the 32.65 GB library residency.
- **IDs:** none directly, but it buys the headroom that A and the prefilter relaxation need.

### C. tcmalloc by default
Measured: peak 189 -> 106 GB (-44%), identical IDs, +21% wall on the old build.
- **Memory:** large, free, already validated.
- **IDs:** none.
- **Cost:** wall-clock. Worth re-measuring on the current 93 GB build, where the fragmentation
  share may differ.

### D. Feature representation
32.2 M separate `MetaInfo` heap allocations, one per feature and subordinate, each holding a
string-keyed `flat_map`. Replacing string-keyed metadata with a fixed sub-score vector would remove
most of the allocation count.
- **Memory:** unquantified, but the allocation count is the fragmentation driver.
- **IDs:** none directly.
- **Risk:** touches the OpenMS `Feature` interface -- clean-room constrained, so ODIA would need its
  own scored-feature type.

### E. Stream features out during extraction
Currently 2.07 M features are retained to the end. Scoring needs them all only because the
classifier trains globally.
- Subsumed by A if batching lands: a batch's features can be scored and discarded.

---

## 4. Revised priority, folding memory into the ID plan

| # | action | IDs | memory | note |
|---|---|---|---|---|
| 1 | finish the x1.2 window ladder (720, 864) | + | neutral | queued; the 600 s default may be too low |
| 2 | delete the three calibration gates | +122 measured | none | free |
| 3 | prefilter relaxation measurement (`pf_frag3`, `pf_peaks3k`) | up to +923 | **worse** | running; quantifies the trade |
| 4 | `CompactLibrary` on the load path | none | **-25 GB** | buys the headroom for 3 and 6 |
| 5 | `pRT` + sqrt-compressed `pdRT` scoring features | targets 1,261 | neutral | |
| 6 | **library batching** | +923 and removes the ceiling | **bounded** | the structural fix; retires the trade entirely |
| 7 | tcmalloc default | none | -44% | re-measure on the current build |
| 8 | empirical library export + re-search | provenance + prefilter both | +1 pass | |

The ordering is deliberate: 4 before 3 and 6, because relaxing the prefilter or batching without
first cutting library residency runs into the memory wall that motivated the prefilter in the first
place.

---

## 5. What would falsify this framing

- **The prefilter relaxation gains few IDs.** Then the 923 are not recoverable in practice (they may
  be unidentifiable at any sensible score threshold), and the memory/ID coupling is weaker than
  claimed. `pf_frag3` and `pf_peaks3k` test exactly this.
- **`CompactLibrary` on the load path does not reduce peak RSS.** Peak occurs during
  `write_parquet_bundle`, not `library_load`, so cutting library residency may move the floor
  without moving the peak.
- **Batching changes the FDR.** If per-batch score distributions differ enough that a global
  q-value over the union is not comparable to today's, the ID counts are not comparable either and
  the +923 claim needs restating.
