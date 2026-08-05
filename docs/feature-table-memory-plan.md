# Replacing OpenMS `Feature` with a flat feature table

## The number, and where it comes from

Releasing pass-1's `FeatureMap` drops `in_use` from 29.64 GB to 9.09 GB, so the features are
**20.55 GB of live data**. (RSS falls 38.14 GB across the same release, but the difference is
fragmented pages that `malloc_trim` returns at the same moment — an earlier note in this repo
quoted the 35.5 GB RSS figure as the feature size and that was wrong.)

The same content, written out, is **0.95 GB of parquet across three tables**:

| table | rows | cols | on disk |
|---|---|---|---|
| `features.parquet` | 2,070,355 | 65 (62 double, 3 int64) | 0.53 GB |
| `feature_transition.parquet` | 21,896,391 | 44 (41 double, 3 int64) | 0.40 GB |
| `feature_precursor.parquet` | 1,226,331 | 5 | 0.02 GB |

A 21× representation penalty. It is not mysterious — `odia_features.h` already records the
component walk: per precursor, 4.89 `Feature` objects at 296 B, 71.3 subordinates at 296 B,
269 meta values in string-keyed flat_maps, and **76 separate `MetaInfo` heap allocations**.
`MetaInfoInterface` is an 8-byte pointer to a heap `MetaInfo` holding a `flat_map<UInt,
DataValue>`, so the meta values are invisible to `sizeof()` and are the larger term.

## Target

Flat columnar, one contiguous block per column, computed from the real schemas:

| | double throughout | float32 score columns |
|---|---|---|
| features | 1.00 GB | 0.52 GB |
| feature_transition | 7.18 GB | 3.83 GB |
| feature_precursor | 0.04 GB | 0.04 GB |
| **total** | **8.22 GB** | **4.4 GB** |

Against 20.55 GB today: **2.5× / 4.7×**. The second column also removes ~157M allocations
(2.07M features × 76), which matters beyond its own bytes — the allocator's retained pool is
the single largest item at peak, and it is fed by exactly this kind of small short-lived churn.

## Phasing, and why it is in this order

**The precision question is deferred on purpose.** A float32 spectrum store — 0.059 ppm worst
error against a 10 ppm window, i.e. two orders of margin — nevertheless came back 22 peptides
short (6,002 vs 6,024) on a full run, and that is *still unexplained*. Narrowing 41 more
columns to float32 on the same reasoning, in the same change as a representation swap, would
make an ID regression impossible to attribute. So:

### Phase A — flat, still `double` (saves ~12.3 GB, zero precision change)

1. `FeatureTable` (exists, unit-tested, never wired) covers `features.parquet`'s shape:
   ids + RT/width/IM + one contiguous score block, column names stored once.
2. Add `TransitionTable` for the 21.9M subordinate rows — the larger table, and the one with
   no existing counterpart.
3. Convert at the point features are produced and drop the OpenMS objects immediately.
   **Pass 1 is the place to start**: its `FeatureMap` is never written to parquet (it feeds
   anchor selection only, `for_anchors=true`) and it holds the global peak. Pass 1 alone is
   the whole 20.55 GB at the moment RSS is highest.

Because no value changes, IDs must be **bit-identical**. That is the acceptance test, and it
is a real one — the ordinary path is deterministic (two runs shared 6,089 peptides, 0 churn),
so any difference at all is a defect, not noise.

### Phase B — the ODIA-owned parquet writer (unblocks pass 2)

Pass 2's features *are* written, so dropping the OpenMS objects there needs a writer.
`OpenSwathOSWParquetWriter` is 1,337 lines and takes `const FeatureMap&`; OpenMS may not be
modified, so this is an ODIA-side build: 65 + 44 + 5 columns through Arrow builders, plus the
`clearSignBit()` id convention (`OpenSwathOSWParquetWriter.cpp:603`) that the join depends on
— we already shipped a bug where our writer used a different one and the join silently
produced zero matches.

Acceptance: the emitted bundle must be **byte-identical** to the current writer's on the same
run, or differ only in ways enumerated and justified.

### Phase C — narrow the score columns (saves a further ~3.8 GB), measured not assumed

Only after A and B are in and reproducing IDs exactly. Change precision alone, nothing else,
and compare identifications directly. If it costs peptides, keep double for the affected
columns — 8.22 GB was already the bulk of the win.

## What this does not fix

The peak is currently `64.23 GB` during `setup/mass_calibration`, with `extract_pass1_wide`
reaching 63.97 GB. Removing 12–16 GB of features leaves the extraction working set, the
compact library (1.54 GB after prefilter, 18.13 GB at load) and glibc's retained pool. This
is the largest single item, not the last one.
