# Thread-pool chromatogram pipeline — design, adversarial review, and why most of it was dropped

A design (D1–D6) was written for a thread-pool chromatogram pipeline and sent to two independent
reviewers (kimi, codex). **The review killed most of it.** This records what was proposed, what the
reviewers found, what my own arithmetic then found that both reviewers missed, and what was actually
built.

## What was proposed

- **D1** flat work list of (window, batch) pairs, sorted longest-first (LPT)
- **D2** reference-counted SWATH-window residency with a semaphore cap
- **D3** two-stage pipeline: extract pool → bounded queue → score pool
- **D4** per-worker pooled `ChromatogramStore`
- **D5** per-worker feature buffers to kill the `osw_write_out` critical section
- **D6** memory sizing from `n_transitions × (24 + 4 × points)`

## What the reviewers found

### D3 is wrong, and D3 conflicts with D4 (kimi)

> "D4 says each worker owns a store reused across units — but in a pipeline, the store must be handed
> to the score worker, so **D4 conflicts with D3**."

A store cannot be both worker-owned-and-reused and passed through a queue. Beyond that, kimi's
argument against the pipeline is that a flat pool gives the *same* explicit memory bound
(`workers × batch_bytes`) with no queue, no handoff and no E:S ratio to tune.

Worse, on Q8: **if the phase is bandwidth-bound, D3 is actively harmful.**

> "the store is written by E, read by S — cross-core transfer, extra LLC misses... with flat pool,
> same core writes and reads → cache-hot. So if bandwidth-bound, the flat pool wins again, and D3 is
> actively harmful."

### D2 has a race (codex)

> "D2's 'load if outstanding was 0' conflates lifecycle with demand. **A counter of remaining tasks
> cannot safely be both the residency state and the active-user count.**"

Correct, and fatal as written. kimi independently found two more problems with D2: the MS1 map is a
*second* resource every unit needs, making it a two-permit acquisition and a deadlock risk; and LPT
ordering actively fights window locality, because units sorted by size come from many different
windows, maximising the number resident and making workers block on permits while work exists — the
very idle-core problem the design was meant to fix.

### The priority is wrong (kimi, Q7)

> "extraction→0 = 1.90x; **tail→0 = 2.11x**... ship the memory fixes (already built/tested), do the
> tail (bigger Amdahl slice, lower risk), and only then revisit scheduling if occupancy is proven to
> be imbalance not bandwidth."

### The already-built fixes ARE the bandwidth fixes (kimi, Q8)

> "With doubling realloc, ~2.7M reallocations/batch copying ~2x data → several hundred GB of memcpy
> traffic. **That's the bandwidth hog — and D4/reserve/compact store directly reduce it**
> (16→4 B = 4x less traffic, realloc copies gone)."

### Numerics (kimi, Q5)

The float32 risk is **not** bias but *discontinuity*: argmax/apex selection, MI histogram bin edges,
and ties flip a discrete choice. And accumulation must stay in double — a float32 sum over 2392
points with apex 1e6 and tail 1e2 loses the tail, up to ~1e-4 relative. Recommended acceptance test:
max |Δscore| < 1e-4 and **zero** feature-set changes at q=0.01.

### `reserve` (kimi, Q6)

Safe: glibc serves a 38 KB request from bins at 16-byte granularity, so the 42% really is recovered.
Under jemalloc/tcmalloc the size classes are coarser (~1.25× steps) and ~20% would survive. And
reserve is never *semantically* wrong — a wrong count is only suboptimal, provided the code uses
`reserve` and not `resize`. (It does.)

## What BOTH reviewers missed — the premise of D1 is false

D1 assumes there are many (window, batch) units to schedule. There are not:

```
423,079 precursors over 150 windows = 2,821 per window (mean)
batch_size = 10,000  ->  batches per window = 1
per-window spread: min 128 -> 1 batch, median 1,278 -> 1 batch, max 13,304 -> 2 batches
```

**Essentially every window has exactly one batch.** So:

- D1's "flat list of several hundred to a few thousand units" does not exist — it is 150 units, the
  same 150 the current code already has.
- kimi's "strongest adversarial point" — that fixing the `threads_outer_loop_ = -1` division would
  give nested batch parallelism and make D2 unnecessary — is also wrong *for this workload*. Nested
  OpenMP **is** compiled in (`MT_ENABLE_NESTED_OPENMP:BOOL=ON`), so the team-size bug is real, but
  fixing it would parallelise a loop with one iteration.

The actual prerequisite for finer parallelism is to **shrink the batch**, which is already an
exposed option (`-batchSize`), and which also bounds the chromatogram term:

```
-batchSize 10,000 -> ~150 units      (today)
-batchSize  1,000 -> ~300 units
-batchSize    500 -> ~450 units
```

## Verdict

| | status | why |
|---|---|---|
| D1 flat list + LPT | **dropped for now** | premise false: ~1 batch/window. Needs `-batchSize` shrunk first |
| D2 residency refcount | **dropped** | codex: lifecycle/demand race. kimi: MS1 two-permit deadlock, LPT fights locality |
| D3 two-stage pipeline | **dropped** | kimi: conflicts with D4; actively harmful if bandwidth-bound |
| D4 pooled store | **superseded** | exact `reserve` captures the allocation win in one place; pooling adds the rest, second-order |
| D5 kill critical section | **kept, not built** | still valid, unmeasured |
| D6 sizing model | **revised** | missing: resident window spectra (kimi calls it the dominant term), MS1 chromatograms, scoring temporaries, feature buffers |

**What was actually built instead** — the three things the review agreed were the real wins, all
measured:

1. `chrom_list.clear()` after `return_chromatogram` — removes a full duplicate of every chromatogram
   held across the whole scoring phase. Under A/B against the 366.8 GB baseline.
2. Exact `reserve()` in `ChromatogramExtractorAlgorithm` — kills 42% capacity slack (~156 GB at this
   shape) and ~2.7M reallocations per batch. All 15 OpenMS extractor tests pass unchanged.
3. `ChromatogramStore` (4.09×, 225 → 55 GB) — built and tested, **not yet integrated**.

The lesson worth keeping: the elaborate design was the wrong instinct. Two of its six parts were
unsafe, one was self-contradictory, one rested on a false premise about the data, and the parts that
actually mattered were three localised changes that needed no new architecture at all.
