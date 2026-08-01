# GBT parallelisation — design

Measured cost model of one `GBT::fit` (500k training rows, 29 features, 120 trees, depth 4, 64 bins):

| stage | ops | share |
|---|---:|---:|
| **histogram build** | 6.96e9 | **92.6%** |
| row routing | 0.24e9 | 3.2% |
| prediction update | 0.24e9 | 3.2% |
| grad/hess | 0.06e9 | 0.8% |
| split finding | 0.01e9 | 0.2% |

So there is exactly one thing to parallelise. Everything else is rounding error, and parallelising
it anyway is free only because it is on the same row axis.

## P1 — Histogram building: chunked row partition, fixed-order reduction

The histogram is `HG[node][feature][bin] += grad[row]` over all rows — a reduction, so threads
cannot write it directly. Standard fix is one private buffer per thread, reduced at the end.

**But a per-THREAD buffer makes the result depend on the thread count**, because floating-point
addition is not associative and the reduction order changes with the number of partials. These
scores produce q-values; a run that gives different answers on a different machine is not
acceptable.

So the partition is by a **fixed number of chunks, independent of the thread count**:

```
NCHUNK = min(32, max(1, n_rows/1000))          // depends only on the DATA
chunk c covers rows [c*n/NCHUNK, (c+1)*n/NCHUNK)
threads take chunks dynamically; each chunk sums into its OWN buffer
reduce buffers in chunk order 0..NCHUNK-1     // fixed order => bit-identical
```

Threads may execute chunks in any order and the answer is unchanged, because each chunk's partial
sum is independent of scheduling and the final combination order is fixed by index.

Buffer cost: `NCHUNK × nodes_at_level × features × bins × 2 doubles`
= 32 × 16 × 29 × 64 × 16 B ≈ **15 MB** — negligible against the 1.4 GB the phase already uses.

## P2 — Flatten the binned matrix

`std::vector<std::vector<uint8_t>>` is one heap allocation per row, so the inner loop over features
chases a pointer per row and the bins of one row are not adjacent to the next row's. A flat
`vector<uint8_t>` with `row*stride + f` makes the whole histogram pass a contiguous sweep. Same
arithmetic, dramatically better cache behaviour — and it is the 92.6% loop.

## P3 — Nesting

`GBT::fit` is called from inside the fold loop, which is already `#pragma omp parallel for`. By
default a nested parallel region gets a team of one, so the inner pragmas would do nothing.

`GBTParams::n_threads` is set by the caller to `max(1, omp_get_max_threads() / n_folds)` and passed
to `num_threads(...)`, with `omp_set_max_active_levels(2)` enabled once. Total width is then
`n_folds × n_threads` with no oversubscription, instead of `n_folds` alone.

## What is NOT parallelised, and why

**Split finding** (0.2%): it is a loop over nodes × features × bins with a `best_gain` comparison.
Parallelising a max-reduction with a tie rule risks changing which (feature, bin) wins on ties, and
the tie rule is what makes the fit reproducible. 0.2% is not worth that risk.

**Boosting rounds**: inherently sequential — tree `t` fits the residual left by tree `t-1`.

**Bin edges**: per feature and independent, but it runs once per fit, not per tree.
