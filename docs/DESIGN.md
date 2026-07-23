# Design

Working architecture for targeted DIA extraction. Revised after an independent
adversarial review (two external reviewers) that invalidated three load-bearing
claims of the first draft. Those corrections are recorded here rather than
quietly dropped, because they determine what gets built.

## How the incumbents differ

| | **OpenSWATH** | **DIA-NN** |
|---|---|---|
| Paradigm | Peptide-centric | Peptide-centric + spectrum-centric arbitration |
| Intermediate | Materialises full-range XIC chromatograms as first-class objects | Extracts in-memory, locally around a predicted RT |
| Scores | ~40 | 69 |
| Classifier | LDA / mProphet, or PyProphet / Percolator | Linear classifier for peak choice → DNN ensemble for q-value |
| Interference | Not explicitly modelled | RT-collision arbitration (ID) + reference-fragment subtraction (quant) |
| Hardware | CPU, OpenMP | CPU, OpenMP — **no GPU** |

The decisive difference is the **intermediate**. OpenSWATH's
extract-everything-then-score design makes the chromatogram the unit of work,
costing memory and IO proportional to *library size × RT range*. DIA-NN makes
the *candidate* the unit of work. We follow the latter.

## The kernel

> For each candidate precursor, gather its *F* fragment m/z values across *T*
> retention-time points within the containing isolation window, forming an
> **F × T trace matrix**; reduce to a score vector; classify.

## Three corrections from review

**1. Routing is not filtering.** The first draft proposed rejecting candidates
"by arithmetic on mass, RT and IM". This is close to vacuous: the isolation
window is a *partition* — every candidate belongs to some window — and a
predicted RT with a tolerance merely assigns a candidate to a cell. Neither
discards anything. Summed over cells, essentially all candidates still reach the
kernel.

The only mechanisms that genuinely reject are **ion mobility** and a
**fragment-presence prefilter**. The prefilter was unspecified, which means the
algorithm was unspecified. It is now the central design problem.

**2. This is the GPU's worst case, not its best.** Reducing an F × T matrix
looks like a dense batched operation, but the *gather* that fills it is not:
matching a fragment m/z within tolerance against an m/z-sorted scan is a binary
search, giving on the order of thousands of **dependent** random accesses per
candidate. Effective arithmetic intensity lands near 0.02–1 op/byte — far below
an H100's compute balance. This is latency-bound sparse lookup. A many-core CPU
with large aggregate cache may simply win.

**3. Resident indices do not fit.** A one-hour Orbitrap Astral run is on the
order of 10⁹–10¹⁰ centroided peaks — tens to hundreds of GB. timsTOF diaPASEF
adds an ion-mobility dimension and is worse. Streaming by RT block is mandatory.

### The consequent design decision

Store each (window, RT-bin) scan as a **direct-addressed m/z bin array** at the
extraction tolerance, so a fragment lookup is `idx = (mz - mz0) * inv_bin` —
O(1), coalesced, no dependent chain. The same structure supplies the missing
prefilter for free: an **occupancy bitmap** per cell lets a candidate be rejected
by a popcount of its fragment mask against the observed mask, before any gather.

This converts the hot path from pointer-chasing into bounded regular access. It
is the core technical bet of the project, and it must be **validated by
measurement on real data before any kernel is written** — specifically, whether
the bitmap rejects enough candidates to matter.

## Pipeline

```
A. INGEST         mzPeak / mzML → per-window, RT-binned, direct-addressed
                  m/z bins + occupancy bitmaps          [streamed, not resident]
B. CALIBRATE      mass calibration, RT and IM alignment
                  (must precede C — the filter depends on calibrated tolerances)
C. CANDIDATES     library, or FASTA → digest → predicted RT/IM/fragments
                  streamed, never materialised
D. PREFILTER      IM window + occupancy-bitmap popcount     ← the real filter
E. EXTRACT+SCORE  gather F×T → score vector                 ← the hot kernel
F. CLASSIFY       score vector → NN (ONNX Runtime)
G. ARBITRATE+FDR  RT-collision and nested-family resolution; q-values
H. QUANTIFY       reference-fragment interference subtraction; normalisation
```

MS1 extraction and precursor-level scoring run alongside E; DIA-NN leans on MS1
heavily and omitting it forfeits a large share of the available evidence.

## Parallelism

| Level | Unit | Notes |
|---|---|---|
| 1 | run / raw file | trivially parallel |
| 2 | isolation window × RT block | independent m/z ranges, no sharing |
| 3 | candidate precursor | the level that matters |
| 4 | fragment × RT point | SIMD lane |

Structure-of-Arrays throughout, so level 4 is coalesced and level 3 is a stride.
Kernels are written as plain data-parallel functions over SoA spans — no
abstraction framework. If a GPU port is ever justified by profiling, the same
source compiles for it.

## Non-functional requirements

- **Must run on a 16-core / 128 GB workstation with no GPU.** An open tool that
  requires a 200-core machine cannot be used or reproduced by the labs it is
  meant to serve. Large hardware is for development and benchmarking only.
- Bounded, predictable memory. Streaming, not residency.
- Every performance and FDR claim reproducible from public data.

## Success criterion

Not "more identifications than DIA-NN" — with a search space this large and a
compromised decoy model, more identifications is the expected symptom of broken
FDR.

Split by path, because review showed a single criterion conflates two problems:

- **Library path** — FDR calibration *measured* by allotype entrapment, reported
  together with the known biases of that instrument, and never tuned against
  (tuning until a benchmark reads 1% is overfitting the benchmark).
- **FASTA / non-specific path** — entrapment over a ~10⁴-peptide library says
  nothing about calibration in a ~10⁸-candidate space. This needs a separate
  null model. **Unsolved, and tracked as the project's central research
  question.** No calibration claim will be made here until it is solved.
