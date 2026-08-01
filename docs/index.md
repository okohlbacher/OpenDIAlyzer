# Documentation index

Both external reviews flagged the same thing: 22 documents, five overlapping "review" files, three
competing plan/status files, and no notion of which is current. This is the table of contents.

## Start here

| Document | What it is |
|---|---|
| [`../README.md`](../README.md) | Why the project exists, and what it is not |
| [`OpenDIAlyzer-manual.md`](OpenDIAlyzer-manual.md) | **The executable**: synopsis, inputs/outputs, every option, operational notes |
| [`DESIGN.md`](DESIGN.md) | Architecture of the extraction/scoring engine |
| [`../vendored-patches/README.md`](../vendored-patches/README.md) | The OpenMS delta — base commit, every patch, how to regenerate and verify |

## Current — measured, and still true

| Document | Subject |
|---|---|
| [`OpenDIAlyzer-three-tool-report.md`](OpenDIAlyzer-three-tool-report.md) | ODIA vs DIA-NN 2.0 vs DIA-NN 1.7.12 vs OpenSwathWorkflow |
| [`OpenDIAlyzer-memory-review.md`](OpenDIAlyzer-memory-review.md) | Where peak RSS goes, with a fitted model |
| [`OpenDIAlyzer-scoring-fdr-backlog.md`](OpenDIAlyzer-scoring-fdr-backlog.md) | Scoring/FDR state and what is not implemented |
| [`OpenDIAlyzer-lda-vs-pyprophet.md`](OpenDIAlyzer-lda-vs-pyprophet.md) | Line-by-line comparison of the semi-supervised classifier |
| [`OpenDIAlyzer-chromatogram-speedup-review.md`](OpenDIAlyzer-chromatogram-speedup-review.md) | Extraction/scoring cost, ranked |
| [`OVERNIGHT-LOG.md`](OVERNIGHT-LOG.md) | Running log of iterations, with the failures kept in |
| [`DIA-NN-licence-compliance.md`](DIA-NN-licence-compliance.md) | Licence position; **read before touching anything DIA-NN-adjacent** |
| [`OpenDIALibGen.md`](OpenDIALibGen.md) | The predicted-library generator (second executable) |

## Reference

| Document | Subject |
|---|---|
| [`OpenDIAlyzer-requirements.md`](OpenDIAlyzer-requirements.md) | Requirements |
| [`OpenDIAlyzer-reference-datasets.md`](OpenDIAlyzer-reference-datasets.md) | Benchmark datasets and provenance |
| [`OpenDIAlyzer-mzpeak-backlog.md`](OpenDIAlyzer-mzpeak-backlog.md) | mzPeak integration backlog |

## Superseded — kept for provenance, do not plan from these

| Document | Superseded by |
|---|---|
| `OpenDIAlyzer-merger-plan.md` | `DESIGN.md` |
| `OpenDIAlyzer-plan-synthesized.md` | `DESIGN.md` + the backlogs |
| `OpenDIAlyzer-status-assessment.md` | `OpenDIAlyzer-three-tool-report.md` |
| `OpenDIAlyzer-pyprophet-integration.md` | abandoned: pyprophet was dropped after three runs (5h23m, 8h33m, 2h50m with 2.5 h of no output). Its role is filled by ODIA's own scorer |
| `OpenDIAlyzer-parallel-efficiency.md` | measured against 224 cores on a **shared, oversubscribed** node; the occupancy figures overstate the deficit. Use CPU-seconds and peak RSS instead |
| `OpenDIAlyzer-calibration-vs-diann.md`, `OpenDIAlyzer-mz-rt-calibration-review.md` | partly overtaken; the RT-direction and MS1-duplication defects they discuss are fixed |
| `OpenDIAlyzer-perf-review.md`, `OpenDIAlyzer-parallel-streaming-review.md`, `OpenDIAlyzer-chrom-pool-review.md` | folded into `OpenDIAlyzer-memory-review.md` and the speedup review |

## A caveat that applies to every performance number here

The benchmark node is **shared** — 18 users, load average 180–233 on 224 cores. Peak RSS and total
CPU-seconds are contention-independent and comparable. Wall clock and "average cores" measure how
many cores a run *won*, not how well it scales. `experiments/bench_run_fair.sh` records load on both
sides of a run so later numbers carry their own context.
