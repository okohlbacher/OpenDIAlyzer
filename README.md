# OpenDIAlyzer

Open error control for **non-specific DIA search**, with a focus on
**immunopeptidomics** — and, if the method earns it, an extraction engine to go
with it.

> **Status: early. No usable release.** The architecture has been through two
> rounds of adversarial review, which changed it substantially. See
> [docs/DESIGN.md](docs/DESIGN.md).

## Why

Immunopeptidomics breaks the assumptions every mainstream DIA tool is built on.
HLA class I peptides are non-tryptic, 8–12 residues, frequently singly charged,
and drawn from a search space two to three orders of magnitude larger than a
tryptic one. Worse, they come in **nested families** (`SLYNTVATL`,
`SLYNTVATLY`, `LYNTVATL`) that share most of their fragment ions — which
undermines the decoy model that target–decoy FDR depends on.

A published four-tool DIA immunopeptidomics benchmark
([Pandey et al., MCP 2023](https://doi.org/10.1016/j.mcpro.2023.100515))
reported that **every tool's real false-discovery rate exceeded its nominal 1%**:

| Tool | Nominal | Measured (HLA-B\*57:01) | Measured (HLA-A\*02:01) |
|---|---|---|---|
| DIA-NN | 1% | 1.45% | 3.60% |
| Skyline | 1% | 1.25% | 4.44% |
| Spectronaut | 1% | 2.03% | 2.69% |
| PEAKS | 1% | 1.75% | 4.25% |

Only 33–34% of peptides were found by all four tools, and the authors
recommended running at least two pipelines and intersecting them.

> **Treat these numbers with caution — we do.** That study used **DIA-NN 1.8.0**,
> five major versions before the 2.0 rewrite that added tunable decoys and
> non-specific digest support aimed explicitly at immunopeptidomics. It computed
> FDR by counting cross-allotype identifications as false positives, but the
> methods do not state whether isoleucine/leucine ambiguity was controlled — and
> DIA-NN's own documentation calls I/L equivalence *"essential when benchmarking
> using entrapment databases"*. Nor do they address peptides that genuinely bind
> both allotypes, which would be counted as false positives while being correct.
> Both effects inflate the apparent FDR.
>
> So this table motivates the question; it does not settle it. Re-measuring it
> on a current DIA-NN with those confounds controlled is the project's first
> piece of work, and it may well shrink the gap.

Separately, DIA-NN's own documentation states that InfinDIA — the mode you must
use for non-specific search — produces q-values that "may deviate from
externally controlled FDR estimates by a factor of up to approximately 1.5–2×
in either direction."

**What this project does not yet know**, and is measuring first: whether that
error changes any real conclusion. Adversarial review made the fair point that
the field's *de facto* error control is the binding-affinity post-filter
(NetMHCpan %rank), which discards most reported identifications regardless of
q-value. If a 2× spectral-FDR error washes out downstream, the premise is weak
and we should say so.

## Approach: statistics first, engine second

Two independent reviews concluded that building another DIA engine is not
justified — DIA-NN owns throughput, OpenSWATH is open and validated. So:

1. **Measure what the error costs.** Run DIA-NN non-specific single-run and
   quantify the deviation's real impact. This can end the project.
2. **Test whether family-level FDR is even possible.** If the fragment-sharing
   graph over non-specific human peptides percolates into one giant component,
   family-level FDR is proteome-level FDR and the idea dies. Pure computation.
3. **Declare the error unit.** Per precursor, peptidoform, family, allele, run,
   patient? A q-value without a stated discovery unit is not statistics.
4. **Build the method as a post-processor** over DIA-NN and OpenSWATH outputs.
   The statistics does not need a new engine.
5. **Build extraction only where existing engines cannot emit what the method
   needs.**

Speed is explicitly *not* a goal: DIA-NN searches the whole non-specific human
space in roughly 45 seconds on 64 CPU cores. There is no performance crisis.

## Design principles

1. **Correctness before speed.** A fast tool with mis-calibrated FDR
   manufactures false biology.
2. **Runs on a workstation.** 16 cores, 128 GB, no GPU. Large machines are for
   development, not a condition of use.
3. **Measured, not asserted.** Every performance and FDR claim is reproducible
   from public data. Selectivity claims must come with sensitivity.
4. **Validate, never tune.** Tuning until a benchmark reads 1% is overfitting
   the benchmark.

## Relationship to other software

**OpenDIAlyzer contains no DIA-NN code.** None was available to copy — DIA-NN's
source is not published. See [CLEAN-ROOM.md](CLEAN-ROOM.md).

- **[OpenMS](https://github.com/OpenMS/OpenMS)** (BSD-3) — a prerequisite, not a
  host. We link against an installation and live outside its source tree, and
  never push to it.
- **[mzPeak](https://github.com/OpenMS/mzpeak)** — reader/writer for the
  HUPO-PSI mzPeak format.
- **OpenSWATH** (in OpenMS) — the open reference implementation of targeted DIA
  extraction; our comparison baseline, and the source of the validated score
  definitions and test vectors we use as a correctness oracle.

## Build

```bash
cmake -B build -DOpenMS_DIR=/path/to/openms/lib/cmake/OpenMS
cmake --build build -j
```

## Licence

BSD-3-Clause. See [LICENSE](LICENSE).
