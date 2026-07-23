# OpenDIAlyzer

Open, targeted extraction of DIA mass traces for peptide identification and
quantification — with a focus on **immunopeptidomics**.

> **Status: early design.** No usable release yet. The architecture is being
> established and adversarially reviewed before implementation. See
> [docs/DESIGN.md](docs/DESIGN.md).

## Why

Immunopeptidomics breaks the assumptions every mainstream DIA tool is built on.
HLA class I peptides are non-tryptic, 8–12 residues, frequently singly charged,
and drawn from a non-specific search space two to three orders of magnitude
larger than a tryptic one. Worse, HLA peptides come in **nested families**
(`SLYNTVATL`, `SLYNTVATLY`, `LYNTVATL`) that share most of their fragment ions —
which compromises the decoy model that target-decoy FDR depends on.

The consequences are measurable. In a published four-tool DIA immunopeptidomics
benchmark ([Pandey et al., MCP 2023](https://doi.org/10.1016/j.mcpro.2023.100515)),
**every tool's real false-discovery rate exceeded its nominal 1%**:

| Tool | Nominal FDR | Measured FDR (HLA-B\*57:01) | Measured FDR (HLA-A\*02:01) |
|---|---|---|---|
| DIA-NN | 1% | 1.45% | 3.60% |
| Skyline | 1% | 1.25% | 4.44% |
| Spectronaut | 1% | 2.03% | 2.69% |
| PEAKS | 1% | 1.75% | 4.25% |

Only 33–34% of peptides were found by all four tools, and the authors' own
recommendation was to run at least two pipelines and intersect them.

**So the goal here is not "more identifications".** With a search space this
large and a compromised null model, more identifications is the expected symptom
of broken FDR. The goal is *trustworthy* identifications at competitive depth.

## Design principles

1. **Correctness before speed.** A fast tool with mis-calibrated FDR is worse
   than useless in this field — it manufactures false biology.
2. **Runs on a workstation.** 16 cores, 128 GB, no GPU must be enough for a
   standard immunopeptidomics DIA experiment. Large machines are for developing
   and benchmarking, not a requirement of use.
3. **Measured, not asserted.** Every performance and FDR claim is backed by a
   reproducible benchmark on public data.
4. **OpenMS at the boundaries.** File formats, chemistry and validation come
   from OpenMS; the extraction hot path uses our own data structures.

## Relationship to other software

**OpenDIAlyzer contains no DIA-NN code.** None was available to copy — DIA-NN's
source is not published. See [CLEAN-ROOM.md](CLEAN-ROOM.md).

- **[OpenMS](https://github.com/OpenMS/OpenMS)** (BSD-3) — a prerequisite, not a
  host. OpenDIAlyzer links against an OpenMS installation and lives outside the
  OpenMS source tree. We never push to it.
- **[mzPeak](https://github.com/OpenMS/mzpeak)** — reader/writer for the
  HUPO-PSI mzPeak format, a separate library from OpenMS.
- **OpenSWATH** (part of OpenMS) — the open reference implementation of targeted
  DIA extraction, and our primary comparison baseline.

## Build

```bash
cmake -B build -DOpenMS_DIR=/path/to/openms/lib/cmake/OpenMS
cmake --build build -j
```

## Licence

BSD-3-Clause. See [LICENSE](LICENSE).
