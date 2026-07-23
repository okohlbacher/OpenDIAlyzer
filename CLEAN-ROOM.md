# Clean-room provenance

OpenDIAlyzer contains **no DIA-NN code**. This is a deliberate choice, and this
document records the facts behind it accurately — including a correction to an
earlier version of this file.

## Correction (2026-07-23)

An earlier version of this document stated that DIA-NN's source "was never
available to copy." **That was wrong**, and the error is instructive: it was
based on the repository's *working tree* without checking its *history*.

`src/diann.cpp` (11 260 lines) **was published** and was later removed:

| | |
|---|---|
| Last content commit | `cff0408`, 2020-10-01 |
| Deleted in | `f2c8a78`, 2021-07-08, "Delete diann.cpp" |
| Era | DIA-NN 1.7.x — the version described by the Nature Methods 2020 paper |
| Licence then in force | **Creative Commons Attribution 4.0 International** |

From `LICENSE.txt` at that commit, verbatim:

> "The file /src/diann.cpp is licensed under the Creative Commons Attribution
> 4.0 International License."

CC BY 4.0 is permissive and, for copies already distributed, irrevocable.
Deleting the file ends the author's future distribution; it does not retract the
grant on the version that was published.

**So we could lawfully copy that code. We choose not to.**

## Current state of the upstream repository

As of 2026-07-23 (`github.com/vdemichev/DiaNN`, commit `125b297`), the working
tree does not contain DIA-NN's source. Of 587 tracked files, the only C/C++
outside vendored third-party libraries is `src/cpu_info.h`. The MSVC project
still references `..\src\diann.cpp`, which is absent from the tree.

Current DIA-NN (2.x) is distributed as proprietary binaries in two tiers —
"Enterprise" (via Aptila Biotech / Thermo Fisher) and "Academia" (free,
feature-limited). **The 2.x algorithms — InfinDIA, QuantUMS, tunable decoy
models, proteoform confidence — have never been published in any form.**

## Our position

**We use the CC BY 4.0 source as documentation, not as a code source.**

- We read it to understand the published 2020 algorithm, as one reads a paper.
- We implement independently, under BSD-3.
- We cite Demichev et al. for any idea traceable to it.
- **No DIA-NN code, headers, model weights or data files are copied, adapted,
  translated or linked into this project.**

The reasoning is practical as much as legal. CC BY is a content licence;
Creative Commons itself advises against using it for software (no patent grant,
no source/object distinction). Copying blocks of it would be lawful but would
pull CC BY terms into an otherwise BSD-3 codebase and complicate downstream
reuse — for no benefit that reading it does not already provide.

The extracted file is kept outside the repository (in a gitignored working
directory), with its contemporaneous `LICENSE.txt` alongside it so provenance
travels with the file.

## What OpenDIAlyzer is derived from

- **Published literature.** Demichev et al., *Nat Methods* 17:41–44 (2020), and
  the CC-BY-4.0 preprint bioRxiv 10.1101/282699.
- **The CC BY 4.0 `diann.cpp` (1.7.x)** — read as documentation, with
  attribution, never copied.
- **Public documentation.** The DIA-NN `README.md` and GUI parameter surface,
  for interface and semantics.
- **OpenSWATH** (OpenMS, BSD-3) — an open, independently developed targeted DIA
  implementation, and the source of validated score definitions and test vectors
  we use as a correctness oracle.
- **Independent design.** All algorithms here are our own.

Ideas, algorithms and mathematical methods described in published work are not
protected by copyright. Reimplementing a published algorithm from its
description is lawful and is the normal mechanism of scientific reproduction.

## Benchmarking

DIA-NN 2.x is used **only as an unmodified black-box baseline**, executed on its
own published binaries. It is not disassembled, decompiled, or reverse
engineered.

> **Open item:** DIA-NN 2.x's `LICENSE.txt` is referenced by its release notes
> but is **not shipped** in `DIA-NN-2.6.1-Academia-Linux.zip`; upstream issue
> [#1311](https://github.com/vdemichev/DiaNN/issues/1311) confirms it exists
> only in the Windows installer's output. Until those terms are obtained and
> read, no published comparative benchmark of 2.x should be considered cleared.
> Note this is a separate question from the 1.7.x CC BY grant above.

## Third-party components we *do* use

| Component | Licence | Use |
|---|---|---|
| OpenMS | BSD-3-Clause | Core MS data structures, IO, OpenSWATH reference |
| mzPeak (C++) | see upstream | mzPeak reader/writer |
| ONNX Runtime | MIT | Neural-network inference |

Neither OpenMS nor DIA-NN is modified by this project, and this project never
pushes to either repository.
