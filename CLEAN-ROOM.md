# Clean-room provenance

OpenDIAlyzer contains **no DIA-NN code**, and none was ever available to copy.

## Finding of fact

As of 2026-07-23, the DIA-NN repository (`github.com/vdemichev/DiaNN`, commit
`125b297`) **does not publish DIA-NN's source code**. Of 587 tracked files, the
only C/C++ outside vendored third-party libraries is `src/cpu_info.h` (321
lines of CPU feature detection, itself third-party). The MSVC project file
`diann/diann.vcxproj` references `..\src\diann.cpp` — **that file is not present
in the repository.**

What the repository does contain: a C# WinForms GUI, `WiffReader.cpp`, MSVC and
installer project files, `README.md`, and vendored copies of Eigen, MSToolkit,
MiniDNN and cranium.

DIA-NN is distributed as proprietary binaries under two tiers — "Enterprise"
(commercial, via Aptila Biotech / Thermo Fisher) and "Academia" (free,
feature-limited). Its licence terms are not distributed with the Linux release.

## What this means

Clean-room separation here is **structural, not procedural**. There is no
DIA-NN source to have been contaminated by. No DIA-NN code, headers, model
weights, or data files are copied, adapted, translated, or linked.

## What OpenDIAlyzer *is* derived from

- **Published literature.** Demichev et al., *Nat Methods* 17:41–44 (2020), and
  the CC-BY-4.0 preprint bioRxiv 10.1101/282699 — algorithmic descriptions
  intended by their authors for public dissemination.
- **Public documentation.** The DIA-NN `README.md` and GUI parameter surface,
  used to understand *interface and semantics*, not implementation.
- **OpenSWATH** (OpenMS, BSD-3) — an open, independently developed targeted DIA
  implementation.
- **Independent design.** All algorithms here are our own, built on OpenMS.

Ideas, algorithms and mathematical methods described in published papers are not
protected by copyright. Reimplementing a published algorithm from its
description is lawful and is the normal mechanism of scientific reproduction.

## Benchmarking

DIA-NN is used **only as an unmodified black-box baseline**, executed on its own
published binaries. It is not disassembled, decompiled, or reverse-engineered.

> **Open item:** DIA-NN's `LICENSE.txt` is referenced by its release notes but
> is **not shipped** in `DIA-NN-2.6.1-Academia-Linux.zip`. Until those terms are
> obtained and read, no published comparative benchmark should be considered
> cleared. See `vault/10-DIA-NN/DIA-NN licence terms are unverified.md`.

## Third-party components we *do* use

| Component | Licence | Use |
|---|---|---|
| OpenMS | BSD-3-Clause | Core MS data structures, IO, OpenSWATH reference |
| mzPeak (C++) | see upstream | mzPeak reader/writer |
| ONNX Runtime | MIT | Neural-network inference |

Neither OpenMS nor DIA-NN is modified by this project, and this project never
pushes to either repository.
