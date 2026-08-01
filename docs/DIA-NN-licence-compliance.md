# DIA-NN licence — what it means for this project

Source: <https://github.com/vdemichev/diann>. Version in use: **DIA-NN 2.0 Academia**,
installed at `/ceph/ibmi/abi/oliver/opt/diann/diann-2.0`. Licence text supplied by Oliver
2026-07-29; this note records how it constrains our use. Not legal advice.

## Clauses that bind us

**§3 — academic research or education only; no commercial or for-profit use.**
OpenDIAlyzer is open-source academic work, so the purpose is permitted.

**§3 — "Usage of DIA-NN in the cloud or on any computer setups that are not in your immediate
possession is only permitted under Collaborative use."**
**RESOLVED 2026-07-29 (Oliver):** the IBMI cluster is in our immediate possession, so this is
not cloud use and Collaborative use is not required. No further action; DIA-NN may run on the
IBMI nodes, including from the shared install at `/ceph/ibmi/abi/oliver/`.

**§6 — attribution is mandatory.** Any shared, distributed or published material produced using
DIA-NN must state so. Use "Processed using DIA-NN" in papers, preprints, READMEs, figures and
benchmark tables that include a DIA-NN arm.

**§3 — no decompiling, reverse engineering, deriving source, or creating derivative works.**
Reinforces the existing clean-room posture (`CLEAN-ROOM.md`): benchmarking DIA-NN's OUTPUT is
fine; inferring its algorithms from its behaviour to reimplement them is not. Our reconstruction
of DIA-NN's algorithm surface came from its **publications and documentation**, which is
legitimate, and must stay that way.

**§8 — Experimental functionality.** Any option NOT referenced in the DIA-NN README is
Experimental, as are `--extract` and `--mgf` explicitly. Publishing results derived from
Experimental functionality before the corresponding peer-reviewed paper appears requires
collaboration and co-authorship with the Demichev laboratory.

## What we have actually run

    diann --f <mzML> --lib <library.tsv> --threads N --qvalue 0.01 \
          --out <out> --temp <dir> --no-prot-inf

All of these are documented README options. **No Experimental functionality used**, so §8 does
not currently apply. Keep it that way: before adding any DIA-NN flag, check it appears in
<https://github.com/vdemichev/DiaNN/blob/master/README.md>, and record the check.

## Standing rules

1. Every benchmark table or figure with a DIA-NN arm carries "Processed using DIA-NN".
2. No DIA-NN flag is used without confirming it is in the README (§8).
3. DIA-NN is used as a black-box reference only — never to derive its internals (§3).
4. No commercial or for-profit use, and no redistribution of the binary (§5).
5. §3 possession: settled — the IBMI cluster is in our immediate possession (2026-07-29).
