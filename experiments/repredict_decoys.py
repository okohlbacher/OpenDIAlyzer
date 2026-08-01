#!/usr/bin/env python3
"""P2: give decoys their OWN predicted MS2 intensities instead of the target's.

OpenSwathDecoyGenerator shuffles the sequence and computes decoy fragment m/z for
that shuffled sequence, but COPIES the target's intensity vector -- so a decoy
inheriting a peaky target spectrum scores as high as the target. Fix: re-predict
each decoy's intensities from its OWN (shuffled) sequence via OpenDIALibGen, and
overwrite only LibraryIntensity (m/z, RT, IM untouched -- codex #7).

    repredict_decoys.py <in_library.tsv> <out_library.tsv> <odia_bin> <model_dir>

Header-driven column lookup: the input library is 29-col OSW; OpenDIALibGen -raw
output is a different (18-col) layout, so both are resolved by column NAME.
"""
import sys, csv, subprocess, os, tempfile, re

IN, OUT, ODIA, MODELDIR = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
CANON = set("ACDEFGHIKLMNPQRSTVWY")

def cols(header):
    return {n: i for i, n in enumerate(header)}

# --- pass 1: collect decoy precursors (unique modseq+charge) from input lib ---
with open(IN) as f:
    r = csv.reader(f, delimiter="\t"); H = cols(next(r))
    iMOD, iPCH, iDEC = H["ModifiedPeptideSequence"], H["PrecursorCharge"], H["Decoy"]
    iFT, iFS, iFC, iIN = H["FragmentType"], H["FragmentSeriesNumber"], H["ProductCharge"], H["LibraryIntensity"]
    decoy_prec = {}
    for c in r:
        if len(c) > iDEC and c[iDEC] == "1":
            decoy_prec[(c[iMOD], c[iPCH])] = None
print(f"decoy precursors to re-predict: {len(decoy_prec):,}", file=sys.stderr)

# --- write Stage-1 peptide list (mod-aware; skip non-canonical) ---
tmp = tempfile.mkdtemp(dir=os.path.dirname(OUT) or ".")
pep_tsv = os.path.join(tmp, "decoy_peps.tsv")
nskip = 0
with open(pep_tsv, "w") as o:
    o.write("PeptideSequence\tModifiedPeptideSequence\tPrecursorCharge\tProteinId\n")
    for (ms, ch) in decoy_prec:
        strip = re.sub(r"\([^)]*\)", "", ms).replace(".", "")
        if len(strip) < 2 or any(a not in CANON for a in strip):
            nskip += 1; continue
        o.write(f"{strip}\t{ms}\t{ch}\tDECOY\n")
print(f"skipped {nskip} non-canonical decoy sequences", file=sys.stderr)

# --- run OpenDIALibGen Stage-1 (-raw) ---
pred_tsv = os.path.join(tmp, "decoy_pred.tsv")
env = dict(os.environ, LD_LIBRARY_PATH="/scratch/kohlbach/opendialyzer/build/openms-onnx/lib")
res = subprocess.run([ODIA, "-in", pep_tsv, "-out", pred_tsv, "-model_dir", MODELDIR,
                      "-threads", "96", "-raw"], env=env,
                     stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
if res.returncode != 0:
    sys.stderr.write("OpenDIALibGen failed:\n" + res.stderr.decode()[-800:] + "\n"); sys.exit(1)

# --- index predicted intensities by (modseq,charge) -> {(type,series,fz): int} ---
pred = {}
with open(pred_tsv) as f:
    r = csv.reader(f, delimiter="\t"); P = cols(next(r))
    pM, pP, pT, pS, pC, pI = (P["ModifiedPeptideSequence"], P["PrecursorCharge"],
                              P["FragmentType"], P["FragmentSeriesNumber"],
                              P["ProductCharge"], P["LibraryIntensity"])
    for c in r:
        if len(c) <= max(pM, pP, pT, pS, pC, pI): continue
        pred.setdefault((c[pM], c[pP]), {})[(c[pT].lower()[:1], c[pS], c[pC])] = c[pI]
print(f"re-predicted {len(pred):,} decoy precursors", file=sys.stderr)

# --- pass 2: rewrite decoy intensities ---
nrepl = norig = 0
with open(IN) as f, open(OUT, "w") as o:
    r = csv.reader(f, delimiter="\t"); w = csv.writer(o, delimiter="\t", lineterminator="\n")
    w.writerow(next(r))
    for c in r:
        if len(c) > iDEC and c[iDEC] == "1":
            p = pred.get((c[iMOD], c[iPCH]))
            fk = (c[iFT].lower()[:1], c[iFS], c[iFC])
            if p and fk in p: c[iIN] = p[fk]; nrepl += 1
            else: norig += 1
        w.writerow(c)
print(f"decoy transitions: {nrepl:,} re-predicted, {norig:,} kept original", file=sys.stderr)
