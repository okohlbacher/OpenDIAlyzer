#!/usr/bin/env python3
"""ROOT-CAUSE GATE (codex/vibe): give OpenSWATH our library with OBSERVED RT and
see if IDs jump. Rewrites NormalizedRetentionTime to DIA-NN's observed run-RT
(report.parquet) for precursors DIA-NN found; predicted-mapped-to-run-scale for
the rest. If IDs rise sharply, RT accuracy is the dominant lever and the fix is a
two-pass calibration, not a new engine. Usage:
    rt_gate.py <first_pass.osw> <diann_report.parquet> <in_lib.tsv> <out_lib.tsv>
"""
import sqlite3, sys, csv
import numpy as np, pandas as pd
csv.field_size_limit(1 << 24)
FIT_OSW, DIANN, IN_LIB, OUT_LIB = sys.argv[1:5]

# transform: library 0-1 RT -> observed run-seconds (from first-pass confident IDs)
c = sqlite3.connect(FIT_OSW)
rows = c.execute("""SELECT p.LIBRARY_RT, f.EXP_RT FROM SCORE_MS2 s
    JOIN FEATURE f ON s.FEATURE_ID=f.ID JOIN PRECURSOR p ON f.PRECURSOR_ID=p.ID
    WHERE s.QVALUE<0.01 AND p.DECOY=0 AND p.LIBRARY_RT IS NOT NULL""").fetchall()
lib = np.array([r[0] for r in rows]); obs = np.array([r[1] for r in rows])
o = np.argsort(lib); lib, obs = lib[o], obs[o]
edges = np.quantile(lib, np.linspace(0, 1, 61)); edges[-1] += 1e-9
idx = np.clip(np.digitize(lib, edges) - 1, 0, 59)
cx, cy = [], []
for b in range(60):
    m = idx == b
    if m.sum() >= 5: cx.append(np.median(lib[m])); cy.append(np.median(obs[m]))
cx, cy = np.array(cx), np.array(cy)
pred_to_s = lambda x: float(np.interp(x, cx, cy))

# DIA-NN observed run-RT (minutes -> seconds), keyed (stripped, charge)
d = pd.read_parquet(DIANN, columns=["Stripped.Sequence","Precursor.Charge","RT","Q.Value"])
d = d[d["Q.Value"] < 0.01]
obs_rt = {(s, int(z)): float(rt) * 60.0
          for s, z, rt in zip(d["Stripped.Sequence"], d["Precursor.Charge"], d["RT"])}
print(f"observed-RT anchors from DIA-NN: {len(obs_rt):,}", file=sys.stderr)

nobs = npred = 0
with open(IN_LIB) as fi, open(OUT_LIB, "w", newline="") as fo:
    r = csv.reader(fi, delimiter="\t"); H = {n: i for i, n in enumerate(next(r))}
    w = csv.writer(fo, delimiter="\t", lineterminator="\n"); w.writerow(list(H.keys()))
    iSEQ, iCH, iRT = H["PeptideSequence"], H["PrecursorCharge"], H["NormalizedRetentionTime"]
    for row in r:
        if len(row) > iRT:
            key = (row[iSEQ], int(row[iCH])) if row[iCH].isdigit() else None
            if key in obs_rt:
                row[iRT] = f"{obs_rt[key]:.4f}"; nobs += 1
            else:
                try: row[iRT] = f"{pred_to_s(float(row[iRT])):.4f}"; npred += 1
                except ValueError: pass
        w.writerow(row)
print(f"rewrote RT: {nobs:,} rows observed, {npred:,} rows predicted-to-run-scale -> {OUT_LIB}",
      file=sys.stderr)
