#!/usr/bin/env python3
"""Iterative RT recalibration (DIA-NN-style), step 1 of the recalibration engine.

OpenSWATH calibrates RT once, from ~100 CiRT anchor peptides. This re-fits the
library RT -> observed RT transform using ALL first-pass confident IDs (thousands),
which covers the gradient far better, then rewrites the library's RT into run-scale
so a second search can use a tight window (less interference -> more IDs).

    recalibrate_rt.py <first_pass.osw> <in_library.tsv> <out_library.tsv>

Robust binned-median fit (isotonic-ish), applied to every precursor by interpolation.
Reports the recalibrated residual vs the first-pass residual.
"""
import sqlite3, sys, csv, bisect
import numpy as np

FIT_OSW, IN_LIB, OUT_LIB = sys.argv[1], sys.argv[2], sys.argv[3]
csv.field_size_limit(1 << 24)

# --- gather (library_rt, observed_rt) from first-pass confident IDs ---
c = sqlite3.connect(FIT_OSW)
rows = c.execute("""SELECT p.LIBRARY_RT, f.EXP_RT FROM SCORE_MS2 s
    JOIN FEATURE f ON s.FEATURE_ID=f.ID JOIN PRECURSOR p ON f.PRECURSOR_ID=p.ID
    WHERE s.QVALUE<0.01 AND p.DECOY=0 AND p.LIBRARY_RT IS NOT NULL""").fetchall()
lib = np.array([r[0] for r in rows], float); obs = np.array([r[1] for r in rows], float)
order = np.argsort(lib); lib, obs = lib[order], obs[order]
print(f"fit from {len(lib):,} confident IDs", file=sys.stderr)

# --- robust binned-median transform library_rt -> observed_rt ---
nb = 60
edges = np.quantile(lib, np.linspace(0, 1, nb + 1)); edges[-1] += 1e-9
idx = np.clip(np.digitize(lib, edges) - 1, 0, nb - 1)
cx, cy = [], []
for b in range(nb):
    m = idx == b
    if m.sum() >= 5:
        cx.append(np.median(lib[m])); cy.append(np.median(obs[m]))
cx, cy = np.array(cx), np.array(cy)
def transform(x): return np.interp(x, cx, cy)

resid_before = obs - np.interp(lib, [lib.min(), lib.max()],
                               [obs.min(), obs.max()])  # ~ the crude linear map
resid_after = obs - transform(lib)
print(f"residual sd: crude-linear {resid_before.std():.1f}s  ->  recalibrated {resid_after.std():.1f}s",
      file=sys.stderr)
print(f"recalibrated 95%|resid| = {np.percentile(np.abs(resid_after),95):.1f}s", file=sys.stderr)

# --- rewrite library RT (col NormalizedRetentionTime) to run-scale refined RT ---
with open(IN_LIB) as fi, open(OUT_LIB, "w", newline="") as fo:
    r = csv.reader(fi, delimiter="\t"); H = {n: i for i, n in enumerate(next(r))}
    w = csv.writer(fo, delimiter="\t", lineterminator="\n")
    w.writerow(list(H.keys())); rt = H["NormalizedRetentionTime"]
    n = 0
    for row in r:
        if len(row) > rt:
            try: row[rt] = f"{float(transform(float(row[rt]))):.4f}"
            except ValueError: pass
        w.writerow(row); n += 1
print(f"rewrote {n:,} rows -> {OUT_LIB}", file=sys.stderr)
