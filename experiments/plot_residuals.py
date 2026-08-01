#!/usr/bin/env python3
"""dRT(RT) and dmz(mz) for our predicted library vs DIA-NN, from high-confidence IDs.
Shows WHY OpenSWATH+our-library underperforms: our single-calibration RT residual
vs DIA-NN's iteratively-recalibrated one. Usage: plot_residuals.py <out.png>"""
import sqlite3, sys
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = sys.argv[1]
OURS = "/scratch/kohlbach/opendialyzer/bench/iter/P1P2_FULL/features.osw"
DIANN = "/scratch/kohlbach/opendialyzer/bench/agxt_arms/diann_pred/report.parquet"

# --- ours: DELTA_RT (residual, s) vs EXP_RT (s); mass dev vs precursor mz ---
c = sqlite3.connect(OURS)
o = pd.read_sql("""SELECT f.EXP_RT rt, f.DELTA_RT drt, p.PRECURSOR_MZ mz, m.VAR_MASSDEV_SCORE mdev
    FROM SCORE_MS2 s JOIN FEATURE f ON s.FEATURE_ID=f.ID
    JOIN PRECURSOR p ON f.PRECURSOR_ID=p.ID
    JOIN FEATURE_MS2 m ON m.FEATURE_ID=f.ID
    WHERE s.QVALUE<0.01 AND p.DECOY=0""", c)
print(f"ours: {len(o):,} IDs  RT residual sd={o.drt.std():.1f}s  95%={np.percentile(o.drt.abs(),95):.1f}s")

# --- DIA-NN: (RT-Predicted.RT) minutes -> s ; Ms1.Apex.Mz.Delta ---
d = pd.read_parquet(DIANN, columns=["RT","Predicted.RT","Precursor.Mz","Ms1.Apex.Mz.Delta","Q.Value"])
d = d[d["Q.Value"] < 0.01]
d["drt"] = (d["RT"] - d["Predicted.RT"]) * 60.0   # minutes -> seconds
d["rt_s"] = d["RT"] * 60.0
print(f"DIA-NN: {len(d):,} IDs  RT residual sd={d.drt.std():.1f}s  95%={np.percentile(d.drt.abs(),95):.1f}s")

fig, ax = plt.subplots(1, 2, figsize=(13, 5.2))
# Panel A: dRT vs RT
ax[0].scatter(o.rt, o.drt, s=3, alpha=0.15, color="#d1495b",
              label=f"ours (sd {o.drt.std():.0f}s)")
ax[0].scatter(d.rt_s, d.drt, s=3, alpha=0.15, color="#2e86ab",
              label=f"DIA-NN (sd {d.drt.std():.0f}s)")
ax[0].axhline(0, color="k", lw=.6)
ax[0].set_xlabel("observed RT (s)"); ax[0].set_ylabel("RT residual: predicted - observed (s)")
ax[0].set_title("dRT(RT): RT prediction error at high-confidence IDs")
ax[0].set_ylim(-800, 800); ax[0].legend(markerscale=4)

# Panel B: dmz vs mz.  Ours: library m/z is computed exact (transition-level
# median |d|=3e-5), so the informative mass error is DIA-NN's observed MS1 delta.
ppm = d["Ms1.Apex.Mz.Delta"] / d["Precursor.Mz"] * 1e6
ax[1].scatter(d["Precursor.Mz"], ppm, s=3, alpha=0.15, color="#2e86ab",
              label=f"DIA-NN MS1 (sd {ppm.std():.1f} ppm)")
ax[1].axhline(0, color="k", lw=.6)
ax[1].set_xlabel("precursor m/z"); ax[1].set_ylabel("MS1 mass error (ppm)")
ax[1].set_title("dmz(mz): our library m/z is exact; DIA-NN observed MS1 error shown")
ax[1].set_ylim(-20, 20); ax[1].legend(markerscale=4)
fig.tight_layout(); fig.savefig(OUT, dpi=130)
print("wrote", OUT)
