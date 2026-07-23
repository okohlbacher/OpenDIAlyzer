#!/usr/bin/env python3
"""Add competitive features to an OpenSWATH .osw, DIA-NN style.

OpenSWATH scores every candidate peak group in isolation: how well does this
precursor's trace look, on its own terms. DIA-NN 1.7.x additionally scores each
candidate *relative to its rivals at the same locus* and feeds that to the
classifier as ordinary inputs (`pBestCorrDelta`, `pTotCorrSum` in diann.cpp).
That is cheap, and it is information OpenSWATH's feature vector simply does not
contain: two candidates with identical isolated scores are not equally credible
if one of them is the best explanation of its retention-time region and the
other is the fourth-best.

Idea from Demichev et al., Nat Methods 17:41-44 (2020); diann.cpp @ cff0408,
CC BY 4.0. Implemented independently -- no code copied.

Adds four VAR_ columns to FEATURE_MS2, which PyProphet picks up automatically:

  VAR_ODIA_BEST_DELTA  own score minus the best *rival* score at this locus
                       (positive = this is the best explanation here)
  VAR_ODIA_TOT_RATIO   log(own / total at locus) -- share of the local evidence
  VAR_ODIA_RANK        fractional rank among competitors, 0 = best
  VAR_ODIA_N_COMPET    log1p(number of competitors)

CRITICAL: targets and decoys compete in the same pool. Letting only targets
compete would penalise targets and leave decoys untouched, inflating the
target-decoy separation and destroying FDR validity. The competition must be
blind to decoy status -- that is what keeps target-decoy symmetry.

    odia_rescore.py <features.osw> [--score VAR_XCORR_SHAPE] [--rt-tol 15]
"""
import sqlite3
import sys

import numpy as np
import pandas as pd

# Precursors only compete if the instrument could have co-isolated them, i.e.
# same DIA isolation window. 24 Th matches this acquisition; see odia-info.
SWATH_WIDTH = 24.0


def competitive_features(df, rt_tol=15.0):
    """df: feature_id, mz, rt, score, decoy. Returns df + 4 competitive columns.

    Locus = same isolation window and within rt_tol seconds. rt_tol should be
    about the chromatographic peak half-width, so that genuinely co-eluting
    candidates compete and unrelated ones do not.
    """
    out = df.reset_index(drop=True).copy()
    out["_win"] = np.floor(out["mz"] / SWATH_WIDTH).astype(int)

    # Indexed by ROW POSITION IN `out`, not position within a group. Writing at
    # within-group positions makes every window overwrite the same low indices
    # and silently zeroes most of the output -- a single-window test cannot see
    # this, because there the two indexings coincide.
    best = np.zeros(len(out))
    total = np.zeros(len(out))
    rank = np.zeros(len(out))
    ncomp = np.zeros(len(out))

    for _, g in out.groupby("_win", sort=False):
        order = np.argsort(g["rt"].to_numpy())
        rt_s = g["rt"].to_numpy()[order]
        sc_s = g["score"].to_numpy()[order]
        row_s = g.index.to_numpy()[order]          # positions in `out`
        lo = np.searchsorted(rt_s, rt_s - rt_tol, side="left")
        hi = np.searchsorted(rt_s, rt_s + rt_tol, side="right")
        for i, row in enumerate(row_s):
            s, e = lo[i], hi[i]
            window = sc_s[s:e]
            n = len(window) - 1                       # exclude self
            own = sc_s[i]
            if n > 0:
                # best RIVAL: drop one instance of self, not all ties
                rival = np.delete(window, i - s)
                best[row] = own - rival.max()
                rank[row] = (window > own).sum() / n
            total[row] = np.log((abs(own) + 1e-6) / (abs(window.sum()) + 1e-6))
            ncomp[row] = n

    out["VAR_ODIA_BEST_DELTA"] = best
    out["VAR_ODIA_TOT_RATIO"] = total
    out["VAR_ODIA_RANK"] = rank
    out["VAR_ODIA_N_COMPET"] = np.log1p(ncomp)
    return out.drop(columns=["_win"])


def selftest():
    # Candidates 1-4 in one isolation window; 5-7 in a DIFFERENT window.
    # The second window is essential: with only one window, a bug that indexes
    # by within-group position instead of global row position is invisible.
    df = pd.DataFrame({
        "feature_id": [1, 2, 3, 4, 5, 6, 7],
        "mz": [500.0, 500.0, 500.0, 500.0, 900.0, 900.0, 900.0],
        "rt": [100.0, 105.0, 110.0, 5000.0, 200.0, 205.0, 9000.0],
        "score": [0.9, 0.5, 0.3, 0.4, 0.7, 0.2, 0.6],
        "decoy": [0, 0, 1, 0, 0, 1, 0],
    })
    r = competitive_features(df, rt_tol=15.0).set_index("feature_id")

    # Second window must be computed independently and NOT overwritten.
    assert abs(r.loc[5, "VAR_ODIA_BEST_DELTA"] - 0.5) < 1e-9, r.loc[5]
    assert abs(r.loc[6, "VAR_ODIA_BEST_DELTA"] - (-0.5)) < 1e-9, r.loc[6]
    assert r.loc[7, "VAR_ODIA_N_COMPET"] == 0.0      # isolated in RT
    assert r.loc[5, "VAR_ODIA_N_COMPET"] > 0
    # a peptide in window A must never compete with one in window B
    assert r.loc[1, "VAR_ODIA_N_COMPET"] == np.log1p(2)
    # best candidate beats its best rival by 0.9-0.5
    assert abs(r.loc[1, "VAR_ODIA_BEST_DELTA"] - 0.4) < 1e-9, r.loc[1]
    # second-best is behind the leader
    assert abs(r.loc[2, "VAR_ODIA_BEST_DELTA"] - (-0.4)) < 1e-9
    assert r.loc[1, "VAR_ODIA_RANK"] == 0.0            # nothing beats it
    assert abs(r.loc[3, "VAR_ODIA_RANK"] - 1.0) < 1e-9  # worst of three
    # the isolated one has no competitors
    assert r.loc[4, "VAR_ODIA_BEST_DELTA"] == 0.0
    assert r.loc[4, "VAR_ODIA_N_COMPET"] == 0.0
    assert r.loc[1, "VAR_ODIA_N_COMPET"] > 0
    # decoys must participate: feature 3 is a decoy and still got a rank
    assert r.loc[3, "VAR_ODIA_N_COMPET"] > 0
    # share of local evidence is negative (own < total) and ordered
    assert r.loc[1, "VAR_ODIA_TOT_RATIO"] > r.loc[2, "VAR_ODIA_TOT_RATIO"]
    print("selftest OK")


def main():
    if len(sys.argv) == 2 and sys.argv[1] == "--selftest":
        return selftest()
    if len(sys.argv) < 2:
        raise SystemExit(__doc__)
    osw = sys.argv[1]
    score = sys.argv[3] if len(sys.argv) > 3 and sys.argv[2] == "--score" else "VAR_XCORR_SHAPE"
    rt_tol = float(sys.argv[5]) if len(sys.argv) > 5 and sys.argv[4] == "--rt-tol" else 15.0

    con = sqlite3.connect(osw)
    df = pd.read_sql(f"""
        SELECT f.ID feature_id, pr.PRECURSOR_MZ mz, f.EXP_RT rt,
               ms2.{score} score, pr.DECOY decoy
        FROM FEATURE f
        JOIN FEATURE_MS2 ms2 ON ms2.FEATURE_ID = f.ID
        JOIN PRECURSOR pr ON pr.ID = f.PRECURSOR_ID
    """, con)
    print(f"{len(df):,} features, score={score}, rt_tol={rt_tol}s")

    r = competitive_features(df, rt_tol)
    cols = ["VAR_ODIA_BEST_DELTA", "VAR_ODIA_TOT_RATIO",
            "VAR_ODIA_RANK", "VAR_ODIA_N_COMPET"]

    existing = {c[1] for c in con.execute("PRAGMA table_info(FEATURE_MS2)")}
    for c in cols:
        if c not in existing:
            con.execute(f"ALTER TABLE FEATURE_MS2 ADD COLUMN {c} REAL")
    con.executemany(
        f"UPDATE FEATURE_MS2 SET {'=?,'.join(cols)}=? WHERE FEATURE_ID=?",
        r[cols + ["feature_id"]].itertuples(index=False, name=None))
    con.commit()

    print(f"median competitors per feature: {np.expm1(r['VAR_ODIA_N_COMPET']).median():.1f}")
    for lbl, d in (("target", r[r.decoy == 0]), ("decoy", r[r.decoy == 1])):
        print(f"  {lbl:7s} best_delta median {d['VAR_ODIA_BEST_DELTA'].median():+.4f}"
              f"   rank median {d['VAR_ODIA_RANK'].median():.3f}")
    print(f"wrote {len(cols)} columns to {osw}")


if __name__ == "__main__":
    main()
