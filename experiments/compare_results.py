#!/usr/bin/env python3
"""Compare DIA-NN and OpenSWATH identifications on the same run and library.

Reports counts, overlap, and the length/charge breakdowns that matter for
immunopeptidomics -- a tool can win on raw count while losing the short and
singly-charged ligands that are the point of the experiment, so counts alone are
not a verdict.

    compare_results.py <diann_report.parquet> <openswath_scored.osw> [qvalue]

The .osw must have been scored by PyProphet; a raw OpenSwathWorkflow output has
features but no q-values and cannot be compared at a fixed FDR.
"""
import sqlite3
import sys

import pandas as pd

QCOL_CANDIDATES = ["Q.Value", "Global.Q.Value", "PG.Q.Value"]
SEQCOL_CANDIDATES = ["Stripped.Sequence", "Peptide.Sequence", "PeptideSequence"]


def pick(df, names, what):
    for n in names:
        if n in df.columns:
            return n
    raise SystemExit(f"error: no {what} column; have {sorted(df.columns)[:25]}")


def load_diann(path, q):
    d = pd.read_parquet(path)
    seq, qcol = pick(d, SEQCOL_CANDIDATES, "sequence"), pick(d, QCOL_CANDIDATES, "q-value")
    d = d[d[qcol] <= q]
    charge = "Precursor.Charge" if "Precursor.Charge" in d else None
    return pd.DataFrame({
        "seq": d[seq],
        "charge": d[charge] if charge else pd.NA,
    }).drop_duplicates()


def load_osw(path, q):
    # PyProphet writes SCORE_MS2; join up to the peptide via the mapping table.
    sql = """
        SELECT DISTINCT p.UNMODIFIED_SEQUENCE AS seq, pr.CHARGE AS charge
        FROM SCORE_MS2 s
        JOIN FEATURE f  ON f.ID = s.FEATURE_ID
        JOIN PRECURSOR pr ON pr.ID = f.PRECURSOR_ID
        JOIN PRECURSOR_PEPTIDE_MAPPING m ON m.PRECURSOR_ID = pr.ID
        JOIN PEPTIDE p ON p.ID = m.PEPTIDE_ID
        WHERE s.QVALUE <= ? AND pr.DECOY = 0
    """
    with sqlite3.connect(path) as con:
        try:
            return pd.read_sql(sql, con, params=(q,))
        except pd.errors.DatabaseError as e:
            raise SystemExit(f"error: {e}\nWas the .osw scored by PyProphet?")


def breakdown(df, label):
    n = len(df)
    print(f"\n{label}: {n:,} precursors, {df['seq'].nunique():,} peptides")
    if not n:
        return
    lens = df["seq"].str.len().value_counts().sort_index()
    print("  length:", "  ".join(f"{k}:{v}" for k, v in lens.items()))
    if df["charge"].notna().any():
        ch = df["charge"].value_counts().sort_index()
        print("  charge:", "  ".join(f"{k}+:{v}" for k, v in ch.items()))


def selftest():
    a = pd.DataFrame({"seq": ["SLYNTVATL", "KIFGSLAFL", "GILGFVFTL"], "charge": [2, 2, 1]})
    b = pd.DataFrame({"seq": ["SLYNTVATL", "GILGFVFTL", "NLVPMVATV"], "charge": [2, 1, 2]})
    sa, sb = set(a["seq"]), set(b["seq"])
    shared = sa & sb
    assert len(shared) == 2, shared
    assert abs(len(shared) / len(sa | sb) - 0.5) < 1e-9, "Jaccard 2/4"
    assert len(sa - sb) == 1 and len(sb - sa) == 1
    breakdown(a, "selftest-A")
    print("selftest OK")


def main():
    if len(sys.argv) == 2 and sys.argv[1] == "--selftest":
        return selftest()
    if len(sys.argv) < 3:
        raise SystemExit(__doc__)
    q = float(sys.argv[3]) if len(sys.argv) > 3 else 0.01

    dn, ow = load_diann(sys.argv[1], q), load_osw(sys.argv[2], q)
    breakdown(dn, "DIA-NN")
    breakdown(ow, "OpenSWATH")

    sa, sb = set(dn["seq"]), set(ow["seq"])
    both, union = sa & sb, sa | sb
    print(f"\noverlap (peptide sequences), q <= {q}")
    print(f"  both        {len(both):,}")
    print(f"  DIA-NN only {len(sa - sb):,}")
    print(f"  OpenSWATH   {len(sb - sa):,} only")
    print(f"  Jaccard     {len(both) / len(union):.1%}" if union else "  (empty)")
    # The MCP 2023 four-tool benchmark found only 33-34% agreement, so low
    # overlap is the expected result, not evidence that something is broken.
    print("\n  reference: MCP 2023 found 33-34% agreement across four tools")


if __name__ == "__main__":
    main()
