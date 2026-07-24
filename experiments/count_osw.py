#!/usr/bin/env python3
# Distinct target IDs at 1% MS2 FDR from an OpenSWATH .osw (post pyprophet score).
# Precursor FDR is the clean comparison to DIA-NN's Precursors.Identified; peptide
# and protein are distinct-counts among 1% precursors (no separate peptide/protein
# FDR was run), reported as proxies.
import sqlite3, sys

def counts(fp):
    c = sqlite3.connect(fp)
    q = """SELECT f.PRECURSOR_ID FROM SCORE_MS2 s
           JOIN FEATURE f ON s.FEATURE_ID=f.ID
           JOIN PRECURSOR p ON f.PRECURSOR_ID=p.ID
           WHERE s.QVALUE<0.01 AND p.DECOY=0"""
    prec = set(r[0] for r in c.execute(q))
    if not prec:
        return dict(prec=0, pep=0, prot=0)
    ids = ",".join(str(x) for x in prec)
    pep = c.execute(f"""SELECT COUNT(DISTINCT pep.MODIFIED_SEQUENCE)
        FROM PRECURSOR_PEPTIDE_MAPPING m JOIN PEPTIDE pep ON m.PEPTIDE_ID=pep.ID
        WHERE m.PRECURSOR_ID IN ({ids})""").fetchone()[0]
    try:
        prot = c.execute(f"""SELECT COUNT(DISTINCT ppm.PROTEIN_ID)
            FROM PRECURSOR_PEPTIDE_MAPPING m
            JOIN PEPTIDE_PROTEIN_MAPPING ppm ON m.PEPTIDE_ID=ppm.PEPTIDE_ID
            WHERE m.PRECURSOR_ID IN ({ids})""").fetchone()[0]
    except Exception:
        prot = -1
    return dict(prec=len(prec), pep=pep, prot=prot)

for fp in sys.argv[1:]:
    try:
        r = counts(fp)
        print(f"{fp}\tprec={r['prec']}\tpep={r['pep']}\tprot={r['prot']}")
    except Exception as e:
        print(f"{fp}\tERR {e}")
