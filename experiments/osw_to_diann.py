#!/usr/bin/env python3
"""Convert our OSW-TSV library (targets only) to a DIA-NN TSV library, so DIA-NN's
engine can search OUR predicted library. This isolates library vs engine: if
DIA-NN + our library gives ~40k it was OpenSWATH's fault; if ~5k it's our library.

Header-driven. Modified.Sequence: strip our leading N-term '.' marker (DIA-NN puts
the N-term mod as (UniMod:1) at the start with no dot). Let DIA-NN make its own
decoys, so emit targets only.  Usage: osw_to_diann.py <in.tsv> <out.tsv>
"""
import sys, csv
IN, OUT = sys.argv[1], sys.argv[2]
csv.field_size_limit(1 << 24)
with open(IN) as f, open(OUT, "w", newline="") as o:
    r = csv.reader(f, delimiter="\t"); H = {n: i for i, n in enumerate(next(r))}
    w = csv.writer(o, delimiter="\t", lineterminator="\n")
    w.writerow(["ModifiedPeptide", "StrippedPeptide", "PrecursorCharge", "PrecursorMz",
                "Tr_recalibrated", "IonMobility", "ProductMz", "LibraryIntensity",
                "FragmentType", "FragmentCharge", "FragmentSeriesNumber", "ProteinId"])
    g = lambda c, n: c[H[n]]
    n = 0
    for c in r:
        if len(c) <= H["Decoy"] or c[H["Decoy"]] != "0":  # targets only
            continue
        mod = g(c, "ModifiedPeptideSequence")
        if mod.startswith("."):
            mod = mod[1:]                      # drop N-term dot marker
        w.writerow([mod, g(c, "PeptideSequence"), g(c, "PrecursorCharge"),
                    g(c, "PrecursorMz"), g(c, "NormalizedRetentionTime"),
                    g(c, "PrecursorIonMobility"), g(c, "ProductMz"),
                    g(c, "LibraryIntensity"), g(c, "FragmentType"),
                    g(c, "ProductCharge"), g(c, "FragmentSeriesNumber"),
                    g(c, "ProteinId")])
        n += 1
print(f"wrote {n:,} transition rows", file=sys.stderr)
