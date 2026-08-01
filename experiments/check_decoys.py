#!/usr/bin/env python3
"""Is our target/decoy pair exchangeable, or did the decoy inherit target values?

OpenSwathDecoyGenerator shuffles the SEQUENCE but may reassign the target's
fragment intensities / RT / IM rather than predicting them for the shuffled
peptide. If so, PyProphet can learn decoy-construction artefacts instead of
true-vs-false, and the nominal 1% FDR is not trustworthy.

    check_decoys.py <library.tsv>
"""
import sys, collections

def main(fp):
    # 29-col OpenSWATH TSV: 6 RT, 7 PeptideSequence, 8 ModSeq, 5 LibraryIntensity,
    # 22 PrecursorIonMobility, 23 TransitionGroupId, 25 Decoy
    tgt, dec = {}, {}
    with open(fp) as f:
        hdr = f.readline().rstrip("\n").split("\t")
        ix = {n: i for i, n in enumerate(hdr)}
        for line in f:
            c = line.rstrip("\n").split("\t")
            if len(c) <= ix["Decoy"]:
                continue
            key = c[ix["TransitionGroupId"]]
            rec = (c[ix["NormalizedRetentionTime"]], c[ix["PrecursorIonMobility"]],
                   c[ix["PeptideSequence"]], c[ix["LibraryIntensity"]])
            (dec if c[ix["Decoy"]] == "1" else tgt)[key] = rec

    print(f"target groups: {len(tgt):,}   decoy groups: {len(dec):,}")

    # pair by stripping a DECOY_ prefix/suffix if present
    paired = 0
    same_rt = same_im = same_seq = 0
    for dk, dv in dec.items():
        tk = dk.replace("DECOY_", "").replace("_DECOY", "")
        tv = tgt.get(tk)
        if not tv:
            continue
        paired += 1
        if dv[0] == tv[0]: same_rt += 1
        if dv[1] == tv[1]: same_im += 1
        if dv[2] == tv[2]: same_seq += 1
    print(f"paired target/decoy groups: {paired:,}")
    if paired:
        print(f"  identical RT   : {same_rt:,} ({100*same_rt/paired:.1f}%)  <- copied if ~100%")
        print(f"  identical 1/K0 : {same_im:,} ({100*same_im/paired:.1f}%)  <- copied if ~100%")
        print(f"  identical seq  : {same_seq:,} ({100*same_seq/paired:.1f}%)  <- should be ~0%")
    else:
        # fall back to distribution comparison
        import statistics as st
        for nm, i in (("RT", 0), ("1/K0", 1)):
            t = [float(v[i]) for v in tgt.values() if v[i] not in ("", "NA")]
            d = [float(v[i]) for v in dec.values() if v[i] not in ("", "NA")]
            if t and d:
                print(f"  {nm}: target mean={st.mean(t):.4f} sd={st.pstdev(t):.4f} | "
                      f"decoy mean={st.mean(d):.4f} sd={st.pstdev(d):.4f}")
        # exact-value overlap is the tell-tale for copied properties
        tv = collections.Counter(v[0] for v in tgt.values())
        dv = collections.Counter(v[0] for v in dec.values())
        shared = sum((tv & dv).values())
        print(f"  RT exact-value overlap: {shared:,} of {sum(dv.values()):,} decoys "
              f"({100*shared/max(1,sum(dv.values())):.1f}%)  <- high => RT copied")

if __name__ == "__main__":
    main(sys.argv[1])
