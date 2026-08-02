#!/usr/bin/env python3
"""Where do the missing IDs go?

Compares ODIA's identifications against a reference ID list (SEQUENCE_CHARGE per line) and, for
every reference precursor ODIA does NOT report at q<0.01, asks what happened to it. That is the
question that separates the competing explanations:

  in library?  extracted?  best q     ->  diagnosis
  no           --          --             LIBRARY: not searchable, cannot be an ODIA defect
  yes          no          --             DETECTION: never extracted -- RT/mz window or prefilter
  yes          yes         0.01-0.1       FDR/SCORING: found, ranked, just under the cutoff
  yes          yes         >0.5           DISCRIMINATION: extracted but scored like noise

usage: id_gap.py <run.oswpq> <reference_ids.txt>
"""
import sys, glob, zipfile, tempfile, collections

try:
    import pyarrow.parquet as pq, pyarrow.compute as pc
except ImportError:
    sys.exit("needs pyarrow (use /scratch/kohlbach/mzpenv/bin/python)")


def load(bundle):
    d = tempfile.mkdtemp(prefix="idgap_")
    with zipfile.ZipFile(bundle) as z:
        for n in z.namelist():
            if n.endswith(("features.parquet", "score_ms2.parquet", "precursors.parquet")):
                z.extract(n, d)
    feat = pq.read_table(glob.glob(f"{d}/runs/*/features.parquet")[0],
                         columns=["feature_id", "precursor_id", "exp_rt"])
    sc = pq.read_table(glob.glob(f"{d}/runs/*/score_ms2.parquet")[0],
                       columns=["feature_id", "qvalue", "score"])
    prec = pq.read_table(f"{d}/library/precursors.parquet",
                         columns=["precursor_id", "precursor_mz", "charge", "library_rt",
                                  "decoy", "unmodified_sequence", "modified_sequence"])
    return feat.join(sc, keys="feature_id"), prec


def main(bundle, reflist):
    fs, prec = load(bundle)
    tgt = prec.filter(pc.invert(pc.cast(prec["decoy"], "bool")))

    # library key -> precursor_id, and the reverse
    keys, pids = [], tgt["precursor_id"].to_pylist()
    for seq, ch in zip(tgt["unmodified_sequence"].to_pylist(), tgt["charge"].to_pylist()):
        keys.append(f"{seq}_{ch}")
    lib_key = dict(zip(keys, pids))
    pid_meta = {p: (mz, ch, rt) for p, mz, ch, rt in
                zip(pids, tgt["precursor_mz"].to_pylist(), tgt["charge"].to_pylist(),
                    tgt["library_rt"].to_pylist())}

    # best q-value per precursor across its features
    best_q, best_s = {}, {}
    for p, q, s in zip(fs["precursor_id"].to_pylist(), fs["qvalue"].to_pylist(),
                       fs["score"].to_pylist()):
        if q is None:                      # a feature with no score is not evidence either way
            continue
        if p not in best_q or q < best_q[p]:
            best_q[p], best_s[p] = q, s

    odia = {p for p, q in best_q.items() if q < 0.01 and p in pid_meta}
    ref = {l.strip() for l in open(reflist) if l.strip()}
    ref_pid = {lib_key[k] for k in ref if k in lib_key}

    print(f"reference list      : {len(ref)} precursors")
    print(f"  present in library: {len(ref_pid)}  ({len(ref) - len(ref_pid)} absent -> not searchable)")
    print(f"ODIA q<0.01 targets : {len(odia)}")
    print(f"  overlap with ref  : {len(odia & ref_pid)}")
    print(f"  ODIA-only         : {len(odia - ref_pid)}")
    print(f"  reference-only    : {len(ref_pid - odia)}   <-- the gap to explain")
    print()

    buckets = collections.Counter()
    detail = collections.defaultdict(list)
    for p in ref_pid - odia:
        if p not in best_q:
            buckets["DETECTION: never extracted"] += 1
            detail["never"].append(p)
        else:
            q = best_q[p]
            b = ("FDR: 0.01 <= q < 0.05" if q < 0.05 else
                 "FDR: 0.05 <= q < 0.20" if q < 0.20 else
                 "DISCRIMINATION: 0.20 <= q < 0.50" if q < 0.50 else
                 "DISCRIMINATION: q >= 0.50")
            buckets[b] += 1
            detail[b].append((p, q, best_s.get(p)))
    tot = max(1, len(ref_pid - odia))
    for b, n in sorted(buckets.items(), key=lambda kv: -kv[1]):
        print(f"  {n:6d}  {100*n/tot:5.1f}%  {b}")

    # characterise the never-extracted ones -- are they odd in mz, charge or RT?
    if detail["never"]:
        mz = [pid_meta[p][0] for p in detail["never"] if p in pid_meta]
        ch = collections.Counter(pid_meta[p][1] for p in detail["never"] if p in pid_meta)
        rt = sorted(pid_meta[p][2] for p in detail["never"] if p in pid_meta)
        found_mz = sorted(pid_meta[p][0] for p in odia if p in pid_meta)
        if mz and found_mz:
            mz.sort()
            print()
            print(f"  never-extracted: mz median {mz[len(mz)//2]:.1f} "
                  f"(identified median {found_mz[len(found_mz)//2]:.1f}), "
                  f"charges {dict(ch.most_common(4))}, "
                  f"library_rt {rt[0]:.3f}..{rt[-1]:.3f}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
