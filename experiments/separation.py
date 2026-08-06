#!/usr/bin/env python3
"""How well does the classifier separate targets from decoys, and where do the missed refs sit?

The gap is entirely scoring-side, so the question is whether the SCORE has the information and
the threshold is just tight (near-misses cluster near the cutoff, decoys well below), or whether
the score does not carry the information at all (missed refs sit inside the decoy distribution).
Those imply completely different fixes: a better-calibrated FDR vs better evidence.
"""
import sys, glob, zipfile, tempfile, collections, statistics as st
import pyarrow.parquet as pq

bundle, reflist = sys.argv[1], sys.argv[2]
d = tempfile.mkdtemp(prefix="sep_")
with zipfile.ZipFile(bundle) as z:
    for n in z.namelist():
        if n.endswith(("features.parquet", "score_ms2.parquet", "precursors.parquet")):
            z.extract(n, d)
feat = pq.read_table(glob.glob(f"{d}/runs/*/features.parquet")[0],
                     columns=["feature_id", "precursor_id"]).to_pydict()
sc = pq.read_table(glob.glob(f"{d}/runs/*/score_ms2.parquet")[0],
                   columns=["feature_id", "qvalue", "score"]).to_pydict()
prec = pq.read_table(f"{d}/library/precursors.parquet",
                     columns=["precursor_id", "charge", "decoy", "unmodified_sequence"]).to_pydict()

want = set()
for line in open(reflist):
    line = line.strip()
    if not line: continue
    s, _, z2 = line.rpartition('_'); want.add((s, z2))

pid_meta = {p: (dc, sq, ch) for p, dc, sq, ch in
            zip(prec["precursor_id"], prec["decoy"], prec["unmodified_sequence"], prec["charge"])}
score_of = dict(zip(sc["feature_id"], sc["score"]))
q_of = dict(zip(sc["feature_id"], sc["qvalue"]))

# best score per precursor
best = {}
for fid, pid in zip(feat["feature_id"], feat["precursor_id"]):
    s = score_of.get(fid)
    if s is None: continue
    if pid not in best or s > best[pid][0]: best[pid] = (s, q_of.get(fid, 1.0))

tgt, dec, ref_hit, ref_miss = [], [], [], []
for pid, (s, q) in best.items():
    m = pid_meta.get(pid)
    if m is None: continue
    dc, sq, ch = m
    if dc: dec.append(s); continue
    tgt.append(s)
    if (sq, str(ch)) in want:
        (ref_hit if q < 0.01 else ref_miss).append(s)

def qs(v, p): 
    v = sorted(v); return v[int(p*(len(v)-1))] if v else float('nan')
print(f"precursor-level best scores:  targets {len(tgt):,}  decoys {len(dec):,}")
print(f"  reference FOUND (q<0.01): {len(ref_hit):,}     reference MISSED: {len(ref_miss):,}")
print()
print(f"{'stratum':<22} {'n':>8} {'p10':>9} {'median':>9} {'p90':>9}")
for nm, v in [("decoys", dec), ("all targets", tgt),
              ("reference FOUND", ref_hit), ("reference MISSED", ref_miss)]:
    print(f"  {nm:<20} {len(v):>8,} {qs(v,.10):>9.3f} {qs(v,.50):>9.3f} {qs(v,.90):>9.3f}")
print()
d25, d50, d75, d99 = qs(dec,.25), qs(dec,.50), qs(dec,.75), qs(dec,.99)
print(f"decoy IQR = [{d25:.3f}, {d75:.3f}], decoy p99 = {d99:.3f}")
for nm, v in [("all targets", tgt), ("reference MISSED", ref_miss), ("reference FOUND", ref_hit)]:
    inside = sum(1 for s in v if d25 <= s <= d75)
    above  = sum(1 for s in v if s > d99)
    print(f"  {nm:<20} inside decoy IQR {100*inside/len(v):>5.1f}%   above decoy p99 {100*above/len(v):>5.1f}%")
