#!/usr/bin/env python3
"""Are the threshold-setting decoys sitting on top of real peptides?

70 decoys set the 1% cutoff.  If a high-scoring decoy is an INTERFERENCE artefact, it should
co-locate with a confidently identified target: same isolation window (similar precursor m/z) and
overlapping elution.  A shuffled decoy has different fragment m/z from its own target, but it can
still collect signal from whatever real peptide happens to co-elute in the same window.

Test: for each top decoy, is there a confident target (q<0.01) within +-dRT and in the same
isolation window?  Compare against decoys drawn from the BULK of the decoy distribution -- if the
top decoys co-locate more often than bulk decoys, interference is implicated.
"""
import sys, glob, zipfile, tempfile, bisect, random
import pyarrow.parquet as pq
bundle = sys.argv[1]
d = tempfile.mkdtemp(prefix="td_")
with zipfile.ZipFile(bundle) as z:
    for n in z.namelist():
        if n.endswith(("features.parquet","score_ms2.parquet","precursors.parquet")):
            z.extract(n, d)
feat = pq.read_table(glob.glob(f"{d}/runs/*/features.parquet")[0],
                     columns=["feature_id","precursor_id","exp_rt"]).to_pydict()
sc = pq.read_table(glob.glob(f"{d}/runs/*/score_ms2.parquet")[0],
                   columns=["feature_id","qvalue","score"]).to_pydict()
prec = pq.read_table(f"{d}/library/precursors.parquet",
                     columns=["precursor_id","precursor_mz","decoy"]).to_pydict()
mz_of = dict(zip(prec["precursor_id"], prec["precursor_mz"]))
dc_of = dict(zip(prec["precursor_id"], prec["decoy"]))
s_of = dict(zip(sc["feature_id"], sc["score"])); q_of = dict(zip(sc["feature_id"], sc["qvalue"]))

best = {}   # pid -> (score, q, rt)
for fid, pid, rt in zip(feat["feature_id"], feat["precursor_id"], feat["exp_rt"]):
    s = s_of.get(fid)
    if s is None: continue
    if pid not in best or s > best[pid][0]: best[pid] = (s, q_of.get(fid,1.0), rt)

conf = sorted((rt, mz_of.get(pid, -1)) for pid,(s,q,rt) in best.items()
              if not dc_of.get(pid,0) and q < 0.01)
conf_rt = [c[0] for c in conf]
decs = sorted(((s,rt,mz_of.get(pid,-1)) for pid,(s,q,rt) in best.items() if dc_of.get(pid,0)),
              reverse=True)
print(f"confident targets: {len(conf):,}   decoys: {len(decs):,}")

def colocated(rt, mz, dRT=20.0, dMZ=12.5):
    lo = bisect.bisect_left(conf_rt, rt-dRT); hi = bisect.bisect_right(conf_rt, rt+dRT)
    return sum(1 for k in range(lo,hi) if abs(conf[k][1]-mz) <= dMZ)

random.seed(5)
groups = [("top 70 decoys (set the cutoff)", decs[:70]),
          ("decoys ranked 71-500",           decs[70:500]),
          ("decoys ranked 501-2000",         decs[500:2000]),
          ("random 2000 from the bulk",      random.sample(decs[5000:], 2000))]
print(f"\n{'group':<34} {'n':>6} {'median co-located':>18} {'frac with >=1':>14}")
for nm, g in groups:
    counts = [colocated(rt,mz) for _,rt,mz in g]
    counts.sort()
    med = counts[len(counts)//2] if counts else 0
    frac = 100*sum(1 for c in counts if c>=1)/len(counts) if counts else 0
    print(f"  {nm:<32} {len(g):>6,} {med:>18} {frac:>13.1f}%")
