#!/usr/bin/env python3
"""What sets the q<0.01 threshold -- the target signal, or the extreme decoy tail?

FDR at score t is roughly #decoys>t / #targets>t.  If only a few dozen decoys sit above the
operating threshold out of ~200k, then those few decoys ARE the threshold, and suppressing them
is far more leveraged than lifting thousands of targets.
"""
import sys, glob, zipfile, tempfile
import pyarrow.parquet as pq
bundle, reflist = sys.argv[1], sys.argv[2]
d = tempfile.mkdtemp(prefix="tail_")
with zipfile.ZipFile(bundle) as z:
    for n in z.namelist():
        if n.endswith(("features.parquet", "score_ms2.parquet", "precursors.parquet")):
            z.extract(n, d)
feat = pq.read_table(glob.glob(f"{d}/runs/*/features.parquet")[0],
                     columns=["feature_id","precursor_id"]).to_pydict()
sc = pq.read_table(glob.glob(f"{d}/runs/*/score_ms2.parquet")[0],
                   columns=["feature_id","qvalue","score"]).to_pydict()
prec = pq.read_table(f"{d}/library/precursors.parquet",
                     columns=["precursor_id","charge","decoy","unmodified_sequence"]).to_pydict()
meta = {p:(dc,sq,ch) for p,dc,sq,ch in zip(prec["precursor_id"],prec["decoy"],
                                           prec["unmodified_sequence"],prec["charge"])}
s_of = dict(zip(sc["feature_id"], sc["score"])); q_of = dict(zip(sc["feature_id"], sc["qvalue"]))
want=set()
for l in open(reflist):
    l=l.strip()
    if l: s,_,z2=l.rpartition('_'); want.add((s,z2))
best={}
for fid,pid in zip(feat["feature_id"],feat["precursor_id"]):
    s=s_of.get(fid)
    if s is None: continue
    if pid not in best or s>best[pid][0]: best[pid]=(s,q_of.get(fid,1.0))
T=[];D=[];MISS=[]
for pid,(s,q) in best.items():
    m=meta.get(pid)
    if not m: continue
    dc,sq,ch=m
    if dc: D.append(s)
    else:
        T.append(s)
        if (sq,str(ch)) in want and q>=0.01: MISS.append(s)
D.sort(); T.sort(); MISS.sort()
# operating threshold: smallest t where #D>t / #T>t <= 0.01
import bisect
def above(v,t): return len(v)-bisect.bisect_right(v,t)
lo,hi=min(D+T),max(D+T); thr=hi
for _ in range(60):
    mid=(lo+hi)/2
    fdr = above(D,mid)/max(1,above(T,mid))
    if fdr<=0.01: thr=mid; hi=mid
    else: lo=mid
print(f"operating threshold (FDR<=1%): score = {thr:.3f}")
print(f"  targets above: {above(T,thr):,}   decoys above: {above(D,thr):,}   of {len(D):,} decoys")
print(f"  => the threshold is set by {above(D,thr)} decoys, {100*above(D,thr)/len(D):.4f}% of the decoy pool")
print()
print(f"missed reference precursors: {len(MISS):,}")
for frac in (0.25,0.5,0.75,0.9):
    print(f"   p{int(frac*100):<3} score {MISS[int(frac*(len(MISS)-1))]:.3f}")
print()
print("if the decoy tail were suppressed, how many missed refs clear a LOWER threshold?")
for t in (thr, thr*0.9, thr*0.75, thr*0.5, 3.966):
    print(f"   threshold {t:>6.3f}: missed refs recovered {above(MISS,t):>5,}  "
          f"(decoys above: {above(D,t):>6,})")
