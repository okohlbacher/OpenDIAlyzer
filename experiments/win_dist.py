#!/usr/bin/env python3
"""RT-error and peak-width vs RT, restricted to KNOWN-REAL precursors.

v1 was unusable for the error distribution: |delta_rt| p99 came out at 427.4 s in every bin against
a pass-2 half-width of 432 s, i.e. the distribution is CENSORED by the extraction window and
"2*p99 + width" simply reproduces the window it was measured in. It also pooled ~4.9 candidate
peaks per precursor, nearly all noise.

Fix: restrict to precursors DIA-NN identified (absent from this pipeline's own scoring, so not
circular), and take each precursor's peak CLOSEST to its prediction as the putative real one. Bin by
EXPERIMENTAL RT in seconds, which is the axis a window is actually defined on.
"""
import io, sys, zipfile
import pyarrow.parquet as pq

bundle, ref_path = sys.argv[1], sys.argv[2]
ref = set(open(ref_path).read().split())
z = zipfile.ZipFile(bundle)
fp = [n for n in z.namelist() if n.endswith("features.parquet")][0]

p = pq.read_table(io.BytesIO(z.read("library/precursors.parquet")),
                  columns=["precursor_id","library_rt","decoy","unmodified_sequence","charge"]).to_pydict()
keep, librt = {}, {}
for pid, lrt, dec, seq, ch in zip(p["precursor_id"], p["library_rt"], p["decoy"],
                                  p["unmodified_sequence"], p["charge"]):
    if dec: continue
    librt[pid] = lrt
    if f"{seq}_{ch}" in ref: keep[pid] = lrt
print(f"library targets {len(librt):,}; of those DIA-NN-confirmed {len(keep):,}")

f = pq.read_table(io.BytesIO(z.read(fp)),
                  columns=["precursor_id","exp_rt","delta_rt","left_width","right_width"]).to_pydict()

best = {}
for pid, ert, drt, lw, rw in zip(f["precursor_id"], f["exp_rt"], f["delta_rt"],
                                 f["left_width"], f["right_width"]):
    if pid not in keep or ert is None or drt is None: continue
    cur = best.get(pid)
    if cur is None or abs(drt) < abs(cur[1]): best[pid] = (ert, drt, lw, rw)
print(f"confirmed precursors with a feature: {len(best):,}")

rows = sorted(((ert, drt, (rw-lw) if lw is not None and rw is not None else None)
               for ert, drt, lw, rw in best.values()), key=lambda r: r[0])
def pct(v,q):
    if not v: return float("nan")
    v=sorted(v); return v[min(len(v)-1,int(q*len(v)))]

NB, n = 10, len(rows)
print()
print(f"{'exp RT (s)':>18} {'n':>7} {'|dRT| med':>10} {'p95':>8} {'p99':>8} {'max':>8} "
      f"{'width med':>10} {'p95':>8}")
print("-"*84)
for b in range(NB):
    c = rows[n*b//NB : n*(b+1)//NB]
    if not c: continue
    d=[abs(r[1]) for r in c]; w=[r[2] for r in c if r[2] and r[2]>0]
    print(f"{c[0][0]:8.0f}-{c[-1][0]:<9.0f} {len(c):>7,} {pct(d,.5):>10.1f} {pct(d,.95):>8.1f} "
          f"{pct(d,.99):>8.1f} {max(d):>8.1f} {pct(w,.5):>10.1f} {pct(w,.95):>8.1f}")
d=[abs(r[1]) for r in rows]; w=[r[2] for r in rows if r[2] and r[2]>0]
print("-"*84)
print(f"{'ALL':>18} {len(rows):>7,} {pct(d,.5):>10.1f} {pct(d,.95):>8.1f} {pct(d,.99):>8.1f} "
      f"{max(d):>8.1f} {pct(w,.5):>10.1f} {pct(w,.95):>8.1f}")
print()
half = 432.0
print(f"censoring check: {100*sum(1 for x in d if x > half-10)/len(d):.2f}% of |dRT| sit within "
      f"10 s of the {half:.0f} s half-window (high % => still censored, distribution not trustworthy)")
print(f"suggested window = 2*p99(|dRT|) + p95(width) = {2*pct(d,.99)+pct(w,.95):.0f} s  (current 864)")
