#!/usr/bin/env python3
"""Falsification test for "the reference losses are weaker than the average candidate".

If reference precursors that PASS (depth>=4) are ALSO weaker than the average depth>=4 candidate,
then "weaker" is a property of REAL PEPTIDES in a pool dominated by chance matches -- not
something specific to the losses, and not evidence that intensity ranking is inverted for the
recovery problem.  If instead reference passers are STRONGER than average while reference losses
are WEAKER, the depth-3 pool is genuinely different and the original reading stands.
"""
import re, sys, statistics as st
dump, ref = sys.argv[1], sys.argv[2]
strip = lambda s: re.sub(r'\(UniMod:\d+\)', '', s).replace('.', '').replace('_', '')
want = set()
for line in open(ref):
    line = line.strip()
    if not line: continue
    seq, _, z = line.rpartition('_'); want.add((strip(seq), z))

best = {}
with open(dump) as fh:
    hdr = fh.readline().rstrip('\n').split('\t'); c = {n: i for i, n in enumerate(hdr)}
    for line in fh:
        f = line.rstrip('\n').split('\t')
        if f[c['decoy']] != '0': continue
        seq, _, z = f[c['id']].rpartition('_'); k = (strip(seq), z)
        depth = int(f[c['ms2_best_fragment_hits']])
        cur = best.get(k)
        if cur is None or depth > cur[0]:
            best[k] = (depth, int(f[c['ms2_qualifying_spectra']]),
                       float(f[c['ms2_sum_intensity']]), float(f[c['ms2_max_intensity']]))

def med(ks, i): 
    v = sorted(best[k][i] for k in ks)
    return v[len(v)//2] if v else float('nan')

print(f"{'stratum':<28} {'n':>9} {'recur':>9} {'sum_int':>13} {'max_int':>12}")
for name, lo, hi in [("depth == 3", 3, 3), ("depth == 4", 4, 4),
                     ("depth == 5", 5, 5), ("depth == 6", 6, 6)]:
    pool = [k for k in best if lo <= best[k][0] <= hi]
    refs = [k for k in pool if k in want]
    if not refs: continue
    print(f"  {name:<26} pool {len(pool):>9,} {med(pool,1):>8.0f} {med(pool,2):>13,.0f} {med(pool,3):>12,.0f}")
    print(f"  {'':<26} ref  {len(refs):>9,} {med(refs,1):>8.0f} {med(refs,2):>13,.0f} {med(refs,3):>12,.0f}")
    r = med(refs,2)/med(pool,2) if med(pool,2) else float('nan')
    print(f"  {'':<26} ratio(sum_int) ref/pool = {r:.2f}")
