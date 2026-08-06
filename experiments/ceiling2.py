#!/usr/bin/env python3
"""Ceiling estimate, at the PRECURSOR level.

Correcting a flaw in the first pass: it ranked depth-3 ROWS and counted any reference precursor
with a depth-3 row.  But a precursor SURVIVES if any of its modified forms passes, so a precursor
whose best form is depth 5 is not a loss even though it also has a depth-3 row.  Aggregate to the
precursor (max depth over its forms) first, then ask about the ones whose BEST is depth 3.
"""
import re, sys, collections
dump, ref = sys.argv[1], sys.argv[2]
strip = lambda s: re.sub(r'\(UniMod:\d+\)', '', s).replace('.', '').replace('_', '')

want = set()
for line in open(ref):
    line = line.strip()
    if not line: continue
    seq, _, z = line.rpartition('_')
    want.add((strip(seq), z))

# aggregate to precursor: best depth, and the evidence of the row achieving it
best = {}
with open(dump) as fh:
    hdr = fh.readline().rstrip('\n').split('\t'); c = {n: i for i, n in enumerate(hdr)}
    for line in fh:
        f = line.rstrip('\n').split('\t')
        if f[c['decoy']] != '0': continue
        seq, _, z = f[c['id']].rpartition('_')
        k = (strip(seq), z)
        depth = int(f[c['ms2_best_fragment_hits']])
        cur = best.get(k)
        if cur is None or depth > cur[0]:
            best[k] = (depth, int(f[c['ms2_qualifying_spectra']]), int(f[c['ms2_hit_count']]),
                       float(f[c['ms2_sum_intensity']]), float(f[c['ms2_max_intensity']]),
                       float(f[c['ms1_max_intensity']]))

d3 = {k: v for k, v in best.items() if v[0] == 3}
refs3 = {k for k in d3 if k in want}
print(f"target precursors, best depth == 3: {len(d3):,}")
print(f"   of which reference losses:       {len(refs3)}")
BUDGET = 91779
print(f"budget to reallocate: {BUDGET:,} ({100*BUDGET/len(d3):.1f}% of the pool)")
base = len(refs3)*BUDGET/len(d3)
print(f"random-selection expectation: {base:.0f}\n")

keys = list(d3)
def evaluate(name, idx):
    order = sorted(keys, key=lambda k: d3[k][idx], reverse=True)
    got = sum(1 for k in order[:BUDGET] if k in refs3)
    print(f"  {name:<34} captured {got:>4} / {len(refs3)}  ({100*got/len(refs3):>5.1f}%)   "
          f"{got/base:>5.2f}x random")
print("ranking the best-depth-3 precursors by evidence already in the dump:")
evaluate("qualifying spectra (recurrence)", 1)
evaluate("ms2 hit count",                   2)
evaluate("ms2 summed intensity",            3)
evaluate("ms2 max intensity",               4)
evaluate("ms1 max intensity",               5)

# Distribution comparison: are the reference losses systematically weak or strong?
import statistics as st
for nm, idx in [("qualifying spectra", 1), ("ms2 sum intensity", 3), ("ms2 max intensity", 4)]:
    r = sorted(d3[k][idx] for k in refs3)
    a = sorted(d3[k][idx] for k in keys)
    med_r = r[len(r)//2]; med_a = a[len(a)//2]
    print(f"    {nm:<20} median: reference {med_r:>12.1f}   all depth-3 {med_a:>12.1f}")
