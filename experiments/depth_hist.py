#!/usr/bin/env python3
"""Depth histogram of the reference precursors, from a threshold-1 evidence dump.

The question this answers is the GO/NO-GO for an evidence-richer prefilter: of the reference
precursors the DEFAULT rule (>=4 of top-6 in one spectrum) deletes, how many actually had
fragment evidence?  A dot product and a mass-error spread are undefined without matched
fragments, so losses at depth 0-1 cannot be rescued by better scoring of the evidence -- they
would need a different mechanism entirely.

Keyed on STRIPPED sequence + charge: the dump carries modified sequences and the reference list
does not, and matching on sequence alone conflates charge states (pf_lost.py records that bug).
Several modified forms can map to one reference key, so a key takes the MAX depth over its rows
-- the best chance that precursor had.
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

best, best_spec, best_ms1 = {}, {}, {}
with open(dump) as fh:
    hdr = fh.readline().rstrip('\n').split('\t')
    c = {n: i for i, n in enumerate(hdr)}
    for line in fh:
        f = line.rstrip('\n').split('\t')
        if f[c['decoy']] != '0': continue                      # targets only
        seq, _, z = f[c['id']].rpartition('_')
        k = (strip(seq), z)
        if k not in want: continue
        d = int(f[c['ms2_best_fragment_hits']])
        if d > best.get(k, -1):
            best[k] = d
            best_spec[k] = int(f[c['ms2_qualifying_spectra']])
            best_ms1[k] = int(f[c['ms1_hit_count']])

print(f"reference precursors:            {len(want)}")
print(f"  present in the dump as targets: {len(best)}")
print(f"  absent from the dump entirely:  {len(want) - len(best)}   (never a candidate at all)")
print()
hist = collections.Counter(best.values())
tot = len(best)
print("depth = ms2_best_fragment_hits (best over modified forms), TARGETS in the reference set")
print(f"{'depth':>6} {'n':>7} {'share':>8}   {'cumulative':>10}")
cum = 0
for d in sorted(hist):
    cum += hist[d]
    print(f"{d:>6} {hist[d]:>7} {100*hist[d]/tot:>7.1f}% {100*cum/tot:>10.1f}%")
print()
lost = {k: v for k, v in best.items() if v < 4}
print(f"DELETED by the default rule (depth < 4): {len(lost)}")
if lost:
    lh = collections.Counter(lost.values())
    for d in sorted(lh):
        print(f"   depth {d}: {lh[d]:>5}  ({100*lh[d]/len(lost):.1f}% of the losses)")
    rescuable = sum(n for d, n in lh.items() if d >= 2)
    print()
    print(f"   depth >= 2 (a dot product and mass-error spread are DEFINED): "
          f"{rescuable}  ({100*rescuable/len(lost):.1f}%)")
    print(f"   depth <  2 (no usable fragment evidence):                     "
          f"{len(lost)-rescuable}  ({100*(len(lost)-rescuable)/len(lost):.1f}%)")
    ms1 = sum(1 for k in lost if best_ms1[k] > 0)
    print(f"   of the losses, with MS1 evidence despite failing MS2:         {ms1}")
# the single-spectrum axis, separately
deep_1spec = sum(1 for k, v in best.items() if v >= 4 and best_spec[k] <= 1)
print()
print(f"depth >= 4 but only ONE qualifying spectrum: {deep_1spec}")
