#!/usr/bin/env python3
"""Population at each hit depth, targets vs decoys.

The go/no-go said 80% of the recoverable losses sit at depth 3.  This asks the harder question:
how CROWDED is depth 3?  We must promote ~766 depth-3 targets without growing the budget, so the
achievable gain is bounded by how well anything can separate them from the depth-3 population as
a whole.  Enrichment = targets/decoys at a depth; 1.0 means the current criterion cannot tell a
real peptide from a shuffled one there, which is precisely where new evidence has to do the work.
"""
import sys, collections
dump = sys.argv[1]
t = collections.Counter(); d = collections.Counter()
with open(dump) as fh:
    hdr = fh.readline().rstrip('\n').split('\t')
    c = {n: i for i, n in enumerate(hdr)}
    for line in fh:
        f = line.rstrip('\n').split('\t')
        depth = int(f[c['ms2_best_fragment_hits']])
        (d if f[c['decoy']] == '1' else t)[depth] += 1
print(f"{'depth':>6} {'targets':>12} {'decoys':>12} {'enrich':>8}   {'cum targets >= depth':>20}")
alld = sorted(set(t) | set(d))
for dep in alld:
    ct = sum(t[x] for x in alld if x >= dep)
    e = t[dep] / d[dep] if d[dep] else float('inf')
    print(f"{dep:>6} {t[dep]:>12,} {d[dep]:>12,} {e:>8.3f} {ct:>20,}")
print()
tt, td = sum(t.values()), sum(d.values())
print(f"total targets {tt:,}  decoys {td:,}")
p4t = sum(t[x] for x in alld if x >= 4); p4d = sum(d[x] for x in alld if x >= 4)
p3t = sum(t[x] for x in alld if x >= 3); p3d = sum(d[x] for x in alld if x >= 3)
print(f"depth>=4 (today's rule): targets {p4t:,}  decoys {p4d:,}  enrich {p4t/max(1,p4d):.3f}")
print(f"depth>=3               : targets {p3t:,}  decoys {p3d:,}  enrich {p3t/max(1,p3d):.3f}")
print()
print(f"To keep the budget fixed, admitting depth-3 targets means dropping the same number from")
print(f"depth>=4.  Depth-3 pool to choose 766 from: {t[3]+d[3]:,} candidates "
      f"({t[3]:,} target / {d[3]:,} decoy).")
