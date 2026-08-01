#!/usr/bin/env python3
"""Keep 1/N of transition groups (target + its paired decoy) for a fast benchmark.
Pairs are kept together by hashing the DECOY_-stripped TransitionGroupId, plus the
CiRT anchors (needed for RT calibration). Usage: subset_library.py <in> <out> [N]"""
import sys, csv, zlib
IN, OUT = sys.argv[1], sys.argv[2]
N = int(sys.argv[3]) if len(sys.argv) > 3 else 8
TG, SEQ, DECOY = 22, 6, 24
CIRT = "/home/kohlbach/openms3/share/OpenMS/CHEMISTRY/cirtkit.tsv"
anchors = set()
try:
    with open(CIRT) as f:
        r = csv.reader(f, delimiter="\t"); next(r)
        for c in r:
            if len(c) > 6: anchors.add(c[6])
except Exception: pass

kept_groups = 0
with open(IN) as fi, open(OUT, "w") as fo:
    r = csv.reader(fi, delimiter="\t"); w = csv.writer(fo, delimiter="\t", lineterminator="\n")
    w.writerow(next(r))
    keep = {}
    for c in r:
        if len(c) <= DECOY: continue
        base = c[TG].replace("DECOY_", "")
        if base not in keep:
            h = zlib.crc32(base.encode()) % N
            keep[base] = (h == 0) or (c[SEQ] in anchors)
            if keep[base]: kept_groups += 1
        if keep[base]:
            w.writerow(c)
print(f"kept {kept_groups:,} groups (1/{N} + CiRT anchors)", file=sys.stderr)
