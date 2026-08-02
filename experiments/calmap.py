#!/usr/bin/env python3
"""Reference calibration map from a run's own confident identifications.

Builds library_rt -> exp_rt from q<0.01 rank-1 target features, then reports the residual
around that map, globally and binned by observed RT.

The residuals are CENSORED at the extraction window: peptides whose true RT deviation exceeds
the window were never extracted, so they cannot appear. A map built at a 600 s window reports a
p99 that says nothing about what a 1435 s window would have found. Comparing two runs at
different windows is how that bias is measured rather than assumed -- see --compare.

usage: calmap.py <run.oswpq> [more.oswpq ...]
"""
import sys, os, json, glob, zipfile, tempfile, statistics as st

try:
    import pyarrow.parquet as pq, pyarrow.compute as pc
except ImportError:
    sys.exit("needs pyarrow (use /scratch/kohlbach/mzpenv/bin/python)")


def load(bundle):
    d = tempfile.mkdtemp(prefix="calmap_")
    with zipfile.ZipFile(bundle) as z:
        for n in z.namelist():
            if n.endswith(("features.parquet", "score_ms2.parquet", "precursors.parquet")):
                z.extract(n, d)
    feat = pq.read_table(glob.glob(f"{d}/runs/*/features.parquet")[0],
                         columns=["feature_id", "precursor_id", "exp_rt", "delta_rt"])
    sc = pq.read_table(glob.glob(f"{d}/runs/*/score_ms2.parquet")[0],
                       columns=["feature_id", "qvalue", "rank"])
    prec = pq.read_table(f"{d}/library/precursors.parquet",
                         columns=["precursor_id", "library_rt", "decoy"])
    t = feat.join(sc, keys="feature_id").join(prec, keys="precursor_id")
    m = pc.and_(pc.and_(pc.less(t["qvalue"], 0.01), pc.equal(t["rank"], 1)),
                pc.invert(pc.cast(t["decoy"], "bool")))
    return t.filter(m)


def refmap(pairs, nknots=40):
    """Binned-median curve through (library_rt, exp_rt). Equal-count bins, so knot density
    follows the data rather than the axis."""
    pairs = sorted(pairs)
    n, ref = len(pairs), []
    for b in range(nknots):
        ch = pairs[b * n // nknots:(b + 1) * n // nknots]
        if ch:
            ref.append((st.median([p[0] for p in ch]), st.median([p[1] for p in ch])))
    return ref


def interp(ref, x):
    if x <= ref[0][0]:
        return ref[0][1]
    if x >= ref[-1][0]:
        return ref[-1][1]
    for i in range(len(ref) - 1):
        if ref[i][0] <= x <= ref[i + 1][0]:
            f = (x - ref[i][0]) / max(1e-12, ref[i + 1][0] - ref[i][0])
            return ref[i][1] + f * (ref[i + 1][1] - ref[i][1])
    return ref[-1][1]


def qs(v, *ps):
    v = sorted(v)
    return [v[min(len(v) - 1, int(p * len(v)))] for p in ps]


def report(bundle):
    r = load(bundle)
    lib = r["library_rt"].to_pylist()
    exp = r["exp_rt"].to_pylist()
    drt = [abs(x) for x in r["delta_rt"].to_pylist()]
    ref = refmap(list(zip(lib, exp)))
    res = [abs(e - interp(ref, l)) for l, e in zip(lib, exp)]
    m, p95, p99 = qs(res, .5, .95, .99)
    dm, d95, d99 = qs(drt, .5, .95, .99)
    name = os.path.basename(bundle).replace(".oswpq", "")
    print(f"{name:20s} n={len(lib):6d}  |res| med {m:6.1f} p95 {p95:6.1f} p99 {p99:6.1f}"
          f"   |delta_rt| med {dm:6.1f} p95 {d95:6.1f} p99 {d99:6.1f}")

    lo, hi, nb = min(exp), max(exp), 8
    bins = [[] for _ in range(nb)]
    for e, d in zip(exp, drt):
        bins[min(nb - 1, int((e - lo) / (hi - lo) * nb))].append(d)
    prof = [round(qs(b, .5)[0], 1) if len(b) >= 5 else None for b in bins]
    print(f"{'':20s} |delta_rt| median by RT octile: {prof}")
    return {"name": name, "n": len(lib), "res_p99": p99, "drt_med": dm, "profile": prof,
            "exp_min": lo, "exp_max": hi}


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    out = [report(b) for b in sys.argv[1:]]
    if len(out) > 1:
        print()
        print("censoring check -- a map built at a NARROW window cannot see what a WIDE one finds:")
        base = max(out, key=lambda o: o["n"])
        for o in out:
            if o is base:
                continue
            print(f"  {o['name']} has {base['n'] - o['n']:+d} fewer IDs than {base['name']}; "
                  f"its p99 residual {o['res_p99']:.1f} s vs {base['res_p99']:.1f} s")
    json.dump(out, open("calmap_summary.json", "w"), indent=1)
