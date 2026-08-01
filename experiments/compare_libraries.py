#!/usr/bin/env python3
"""Transition-by-transition comparison: our predicted library vs a DIA-NN library.

Answers: where do we deviate? Do we miss transitions? Are we systematically off
in RT / IM / precursor m/z / fixed-mod fragment m/z? Is our MS2 INTENSITY
prediction the problem (the suspected cause of the target/decoy collapse)?

    compare_libraries.py <our_lib.tsv(OSW)> <diann_lib.parquet> [n_detail]

Join on (stripped sequence, charge). For shared precursors, match fragments by
(type, series, fragment-charge) and compare m/z + normalised intensity. Reports
aggregate distributions + a random detailed dump of n_detail precursors.
"""
import sys, csv, random, math
import pandas as pd

random.seed(1)  # deterministic subset (Math.random-free requirement is JS; here fine)

OUR = sys.argv[1]
DIANN = sys.argv[2]
N_DETAIL = int(sys.argv[3]) if len(sys.argv) > 3 else 8

# OSW TSV 0-indexed cols
PMZ, PRODMZ, PCH, PRODCH, INT, RT, SEQ, MODSEQ = 0, 1, 2, 3, 4, 5, 6, 7
FTYPE, FSER, IM, DECOY = 17, 18, 21, 24

# ---- DIA-NN library (measured/predicted) -> per-precursor fragment tables ----
d = pd.read_parquet(DIANN, columns=[
    "Stripped.Sequence", "Precursor.Charge", "Product.Mz", "Relative.Intensity",
    "Fragment.Type", "Fragment.Charge", "Fragment.Series.Number",
    "Fragment.Loss.Type", "RT", "IM", "Precursor.Mz", "Decoy"])
d = d[d["Decoy"] == 0]
# fragment key: type/series/charge, only b/y no loss (match our library's content)
d = d[d["Fragment.Loss.Type"].isin(["noloss", "None", ""]) | d["Fragment.Loss.Type"].isna()]
dn = {}
for (seq, ch), g in d.groupby(["Stripped.Sequence", "Precursor.Charge"]):
    frags = {}
    for _, r in g.iterrows():
        frags[(str(r["Fragment.Type"]).lower()[:1], int(r["Fragment.Series.Number"]), int(r["Fragment.Charge"]))] = \
            (float(r["Product.Mz"]), float(r["Relative.Intensity"]))
    dn[(seq, int(ch))] = dict(rt=float(g["RT"].iloc[0]), im=float(g["IM"].iloc[0]),
                              pmz=float(g["Precursor.Mz"].iloc[0]), frags=frags)
print(f"DIA-NN precursors: {len(dn):,}")

# ---- our library (stream; keep only shared precursors) ----
ours = {}
with open(OUR) as f:
    r = csv.reader(f, delimiter="\t"); next(r)
    for c in r:
        if len(c) <= DECOY or c[DECOY] != "0":
            continue
        key = (c[SEQ], int(c[PCH]))
        if key not in dn:
            continue
        o = ours.setdefault(key, dict(rt=float(c[RT]), im=float(c[IM]),
                                      pmz=float(c[PMZ]), frags={}))
        o["frags"][(c[FTYPE].lower()[:1], int(c[FSER]), int(c[PRODCH]))] = \
            (float(c[PRODMZ]), float(c[INT]))
shared = list(ours.keys())
print(f"our precursors total-in-shared: {len(shared):,}  (DIA-NN {len(dn):,})")

# ---- aggregate stats over all shared precursors ----
import statistics as st
def svec(fr, keys):
    v = [fr[k][1] for k in keys]; n = math.sqrt(sum(x*x for x in v)) or 1.0
    return [x/n for x in v]

n_ours_f, n_dn_f, n_match_f, n_miss_f, n_extra_f = [], [], [], [], []
mz_absdiff, spectral_angle, pearson = [], [], []
rt_pairs, im_pairs, pmz_absdiff = [], [], []

for key in shared:
    o, dd = ours[key], dn[key]
    of, df = o["frags"], dd["frags"]
    ok, dk = set(of), set(df)
    both = ok & dk
    n_ours_f.append(len(ok)); n_dn_f.append(len(dk))
    n_match_f.append(len(both)); n_miss_f.append(len(dk - ok)); n_extra_f.append(len(ok - dk))
    for k in both:
        mz_absdiff.append(abs(of[k][0] - df[k][0]))
    if len(both) >= 3:
        bk = sorted(both)
        ov, dv = svec(of, bk), svec(df, bk)
        cos = sum(a*b for a, b in zip(ov, dv))
        spectral_angle.append(1 - 2*math.acos(min(1.0, cos))/math.pi)  # 1=identical
        mo, md = st.mean(ov), st.mean(dv)
        num = sum((a-mo)*(b-md) for a, b in zip(ov, dv))
        den = math.sqrt(sum((a-mo)**2 for a in ov)*sum((b-md)**2 for b in dv)) or 1.0
        pearson.append(num/den)
    rt_pairs.append((o["rt"], dd["rt"])); im_pairs.append((o["im"], dd["im"]))
    pmz_absdiff.append(abs(o["pmz"] - dd["pmz"]))

def pct(x, p):
    xs = sorted(x); return xs[min(len(xs)-1, int(p*len(xs)))] if xs else float("nan")
print("\n=== FRAGMENT COVERAGE (per shared precursor) ===")
print(f"  our frags   median {pct(n_ours_f,.5):.0f}   DIA-NN median {pct(n_dn_f,.5):.0f}")
print(f"  matched     median {pct(n_match_f,.5):.0f}")
print(f"  WE MISS (DIA-NN has, we don't)  median {pct(n_miss_f,.5):.0f}  mean {st.mean(n_miss_f):.2f}")
print(f"  we extra    median {pct(n_extra_f,.5):.0f}")
print("\n=== FRAGMENT m/z (shared frags; tests chemistry + fixed mods) ===")
print(f"  |dmz| median {pct(mz_absdiff,.5):.5f}  p95 {pct(mz_absdiff,.95):.5f}  max {max(mz_absdiff):.4f}")
print("\n=== MS2 INTENSITY agreement (the suspected blocker) ===")
print(f"  spectral-angle-sim  median {pct(spectral_angle,.5):.3f}  (1=identical, ~0=orthogonal)  n={len(spectral_angle):,}")
print(f"  Pearson r           median {pct(pearson,.5):.3f}")
print(f"  frac precursors r<0.5  {sum(1 for x in pearson if x<0.5)/len(pearson):.1%}")
print("\n=== PRECURSOR m/z ===")
print(f"  |dmz| median {pct(pmz_absdiff,.5):.5f}  p95 {pct(pmz_absdiff,.95):.5f}")
print("\n=== RT relationship (our NormalizedRT vs DIA-NN RT) ===")
ox = [a for a,_ in rt_pairs]; oy = [b for _,b in rt_pairs]
mx, my = st.mean(ox), st.mean(oy)
num = sum((a-mx)*(b-my) for a,b in rt_pairs); den = math.sqrt(sum((a-mx)**2 for a in ox)*sum((b-my)**2 for b in oy)) or 1
slope = num/(sum((a-mx)**2 for a in ox) or 1)
print(f"  Pearson r {num/den:.4f}   fit DIANN_rt = {slope:.2f}*ourRT + {my-slope*mx:.2f}")
print(f"  our RT range {min(ox):.3f}..{max(ox):.3f}   DIA-NN RT range {min(oy):.1f}..{max(oy):.1f}")
ix = [a for a,_ in im_pairs]; iy = [b for _,b in im_pairs]
mix, miy = st.mean(ix), st.mean(iy)
inum = sum((a-mix)*(b-miy) for a,b in im_pairs); iden = math.sqrt(sum((a-mix)**2 for a in ix)*sum((b-miy)**2 for b in iy)) or 1
print(f"\n=== IM relationship ===\n  Pearson r {inum/iden:.4f}   our range {min(ix):.3f}..{max(ix):.3f}  DIA-NN {min(iy):.3f}..{max(iy):.3f}")

# ---- detailed transition dump for a random subset ----
print(f"\n=== DETAILED transition-by-transition ({N_DETAIL} random shared precursors) ===")
for key in random.sample(shared, min(N_DETAIL, len(shared))):
    o, dd = ours[key], dn[key]
    print(f"\n{key[0]}/{key[1]}+  ourRT={o['rt']:.3f} diannRT={dd['rt']:.1f}  ourIM={o['im']:.3f} diannIM={dd['im']:.3f}  ourPmz={o['pmz']:.4f} diannPmz={dd['pmz']:.4f}")
    allk = sorted(set(o["frags"]) | set(dd["frags"]))
    print("   frag        our_mz   dn_mz    our_int  dn_int")
    for k in allk:
        om, oi = o["frags"].get(k, (float("nan"), float("nan")))
        dm, di = dd["frags"].get(k, (float("nan"), float("nan")))
        tag = "" if (k in o["frags"] and k in dd["frags"]) else ("  <-WE MISS" if k not in o["frags"] else "  <-we-extra")
        print(f"   {k[0]}{k[1]}^{k[2]}     {om:8.3f} {dm:8.3f}  {oi:7.4f} {di:7.4f}{tag}")
