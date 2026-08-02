#!/usr/bin/env python3
"""Which prefilter rule removes the reference IDs we lose?

Takes the per-candidate evidence dump (-prefilter_out) and a reference ID list, isolates the
reference precursors that did NOT survive, and attributes each to the rule that removed it.

The rules, with their defaults:

    ms2_top_transitions_per_precursor = 6     index only the top-6 fragments by PREDICTED intensity
    ms2_min_fragment_hits             = 4     require 4 of those 6...
    ms2_top_peaks_per_spectrum        = 1000  ...within the top 1000 peaks...
    ms2_min_qualifying_spectra        = 1     ...of a SINGLE spectrum
    prefilter_mz_extraction_window    = 10    at 10 ppm

so a candidate survives iff ms2_best_fragment_hits >= 4 in at least one spectrum.

Attribution matters because the two obvious relaxations both made identifications WORSE
(6930 -> 5580 at 3-of-6; 6567 at top-3000 peaks), so a blanket loosening is not the answer. What
distinguishes the recoverable cases:

  hits == 3        : just under. But 3-of-6 was tested and lost 1350 IDs, so keeping THESE
                     specifically would need something narrower than lowering the threshold.
  hits in 1..2     : far under on fragment evidence.
  hits == 0        : no fragment evidence at all -- the prefilter is right and the reference tool
                     is finding them by some other route.
  hits >= 4 but    : the SINGLE-SPECTRUM requirement is what removed them, not the hit count. This
  not survived       is a different axis and relaxing it does not admit the same noise.
  ms1_hit_count>0  : MS1 evidence exists where MS2 does not -- a rescue route that does not touch
  with low MS2       the MS2 noise floor at all.

usage: pf_lost.py <pf_evidence.tsv> <reference_ids.txt>
"""
import sys, collections


def main(ev_path, ref_path):
    ref = {l.strip() for l in open(ref_path) if l.strip()}
    print(f"reference list: {len(ref)} precursors")

    rows = []
    with open(ev_path) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        idx = {n: i for i, n in enumerate(head)}
        need = ["sequence", "decoy", "survived", "ms2_best_fragment_hits", "ms2_hit_count",
                "ms2_qualifying_spectra", "ms1_hit_count"]
        missing = [n for n in need if n not in idx]
        if missing:
            sys.exit(f"evidence dump lacks columns {missing}; rebuild with the extended -prefilter_out")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < len(head):
                continue
            rows.append(f)
    print(f"evidence rows: {len(rows)}")

    # The dump keys on compound id; the reference list is SEQUENCE_CHARGE. The evidence carries the
    # sequence but not the charge, so match on sequence and report the ambiguity rather than hiding
    # it -- a sequence with two charge states counts once here.
    ref_seq = {k.rsplit("_", 1)[0] for k in ref}
    print(f"reference distinct sequences: {len(ref_seq)}")

    lost, kept, buckets = 0, 0, collections.Counter()
    ms1_rescue, spectra_rule, examples = 0, 0, []
    for f in rows:
        if f[idx["decoy"]] == "1":
            continue
        if f[idx["sequence"]] not in ref_seq:
            continue
        if f[idx["survived"]] == "1":
            kept += 1
            continue
        lost += 1
        hits = int(f[idx["ms2_best_fragment_hits"]] or 0)
        tot = int(f[idx["ms2_hit_count"]] or 0)
        qual = int(f[idx["ms2_qualifying_spectra"]] or 0)
        ms1 = int(f[idx["ms1_hit_count"]] or 0)

        if hits >= 4:
            buckets["SINGLE-SPECTRUM rule (>=4 hits but not qualifying)"] += 1
            spectra_rule += 1
        elif hits == 3:
            buckets["hits == 3 (one short of the threshold)"] += 1
        elif hits in (1, 2):
            buckets[f"hits == {hits} (far under)"] += 1
        else:
            buckets["hits == 0 (no fragment evidence at all)"] += 1
        if hits < 4 and ms1 > 0:
            ms1_rescue += 1
        if len(examples) < 5:
            examples.append((f[idx["sequence"]][:22], hits, tot, qual, ms1))

    print()
    print(f"reference precursors in the evidence dump: kept {kept}, LOST {lost}")
    if lost == 0:
        print("nothing lost -- check that the dump and the reference list describe the same run")
        return
    print()
    print("attribution of the losses:")
    for b, n in sorted(buckets.items(), key=lambda kv: -kv[1]):
        print(f"  {n:6d}  {100.0*n/lost:5.1f}%  {b}")
    print()
    print(f"  {ms1_rescue:6d}  {100.0*ms1_rescue/lost:5.1f}%  of the losses have MS1 evidence "
          f"(ms1_hit_count > 0) despite failing MS2")
    print()
    print("examples (sequence, best_hits, total_hits, qualifying_spectra, ms1_hits):")
    for e in examples:
        print(f"  {e[0]:24s} {e[1]:3d} {e[2]:4d} {e[3]:3d} {e[4]:4d}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
