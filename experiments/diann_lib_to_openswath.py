#!/usr/bin/env python3
"""Convert a DIA-NN spectral library to the OpenSWATH transition TSV format.

Both tools must search the *same* library, otherwise the comparison measures
library quality rather than extraction algorithm. DIA-NN's library columns map
onto OpenSWATH's almost one-to-one -- this is mostly a rename plus a couple of
derived columns.

Target column names come from TransitionTSVFile::header_names_ in OpenMS 3.6.

    python3 diann_lib_to_openswath.py lib.parquet lib_osw.tsv
"""
import sys

import pandas as pd

# DIA-NN column -> OpenSWATH column. First DIA-NN name present wins, so
# alternatives can be listed for fields DIA-NN has renamed across versions.
RENAME = {
    "PrecursorMz": ["PrecursorMz"],
    "ProductMz": ["ProductMz"],
    "PrecursorCharge": ["PrecursorCharge"],
    "ProductCharge": ["FragmentCharge", "ProductCharge"],
    "LibraryIntensity": ["LibraryIntensity", "RelativeIntensity"],
    "NormalizedRetentionTime": ["Tr_recalibrated", "RetentionTime", "iRT"],
    "PeptideSequence": ["PeptideSequence", "StrippedPeptide"],
    "ModifiedPeptideSequence": ["ModifiedPeptide", "FullUniModPeptideName",
                                "ModifiedPeptideSequence"],
    "ProteinId": ["ProteinGroup", "ProteinName", "ProteinId"],
    "GeneName": ["Genes", "GeneName"],
    "FragmentType": ["FragmentType"],
    "FragmentSeriesNumber": ["FragmentSeriesNumber", "FragmentNumber"],
    "TransitionGroupId": ["transition_group_id", "TransitionGroupId"],
    "TransitionId": ["transition_name", "TransitionId"],
    "Decoy": ["decoy", "Decoy"],
    "PrecursorIonMobility": ["IonMobility", "PrecursorIonMobility"],
}

REQUIRED = ["PrecursorMz", "ProductMz", "PrecursorCharge", "LibraryIntensity",
            "NormalizedRetentionTime", "PeptideSequence", "ProteinId"]


def convert(df):
    out = pd.DataFrame()
    for target, sources in RENAME.items():
        for s in sources:
            if s in df.columns:
                out[target] = df[s]
                break

    missing = [c for c in REQUIRED if c not in out.columns]
    if missing:
        raise SystemExit(f"error: library lacks required columns: {missing}\n"
                         f"available: {sorted(df.columns)}")

    # OpenSWATH needs stable group/transition ids; DIA-NN omits them when the
    # library is predicted rather than empirical.
    if "TransitionGroupId" not in out.columns:
        out["TransitionGroupId"] = (out["ModifiedPeptideSequence"].astype(str)
                                    + "_" + out["PrecursorCharge"].astype(str))
    if "TransitionId" not in out.columns:
        out["TransitionId"] = (out["TransitionGroupId"].astype(str) + "_"
                               + out.groupby("TransitionGroupId").cumcount().astype(str))
    if "Decoy" not in out.columns:
        out["Decoy"] = 0

    # OpenSWATH treats every transition as detecting unless told otherwise.
    out["DetectingTransition"] = 1
    out["IdentifyingTransition"] = 0
    out["QuantifyingTransition"] = 1
    if "ModifiedPeptideSequence" not in out.columns:
        out["ModifiedPeptideSequence"] = out["PeptideSequence"]
    return out


def selftest():
    df = pd.DataFrame({
        "PrecursorMz": [500.1, 500.1], "ProductMz": [301.2, 402.3],
        "PrecursorCharge": [2, 2], "FragmentCharge": [1, 1],
        "LibraryIntensity": [1000.0, 500.0], "Tr_recalibrated": [30.5, 30.5],
        "PeptideSequence": ["SLYNTVATL"] * 2, "ModifiedPeptide": ["SLYNTVATL"] * 2,
        "ProteinGroup": ["P12345"] * 2, "FragmentType": ["y", "b"],
        "FragmentSeriesNumber": [3, 4],
    })
    out = convert(df)
    assert list(out["ProductCharge"]) == [1, 1]
    assert list(out["NormalizedRetentionTime"]) == [30.5, 30.5]
    assert out["TransitionGroupId"].nunique() == 1, "same precursor -> one group"
    assert out["TransitionId"].nunique() == 2, "distinct transitions"
    assert (out["Decoy"] == 0).all()
    assert all(c in out.columns for c in REQUIRED)
    print("selftest OK")


def main():
    if len(sys.argv) == 2 and sys.argv[1] == "--selftest":
        return selftest()
    if len(sys.argv) != 3:
        raise SystemExit("usage: diann_lib_to_openswath.py <lib.parquet|tsv> <out.tsv>")

    src, dst = sys.argv[1], sys.argv[2]
    df = (pd.read_parquet(src) if src.endswith((".parquet", ".pq"))
          else pd.read_csv(src, sep="\t"))
    out = convert(df)
    out.to_csv(dst, sep="\t", index=False)
    print(f"{len(df):,} rows -> {dst}")
    print(f"{out['TransitionGroupId'].nunique():,} precursors, "
          f"{len(out) / max(1, out['TransitionGroupId'].nunique()):.1f} transitions each")


if __name__ == "__main__":
    main()
