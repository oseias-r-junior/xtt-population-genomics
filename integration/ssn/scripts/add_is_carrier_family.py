#!/usr/bin/env python3
"""
add_is_carrier_family.py -- augments ssn_nodes_final.tsv with an
IS_carrier_family column.

Context: hosts (PPH/GI/PLS nodes) that contain an internal IS element
already carry a boolean IS_carrier_flag, but not WHICH IS family is
inside them (that information was only available on the IS node itself,
and most of those IS nodes were later removed from the network as
redundant, since the info is now represented via the flag/ring instead).

This script re-derives, for each flagged host, the family of the IS
instance(s) overlapping its coordinates, using the IS source matrix
(GBK-curated, with Family/Genome/Contig/Start/End columns) -- the same
one used upstream in build_ssn_network.py. If a host contains IS from
more than one family (should be rare), the DOMINANT (most frequent)
family is kept, and the count of distinct families found is also
recorded, so this simplification is auditable rather than silent.

Source file: mobilome/insertion_sequences/matrices/is_gbk_confirmed_instances.tsv
(columns: Genome, Contig, Family, isBegin, isEnd, isLen, GBK_Confirmed,
GBK_file_found) -- the GBK-curated IS instance table already used
upstream (10,140 -> 8,698 after curation). Only rows with
GBK_Confirmed truthy are used here, matching that curation.

Usage
-----
python add_is_carrier_family.py
"""
from collections import Counter
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[3]
NODES_FILE = PROJECT_ROOT / "integration" / "ssn" / "matrices" / "ssn_nodes_final.tsv"
IS_SOURCE_FILE = (PROJECT_ROOT / "mobilome" / "insertion_sequences" / "matrices"
                   / "is_gbk_confirmed_instances.tsv")


def load_is_source():
    if not IS_SOURCE_FILE.exists():
        print(f"[FATAL] Expected IS source file not found: {IS_SOURCE_FILE}")
        return None
    df = pd.read_csv(IS_SOURCE_FILE, sep="\t")
    df = df.rename(columns={"isBegin": "Start", "isEnd": "End"})
    df["Contig"] = df["Contig"].astype(str)
    if "GBK_Confirmed" in df.columns:
        before = len(df)
        truthy = df["GBK_Confirmed"].astype(str).str.strip().str.lower().isin(
            {"true", "1", "yes", "y", "t"})
        df = df[truthy].copy()
        print(f"[INFO] IS source rows: {before} total, {len(df)} GBK_Confirmed")
    else:
        print(f"[INFO] IS source rows: {len(df)} (no GBK_Confirmed column found, using all)")
    return df


def main():
    is_df = load_is_source()
    if is_df is None:
        return
    print(f"[INFO] Using IS source file: {IS_SOURCE_FILE}")

    nodes = pd.read_csv(NODES_FILE, sep="\t")
    if "IS_carrier_flag" not in nodes.columns:
        print("[FATAL] ssn_nodes_final.tsv has no IS_carrier_flag column -- run build_ssn_network.py first.")
        return
    nodes["Contig"] = nodes["Contig"].astype(str)

    hosts_mask = nodes["Category"].isin(["PPH", "GI", "PLS"]) & (nodes["IS_carrier_flag"] == True)  # noqa: E712
    hosts = nodes[hosts_mask]
    print(f"[INFO] {len(hosts)} host nodes flagged as IS-carrying")

    is_by_genome_contig = {}
    for _, row in is_df.iterrows():
        key = (row["Genome"], str(row["Contig"]))
        fam = row.get("Family", "unknown")
        fam = fam if pd.notna(fam) and str(fam).strip() else "unknown"
        is_by_genome_contig.setdefault(key, []).append((row["Start"], row["End"], fam))

    dominant_family, n_families = {}, {}
    for _, row in hosts.iterrows():
        key = (row["Genome"], str(row["Contig"]))
        h_start, h_end = row["Start"], row["End"]
        fams = []
        for is_start, is_end, fam in is_by_genome_contig.get(key, []):
            overlap = min(h_end, is_end) - max(h_start, is_start)
            if overlap > 0:
                fams.append(fam)
        if fams:
            counts = Counter(fams)
            dominant_family[row["NodeID"]] = counts.most_common(1)[0][0]
            n_families[row["NodeID"]] = len(counts)

    nodes["IS_carrier_family"] = nodes["NodeID"].map(dominant_family).fillna("")
    nodes["IS_carrier_n_families"] = nodes["NodeID"].map(n_families).fillna(0).astype(int)

    matched = (nodes["IS_carrier_family"] != "").sum()
    print(f"[INFO] Hosts with a resolved dominant IS family: {matched} of {hosts_mask.sum()}")
    print("[INFO] Dominant family distribution among IS-carrying hosts:")
    print(nodes.loc[nodes["IS_carrier_family"] != "", "IS_carrier_family"].value_counts().to_string())
    multi = int((nodes["IS_carrier_n_families"] > 1).sum())
    print(f"[INFO] Hosts carrying >1 distinct IS family (simplified to dominant): {multi}")
    unresolved = int(hosts_mask.sum() - matched)
    if unresolved:
        print(f"[WARN] {unresolved} flagged hosts had no coordinate overlap found in {IS_SOURCE_FILE.name} "
              "-- check Genome/Contig naming consistency if this number is large.")

    nodes.to_csv(NODES_FILE, sep="\t", index=False)
    print(f"[DONE] Updated {NODES_FILE} with IS_carrier_family / IS_carrier_n_families columns.")


if __name__ == "__main__":
    main()
