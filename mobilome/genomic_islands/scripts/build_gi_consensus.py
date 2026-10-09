#!/usr/bin/env python3
"""
mobilome/genomic_islands/scripts/build_gi_consensus.py

3-way region consensus for genomic islands: AlienHunter x IslandPath x
PanISLE (recalibrated). Mirrors the prophage consensus methodology
(build_prophage_consensus_regions.py): regions from all 3 tools are pooled
per genome and merged via Union-Find whenever they show reciprocal overlap
>= 0.5 (i.e. the overlap is at least 50% of BOTH regions' lengths). This
avoids any pre-consensus bias toward a pair of tools run first.

Each merged group becomes one consensus genomic island, annotated with:
    - how many / which of the 3 tools support it (1-3)
    - a functional classification, taken from PanISLE (the only one of the
      3 tools that produces a genuine functional call: virulence /
      resistance / symbiotic / metabolic / unknown). If no PanISLE region
      is part of the merged group, classification is "No_PanISLE_Support".
    - an AlienHunter rRNA-operon caution flag, if any AlienHunter region in
      the group was flagged (see parse_alienhunter_islands.py docstring --
      a known AlienHunter false-positive mode, not an automatic exclusion).

NOTE: like the prophage/plasmid consensus scripts in this pipeline, overlap
is computed on (Genome, Start, End) without a separate Contig key. This
matches the established project convention (see build_prophage_consensus_regions.py)
and is safe in practice because these are near-complete assemblies where GI
calls fall within a single dominant contig; it is flagged here for awareness
in case a future genome set includes more fragmented assemblies.
"""

from pathlib import Path
import sys

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parents[2]

sys.path.append(str(PROJECT_ROOT / "scripts"))

from config import *

import pandas as pd

# ======================================================
# INPUT
# ======================================================

ALIENHUNTER_FILE = ALIENHUNTER / "matrices" / "alienhunter_islands.tsv"
ISLANDPATH_FILE = ISLANDPATH_MATRICES / "islandpath_coordinates.tsv"
PANISLE_FILE = PANISLE_MATRICES / "panisle_islands.tsv"

OVERLAP_THRESHOLD = 0.5

# ======================================================
# OUTPUT
# ======================================================

GI_MATRICES.mkdir(parents=True, exist_ok=True)
GI_STATS.mkdir(parents=True, exist_ok=True)

OUT_CONSENSUS = GI_MATRICES / "gi_consensus_regions.tsv"
OUT_GENOME_SUMMARY = GI_STATS / "gi_consensus_genome_summary.tsv"

# ======================================================
# LOAD
# ======================================================

alienhunter = pd.read_csv(ALIENHUNTER_FILE, sep="\t")
islandpath = pd.read_csv(ISLANDPATH_FILE, sep="\t")
panisle = pd.read_csv(PANISLE_FILE, sep="\t")


def reciprocal_overlap(s1, e1, s2, e2):
    ov = min(e1, e2) - max(s1, s2)
    if ov <= 0:
        return 0.0
    len1 = e1 - s1
    len2 = e2 - s2
    return min(ov / len1, ov / len2)


# ======================================================
# POOL ALL REGIONS PER GENOME (tagged by source tool)
# ======================================================

pooled = []

for _, row in alienhunter.iterrows():
    pooled.append({
        "Genome": row["Genome"],
        "Tool": "AlienHunter",
        "Start": int(row["Start"]),
        "End": int(row["End"]),
        "RRNA_Operon_Flag": bool(row["RRNA_Operon_Flag"]),
        "Classification": None,
    })

for _, row in islandpath.iterrows():
    pooled.append({
        "Genome": row["Genome"],
        "Tool": "IslandPath",
        "Start": int(row["Start"]),
        "End": int(row["End"]),
        "RRNA_Operon_Flag": False,
        "Classification": None,
    })

for _, row in panisle.iterrows():
    pooled.append({
        "Genome": row["Genome"],
        "Tool": "PanISLE",
        "Start": int(row["Start"]),
        "End": int(row["End"]),
        "RRNA_Operon_Flag": False,
        "Classification": row.get("Classification", "unknown"),
    })

pooled_df = pd.DataFrame(pooled)

# ======================================================
# UNION-FIND PER GENOME
# ======================================================

class UnionFind:
    def __init__(self, n):
        self.parent = list(range(n))

    def find(self, x):
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, x, y):
        rx, ry = self.find(x), self.find(y)
        if rx != ry:
            self.parent[ry] = rx


records = []

for genome, group in pooled_df.groupby("Genome"):
    group = group.reset_index(drop=True)
    n = len(group)
    uf = UnionFind(n)

    for i in range(n):
        for j in range(i + 1, n):
            ov = reciprocal_overlap(
                group.loc[i, "Start"], group.loc[i, "End"],
                group.loc[j, "Start"], group.loc[j, "End"],
            )
            if ov >= OVERLAP_THRESHOLD:
                uf.union(i, j)

    clusters = {}
    for i in range(n):
        root = uf.find(i)
        clusters.setdefault(root, []).append(i)

    for cluster_idx, indices in enumerate(clusters.values(), start=1):
        sub = group.loc[indices]

        tools = sorted(sub["Tool"].unique().tolist())
        start = int(sub["Start"].min())
        end = int(sub["End"].max())

        panisle_calls = sub[sub["Tool"] == "PanISLE"]["Classification"]
        non_unknown = panisle_calls[panisle_calls != "unknown"]
        if len(non_unknown) > 0:
            classification = non_unknown.mode().iat[0]
        elif len(panisle_calls) > 0:
            classification = "unknown"
        else:
            classification = "No_PanISLE_Support"

        rrna_flag = bool(sub["RRNA_Operon_Flag"].any())

        records.append({
            "Genome": genome,
            "Region": f"GI_{cluster_idx}",
            "Start": start,
            "End": end,
            "Length": end - start,
            "N_Tools_Support": len(tools),
            "Tools_Supporting": ",".join(tools),
            "Classification": classification,
            "AlienHunter_rRNA_Operon_Flag": rrna_flag,
        })

result = pd.DataFrame(records)
# ======================================================
# GI x PPH EXCLUSION (Uceda-Campos convention): "GI regions overlapping to
# prophage regions were not considered." A GI region is dropped if it
# overlaps (any overlap, same Genome) a region from the high-confidence
# (all-tools) prophage consensus set.
# ======================================================
PROPHAGE_HIGH_CONF_FILE = PROJECT_ROOT / "mobilome" / "prophages" / "matrices" / "prophage_high_confidence.tsv"
if PROPHAGE_HIGH_CONF_FILE.exists():
    pph_hc = pd.read_csv(PROPHAGE_HIGH_CONF_FILE, sep="\t")
    pph_by_genome = {g: sub[["Start", "Stop"]].values.tolist() for g, sub in pph_hc.groupby("Genome")}

    def _overlaps_prophage(row):
        for pstart, pstop in pph_by_genome.get(row["Genome"], []):
            if row["Start"] <= pstop and row["End"] >= pstart:
                return True
        return False

    before_n = len(result)
    excl_mask = result.apply(_overlaps_prophage, axis=1)
    excluded_gi = result[excl_mask]
    result = result[~excl_mask].copy()
    print(f"GI regions excluded (overlap with high-confidence prophage): {excl_mask.sum()} of {before_n}")
    if len(excluded_gi):
        print(excluded_gi[["Genome", "Region", "Start", "End"]].to_string(index=False))
else:
    print(f"[WARN] {PROPHAGE_HIGH_CONF_FILE} not found -- skipping GI x PPH exclusion")

result = result.sort_values(["Genome", "Start"]).reset_index(drop=True)
result.to_csv(OUT_CONSENSUS, sep="\t", index=False)

genome_summary = (
    result
    .groupby("Genome")
    .agg(
        Total_GI_regions=("Region", "count"),
        N_supported_by_3_tools=("N_Tools_Support", lambda s: (s == 3).sum()),
        N_supported_by_2_tools=("N_Tools_Support", lambda s: (s == 2).sum()),
        N_supported_by_1_tool=("N_Tools_Support", lambda s: (s == 1).sum()),
        Total_GI_bp=("Length", "sum"),
    )
    .reset_index()
)
genome_summary.to_csv(OUT_GENOME_SUMMARY, sep="\t", index=False)

print(f"Genomes with >=1 GI consensus region: {result['Genome'].nunique()}")
print(f"Total consensus GI regions: {len(result)}")
print()
print("Support level:")
print(result["N_Tools_Support"].value_counts().sort_index().to_string())
print()
print("Classification:")
print(result["Classification"].value_counts().to_string())
print()
print(f"Regions flagged as rRNA-operon-overlap (AlienHunter): {int(result['AlienHunter_rRNA_Operon_Flag'].sum())}")
print()
print(f"Saved: {OUT_CONSENSUS}")
print(f"Saved: {OUT_GENOME_SUMMARY}")
