#!/usr/bin/env python3
"""
mobilome/prophages/scripts/build_prophage_consensus_regions.py

Genuine 3-way, region-level prophage consensus: PHASTEST x VirSorter2 x
geNomad, built from scratch with all three tools competing on equal footing
(no PHASTEST/VirSorter2 "pre-consensus" with geNomad bolted on afterward).

For each genome, region calls from the three tools are clustered by
RECIPROCAL overlap >= 50% (overlap must be >=50% of BOTH region lengths,
i.e. bedtools -f 0.5 -r style) using a union-find over all calls. Each
resulting cluster is one consensus prophage region, annotated with:
    - how many / which tools support it (1-3)
    - the union coordinate extent
    - per-tool provenance fields (PHASTEST completeness/phage, VirSorter2
      group/score, geNomad taxonomy/score)
    - a taxonomic classification (Inovirus / Caudoviricetes / Other_classified
      / No_GeNomad_Support), always derived only from geNomad (the only tool
      of the three with a real taxonomy)

This replaces the old two-step pipeline (build_prophage_consensus_matrix.py
PHASTEST x VirSorter2 genome-count consensus + prophage_taxonomy.py geNomad
bolt-on), which gave PHASTEST/VirSorter2 first-class status and geNomad
second-class status in defining region boundaries.
"""

from pathlib import Path
import sys

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parents[2]

sys.path.append(str(PROJECT_ROOT / "scripts"))
sys.path.append(str(PROJECT_ROOT / "scripts" / "metadata"))

from config import *
from normalize_genome_names import normalize_genome_name

import pandas as pd

# ======================================================
# INPUT
# ======================================================

PHASTEST_REGIONS = PROPHAGE_MATRICES / "prophage_regions.tsv"
VIRSORTER2_REGIONS = VIRSORTER2_MATRICES / "virsorter2_regions.tsv"

# ======================================================
# OUTPUT
# ======================================================

PROPHAGE_MATRICES.mkdir(parents=True, exist_ok=True)
PROPHAGE_STATS.mkdir(parents=True, exist_ok=True)

OUT_REGIONS = PROPHAGE_MATRICES / "prophage_consensus_regions.tsv"
OUT_GENOME_SUMMARY = PROPHAGE_STATS / "prophage_consensus_genome_summary.tsv"

RECIPROCAL_OVERLAP_THRESHOLD = 0.5


# ======================================================
# UNION-FIND
# ======================================================

class UnionFind:
    def __init__(self, n):
        self.parent = list(range(n))

    def find(self, x):
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.parent[ra] = rb


def reciprocal_overlap(s1, e1, s2, e2):
    ov = min(e1, e2) - max(s1, s2)
    if ov <= 0:
        return 0.0
    len1 = e1 - s1
    len2 = e2 - s2
    if len1 <= 0 or len2 <= 0:
        return 0.0
    return min(ov / len1, ov / len2)


def classify(has_genomad_call: bool, taxonomy: str) -> str:
    if not has_genomad_call:
        return "No_GeNomad_Support"
    if not isinstance(taxonomy, str) or taxonomy.strip() == "":
        return "Other_classified"
    if "Inoviridae" in taxonomy:
        return "Inovirus"
    if "Caudoviricetes" in taxonomy:
        return "Caudoviricetes"
    return "Other_classified"


# ======================================================
# LOAD PHASTEST REGIONS
# ======================================================

phastest = pd.read_csv(PHASTEST_REGIONS, sep="\t")
phastest["Genome"] = phastest["Genome"].astype(str).apply(normalize_genome_name)

# ======================================================
# LOAD VIRSORTER2 REGIONS
# ======================================================

virsorter2 = pd.read_csv(VIRSORTER2_REGIONS, sep="\t")
virsorter2["Genome"] = virsorter2["Genome"].astype(str).apply(normalize_genome_name)


# ======================================================
# LOAD GENOMAD REGIONS (per genome, on demand)
# ======================================================

def load_genomad_calls(genome: str) -> pd.DataFrame:
    summary_file = (
        GENOMAD_RAW
        / genome
        / f"{genome}_summary"
        / f"{genome}_virus_summary.tsv"
    )

    cols = ["start", "end", "virus_score", "taxonomy"]

    if not summary_file.exists():
        return pd.DataFrame(columns=cols)

    df = pd.read_csv(summary_file, sep="\t")

    if df.empty:
        return pd.DataFrame(columns=cols)

    coords = df["coordinates"].astype(str).str.split("-", expand=True)
    df["start"] = pd.to_numeric(coords[0], errors="coerce")
    df["end"] = pd.to_numeric(coords[1], errors="coerce")

    missing = df["start"].isna() | df["end"].isna()
    df.loc[missing, "start"] = 1
    df.loc[missing, "end"] = df.loc[missing, "length"]

    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)

    return df[cols]


# ======================================================
# BUILD PER-GENOME CALL LISTS
# ======================================================

all_genomes = sorted(
    set(phastest["Genome"])
    | set(virsorter2["Genome"])
    | {p.name for p in GENOMAD_RAW.iterdir() if p.is_dir()}
)

consensus_records = []
genome_summary_records = []

for genome in all_genomes:

    calls = []  # each: dict(tool, start, stop, **fields)

    for _, row in phastest[phastest["Genome"] == genome].iterrows():
        calls.append({
            "tool": "PHASTEST",
            "start": int(row["Start"]),
            "stop": int(row["Stop"]),
            "PHASTEST_Completeness": row.get("Completeness", ""),
            "PHASTEST_Most_Common_Phage": row.get("Most_Common_Phage", ""),
            "PHASTEST_GC": row.get("GC", None),
        })

    for _, row in virsorter2[virsorter2["Genome"] == genome].iterrows():
        calls.append({
            "tool": "VirSorter2",
            "start": int(row["Start"]),
            "stop": int(row["Stop"]),
            "VirSorter2_Group": row.get("Group", ""),
            "VirSorter2_Score": row.get("Score", None),
            "VirSorter2_Hallmark_cnt": row.get("Hallmark_cnt", None),
        })

    genomad_df = load_genomad_calls(genome)
    for _, row in genomad_df.iterrows():
        calls.append({
            "tool": "geNomad",
            "start": int(row["start"]),
            "stop": int(row["end"]),
            "GeNomad_Taxonomy": row.get("taxonomy", ""),
            "GeNomad_Virus_Score": row.get("virus_score", None),
        })

    if not calls:
        continue

    # ---- cluster calls by reciprocal overlap >= threshold ----
    n = len(calls)
    uf = UnionFind(n)

    for i in range(n):
        for j in range(i + 1, n):
            ro = reciprocal_overlap(
                calls[i]["start"], calls[i]["stop"],
                calls[j]["start"], calls[j]["stop"]
            )
            if ro >= RECIPROCAL_OVERLAP_THRESHOLD:
                uf.union(i, j)

    clusters = {}
    for i in range(n):
        root = uf.find(i)
        clusters.setdefault(root, []).append(calls[i])

    # ---- summarize each cluster into one consensus region ----
    region_list = []
    for cluster_calls in clusters.values():
        start = min(c["start"] for c in cluster_calls)
        stop = max(c["stop"] for c in cluster_calls)
        tools = sorted({c["tool"] for c in cluster_calls})

        phastest_calls = [c for c in cluster_calls if c["tool"] == "PHASTEST"]
        virsorter2_calls = [c for c in cluster_calls if c["tool"] == "VirSorter2"]
        genomad_calls = [c for c in cluster_calls if c["tool"] == "geNomad"]

        best_genomad = None
        if genomad_calls:
            best_genomad = max(
                genomad_calls,
                key=lambda c: (c["GeNomad_Virus_Score"] or 0)
            )

        region_list.append({
            "Genome": genome,
            "Start": start,
            "Stop": stop,
            "Length_bp": stop - start,
            "N_Tools_Support": len(tools),
            "Tools_Supporting": ",".join(tools),
            "Classification": classify(
                bool(genomad_calls),
                best_genomad["GeNomad_Taxonomy"] if best_genomad else ""
            ),
            "GeNomad_Taxonomy": best_genomad["GeNomad_Taxonomy"] if best_genomad else "",
            "GeNomad_Virus_Score": best_genomad["GeNomad_Virus_Score"] if best_genomad else None,
            "PHASTEST_Completeness": phastest_calls[0]["PHASTEST_Completeness"] if phastest_calls else "",
            "PHASTEST_Most_Common_Phage": phastest_calls[0]["PHASTEST_Most_Common_Phage"] if phastest_calls else "",
            "VirSorter2_Group": virsorter2_calls[0]["VirSorter2_Group"] if virsorter2_calls else "",
            "VirSorter2_Score": virsorter2_calls[0]["VirSorter2_Score"] if virsorter2_calls else None,
        })

    region_list.sort(key=lambda r: r["Start"])
    for i, r in enumerate(region_list, start=1):
        r["Region"] = i
        consensus_records.append(r)

    genome_summary_records.append({
        "Genome": genome,
        "Total_consensus_regions": len(region_list),
        "N_supported_by_3_tools": sum(1 for r in region_list if r["N_Tools_Support"] == 3),
        "N_supported_by_2_tools": sum(1 for r in region_list if r["N_Tools_Support"] == 2),
        "N_supported_by_1_tool": sum(1 for r in region_list if r["N_Tools_Support"] == 1),
    })


# ======================================================
# WRITE OUTPUT
# ======================================================

cols_order = [
    "Genome", "Region", "Start", "Stop", "Length_bp",
    "N_Tools_Support", "Tools_Supporting", "Classification",
    "GeNomad_Taxonomy", "GeNomad_Virus_Score",
    "PHASTEST_Completeness", "PHASTEST_Most_Common_Phage",
    "VirSorter2_Group", "VirSorter2_Score",
]

result = pd.DataFrame(consensus_records)[cols_order]
result.to_csv(OUT_REGIONS, sep="\t", index=False)

genome_summary = pd.DataFrame(genome_summary_records)
genome_summary.to_csv(OUT_GENOME_SUMMARY, sep="\t", index=False)

print(f"Genomes with >=1 consensus region: {len(genome_summary)}")
print(f"Total consensus regions: {len(result)}")
print()
print("Support level (N tools agreeing on a region):")
print(result["N_Tools_Support"].value_counts().sort_index().to_string())
print()
print("Classification:")
print(result["Classification"].value_counts().to_string())
print()
print(f"Saved: {OUT_REGIONS}")
print(f"Saved: {OUT_GENOME_SUMMARY}")
