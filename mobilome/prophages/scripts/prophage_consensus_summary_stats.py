#!/usr/bin/env python3
"""
mobilome/prophages/scripts/prophage_consensus_summary_stats.py

Summary stats + plots for the 3-way (PHASTEST x VirSorter2 x geNomad)
region-level prophage consensus (prophage_consensus_regions.tsv).

Replaces the old 2-tool (PHASTEST x VirSorter2) genome-count version.
"""

from pathlib import Path
import sys

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parents[2]

sys.path.append(str(PROJECT_ROOT / "scripts"))

from config import *

import pandas as pd
import matplotlib.pyplot as plt

# ==========================================================
# INPUT
# ==========================================================

REGIONS_FILE = PROPHAGE_MATRICES / "prophage_consensus_regions.tsv"
GENOME_SUMMARY_FILE = PROPHAGE_STATS / "prophage_consensus_genome_summary.tsv"

# ==========================================================
# OUTPUT
# ==========================================================

PROPHAGE_STATS.mkdir(parents=True, exist_ok=True)
PROPHAGE_PLOTS.mkdir(parents=True, exist_ok=True)

OUT_STATS = PROPHAGE_STATS / "prophage_consensus_summary.tsv"
OUT_PREVALENCE = PROPHAGE_PLOTS / "prophage_consensus_regions_prevalence.png"
OUT_SUPPORT = PROPHAGE_PLOTS / "prophage_consensus_support_level.png"
OUT_CLASSIFICATION = PROPHAGE_PLOTS / "prophage_consensus_classification.png"
OUT_TOOLS_COMBO = PROPHAGE_PLOTS / "prophage_consensus_tools_combination.png"

# ==========================================================
# LOAD
# ==========================================================

regions = pd.read_csv(REGIONS_FILE, sep="\t")
genome_summary = pd.read_csv(GENOME_SUMMARY_FILE, sep="\t")

print(f"Genomes with >=1 consensus region: {len(genome_summary)}")
print(f"Total consensus regions: {len(regions)}")

# ==========================================================
# SUMMARY TABLE
# ==========================================================

summary = pd.DataFrame({
    "Metric": [
        "Genomes_with_regions",
        "Total_consensus_regions",
        "Mean_regions_per_genome",
        "Median_regions_per_genome",
        "Max_regions_per_genome",
        "Regions_supported_by_3_tools",
        "Regions_supported_by_2_tools",
        "Regions_supported_by_1_tool",
        "Percent_3_tool_support",
        "Classification_Inovirus",
        "Classification_Caudoviricetes",
        "Classification_Other_classified",
        "Classification_No_GeNomad_Support",
    ],
    "Value": [
        len(genome_summary),
        len(regions),
        genome_summary["Total_consensus_regions"].mean(),
        genome_summary["Total_consensus_regions"].median(),
        genome_summary["Total_consensus_regions"].max(),
        (regions["N_Tools_Support"] == 3).sum(),
        (regions["N_Tools_Support"] == 2).sum(),
        (regions["N_Tools_Support"] == 1).sum(),
        100 * (regions["N_Tools_Support"] == 3).mean(),
        (regions["Classification"] == "Inovirus").sum(),
        (regions["Classification"] == "Caudoviricetes").sum(),
        (regions["Classification"] == "Other_classified").sum(),
        (regions["Classification"] == "No_GeNomad_Support").sum(),
    ]
})

summary.to_csv(OUT_STATS, sep="\t", index=False)

# ==========================================================
# PREVALENCE (regions per genome)
# ==========================================================

counts = (
    genome_summary["Total_consensus_regions"]
    .value_counts()
    .sort_index()
)

plt.figure(figsize=(6, 4))
plt.bar(counts.index.astype(str), counts.values)
plt.xlabel("Consensus prophage regions")
plt.ylabel("Genomes")
plt.title("Prophage regions per genome")
plt.tight_layout()
plt.savefig(OUT_PREVALENCE, dpi=300)
plt.close()

# ==========================================================
# SUPPORT LEVEL (how many tools agree per region)
# ==========================================================

support = (
    regions["N_Tools_Support"]
    .value_counts()
    .sort_index()
)

plt.figure(figsize=(4, 4))
plt.bar(support.index.astype(str), support.values)
plt.xlabel("Tools supporting region")
plt.ylabel("Regions")
plt.title("3-way consensus support level")
plt.tight_layout()
plt.savefig(OUT_SUPPORT, dpi=300)
plt.close()

# ==========================================================
# CLASSIFICATION
# ==========================================================

classification = regions["Classification"].value_counts()

plt.figure(figsize=(6, 4))
plt.bar(classification.index, classification.values)
plt.ylabel("Regions")
plt.title("Prophage taxonomic classification (geNomad)")
plt.xticks(rotation=20, ha="right")
plt.tight_layout()
plt.savefig(OUT_CLASSIFICATION, dpi=300)
plt.close()

# ==========================================================
# TOOLS COMBINATION (which tool(s) detected each region)
# ==========================================================

combo = regions["Tools_Supporting"].value_counts()

plt.figure(figsize=(8, 4))
plt.bar(combo.index, combo.values)
plt.ylabel("Regions")
plt.title("Detecting tool combination")
plt.xticks(rotation=30, ha="right")
plt.tight_layout()
plt.savefig(OUT_TOOLS_COMBO, dpi=300)
plt.close()

# ==========================================================
# REPORT
# ==========================================================

print()
print(f"3-tool support: {(regions['N_Tools_Support'] == 3).sum()} / {len(regions)}")
print(f"Classification breakdown:\n{classification.to_string()}")
print()
print(OUT_STATS)
print(OUT_PREVALENCE)
print(OUT_SUPPORT)
print(OUT_CLASSIFICATION)
print(OUT_TOOLS_COMBO)
