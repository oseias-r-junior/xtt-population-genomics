#!/usr/bin/env python3

from pathlib import Path
import sys

SCRIPT_DIR = Path(__file__).resolve().parent

PROJECT_ROOT = SCRIPT_DIR.parents[2]

sys.path.append(
    str(PROJECT_ROOT / "scripts")
)

from config import *
from normalize_genome_names import normalize_genome_name

import pandas as pd
from scipy.stats import kruskal

# ======================================
# INPUT
# ======================================

REGIONS_FILE = (
    PROPHAGE_MATRICES
    / "prophage_regions.tsv"
)

CLADE_FILE = GENOME_CLADES

# ======================================
# OUTPUT
# ======================================

PROPHAGE_STATS.mkdir(
    parents=True,
    exist_ok=True
)

LENGTH_OUTFILE = (
    PROPHAGE_STATS
    / "prophage_length_by_clade.tsv"
)

GC_OUTFILE = (
    PROPHAGE_STATS
    / "prophage_gc_by_clade.tsv"
)

COMPLETENESS_OUTFILE = (
    PROPHAGE_STATS
    / "prophage_completeness_by_clade.tsv"
)

KRUSKAL_LENGTH = (
    PROPHAGE_STATS
    / "prophage_length_kruskal.tsv"
)

KRUSKAL_GC = (
    PROPHAGE_STATS
    / "prophage_gc_kruskal.tsv"
)

# ======================================
# LOAD
# ======================================

regions = pd.read_csv(
    REGIONS_FILE,
    sep="\t"
)

clades = pd.read_csv(
    CLADE_FILE,
    sep="\t"
)

# ======================================
# NORMALIZE
# ======================================

regions["Genome"] = (
    regions["Genome"]
    .astype(str)
    .apply(normalize_genome_name)
)

clades["Genome"] = (
    clades["Genome"]
    .astype(str)
    .apply(normalize_genome_name)
)

# ======================================
# MERGE
# ======================================

df = regions.merge(
    clades,
    on="Genome",
    how="inner"
)

print(
    f"Merged prophage regions: {len(df)}"
)

# ======================================
# LENGTH STATS
# ======================================

length_stats = (
    df
    .groupby("Clade")
    ["Length_bp"]
    .agg(
        N_Prophages="count",
        Mean_Length_bp="mean",
        Median_Length_bp="median",
        SD_Length_bp="std"
    )
    .reset_index()
)

length_stats.to_csv(
    LENGTH_OUTFILE,
    sep="\t",
    index=False
)

# ======================================
# GC STATS
# ======================================

gc_stats = (
    df
    .groupby("Clade")
    ["GC"]
    .agg(
        Mean_GC="mean",
        Median_GC="median",
        SD_GC="std"
    )
    .reset_index()
)

gc_stats.to_csv(
    GC_OUTFILE,
    sep="\t",
    index=False
)

# ======================================
# COMPLETENESS
# ======================================

comp = pd.crosstab(
    df["Clade"],
    df["Completeness"]
)
print(comp.columns)
print(comp.head())

for col in [
    "intact",
    "questionable",
    "incomplete"
]:

    if col not in comp.columns:

        comp[col] = 0

comp = comp.reset_index()

comp = comp[
    [
        "intact",
        "questionable",
        "incomplete"
    ]
]

comp = comp.reset_index()

comp.columns = [
    "Clade",
    "Intact",
    "Questionable",
    "Incomplete"
]

comp.to_csv(
    COMPLETENESS_OUTFILE,
    sep="\t",
    index=False
)

# ======================================
# KRUSKAL - LENGTH
# ======================================

length_groups = []

for clade in sorted(
    df["Clade"].unique()
):

    length_groups.append(

        df.loc[
            df["Clade"] == clade,
            "Length_bp"
        ]

    )

if len(length_groups) >= 2:

    stat, pvalue = kruskal(
        *length_groups
    )

    pd.DataFrame([{

        "Statistic":
            stat,

        "P_value":
            pvalue

    }]).to_csv(
        KRUSKAL_LENGTH,
        sep="\t",
        index=False
    )

# ======================================
# KRUSKAL - GC
# ======================================

gc_groups = []

for clade in sorted(
    df["Clade"].unique()
):

    gc_groups.append(

        df.loc[
            df["Clade"] == clade,
            "GC"
        ]

    )

if len(gc_groups) >= 2:

    stat, pvalue = kruskal(
        *gc_groups
    )

    pd.DataFrame([{

        "Statistic":
            stat,

        "P_value":
            pvalue

    }]).to_csv(
        KRUSKAL_GC,
        sep="\t",
        index=False
    )

# ======================================
# REPORT
# ======================================

print(
    f"Saved: {LENGTH_OUTFILE}"
)

print(
    f"Saved: {GC_OUTFILE}"
)

print(
    f"Saved: {COMPLETENESS_OUTFILE}"
)

print(
    f"Saved: {KRUSKAL_LENGTH}"
)

print(
    f"Saved: {KRUSKAL_GC}"
)

print("\nFinished.")