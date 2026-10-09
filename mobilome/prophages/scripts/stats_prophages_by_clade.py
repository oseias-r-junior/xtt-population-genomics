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

PROPHAGE_FILE = (
    PROPHAGE_MATRICES
    / "prophage_summary.tsv"
)

CLADE_FILE = GENOME_CLADES

# ======================================
# OUTPUT
# ======================================

PROPHAGE_STATS.mkdir(
    parents=True,
    exist_ok=True
)

SUMMARY_OUTFILE = (
    PROPHAGE_STATS
    / "prophage_clade_summary.tsv"
)

PRESENCE_OUTFILE = (
    PROPHAGE_STATS
    / "prophage_clade_presence.tsv"
)

KRUSKAL_OUTFILE = (
    PROPHAGE_STATS
    / "prophage_kruskal.tsv"
)

# ======================================
# LOAD
# ======================================

prophages = pd.read_csv(
    PROPHAGE_FILE,
    sep="\t"
)

clades = pd.read_csv(
    CLADE_FILE,
    sep="\t"
)

# ======================================
# NORMALIZE
# ======================================

prophages["Genome"] = (
    prophages["Genome"]
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

df = prophages.merge(
    clades,
    on="Genome",
    how="inner"
)

print(
    f"Merged genomes: {len(df)}"
)

# ======================================
# SUMMARY
# ======================================

summary = (
    df
    .groupby("Clade")
    ["Total_Prophages"]
    .agg(
        N_Genomes="count",
        Mean_Prophages="mean",
        Median_Prophages="median",
        SD_Prophages="std"
    )
    .reset_index()
)

summary.to_csv(
    SUMMARY_OUTFILE,
    sep="\t",
    index=False
)

# ======================================
# PRESENCE / ABSENCE
# ======================================

df["Has_Prophage"] = (
    df["Total_Prophages"] > 0
)

presence_rows = []

for clade, subdf in df.groupby("Clade"):

    n_total = len(subdf)

    n_with = (
        subdf["Has_Prophage"]
        .sum()
    )

    n_without = (
        n_total - n_with
    )

    presence_rows.append({

        "Clade":
            clade,

        "N_Genomes":
            n_total,

        "With_Prophage":
            n_with,

        "Without_Prophage":
            n_without,

        "Presence_Rate":
            round(
                n_with / n_total,
                3
            )

    })

presence = pd.DataFrame(
    presence_rows
)

presence.to_csv(
    PRESENCE_OUTFILE,
    sep="\t",
    index=False
)

# ======================================
# KRUSKAL-WALLIS
# ======================================

groups = []

for clade in sorted(
    df["Clade"].unique()
):

    groups.append(

        df.loc[
            df["Clade"] == clade,
            "Total_Prophages"
        ]

    )

stat, pvalue = kruskal(
    *groups
)

kruskal_df = pd.DataFrame([{

    "Statistic":
        stat,

    "P_value":
        pvalue

}])

kruskal_df.to_csv(
    KRUSKAL_OUTFILE,
    sep="\t",
    index=False
)

# ======================================
# REPORT
# ======================================

print(
    f"Saved: {SUMMARY_OUTFILE}"
)

print(
    f"Saved: {PRESENCE_OUTFILE}"
)

print(
    f"Saved: {KRUSKAL_OUTFILE}"
)

print("\nFinished.")