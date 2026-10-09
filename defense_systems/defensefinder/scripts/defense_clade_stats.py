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
import numpy as np

from scipy.stats import kruskal

import matplotlib.pyplot as plt

# ======================================
# INPUT
# ======================================

MATRIX_FILE = (
    DEFENSEFINDER_MATRICES
    / "defensefinder_matrix.tsv"
)

CLADE_FILE = GENOME_CLADES

# ======================================
# OUTPUT
# ======================================

OUTDIR = (
    DEFENSEFINDER_STATS
    / "defense_clade_stats"
)

OUTDIR.mkdir(
    parents=True,
    exist_ok=True
)

MEAN_FILE = (
    OUTDIR
    / "mean_defense_by_clade.tsv"
)

KRUSKAL_FILE = (
    OUTDIR
    / "kruskal_results.tsv"
)

BOXPLOT_FILE = (
    OUTDIR
    / "defense_boxplots_by_clade.png"
)

# ======================================
# LOAD DEFENSE MATRIX
# ======================================

df = pd.read_csv(
    MATRIX_FILE,
    sep="\t"
)

df["Genome"] = (
    df["Genome"]
    .astype(str)
    .apply(normalize_genome_name)
)

# ======================================
# LOAD CLADES
# ======================================

clades = pd.read_csv(
    CLADE_FILE,
    sep="\t",
    comment="#",
    header=None
)

clades.columns = [
    "Genome",
    "Clade"
]

df["Genome"] = (
    df["Genome"]
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

merged = pd.merge(
    df,
    clades,
    on="Genome",
    how="inner"
)

print(
    f"Genomes after merge: {len(merged)}"
)

feature_cols = [
    c
    for c in merged.columns
    if c not in ["Genome", "Clade"]
]

# ======================================
# MEAN BY CLADE
# ======================================

mean_df = (
    merged
    .groupby("Clade")[feature_cols]
    .mean()
    .T
)

mean_df.to_csv(
    MEAN_FILE,
    sep="\t"
)

# ======================================
# KRUSKAL TESTS
# ======================================

results = []

for subtype in feature_cols:

    groups = []

    for clade in sorted(
        merged["Clade"].unique()
    ):

        vals = merged.loc[
            merged["Clade"] == clade,
            subtype
        ]

        groups.append(
            vals.values
        )

    try:

        H, p = kruskal(
            *groups
        )

    except ValueError:

        H = np.nan
        p = np.nan

    results.append({
        "Subtype": subtype,
        "H_statistic": H,
        "P_value": p
    })

kruskal_df = pd.DataFrame(
    results
)

kruskal_df = kruskal_df.sort_values(
    "P_value"
)

kruskal_df.to_csv(
    KRUSKAL_FILE,
    sep="\t",
    index=False
)

# ======================================
# TOTAL DEFENSE SYSTEMS
# ======================================

merged["Total_defense_systems"] = (
    merged[feature_cols]
    .sum(axis=1)
)

# ======================================
# BOXPLOT
# ======================================

clade_order = sorted(
    merged["Clade"].unique()
)

box_data = []

for clade in clade_order:

    vals = merged.loc[
        merged["Clade"] == clade,
        "Total_defense_systems"
    ]

    box_data.append(
        vals
    )

plt.figure(
    figsize=(8,6)
)

plt.boxplot(
    box_data,
    labels=clade_order
)

plt.ylabel(
    "Total defense systems"
)

plt.xlabel(
    "Clade"
)

plt.title(
    "Defense systems per genome by clade"
)

plt.tight_layout()

plt.savefig(
    BOXPLOT_FILE,
    dpi=300
)

plt.close()

# ======================================
# REPORT
# ======================================

print("\nSaved:")

print(MEAN_FILE)

print(KRUSKAL_FILE)

print(BOXPLOT_FILE)