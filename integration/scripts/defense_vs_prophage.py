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

from scipy.stats import (
    pearsonr,
    spearmanr,
    mannwhitneyu
)

import matplotlib.pyplot as plt
import seaborn as sns

# =========================
# INPUT
# =========================

DEFENSE_FILE = (
    INTEGRATION_MATRICES
    / "defensefinder_matrix.tsv"
)

PROPHAGE_FILE = (
    PROPHAGE_MATRICES
    / "prophage_family_presence_absence.tsv"
)

# =========================
# OUTPUT
# =========================

OUT_STATS = (
    INTEGRATION_STATS
    / "defense_vs_prophage.tsv"
)

SCATTER_PLOT = (
    INTEGRATION_PLOTS
    / "defense_vs_prophage_scatter.png"
)


BOXPLOT_OUT = (
    INTEGRATION_PLOTS
    / "defense_vs_prophage_boxplots.png"
)

# =========================
# LOAD DATA
# =========================

defense = pd.read_csv(
    DEFENSE_FILE,
    sep="\t"
)

prophage = pd.read_csv(
    PROPHAGE_FILE,
    sep="\t"
)

defense["Genome"] = (
    defense["Genome"]
    .astype(str)
    .apply(normalize_genome_name)
)

prophage["Genome"] = (
    prophage["Genome"]
    .astype(str)
    .apply(normalize_genome_name)
)

# =========================
# TOTAL COUNTS
# =========================

defense["Total_defense"] = (
    defense.drop(columns=["Genome"])
    .sum(axis=1)
)

prophage["Total_prophages"] = (
    prophage.drop(columns=["Genome"])
    .sum(axis=1)
)

# =========================
# MERGE
# =========================

merged = defense.merge(
    prophage[
        ["Genome", "Total_prophages"]
    ],
    on="Genome"
)

print(
    f"Genomes after merge: {len(merged)}"
)

# =========================
# CORRELATIONS
# =========================

pearson_r, pearson_p = pearsonr(
    merged["Total_defense"],
    merged["Total_prophages"]
)

spearman_rho, spearman_p = spearmanr(
    merged["Total_defense"],
    merged["Total_prophages"]
)

results = [
    {
        "Test": "Pearson",
        "Statistic": pearson_r,
        "P_value": pearson_p
    },
    {
        "Test": "Spearman",
        "Statistic": spearman_rho,
        "P_value": spearman_p
    }
]

# =========================
# INDIVIDUAL SYSTEMS
# =========================

systems = [
    "Gabija",
    "PrrC",
    "CapRel",
    "Retron_XII",
    "Lamassu_Mrr"
]

for system in systems:

    if system not in merged.columns:
        continue

    present = merged.loc[
        merged[system] == 1,
        "Total_prophages"
    ]

    absent = merged.loc[
        merged[system] == 0,
        "Total_prophages"
    ]

    if len(present) == 0 or len(absent) == 0:
        continue

    U, p = mannwhitneyu(
        present,
        absent,
        alternative="two-sided"
    )

    results.append({
        "Test": system,
        "Statistic": U,
        "P_value": p,
        "Mean_present": present.mean(),
        "Mean_absent": absent.mean()
    })

results_df = pd.DataFrame(
    results
)

results_df.to_csv(
    OUT_STATS,
    sep="\t",
    index=False
)

# =========================
# SCATTERPLOT
# =========================

plt.figure(figsize=(7,6))

sns.regplot(
    data=merged,
    x="Total_defense",
    y="Total_prophages",
    scatter_kws={"alpha":0.7}
)

plt.title(
    "Defense systems vs prophage load"
)

plt.tight_layout()

plt.savefig(
    SCATTER_PLOT,
    dpi=300
)

plt.close()

# =========================
# BOXPLOTS
# =========================

plot_rows = []

for system in systems:

    if system not in merged.columns:
        continue

    tmp = merged[
        [system, "Total_prophages"]
    ].copy()

    tmp["System"] = system

    tmp["Presence"] = tmp[system]

    plot_rows.append(
        tmp[
            [
                "System",
                "Presence",
                "Total_prophages"
            ]
        ]
    )

plot_df = pd.concat(
    plot_rows,
    ignore_index=True
)

plt.figure(figsize=(10,6))

sns.boxplot(
    data=plot_df,
    x="System",
    y="Total_prophages",
    hue="Presence"
)

plt.xticks(rotation=45)

plt.tight_layout()

plt.savefig(
    BOXPLOT_OUT,
    dpi=300
)

plt.close()

# =========================
# SUMMARY
# =========================

print("\nSaved:")
print(OUT_STATS)
print(SCATTER_PLOT)
print(BOXPLOT_OUT)

print("\nCorrelation summary:")
print(results_df.head())