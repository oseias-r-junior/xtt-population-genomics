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

# =====================================
# INPUT
# =====================================

DEFENSE_FILE = (
    INTEGRATION_MATRICES
    / "defensefinder_matrix.tsv"
)

IS_FILE = (
    IS_MATRICES
    / "is_family_matrix.tsv"
)

# =====================================
# OUTPUT
# =====================================

OUT_STATS = (
    INTEGRATION_STATS
    / "defense_vs_is_total.tsv"
)

SCATTER_PLOT = (
    INTEGRATION_PLOTS
    / "defense_vs_is_scatter.png"
)

BOXPLOT_OUT = (
    INTEGRATION_PLOTS
    / "defense_vs_is_boxplots.png"
)

# =====================================
# LOAD
# =====================================

defense = pd.read_csv(
    DEFENSE_FILE,
    sep="\t"
)

is_df = pd.read_csv(
    IS_FILE,
    sep="\t"
)

defense["Genome"] = (
    defense["Genome"]
    .astype(str)
    .apply(normalize_genome_name)
)

is_df["Genome"] = (
    is_df["Genome"]
    .astype(str)
    .apply(normalize_genome_name)
)

# =====================================
# TOTAL COUNTS
# =====================================

defense["Total_defense"] = (
    defense
    .drop(columns=["Genome"])
    .sum(axis=1)
)

is_df["Total_IS"] = (
    is_df
    .drop(columns=["Genome"])
    .sum(axis=1)
)

# =====================================
# MERGE
# =====================================

merged = defense.merge(
    is_df[
        ["Genome", "Total_IS"]
    ],
    on="Genome"
)

print(
    f"Genomes after merge: {len(merged)}"
)

# =====================================
# CORRELATION
# =====================================

pearson_r, pearson_p = pearsonr(
    merged["Total_defense"],
    merged["Total_IS"]
)

spearman_rho, spearman_p = spearmanr(
    merged["Total_defense"],
    merged["Total_IS"]
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

# =====================================
# TOP DEFENSE SYSTEMS
# =====================================

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
        "Total_IS"
    ]

    absent = merged.loc[
        merged[system] == 0,
        "Total_IS"
    ]

    if (
        len(present) == 0
        or
        len(absent) == 0
    ):
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

# =====================================
# SCATTER
# =====================================

plt.figure(
    figsize=(7,6)
)

sns.regplot(
    data=merged,
    x="Total_defense",
    y="Total_IS",
    scatter_kws={"alpha":0.7}
)

plt.title(
    "Defense systems vs IS load"
)

plt.tight_layout()

plt.savefig(
    SCATTER_PLOT,
    dpi=300
)

plt.close()

# =====================================
# BOXPLOTS
# =====================================

plot_data = []

for system in systems:

    if system not in merged.columns:
        continue

    tmp = merged[
        [system, "Total_IS"]
    ].copy()

    tmp["System"] = system

    tmp["Presence"] = np.where(
        tmp[system] == 1,
        "Present",
        "Absent"
    )

    plot_data.append(
        tmp[
            ["System",
             "Presence",
             "Total_IS"]
        ]
    )

plot_df = pd.concat(
    plot_data,
    ignore_index=True
)

plt.figure(
    figsize=(10,6)
)

sns.boxplot(
    data=plot_df,
    x="System",
    y="Total_IS",
    hue="Presence"
)

plt.xticks(
    rotation=45,
    ha="right"
)

plt.tight_layout()

plt.savefig(
    BOXPLOT_OUT,
    dpi=300
)

plt.close()

print("\nSaved:")
print(OUT_STATS)
print(SCATTER_PLOT)
print(BOXPLOT_OUT)