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

import matplotlib.pyplot as plt
import seaborn as sns

from scipy.stats import pearsonr

# ======================================
# INPUT
# ======================================

PHASTEST_FILE = (
    PROPHAGE_MATRICES
    / "prophage_summary.tsv"
)

VIRSORTER_FILE = (
    VIRSORTER2_MATRICES
    / "virsorter2_summary_matrix.tsv"
)

# ======================================
# OUTPUT
# ======================================

OUT_STATS = (
    PROPHAGE_STATS
    / "virsorter2_vs_phastest.tsv"
)

OUT_SCATTER = (
    PROPHAGE_PLOTS
    / "virsorter2_vs_phastest_scatter.png"
)

OUT_BOXPLOT = (
    PROPHAGE_PLOTS
    / "virsorter2_vs_phastest_boxplot.png"
)

OUT_STATS.parent.mkdir(
    parents=True,
    exist_ok=True
)

OUT_SCATTER.parent.mkdir(
    parents=True,
    exist_ok=True
)

# ======================================
# LOAD
# ======================================

phastest = pd.read_csv(
    PHASTEST_FILE,
    sep="\t"
)

virsorter = pd.read_csv(
    VIRSORTER_FILE,
    sep="\t"
)

# ======================================
# NORMALIZE
# ======================================

phastest["Genome"] = (
    phastest["Genome"]
    .astype(str)
    .apply(normalize_genome_name)
)

virsorter["Genome"] = (
    virsorter["Genome"]
    .astype(str)
    .apply(normalize_genome_name)
)

# ======================================
# KEEP COLUMNS
# ======================================

phastest = phastest[
    [
        "Genome",
        "Total_Prophages",
        "Intact",
        "Questionable",
        "Incomplete"
    ]
]

virsorter = virsorter[
    [
        "Genome",
        "Total_sequences",
        "Total_length",
        "Mean_length"
    ]
]

# ======================================
# MERGE
# ======================================

df = (
    phastest
    .merge(
        virsorter,
        on="Genome",
        suffixes=(
            "_PHASTEST",
            "_VirSorter2"
        )
    )
)

print(
    f"Genomes: {len(df)}"
)

# ======================================
# AGREEMENT
# ======================================

df["Difference"] = (
    df["Total_sequences"]
    -
    df["Total_Prophages"]
)

df["Agreement"] = (
    df["Difference"]
    == 0
)

# ======================================
# CORRELATION
# ======================================

r, p = pearsonr(
    df["Total_Prophages"],
    df["Total_sequences"]
)

print(
    f"Pearson r = {r:.3f}"
)

print(
    f"P-value = {p:.3g}"
)

# ======================================
# SAVE TABLE
# ======================================

df.to_csv(
    OUT_STATS,
    sep="\t",
    index=False
)

# ======================================
# SCATTER
# ======================================

plt.figure(
    figsize=(6,6)
)

sns.regplot(
    data=df,
    x="Total_Prophages",
    y="Total_sequences",
    scatter_kws={
        "s":70
    }
)

plt.text(
    0.05,
    0.95,
    (
        f"r = {r:.2f}\n"
        f"p = {p:.3g}"
    ),
    transform=plt.gca().transAxes,
    verticalalignment="top",
    bbox=dict(
        facecolor="white",
        alpha=0.8
    )
)

plt.xlabel(
    "PHASTEST prophages"
)

plt.ylabel(
    "VirSorter2 sequences"
)

plt.title(
    "VirSorter2 vs PHASTEST"
)

plt.tight_layout()

plt.savefig(
    OUT_SCATTER,
    dpi=300
)

plt.close()

# ======================================
# BOXPLOT
# ======================================

plot_df = pd.DataFrame({

    "Detector":

        (
            ["PHASTEST"] * len(df)
        )
        +
        (
            ["VirSorter2"] * len(df)
        ),

    "Predicted_sequences":

        pd.concat(
            [
                df["Total_Prophages"],
                df["Total_sequences"]
            ],
            ignore_index=True
        )

})

plt.figure(
    figsize=(5,6)
)

sns.boxplot(
    data=plot_df,
    x="Detector",
    y="Predicted_sequences"
)

sns.stripplot(
    data=plot_df,
    x="Detector",
    y="Predicted_sequences",
    color="black",
    alpha=0.6,
    jitter=True
)

plt.ylabel(
    "Predicted prophages"
)

plt.xlabel("")

plt.title(
    "Predicted prophages per genome"
)

plt.tight_layout()

plt.savefig(
    OUT_BOXPLOT,
    dpi=300
)

plt.close()

# ======================================
# REPORT
# ======================================

print()

print(
    f"Agreement: {df['Agreement'].sum()} / {len(df)} genomes"
)

print(
    f"Mean difference: {df['Difference'].mean():.2f}"
)

print()

print(OUT_STATS)
print(OUT_SCATTER)
print(OUT_BOXPLOT)