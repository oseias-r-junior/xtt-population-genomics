#!/usr/bin/env python3

from pathlib import Path
import sys

SCRIPT_DIR = Path(__file__).resolve().parent

PROJECT_ROOT = SCRIPT_DIR.parents[2]

sys.path.append(
    str(PROJECT_ROOT / "scripts")
)

from config import *

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# =========================
# INPUT
# =========================

INPUT_FILE = (
    DEFENSEFINDER_STATS
    / "defense_clade_stats"
    / "mean_defense_by_clade.tsv"
)

OUTPUT_FILE = (
    DEFENSEFINDER_PLOTS
    / "defense_clade_means_heatmap.png"
)

DEFENSEFINDER_PLOTS.mkdir(
    parents=True,
    exist_ok=True
)

# =========================
# LOAD DATA
# =========================

df = pd.read_csv(
    INPUT_FILE,
    sep="\t",
    index_col=0
)

print(df.shape)

# =========================
# SORT SYSTEMS
# =========================

df = df.loc[
    df.max(axis=1)
      .sort_values(ascending=False)
      .index
]

# =========================
# HEATMAP
# =========================

plt.figure(
    figsize=(8, 12)
)

sns.heatmap(
    df,
    cmap="viridis",
    linewidths=0.5,
    cbar_kws={
        "label": "Mean presence"
    }
)

plt.title(
    "Defense systems prevalence by clade"
)

plt.xlabel(
    "Clade"
)

plt.ylabel(
    "Defense system"
)

plt.tight_layout()

plt.savefig(
    OUTPUT_FILE,
    dpi=300,
    bbox_inches="tight"
)

print("\nSaved:")
print(OUTPUT_FILE)