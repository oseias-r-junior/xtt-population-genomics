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
    / "defense_clustered_heatmap.png"
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
# REMOVE INVARIANT SYSTEMS
# =========================

df = df.loc[
    df.var(axis=1) > 0
]

print(
    f"Variable systems: {df.shape[0]}"
)

# =========================
# CLUSTERED HEATMAP
# =========================

g = sns.clustermap(
    df,
    cmap="viridis",
    linewidths=0.2,
    figsize=(10, 14),
    metric="euclidean",
    method="average",
    cbar_kws={
        "label": "Mean presence"
    }
)

g.fig.suptitle(
    "Defense systems clustered by clade",
    y=1.02
)

plt.savefig(
    OUTPUT_FILE,
    dpi=300,
    bbox_inches="tight"
)

plt.close()

print("\nSaved:")
print(OUTPUT_FILE)