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

from scipy.stats import (
    pearsonr,
    spearmanr
)

import matplotlib.pyplot as plt
import seaborn as sns

# =====================================
# INPUT
# =====================================

CENTRALITY_FILE = (
    INTEGRATION_STATS
    / "cooccurrence"
    / "defense_cooccurrence_centrality.tsv"
)

MOBILOME_FILE = (
    INTEGRATION_STATS
    / "clade_corrected"
    / "defense_vs_mobilome.tsv"
)

# =====================================
# OUTPUT
# =====================================

OUT_STATS = (
    INTEGRATION_STATS
    / "network"
    / "defense_vs_mobilome_centrality.tsv"
)

OUT_PLOT = (
    INTEGRATION_PLOTS
    / "defense_vs_mobilome_centrality.png"
)

OUT_STATS.parent.mkdir(
    parents=True,
    exist_ok=True
)

# =====================================
# LOAD
# =====================================

centrality = pd.read_csv(
    CENTRALITY_FILE,
    sep="\t"
)

mobilome = pd.read_csv(
    MOBILOME_FILE,
    sep="\t"
)

# =====================================
# MERGE
# =====================================

df = centrality.merge(
    mobilome,
    left_on="Subtype",
    right_on="Defense_system"
)

print(
    f"Systems after merge: {len(df)}"
)

# =====================================
# CORRELATIONS
# =====================================

results = []

for metric in [
    "Degree_centrality",
    "Betweenness_centrality",
    "Closeness_centrality",
    "Eigenvector_centrality"
]:

    pearson_r, pearson_p = pearsonr(
        df[metric],
        df["Mobilome_index_effect"]
    )

    spearman_r, spearman_p = spearmanr(
        df[metric],
        df["Mobilome_index_effect"]
    )

    results.append({
        "Metric": metric,
        "Pearson_r": pearson_r,
        "Pearson_p": pearson_p,
        "Spearman_r": spearman_r,
        "Spearman_p": spearman_p
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
# PLOT
# =====================================

fig, axes = plt.subplots(
    2,
    2,
    figsize=(10, 8)
)

metrics = [
    "Degree_centrality",
    "Betweenness_centrality",
    "Closeness_centrality",
    "Eigenvector_centrality"
]

for ax, metric in zip(
    axes.flatten(),
    metrics
):

    sns.regplot(
        data=df,
        x=metric,
        y="Mobilome_index_effect",
        ax=ax,
        scatter_kws={"s": 60}
    )

    r, p = spearmanr(
        df[metric],
        df["Mobilome_index_effect"]
    )

    ax.set_title(
        f"{metric}\n"
        f"Spearman r={r:.2f}, p={p:.3f}"
    )

plt.tight_layout()

plt.savefig(
    OUT_PLOT,
    dpi=300
)

plt.close()

# =====================================
# REPORT
# =====================================

print("\nSaved:")
print(OUT_STATS)
print(OUT_PLOT)