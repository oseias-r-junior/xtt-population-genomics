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

from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

# =========================
# INPUT
# =========================

INPUT_MATRIX = (
    IS_MATRICES
    / "is_family_matrix.tsv"
)

# =========================
# LOAD DATA
# =========================

df = pd.read_csv(INPUT_MATRIX, sep="\t")

families = [c for c in df.columns if c not in ["Genome", "Total_IS"]]

X = df[families]

# =========================
# STANDARDIZE
# =========================

scaler = StandardScaler()

X_scaled = scaler.fit_transform(X)

# =========================
# PCA
# =========================

pca = PCA(n_components=2)

coords = pca.fit_transform(X_scaled)

# =========================
# LOADINGS
# =========================

loadings = pd.DataFrame(
    pca.components_.T,
    columns=["PC1", "PC2"],
    index=families
)

print("\n===== PCA LOADINGS =====\n")

print(loadings)

# Save loadings
IS_STATS.mkdir(
    parents=True,
    exist_ok=True
)

IS_PLOTS.mkdir(
    parents=True,
    exist_ok=True
)

LOADINGS_OUTFILE = (
    IS_STATS
    / "pca_is_loadings.tsv"
)

PLOT_OUTFILE = (
    IS_PLOTS
    / "pca_is_family_loadings.png"
)

# =========================
# PLOT
# =========================

plt.figure(figsize=(8,8))

# Draw arrows
for family in families:

    x = loadings.loc[family, "PC1"]
    y = loadings.loc[family, "PC2"]

    plt.arrow(
        0,
        0,
        x,
        y,
        head_width=0.03,
        length_includes_head=True
    )

    plt.text(
        x * 1.1,
        y * 1.1,
        family,
        fontsize=11
    )

# Axis labels with explained variance
pc1_var = pca.explained_variance_ratio_[0] * 100
pc2_var = pca.explained_variance_ratio_[1] * 100

plt.xlabel(f"PC1 ({pc1_var:.1f}% variance)")
plt.ylabel(f"PC2 ({pc2_var:.1f}% variance)")

plt.title("PCA loadings of IS families")

plt.axhline(0, color='grey', linewidth=0.8)
plt.axvline(0, color='grey', linewidth=0.8)

plt.tight_layout()

loadings.to_csv(
    LOADINGS_OUTFILE,
    sep="\t",
    index=True
)

plt.savefig(
    PLOT_OUTFILE,
    dpi=300,
    bbox_inches="tight"
)

print("\nSaved:")
print("pca_is_family_loadings.png")
print("pca_is_loadings.tsv")