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

from sklearn.decomposition import PCA

import matplotlib.pyplot as plt

# ======================================
# INPUT
# ======================================

MATRIX_FILE = (
    DEFENSEFINDER_MATRICES
    / "defensefinder_matrix.tsv"
)

# ======================================
# OUTPUT
# ======================================

DEFENSEFINDER_PLOTS.mkdir(
    parents=True,
    exist_ok=True
)

COORD_FILE = (
    DEFENSEFINDER_STATS
    / "defense_pca_coordinates.tsv"
)

VAR_FILE = (
    DEFENSEFINDER_STATS
    / "defense_pca_variance.tsv"
)

PLOT_FILE = (
    DEFENSEFINDER_PLOTS
    / "defense_pca.png"
)

# ======================================
# LOAD
# ======================================

df = pd.read_csv(
    MATRIX_FILE,
    sep="\t"
)

genomes = df["Genome"]

X = df.drop(
    columns=["Genome"]
)

# ======================================
# PCA
# ======================================

pca = PCA(
    n_components=2
)

coords = pca.fit_transform(X)

# ======================================
# SAVE COORDS
# ======================================

coord_df = pd.DataFrame({

    "Genome": genomes,

    "PC1": coords[:,0],

    "PC2": coords[:,1]

})

coord_df.to_csv(
    COORD_FILE,
    sep="\t",
    index=False
)

# ======================================
# SAVE VARIANCE
# ======================================

var_df = pd.DataFrame({

    "PC": ["PC1", "PC2"],

    "Variance_explained": (
        pca.explained_variance_ratio_
        * 100
    )

})

var_df.to_csv(
    VAR_FILE,
    sep="\t",
    index=False
)

# ======================================
# PLOT
# ======================================

plt.figure(
    figsize=(7,6)
)

plt.scatter(
    coords[:,0],
    coords[:,1],
    alpha=0.8
)

plt.xlabel(
    f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)"
)

plt.ylabel(
    f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)"
)

plt.title(
    "Defense system PCA"
)

plt.tight_layout()

plt.savefig(
    PLOT_FILE,
    dpi=300
)

plt.close()

print("\nSaved:")
print(COORD_FILE)
print(VAR_FILE)
print(PLOT_FILE)