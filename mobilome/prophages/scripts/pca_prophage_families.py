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
import matplotlib.pyplot as plt

from sklearn.decomposition import PCA

# ======================================
# INPUT
# ======================================

MATRIX_FILE = (
    PROPHAGE_MATRICES
    / "prophage_family_presence_absence.tsv"
)

CLADE_FILE = GENOME_CLADES

# ======================================
# OUTPUT
# ======================================

OUTFIG = (
    PROPHAGE_PLOTS
    / "pca_prophage_families.png"
)

OUTTABLE = (
    PROPHAGE_STATS
    / "pca_prophage_families_coordinates.tsv"
)

PROPHAGE_PLOTS.mkdir(
    parents=True,
    exist_ok=True
)

PROPHAGE_STATS.mkdir(
    parents=True,
    exist_ok=True
)

# ======================================
# LOAD
# ======================================

matrix = pd.read_csv(
    MATRIX_FILE,
    sep="\t"
)

clades = pd.read_csv(
    CLADE_FILE,
    sep="\t"
)

# ======================================
# NORMALIZE
# ======================================

matrix["Genome"] = (
    matrix["Genome"]
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

df = matrix.merge(
    clades,
    on="Genome",
    how="inner"
)

print(
    f"Merged genomes: {len(df)}"
)

# ======================================
# PCA INPUT
# ======================================

X = df.drop(
    columns=[
        "Genome",
        "Clade"
    ]
)

# ======================================
# PCA
# ======================================

pca = PCA(
    n_components=2
)

coords = pca.fit_transform(X)

print(
    pca.explained_variance_ratio_
)

pca_df = pd.DataFrame({

    "Genome":
        df["Genome"],

    "Clade":
        df["Clade"],

    "PC1":
        coords[:, 0],

    "PC2":
        coords[:, 1]

})

# ======================================
# SAVE COORDINATES
# ======================================

pca_df.to_csv(
    OUTTABLE,
    sep="\t",
    index=False
)

print(
    f"Saved: {OUTTABLE}"
)

# ======================================
# PLOT
# ======================================

markers = {

    "K0": "o",
    "K1": "s",
    "K2": "^",
    "Xtc": "D"

}

fig, ax = plt.subplots(
    figsize=(8, 6)
)

for clade in sorted(
    pca_df["Clade"].unique()
):

    subset = pca_df[
        pca_df["Clade"] == clade
    ]

    ax.scatter(

        subset["PC1"],
        subset["PC2"],

        marker=markers.get(
            clade,
            "o"
        ),

        label=clade

    )

pc1_var = (
    pca.explained_variance_ratio_[0]
    * 100
)

pc2_var = (
    pca.explained_variance_ratio_[1]
    * 100
)

ax.set_xlabel(
    f"PC1 ({pc1_var:.1f}%)"
)

ax.set_ylabel(
    f"PC2 ({pc2_var:.1f}%)"
)

ax.set_title(
    "PCA of prophage family profiles"
)

ax.legend(
    title="Clade"
)

plt.tight_layout()

plt.savefig(
    OUTFIG,
    dpi=300,
    bbox_inches="tight"
)

plt.close()

print(
    f"Saved: {OUTFIG}"
)

print("\nFinished.")