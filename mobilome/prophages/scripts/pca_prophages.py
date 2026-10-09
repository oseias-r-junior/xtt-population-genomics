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
import matplotlib.pyplot as plt
import seaborn as sns

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

# ======================================
# INPUT
# ======================================

INPUT_MATRIX = (
    PROPHAGE_MATRICES
    / "prophage_summary.tsv"
)

CLADE_FILE = GENOME_CLADES

# ======================================
# OUTPUT
# ======================================

PROPHAGE_PLOTS.mkdir(
    parents=True,
    exist_ok=True
)

PLOT_OUTFILE = (
    PROPHAGE_PLOTS
    / "pca_prophages.png"
)

# ======================================
# LOAD DATA
# ======================================

df = pd.read_csv(
    INPUT_MATRIX,
    sep="\t"
)

clades = pd.read_csv(
    CLADE_FILE,
    sep="\t"
)

merged = pd.merge(
    df,
    clades,
    on="Genome",
    how="left"
)

# ======================================
# FEATURES
# ======================================

feature_cols = [

    "Total_Prophages",

    "Intact",

    "Questionable",

    "Incomplete"

]

X = merged[
    feature_cols
].fillna(0)


n_unique = X.drop_duplicates().shape[0]

if n_unique < 3:

    print(
        "Not enough variation for PCA."
    )

    sys.exit()

# ======================================
# PCA
# ======================================

scaler = StandardScaler()

X_scaled = scaler.fit_transform(X)

pca = PCA(
    n_components=2
)

coords = pca.fit_transform(
    X_scaled
)

pca_df = pd.DataFrame({

    "Genome":
        merged["Genome"],

    "PC1":
        coords[:, 0],

    "PC2":
        coords[:, 1],

    "Clade":
        merged["Clade"]

})

# ======================================
# PLOT
# ======================================

plt.figure(
    figsize=(8, 6)
)

sns.scatterplot(
    data=pca_df,
    x="PC1",
    y="PC2",
    hue="Clade",
    s=100
)

plt.title(
    "PCA of prophage composition"
)

plt.tight_layout()

plt.savefig(
    PLOT_OUTFILE,
    dpi=300,
    bbox_inches="tight"
)

plt.close()

print(
    f"Saved: {PLOT_OUTFILE}"
)