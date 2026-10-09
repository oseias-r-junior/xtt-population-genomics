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

from sklearn.decomposition import PCA

# ======================================
# INPUT
# ======================================

MATRIX_FILE = (
    PROPHAGE_MATRICES
    / "prophage_family_presence_absence.tsv"
)

# ======================================
# OUTPUT
# ======================================

OUTFILE = (
    PROPHAGE_STATS
    / "prophage_pca_loadings.tsv"
)

# ======================================
# LOAD
# ======================================

df = pd.read_csv(
    MATRIX_FILE,
    sep="\t"
)

# ======================================
# NORMALIZE
# ======================================

df["Genome"] = (
    df["Genome"]
    .astype(str)
    .apply(normalize_genome_name)
)

# ======================================
# PCA INPUT
# ======================================

X = df.drop(
    columns=["Genome"]
)

families = X.columns

# ======================================
# PCA
# ======================================

pca = PCA(
    n_components=2
)

pca.fit(X)

# ======================================
# LOADINGS
# ======================================

loadings = pd.DataFrame(

    pca.components_.T,

    index=families,

    columns=[
        "PC1_loading",
        "PC2_loading"
    ]

)

# ======================================
# ABSOLUTE CONTRIBUTION
# ======================================

loadings["PC1_abs"] = (
    loadings["PC1_loading"]
    .abs()
)

loadings["PC2_abs"] = (
    loadings["PC2_loading"]
    .abs()
)

# ======================================
# SORT
# ======================================

loadings = loadings.sort_values(

    "PC1_abs",

    ascending=False

)

# ======================================
# SAVE
# ======================================

loadings = loadings.reset_index()

loadings.columns = [
    "Family",
    "PC1_loading",
    "PC2_loading",
    "PC1_abs",
    "PC2_abs"
]

loadings.to_csv(

    OUTFILE,
    sep="\t"

)

print(
    f"Saved: {OUTFILE}"
)

print("\nTop PC1 drivers:")

print(
    loadings[
        [
            "PC1_loading",
            "PC2_loading"
        ]
    ]
    .head(10)
)

print("\nExplained variance:")

print(
    pca.explained_variance_ratio_
)

print("\nFinished.")