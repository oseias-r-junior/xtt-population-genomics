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
from sklearn.preprocessing import StandardScaler

# =========================
# INPUTS
# =========================

IS_MATRIX = (
    IS_MATRICES
    / "is_family_matrix.tsv"
)

CLADE_FILE = GENOME_CLADES

# =========================
# LOAD DATA
# =========================

df = pd.read_csv(IS_MATRIX, sep="\t")

clades = pd.read_csv(CLADE_FILE, sep="\t")

df["Genome"] = (
    df["Genome"]
    .astype(str)
    .apply(normalize_genome_name)
)

clades["Genome"] = (
    clades["Genome"]
    .astype(str)
    .apply(normalize_genome_name)
)

df = df.merge(
    clades,
    on="Genome",
    how="inner"
)

print(f"Genomes after merge: {len(df)}")

families = [c for c in df.columns if c not in ["Genome", "Total_IS", "Clade"]]

X = df[families]

# =========================
# PCA
# =========================

scaler = StandardScaler()

X_scaled = scaler.fit_transform(X)

pca = PCA(n_components=2)

coords = pca.fit_transform(X_scaled)

df["PC1"] = coords[:,0]
df["PC2"] = coords[:,1]

# =========================
# PLOT
# =========================

plt.figure(figsize=(8,7))

for clade in sorted(df["Clade"].unique()):

    sub = df[df["Clade"] == clade]

    plt.scatter(
        sub["PC1"],
        sub["PC2"],
        s=80,
        alpha=0.8,
        label=clade
    )

# Labels
plt.xlabel(
    f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}% variance)"
)

plt.ylabel(
    f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}% variance)"
)

plt.title("PCA of IS composition colored by phylogenetic clade")

plt.legend()

plt.tight_layout()

IS_PLOTS.mkdir(
    parents=True,
    exist_ok=True
)

PLOT_OUTFILE = (
    IS_PLOTS
    / "pca_is_clades.png"
)

plt.savefig(
    PLOT_OUTFILE,
    dpi=300,
    bbox_inches="tight"
)

print("\nSaved:")
print("pca_is_clades.png")