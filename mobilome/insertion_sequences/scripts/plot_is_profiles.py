import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

from pathlib import Path
import sys

SCRIPT_DIR = Path(__file__).resolve().parent

PROJECT_ROOT = SCRIPT_DIR.parents[2]

sys.path.append(
    str(PROJECT_ROOT / "scripts")
)

from config import *




# =========================
# Load matrix
# =========================

matrix_file = (
    IS_MATRICES
    / "is_family_matrix.tsv"
)

print(matrix_file)

df = pd.read_csv(
    matrix_file,
    sep="\t"
)

# =========================
# Prepare data
# =========================

families = [
    "IS110",
    "IS1595",
    "IS3",
    "IS4",
    "IS481",
    "IS5",
    "IS605"
]

X = df[families]

# =========================
# HEATMAP
# =========================

plt.figure(figsize=(10, 18))

sns.heatmap(
    X,
    cmap="viridis",
    yticklabels=df["Genome"]
)

plt.title("IS family composition across Xtt genomes")

plt.tight_layout()

IS_PLOTS.mkdir(
    parents=True,
    exist_ok=True
)

outfile = (
    IS_PLOTS
    / "is_family_heatmap.png"
)

plt.savefig(
    outfile,
    dpi=300,
    bbox_inches="tight"
)

print(outfile)

plt.close()

# =========================
# PCA
# =========================

scaler = StandardScaler()

X_scaled = scaler.fit_transform(X)

pca = PCA(n_components=2)

coords = pca.fit_transform(X_scaled)

pca_df = pd.DataFrame({
    "Genome": df["Genome"],
    "PC1": coords[:, 0],
    "PC2": coords[:, 1],
    "Total_IS": df["Total_IS"]
})

plt.figure(figsize=(8, 6))

sns.scatterplot(
    data=pca_df,
    x="PC1",
    y="PC2",
    size="Total_IS",
    sizes=(20, 300)
)

plt.title("PCA of IS family composition")

plt.tight_layout()


IS_PLOTS.mkdir(
    parents=True,
    exist_ok=True
)

outfile = (
    IS_PLOTS
    / "is_family_pca.png"
)

plt.savefig(
    outfile,
    dpi=300,
    bbox_inches="tight"
)

plt.close()

# =========================
# TOTAL IS BARPLOT
# =========================

df_sorted = df.sort_values("Total_IS", ascending=False)

plt.figure(figsize=(18, 6))

sns.barplot(
    data=df_sorted,
    x="Genome",
    y="Total_IS"
)

plt.xticks(rotation=90)

plt.title("Total IS counts across Xtt genomes")

plt.tight_layout()

IS_PLOTS.mkdir(
    parents=True,
    exist_ok=True
)

outfile = (
    IS_PLOTS
    / "is_total_barplot.png"
)

plt.savefig(
    outfile,
    dpi=300,
    bbox_inches="tight"
)

plt.close()

print("Plots generated successfully.")
