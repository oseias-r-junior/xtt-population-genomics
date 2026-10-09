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

from scipy.spatial.distance import pdist, squareform
from scipy.stats import pearsonr
from Bio import Phylo
import matplotlib.pyplot as plt
import random

# =========================
# INPUT FILES
# =========================

IS_MATRIX = (
    IS_MATRICES
    / "is_family_matrix.tsv"
)

TREE_FILE = (
    FINAL_TREES
    / "core_masked_mild.treefile"
)

# =========================
# LOAD IS MATRIX
# =========================

df = pd.read_csv(IS_MATRIX, sep="\t")

df["Genome"] = (
    df["Genome"]
    .astype(str)
    .apply(normalize_genome_name)
)

feature_cols = [
    c for c in df.columns
    if c not in ["Genome", "Total_IS"]
]

df = df.set_index("Genome")

# =========================
# IS DISTANCES
# =========================

is_dist = squareform(
    pdist(df[feature_cols], metric="braycurtis")
)

genomes = list(df.index)

is_dist_df = pd.DataFrame(
    is_dist,
    index=genomes,
    columns=genomes
)

# =========================
# LOAD TREE
# =========================

tree = Phylo.read(TREE_FILE, "newick")

for tip in tree.get_terminals():
    tip.name = normalize_genome_name(tip.name)

tree_taxa = [
    x.name
    for x in tree.get_terminals()
]

common = sorted(
    set(genomes).intersection(tree_taxa)
)

print(f"Common genomes: {len(common)}")

# =========================
# SUBSET IS DIST MATRIX
# =========================

is_dist_df = is_dist_df.loc[
    common,
    common
]

# assert is_dist_df.shape == phylo_dist_df.shape

# assert all(
#     is_dist_df.index
#     == phylo_dist_df.index
# )

# =========================
# PHYLO DISTANCES
# =========================

phylo_dist = np.zeros((len(common), len(common)))

for i, g1 in enumerate(common):
    for j, g2 in enumerate(common):
        phylo_dist[i, j] = tree.distance(g1, g2)

phylo_dist_df = pd.DataFrame(
    phylo_dist,
    index=common,
    columns=common
)

# =========================
# FLATTEN MATRICES
# =========================
assert is_dist_df.shape == phylo_dist_df.shape

assert all(
    is_dist_df.index
    == phylo_dist_df.index
)

is_vals = is_dist_df.values[np.triu_indices_from(is_dist_df, k=1)]

phylo_vals = phylo_dist_df.values[np.triu_indices_from(phylo_dist_df, k=1)]

# =========================
# OBSERVED CORRELATION
# =========================

assert is_dist_df.shape == phylo_dist_df.shape

assert all(
    is_dist_df.index
    == phylo_dist_df.index
)
obs_r, _ = pearsonr(phylo_vals, is_vals)

# =========================
# PERMUTATION TEST
# =========================

n_perm = 999

perm_rs = []

for i in range(n_perm):

    shuffled = np.random.permutation(is_vals)

    r, _ = pearsonr(phylo_vals, shuffled)

    perm_rs.append(r)

perm_rs = np.array(perm_rs)

p_value = np.mean(np.abs(perm_rs) >= np.abs(obs_r))

# =========================
# OUTPUT
# =========================

print("\n===== Mantel-like test =====")

print(f"Observed r = {obs_r:.4f}")

print(f"P-value = {p_value:.5f}")

# =========================
# PLOT
# =========================

plt.figure(figsize=(7,6))

plt.scatter(
    phylo_vals,
    is_vals,
    alpha=0.6
)

plt.xlabel("Phylogenetic distance")

plt.ylabel("IS composition distance")

plt.title("IS composition vs phylogeny")

plt.text(
    0.05,
    0.95,
    f"r = {obs_r:.3f}\np = {p_value:.5f}",
    transform=plt.gca().transAxes,
    verticalalignment='top'
)

plt.tight_layout()

IS_PLOTS.mkdir(
    parents=True,
    exist_ok=True
)

OUTPUT_PLOT = (
    IS_PLOTS
    / "mantel_is_vs_phylogeny.png"
)

plt.savefig(
    OUTPUT_PLOT,
    dpi=300
)


print("\nSaved:")
print("mantel_is_vs_phylogeny.png")