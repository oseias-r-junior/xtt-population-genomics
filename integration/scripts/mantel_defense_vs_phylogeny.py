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

# =========================
# INPUT FILES
# =========================

DEFENSE_MATRIX = (
    INTEGRATION_MATRICES
    / "defensefinder_matrix.tsv"
)

TREE_FILE = (
    FINAL_TREES
    / "core_masked_mild.treefile"
)

# =========================
# LOAD DEFENSE MATRIX
# =========================

df = pd.read_csv(
    DEFENSE_MATRIX,
    sep="\t"
)

df["Genome"] = (
    df["Genome"]
    .astype(str)
    .apply(normalize_genome_name)
)

feature_cols = [
    c for c in df.columns
    if c != "Genome"
]

df = df.set_index("Genome")

# =========================
# DEFENSE DISTANCES
# =========================

defense_dist = squareform(
    pdist(
        df[feature_cols],
        metric="jaccard"
    )
)

genomes = list(df.index)

defense_dist_df = pd.DataFrame(
    defense_dist,
    index=genomes,
    columns=genomes
)

# =========================
# LOAD TREE
# =========================

tree = Phylo.read(
    TREE_FILE,
    "newick"
)

for tip in tree.get_terminals():

    tip.name = normalize_genome_name(
        tip.name
    )

tree_taxa = [
    x.name
    for x in tree.get_terminals()
]

# =========================
# DIAGNOSTICS
# =========================

matrix_only = sorted(
    set(genomes) - set(tree_taxa)
)

tree_only = sorted(
    set(tree_taxa) - set(genomes)
)

print(
    f"\nGenomes only in matrix: {len(matrix_only)}"
)

for g in matrix_only:
    print(g)

print(
    f"\nGenomes only in tree: {len(tree_only)}"
)

for g in tree_only:
    print(g)

# =========================
# COMMON GENOMES
# =========================

common = sorted(
    set(genomes).intersection(tree_taxa)
)

print(
    f"\nCommon genomes: {len(common)}"
)

# =========================
# SUBSET MATRICES
# =========================

defense_dist_df = defense_dist_df.loc[
    common,
    common
]

# =========================
# PHYLO DISTANCES
# =========================

phylo_dist = np.zeros(
    (len(common), len(common))
)

for i, g1 in enumerate(common):

    for j, g2 in enumerate(common):

        phylo_dist[i, j] = tree.distance(
            g1,
            g2
        )

phylo_dist_df = pd.DataFrame(
    phylo_dist,
    index=common,
    columns=common
)

# =========================
# CHECKS
# =========================

assert (
    defense_dist_df.shape
    ==
    phylo_dist_df.shape
)

assert all(
    defense_dist_df.index
    ==
    phylo_dist_df.index
)

# =========================
# FLATTEN MATRICES
# =========================

defense_vals = defense_dist_df.values[
    np.triu_indices_from(
        defense_dist_df,
        k=1
    )
]

phylo_vals = phylo_dist_df.values[
    np.triu_indices_from(
        phylo_dist_df,
        k=1
    )
]

# =========================
# OBSERVED CORRELATION
# =========================

obs_r, _ = pearsonr(
    phylo_vals,
    defense_vals
)

# =========================
# PERMUTATION TEST
# =========================

n_perm = 999

perm_rs = []

for i in range(n_perm):

    shuffled = np.random.permutation(
        defense_vals
    )

    r, _ = pearsonr(
        phylo_vals,
        shuffled
    )

    perm_rs.append(r)

perm_rs = np.array(
    perm_rs
)

p_value = np.mean(
    np.abs(perm_rs)
    >=
    np.abs(obs_r)
)

# =========================
# OUTPUT TABLE
# =========================

INTEGRATION_STATS.mkdir(
    parents=True,
    exist_ok=True
)

RESULT_FILE = (
    INTEGRATION_STATS
    / "defense_phylogeny_association.tsv"
)

pd.DataFrame({
    "Metric": [
        "Pearson_r",
        "P_value",
        "N_genomes"
    ],
    "Value": [
        obs_r,
        p_value,
        len(common)
    ]
}).to_csv(
    RESULT_FILE,
    sep="\t",
    index=False
)

# =========================
# PLOT
# =========================

plt.figure(
    figsize=(7,6)
)

plt.scatter(
    phylo_vals,
    defense_vals,
    alpha=0.6
)

plt.xlabel(
    "Phylogenetic distance"
)

plt.ylabel(
    "Defense composition distance"
)

plt.title(
    "Defense composition vs phylogeny"
)

plt.text(
    0.05,
    0.95,
    f"r = {obs_r:.3f}\np = {p_value:.5f}",
    transform=plt.gca().transAxes,
    verticalalignment="top"
)

plt.tight_layout()

INTEGRATION_PLOTS.mkdir(
    parents=True,
    exist_ok=True
)

PLOT_FILE = (
    INTEGRATION_PLOTS
    / "mantel_defense_vs_phylogeny.png"
)

plt.savefig(
    PLOT_FILE,
    dpi=300
)

plt.close()

# =========================
# REPORT
# =========================

print("\n===== Mantel-like test =====")

print(
    f"Observed r = {obs_r:.4f}"
)

print(
    f"P-value = {p_value:.5f}"
)

print("\nSaved:")

print(PLOT_FILE)

print(RESULT_FILE)