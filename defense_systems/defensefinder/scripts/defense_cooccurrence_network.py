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

import networkx as nx
from networkx.algorithms.community import greedy_modularity_communities
import matplotlib.pyplot as plt

# =====================================
# INPUT
# =====================================

COOC_FILE = (
    DEFENSEFINDER_STATS
    / "cooccurrence"
    / "defense_cooccurrence_significant.tsv"
)

# =====================================
# OUTPUT
# =====================================

OUT_PLOT = (
    DEFENSEFINDER_PLOTS
    / "defense_cooccurrence_network.png"
)

OUT_EDGES = (
    DEFENSEFINDER_STATS
    / "cooccurrence"
    / "defense_cooccurrence_edges.tsv"
)

# =====================================
# LOAD
# =====================================

df = pd.read_csv(
    COOC_FILE,
    sep="\t"
)

# =====================================
# PREPARE EDGE WEIGHTS
# =====================================

finite_max = (
    df.loc[
        np.isfinite(df["Odds_ratio"]),
        "Odds_ratio"
    ]
    .max()
)

# df["Odds_ratio_plot"] = (
#     df["Odds_ratio"]
#     .replace(np.inf, finite_max)
# )

# df["Weight"] = np.log2(
#     df["Odds_ratio_plot"]
# )

# =====================================
# FILTER
# =====================================

df = df[
    (df["FDR"] < 0.01)
    &
    (df["Odds_ratio"] > 2)
].copy()

df["Odds_ratio_plot"] = df["Odds_ratio"].replace(
    np.inf,
    df.loc[
        np.isfinite(df["Odds_ratio"]),
        "Odds_ratio"
    ].max()
)

df["Weight"] = np.log2(
    df["Odds_ratio_plot"]
)

print(
    f"Edges retained: {len(df)}"
)

nodes = (
    set(df["System_A"])
    |
    set(df["System_B"])
)

print(
    f"Nodes retained: {len(nodes)}"
)

df.to_csv(
    OUT_EDGES,
    sep="\t",
    index=False
)

# =====================================
# GRAPH
# =====================================

G = nx.Graph()

for _, row in df.iterrows():

    G.add_edge(
        row["System_A"],
        row["System_B"],
        weight=row["Weight"]
    )

# =====================================
# NODE SIZE
# =====================================

degrees = dict(
    G.degree()
)

node_sizes = [
    500 + degrees[n] * 250
    for n in G.nodes()
]

# =====================================
# EDGE WIDTH
# =====================================

edge_widths = [
    G[u][v]["weight"]
    for u, v in G.edges()
]

# opcional: ampliar visualmente

edge_widths = [
    w * 1.5
    for w in edge_widths
]

# =====================================
# LAYOUT
# =====================================

pos = nx.spring_layout(
    G,
    seed=42,
    k=1.2
)

# =====================================
# PLOT
# =====================================

plt.figure(
    figsize=(12, 10)
)

nx.draw_networkx_nodes(
    G,
    pos,
    node_size=node_sizes,
)

nx.draw_networkx_edges(
    G,
    pos,
    width=edge_widths,
    alpha=0.7
)

# weights = [
#     d["Weight"]
#     for _,_,d in G.edges(data=True)
# ]


nx.draw_networkx_labels(
    G,
    pos,
    font_size=9
)

print(df["Weight"].describe())
print(np.isinf(df["Weight"]).sum())

plt.title(
    "Defense system co-occurrence network"
)

plt.axis("off")

plt.tight_layout()

plt.savefig(
    OUT_PLOT,
    dpi=300,
    bbox_inches="tight"
)

plt.close()

print("\nSaved:")
print(OUT_PLOT)
print(OUT_EDGES)