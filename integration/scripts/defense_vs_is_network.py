#!/usr/bin/env python3

from pathlib import Path
import sys

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parents[2]

sys.path.append(str(PROJECT_ROOT / "scripts"))

from config import *

import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

# =====================================
# INPUT
# =====================================

INFILE = (
    INTEGRATION_STATS
    / "clade_corrected"
    / "defense_vs_is_family_clade_corrected.tsv"
)

# =====================================
# OUTPUT
# =====================================

OUTDIR = (
    INTEGRATION_STATS
    / "networks"
)

OUTDIR.mkdir(
    parents=True,
    exist_ok=True
)

EDGE_FILE = (
    OUTDIR
    / "defense_vs_is_network_edges.tsv"
)

PLOT_FILE = (
    INTEGRATION_PLOTS
    / "defense_vs_is_network.png"
)

# =====================================
# LOAD
# =====================================

df = pd.read_csv(
    INFILE,
    sep="\t"
)

# =====================================
# FILTER
# =====================================

df = df[
    (df["P_value"] < 0.01)
    &
    (df["Effect_size"].abs() > 0.5)
].copy()

print(
    f"Retained associations: {len(df)}"
)

# =====================================
# SAVE EDGES
# =====================================

df.to_csv(
    EDGE_FILE,
    sep="\t",
    index=False
)

# =====================================
# GRAPH
# =====================================

G = nx.Graph()

for _, row in df.iterrows():

    defense = row["Defense_system"]
    is_family = row["IS_family"]

    G.add_node(
        defense,
        bipartite="Defense"
    )

    G.add_node(
        is_family,
        bipartite="IS"
    )

    G.add_edge(
        defense,
        is_family,
        weight=abs(row["Effect_size"]),
        sign=np.sign(row["Effect_size"])
    )

# =====================================
# LAYOUT
# =====================================

defense_nodes = [
    n
    for n, d in G.nodes(data=True)
    if d["bipartite"] == "Defense"
]

is_nodes = [
    n
    for n, d in G.nodes(data=True)
    if d["bipartite"] == "IS"
]

pos = {}

for i, node in enumerate(sorted(defense_nodes)):
    pos[node] = (0, i)

for i, node in enumerate(sorted(is_nodes)):
    pos[node] = (1, i)

# =====================================
# NODE SIZES
# =====================================

degrees = dict(G.degree())

node_sizes = [
    500 + degrees[n] * 250
    for n in G.nodes()
]

# =====================================
# EDGE COLORS
# =====================================

edge_colors = []

for u, v, d in G.edges(data=True):

    if d["sign"] > 0:
        edge_colors.append("darkred")
    else:
        edge_colors.append("steelblue")

# =====================================
# EDGE WIDTHS
# =====================================

edge_widths = [
    d["weight"] * 0.7
    for _, _, d in G.edges(data=True)
]

# =====================================
# PLOT
# =====================================

plt.figure(
    figsize=(12, 10)
)

nx.draw_networkx_nodes(
    G,
    pos,
    node_size=node_sizes
)

nx.draw_networkx_edges(
    G,
    pos,
    edge_color=edge_colors,
    width=edge_widths,
    alpha=0.8
)

nx.draw_networkx_labels(
    G,
    pos,
    font_size=8
)

plt.title(
    "Defense systems vs IS families network"
)

plt.axis("off")

plt.tight_layout()

plt.savefig(
    PLOT_FILE,
    dpi=300,
    bbox_inches="tight"
)

plt.close()

print("\nSaved:")
print(PLOT_FILE)
print(EDGE_FILE)