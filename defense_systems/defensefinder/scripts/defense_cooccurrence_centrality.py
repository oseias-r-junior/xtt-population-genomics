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
import networkx as nx

# =====================================
# INPUT
# =====================================

INPUT_FILE = (
    DEFENSEFINDER_STATS
    / "cooccurrence"
    / "defense_cooccurrence_edges.tsv"
)

# =====================================
# OUTPUT
# =====================================

OUTFILE = (
    DEFENSEFINDER_STATS
    / "cooccurrence"
    / "defense_cooccurrence_centrality.tsv"
)

# =====================================
# LOAD
# =====================================

edges = pd.read_csv(
    INPUT_FILE,
    sep="\t"
)

print(
    f"Edges: {len(edges)}"
)

# =====================================
# BUILD GRAPH
# =====================================

G = nx.Graph()

for _, row in edges.iterrows():

    G.add_edge(
        row["System_A"],
        row["System_B"],
        weight=row["Weight"]
    )

print(
    f"Nodes: {G.number_of_nodes()}"
)

# =====================================
# CENTRALITIES
# =====================================

degree = nx.degree_centrality(G)

betweenness = nx.betweenness_centrality(
    G,
    weight="weight"
)

closeness = nx.closeness_centrality(G)

try:

    eigenvector = nx.eigenvector_centrality(
        G,
        weight="weight",
        max_iter=5000
    )

except Exception:

    eigenvector = {
        n: float("nan")
        for n in G.nodes()
    }

# =====================================
# TABLE
# =====================================

results = []

for node in G.nodes():

    results.append({
        "Subtype": node,
        "Degree_centrality":
            degree.get(node),

        "Betweenness_centrality":
            betweenness.get(node),

        "Closeness_centrality":
            closeness.get(node),

        "Eigenvector_centrality":
            eigenvector.get(node)
    })

results_df = pd.DataFrame(
    results
)

results_df = results_df.sort_values(
    "Eigenvector_centrality",
    ascending=False
)

results_df.to_csv(
    OUTFILE,
    sep="\t",
    index=False
)

print("\nTop hubs:")

print(
    results_df.head(10)
)

print("\nSaved:")
print(OUTFILE)