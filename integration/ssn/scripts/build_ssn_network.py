#!/usr/bin/env python3
"""
integration/ssn/scripts/build_ssn_network.py

Build the MGE Sequence Similarity Network (SSN), mirroring Uceda-Campos et
al. 2022 (Microorganisms) Fig. 5: nodes are individual MGE instances
(prophages, genomic islands, insertion sequences, plasmids), edges are
BLASTN nucleotide similarity hits already filtered to identity >50% /
query coverage >80% (see run_blastn_ssn.slurm + awk filter).

Steps:
    1. Load node metadata (ssn_nodes.tsv, from extract_mge_sequences.py).
    2. Load filtered BLASTN hits (directional; A-vs-B and B-vs-A recorded
       separately) and collapse into undirected edges, keeping the best
       (max) identity/bitscore/coverage per unordered node pair.
    3. Add an "IS_carrier_flag" node attribute to PPH/GI/PLS nodes: True if
       an IS instance's coordinates fall (with any overlap) inside that
       MGE's boundaries, on the same genome and contig.
    4. Assign group labels (PPH-G, GI-G, IS-G, PLS-G) via connected
       components computed SEPARATELY per category (i.e. only counting
       edges between two nodes of the same category) -- matching the
       paper's use of per-class groups even though the full network shown
       includes cross-category edges too.
    5. Export: GraphML (for Cytoscape, matching the paper's own rendering
       tool) + flat nodes/edges TSVs + summary stats.

NOTE: with ~18k nodes, a full matplotlib static render would be an
unreadable hairball. This script exports the network for Cytoscape
(GraphML) rather than attempting an inline static plot.
"""

from pathlib import Path
import sys

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parents[2]

sys.path.append(str(PROJECT_ROOT / "scripts"))

from config import *

import pandas as pd
import networkx as nx

# ======================================================
# INPUT
# ======================================================

NODES_FILE = SSN_MATRICES / "ssn_nodes.tsv"
BLAST_FILE = SSN_BLAST / "all_vs_all_filtered.tsv"

BLAST_COLUMNS = [
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore",
    "qlen", "slen", "qcovhsp",
]

# ======================================================
# OUTPUT
# ======================================================

SSN_MATRICES.mkdir(parents=True, exist_ok=True)

OUT_GRAPHML = SSN_MATRICES / "mge_ssn.graphml"
OUT_NODES = SSN_MATRICES / "ssn_nodes_final.tsv"
OUT_EDGES = SSN_MATRICES / "ssn_edges_final.tsv"

# ======================================================
# 1. LOAD NODES
# ======================================================

nodes_df = pd.read_csv(NODES_FILE, sep="\t")
nodes_df = nodes_df.set_index("NodeID", drop=False)

print(f"Nodes loaded: {len(nodes_df)}")

# ======================================================
# 2. LOAD + COLLAPSE EDGES
# ======================================================

print("Loading BLASTN filtered hits (this may take a minute)...")
hits = pd.read_csv(BLAST_FILE, sep="\t", header=None, names=BLAST_COLUMNS)
print(f"Directional hits loaded: {len(hits)}")

# canonical unordered pair
a = hits["qseqid"].where(hits["qseqid"] < hits["sseqid"], hits["sseqid"])
b = hits["qseqid"].where(hits["qseqid"] < hits["sseqid"], hits["qseqid"])
b = hits["sseqid"].where(hits["qseqid"] < hits["sseqid"], hits["qseqid"])

hits["node_a"] = a
hits["node_b"] = b

edges_df = (
    hits
    .groupby(["node_a", "node_b"], as_index=False)
    .agg(
        pident=("pident", "max"),
        bitscore=("bitscore", "max"),
        qcovhsp=("qcovhsp", "max"),
        length=("length", "max"),
    )
)

print(f"Undirected edges after collapsing: {len(edges_df)}")

# ======================================================
# 3. BUILD GRAPH
# ======================================================

G = nx.Graph()

for node_id, row in nodes_df.iterrows():
    G.add_node(
        node_id,
        Category=row["Category"],
        Genome=row["Genome"],
        Contig=str(row["Contig"]),
        Start=row["Start"],
        End=row["End"],
        Length_bp=int(row["Length_bp"]),
        Classification=str(row["Classification"]),
        N_Tools_Support=str(row["N_Tools_Support"]),
        Tools_Supporting=str(row["Tools_Supporting"]),
        Extra=str(row["Extra"]),
        IS_carrier_flag=False,
        Group="",
    )

# Cross-category edges (e.g. a prophage sequence hitting a free-standing IS
# because the prophage's own sequence contains an embedded IS) are excluded
# here. This information is redundant with IS_carrier_flag (already shown
# as a ring on the host node) and was not part of the original Uceda-Campos
# et al. 2022 SSN design, which never compares sequence similarity across
# MGE categories -- only within each category.
def _category_of(node_id):
    return node_id.split("|")[0]

cross_category_dropped = 0
for _, row in edges_df.iterrows():
    if row["node_a"] in G and row["node_b"] in G:
        if _category_of(row["node_a"]) != _category_of(row["node_b"]):
            cross_category_dropped += 1
            continue
        G.add_edge(
            row["node_a"], row["node_b"],
            pident=float(row["pident"]),
            bitscore=float(row["bitscore"]),
            qcovhsp=float(row["qcovhsp"]),
            alignment_length=int(row["length"]),
        )

print(f"Cross-category edges dropped (redundant with IS_carrier_flag): {cross_category_dropped}")

print(f"Graph: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")

# ======================================================
# 4. IS_CARRIER_FLAG (coordinate overlap, same genome+contig)
# ======================================================

is_nodes = nodes_df[nodes_df["Category"] == "IS"].copy()
is_nodes["Start"] = pd.to_numeric(is_nodes["Start"], errors="coerce")
is_nodes["End"] = pd.to_numeric(is_nodes["End"], errors="coerce")

# index IS instances by (Genome, Contig) for fast lookup
is_by_genome_contig = {}
for node_id, row in is_nodes.iterrows():
    key = (row["Genome"], str(row["Contig"]))
    is_by_genome_contig.setdefault(key, []).append((row["Start"], row["End"]))

carrier_count = 0
for node_id, row in nodes_df.iterrows():
    if row["Category"] == "IS":
        continue
    key = (row["Genome"], str(row["Contig"]))
    is_list = is_by_genome_contig.get(key, [])
    mge_start, mge_end = row["Start"], row["End"]
    for is_start, is_end in is_list:
        if pd.isna(is_start) or pd.isna(is_end):
            continue
        overlap = min(mge_end, is_end) - max(mge_start, is_start)
        if overlap > 0:
            G.nodes[node_id]["IS_carrier_flag"] = True
            carrier_count += 1
            break

print(f"MGE nodes flagged as IS-carrying: {carrier_count}")

# ======================================================
# 4b. EXCLUDE IS INSTANCES CONTAINED WITHIN A PPH/GI/PLS NODE
# ======================================================
# Rationale: if an IS is already represented via the IS_carrier_flag ring on
# its host MGE, keeping it ALSO as an independent IS node in the network is
# redundant (double counting the same biological event) and inflates the IS
# category with instances that provide no new information. Any IS instance
# overlapping (same Genome+Contig, any overlap) a PPH/GI/PLS node already in
# the graph is removed from G entirely (edges touching it are dropped too).
host_nodes = nodes_df[nodes_df["Category"].isin(["PPH", "GI", "PLS"])].copy()
host_by_genome_contig = {}
for node_id, row in host_nodes.iterrows():
    key = (row["Genome"], str(row["Contig"]))
    host_by_genome_contig.setdefault(key, []).append((row["Start"], row["End"]))

is_to_remove = []
for node_id, row in is_nodes.iterrows():
    key = (row["Genome"], str(row["Contig"]))
    host_list = host_by_genome_contig.get(key, [])
    is_start, is_end = row["Start"], row["End"]
    if pd.isna(is_start) or pd.isna(is_end):
        continue
    for h_start, h_end in host_list:
        overlap = min(is_end, h_end) - max(is_start, h_start)
        if overlap > 0:
            is_to_remove.append(node_id)
            break

is_to_remove = [n for n in is_to_remove if n in G]
G.remove_nodes_from(is_to_remove)
print(f"IS nodes excluded (contained within a PPH/GI/PLS node): {len(is_to_remove)} of {len(is_nodes)}")
print(f"Graph after IS exclusion: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")

# ======================================================
# 5. GROUPS (connected components, per-category subgraphs)
# ======================================================

category_prefix = {"PPH": "PPH-G", "GI": "GI-G", "IS": "IS-G", "PLS": "PLS-G"}

for category, prefix in category_prefix.items():
    cat_nodes = [n for n, d in G.nodes(data=True) if d["Category"] == category]
    subG = G.subgraph(cat_nodes)

    components = list(nx.connected_components(subG))
    # sort largest-first for stable, readable numbering
    components.sort(key=len, reverse=True)

    for i, comp in enumerate(components, start=1):
        group_label = f"{prefix}{i}"
        for node_id in comp:
            G.nodes[node_id]["Group"] = group_label

    print(f"{category}: {len(cat_nodes)} nodes -> {len(components)} groups "
          f"(largest: {len(components[0]) if components else 0} nodes)")

# ======================================================
# 6. EXPORT
# ======================================================

nx.write_graphml(G, OUT_GRAPHML)

final_nodes = pd.DataFrame([
    {"NodeID": n, **d} for n, d in G.nodes(data=True)
])
final_nodes.to_csv(OUT_NODES, sep="\t", index=False)

final_edges = pd.DataFrame([
    {"Source": u, "Target": v, **d} for u, v, d in G.edges(data=True)
])
final_edges.to_csv(OUT_EDGES, sep="\t", index=False)

print()
print(f"Saved: {OUT_GRAPHML}")
print(f"Saved: {OUT_NODES}")
print(f"Saved: {OUT_EDGES}")
