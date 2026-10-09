#!/usr/bin/env python3
"""
plot_mge_by_clade.py -- MGE composition per genome, grouped by Xtt clade.

Reproduces the style of Fig. 4 in the Uceda-Campos et al. Xylella paper
("Percentage of mobile genetic elements distributed among the strains..."):
one 100%-stacked vertical bar per genome, showing what fraction of that
genome's MGE complement is prophage / genomic island / insertion sequence /
plasmid. The original figure grouped genomes by assembly completeness
(complete/scaffold/contig); per the PI's explicit instruction, that axis is
replaced here with the phylogenetic clade of X. translucens pv. translucens
(K0 / K1 / K2) -- assembly completeness is NOT used anywhere in this figure.

Data sources (both confirmed by the user):
  - MGE counts per genome/category: integration/ssn/matrices/ssn_nodes_final.tsv
    (the SAME node set used in the SSN network figures -- i.e. AFTER the
    SSN-specific exclusions: IS contained within a PPH/GI/PLS host removed,
    GI overlapping a high-confidence PPH removed. This was an explicit
    methodological choice by the user, for consistency with the other MGE
    figures in the manuscript, not an oversight.)
  - Clade assignment per genome: metadata/genome_clades_master.tsv
    (only the Genome and Clade columns are used -- the rest of that table,
    e.g. Genome_Length/QC columns, is out of scope for this figure per the
    user's instruction not to reproduce it "ipsis literis").

Genomes present in ssn_nodes_final.tsv but absent from the clade table (or
vice-versa) are reported and excluded, rather than silently dropped.

Output:
    integration/clade_correlation/plots/mge_distribution_by_clade.png

Usage
-----
python plot_mge_by_clade.py
"""
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
PROJECT_ROOT = Path(__file__).resolve().parents[3]
NODES_FILE = PROJECT_ROOT / "integration" / "ssn" / "matrices" / "ssn_nodes_final.tsv"
CLADE_FILE = PROJECT_ROOT / "metadata" / "genome_clades_master.tsv"
OUT_DIR = PROJECT_ROOT / "integration" / "clade_correlation" / "plots"
OUT_FILE = OUT_DIR / "mge_distribution_by_clade.png"

# ---------------------------------------------------------------------------
# Visual constants -- reused from the SSN network figures (v9) for
# cross-figure visual consistency in the manuscript.
# ---------------------------------------------------------------------------
CATEGORY_COLOR = {
    "PPH": "#1a1a1a",
    "GI":  "#807dba",
    "IS":  "#fd8d3c",
    "PLS": "#808080",
}
CATEGORY_NAME = {
    "PPH": "Prophages",
    "GI":  "Genomic islands",
    "IS":  "Insertion sequences",
    "PLS": "Plasmids",
}
CATEGORY_ORDER = ["PPH", "GI", "IS", "PLS"]  # stacking order, bottom to top

CLADE_ORDER = ["K0", "K1", "K2"]
CLADE_GAP = 1.5  # extra x-axis gap (in bar-widths) inserted between clade groups


def load_data():
    nodes = pd.read_csv(NODES_FILE, sep="\t")
    clades = pd.read_csv(CLADE_FILE, sep="\t", usecols=["Genome", "Clade"])
    print(f"[INFO] ssn_nodes_final.tsv: {len(nodes)} nodes, "
          f"{nodes['Genome'].nunique()} distinct genomes")
    print(f"[INFO] genome_clades_master.tsv: {len(clades)} genomes with clade assignment "
          f"({clades['Clade'].value_counts().to_dict()})")
    return nodes, clades


def build_counts(nodes: pd.DataFrame, clades: pd.DataFrame) -> pd.DataFrame:
    counts = (
        nodes.groupby(["Genome", "Category"])
        .size()
        .unstack(fill_value=0)
        .reindex(columns=CATEGORY_ORDER, fill_value=0)
    )

    node_genomes = set(counts.index)
    clade_genomes = set(clades["Genome"])

    only_in_nodes = sorted(node_genomes - clade_genomes)
    only_in_clades = sorted(clade_genomes - node_genomes)
    if only_in_nodes:
        print(f"[WARN] {len(only_in_nodes)} genome(s) in ssn_nodes_final.tsv have NO clade "
              f"assignment and will be excluded: {only_in_nodes}")
    if only_in_clades:
        print(f"[WARN] {len(only_in_clades)} genome(s) in genome_clades_master.tsv have NO "
              f"MGE nodes in ssn_nodes_final.tsv (0 across all 4 categories) and will be "
              f"excluded: {only_in_clades}")

    merged = counts.join(clades.set_index("Genome")["Clade"], how="inner")
    merged = merged[merged["Clade"].isin(CLADE_ORDER)]

    dropped_clade_values = set(clades["Clade"].unique()) - set(CLADE_ORDER)
    if dropped_clade_values:
        print(f"[WARN] Clade values present but not in CLADE_ORDER (excluded): "
              f"{dropped_clade_values}")

    total = merged[CATEGORY_ORDER].sum(axis=1)
    zero_total = merged.index[total == 0].tolist()
    if zero_total:
        print(f"[WARN] {len(zero_total)} genome(s) have 0 MGEs across all 4 categories in "
              f"the SSN node set and will be excluded from the plot: {zero_total}")
    merged = merged[total > 0]
    total = total[total > 0]

    for cat in CATEGORY_ORDER:
        merged[f"{cat}_pct"] = merged[cat] / total * 100.0
    merged["Total"] = total

    print(f"[INFO] Final genome count in figure: {len(merged)}")
    print(merged.groupby("Clade").size().reindex(CLADE_ORDER).to_string())

    return merged


def order_genomes(merged: pd.DataFrame) -> list:
    """Group by clade (K0, K1, K2), and within each clade sort by total MGE
    count descending (ties broken alphabetically) -- matches the typical
    left-to-right 'most elements first' convention in this kind of figure."""
    ordered = []
    for clade in CLADE_ORDER:
        sub = merged[merged["Clade"] == clade].sort_values("Total", ascending=False)
        ordered.extend(sub.index.tolist())
    return ordered


def plot(merged: pd.DataFrame, genome_order: list):
    n = len(genome_order)
    fig_width = max(14, n * 0.14)
    fig, ax = plt.subplots(figsize=(fig_width, 8), facecolor="white")
    ax.set_facecolor("white")

    # x positions with extra gaps between clade blocks
    x_positions = []
    x = 0.0
    prev_clade = None
    clade_spans = {}  # clade -> (x_start, x_end)
    for genome in genome_order:
        clade = merged.loc[genome, "Clade"]
        if prev_clade is not None and clade != prev_clade:
            x += CLADE_GAP
        if clade not in clade_spans:
            clade_spans[clade] = [x, x]
        clade_spans[clade][1] = x
        x_positions.append(x)
        prev_clade = clade
        x += 1.0

    bottoms = np.zeros(n)
    for cat in CATEGORY_ORDER:
        vals = merged.loc[genome_order, f"{cat}_pct"].to_numpy()
        ax.bar(x_positions, vals, bottom=bottoms, width=0.85,
               color=CATEGORY_COLOR[cat], edgecolor="white", linewidth=0.15,
               label=CATEGORY_NAME[cat], zorder=2)
        bottoms += vals

    # Clade group labels/brackets above the plot
    y_top = 104
    for clade in CLADE_ORDER:
        if clade not in clade_spans:
            continue
        x0, x1 = clade_spans[clade]
        n_genomes = sum(1 for g in genome_order if merged.loc[g, "Clade"] == clade)
        mid = (x0 + x1) / 2
        ax.plot([x0 - 0.4, x1 + 0.4], [y_top, y_top], color="black", linewidth=1.2,
                clip_on=False, zorder=3)
        ax.plot([x0 - 0.4, x0 - 0.4], [y_top - 1.5, y_top], color="black",
                linewidth=1.2, clip_on=False, zorder=3)
        ax.plot([x1 + 0.4, x1 + 0.4], [y_top - 1.5, y_top], color="black",
                linewidth=1.2, clip_on=False, zorder=3)
        ax.text(mid, y_top + 2.5, f"{clade} (n={n_genomes})", ha="center", va="bottom",
                fontsize=12, fontweight="bold", clip_on=False)

    ax.set_xlim(-1, x_positions[-1] + 1)
    ax.set_ylim(0, 100)
    ax.set_ylabel("Percentage of mobile genetic elements (%)", fontsize=11)
    ax.set_xticks(x_positions)
    ax.set_xticklabels(genome_order, rotation=90, fontsize=5.0)
    ax.tick_params(axis="y", labelsize=10)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)

    ax.set_title(
        "Distribution of mobile genetic elements per genome, grouped by "
        "X. translucens pv. translucens clade (K0/K1/K2)",
        fontsize=13, pad=40,
    )

    legend_elements = [mpatches.Patch(color=CATEGORY_COLOR[c], label=CATEGORY_NAME[c])
                        for c in CATEGORY_ORDER]
    ax.legend(handles=legend_elements, loc="upper center",
              bbox_to_anchor=(0.5, -0.28), ncol=4, fontsize=10, frameon=False)

    fig.text(0.5, 0.005,
              "MGE counts from the SSN network node set (ssn_nodes_final.tsv, post "
              "redundancy-exclusion); each bar sums to 100% of that genome's own "
              "PPH+GI+IS+PLS complement.",
              ha="center", fontsize=7.5, color="#666666", style="italic")

    plt.tight_layout(rect=[0, 0.04, 1, 1])
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_FILE, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"[DONE] Saved: {OUT_FILE} ({OUT_FILE.stat().st_size} bytes)")


def main():
    nodes, clades = load_data()
    merged = build_counts(nodes, clades)
    genome_order = order_genomes(merged)
    plot(merged, genome_order)


if __name__ == "__main__":
    main()
