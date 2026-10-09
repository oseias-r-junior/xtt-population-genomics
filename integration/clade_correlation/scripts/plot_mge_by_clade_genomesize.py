#!/usr/bin/env python3
"""
plot_mge_by_clade_genomesize.py -- MGE genomic footprint per genome, grouped
by Xtt clade.

Companion to plot_mge_by_clade.py, using a different metric per user
feedback: the first version stacked bars by ELEMENT COUNT composition
(each genome's bar summed to 100%), which let insertion sequences (numerous
but short, typically hundreds of bp) visually dominate genomes that are
actually more affected, in genomic real-estate terms, by prophages or
genomic islands (few elements, but each tens of kb). This version instead
stacks bars by % OF GENOME LENGTH occupied by each MGE category:

    pct(genome, category) = sum(Length_bp of that category's elements) /
                             Genome_Length(genome) * 100

Bars are NOT normalized to 100% here -- total bar height is the real
fraction of that genome's sequence that is MGE-derived (in the SSN node
set), so it is directly comparable across genomes and categories.

Data sources (same as plot_mge_by_clade.py):
  - MGE Length_bp per genome/category: integration/ssn/matrices/ssn_nodes_final.tsv
    (same SSN network node set, post redundancy-exclusion -- explicit user choice).
  - Clade + genome length: metadata/genome_clades_master.tsv (Genome, Clade,
    Genome_Length columns only -- not reproducing the full table).

Output:
    integration/clade_correlation/plots/mge_genome_fraction_by_clade.png

Usage
-----
python plot_mge_by_clade_genomesize.py
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
OUT_FILE = OUT_DIR / "mge_genome_fraction_by_clade.png"

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
CLADE_GAP = 1.5


def load_data():
    nodes = pd.read_csv(NODES_FILE, sep="\t")
    clade_meta = pd.read_csv(CLADE_FILE, sep="\t",
                              usecols=["Genome", "Clade", "Genome_Length"])
    print(f"[INFO] ssn_nodes_final.tsv: {len(nodes)} nodes, "
          f"{nodes['Genome'].nunique()} distinct genomes")
    print(f"[INFO] genome_clades_master.tsv: {len(clade_meta)} genomes "
          f"({clade_meta['Clade'].value_counts().to_dict()})")
    return nodes, clade_meta


def build_fractions(nodes: pd.DataFrame, clade_meta: pd.DataFrame) -> pd.DataFrame:
    bp_sums = (
        nodes.groupby(["Genome", "Category"])["Length_bp"]
        .sum()
        .unstack(fill_value=0)
        .reindex(columns=CATEGORY_ORDER, fill_value=0)
    )

    node_genomes = set(bp_sums.index)
    clade_genomes = set(clade_meta["Genome"])
    only_in_nodes = sorted(node_genomes - clade_genomes)
    only_in_clades = sorted(clade_genomes - node_genomes)
    if only_in_nodes:
        print(f"[WARN] {len(only_in_nodes)} genome(s) in ssn_nodes_final.tsv have NO clade/"
              f"length metadata and will be excluded: {only_in_nodes}")
    if only_in_clades:
        print(f"[WARN] {len(only_in_clades)} genome(s) in genome_clades_master.tsv have NO "
              f"MGE nodes in ssn_nodes_final.tsv and will be excluded: {only_in_clades}")

    merged = bp_sums.join(clade_meta.set_index("Genome")[["Clade", "Genome_Length"]],
                           how="inner")
    merged = merged[merged["Clade"].isin(CLADE_ORDER)]

    missing_length = merged["Genome_Length"].isna() | (merged["Genome_Length"] <= 0)
    if missing_length.any():
        print(f"[WARN] {missing_length.sum()} genome(s) have missing/invalid Genome_Length "
              f"and will be excluded: {merged.index[missing_length].tolist()}")
    merged = merged[~missing_length]

    for cat in CATEGORY_ORDER:
        merged[f"{cat}_pct"] = merged[cat] / merged["Genome_Length"] * 100.0
    merged["Total_pct"] = merged[[f"{c}_pct" for c in CATEGORY_ORDER]].sum(axis=1)

    zero_total = merged.index[merged["Total_pct"] == 0].tolist()
    if zero_total:
        print(f"[WARN] {len(zero_total)} genome(s) have 0% MGE genome fraction "
              f"(no elements in SSN node set) and will be excluded: {zero_total}")
    merged = merged[merged["Total_pct"] > 0]

    print(f"[INFO] Final genome count in figure: {len(merged)}")
    print(merged.groupby("Clade").size().reindex(CLADE_ORDER).to_string())
    print(f"[INFO] Total MGE genome-fraction range: "
          f"{merged['Total_pct'].min():.2f}% - {merged['Total_pct'].max():.2f}% "
          f"(mean {merged['Total_pct'].mean():.2f}%)")

    return merged


def order_genomes(merged: pd.DataFrame) -> list:
    ordered = []
    for clade in CLADE_ORDER:
        sub = merged[merged["Clade"] == clade].sort_values("Total_pct", ascending=False)
        ordered.extend(sub.index.tolist())
    return ordered


def plot(merged: pd.DataFrame, genome_order: list):
    n = len(genome_order)
    fig_width = max(14, n * 0.14)
    y_max_data = merged["Total_pct"].max()
    y_max = y_max_data * 1.15  # headroom for clade brackets, added below

    fig, ax = plt.subplots(figsize=(fig_width, 8), facecolor="white")
    ax.set_facecolor("white")

    x_positions = []
    x = 0.0
    prev_clade = None
    clade_spans = {}
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

    y_bracket = y_max_data * 1.06
    y_label = y_max_data * 1.10
    for clade in CLADE_ORDER:
        if clade not in clade_spans:
            continue
        x0, x1 = clade_spans[clade]
        n_genomes = sum(1 for g in genome_order if merged.loc[g, "Clade"] == clade)
        mid = (x0 + x1) / 2
        tick = y_max_data * 0.02
        ax.plot([x0 - 0.4, x1 + 0.4], [y_bracket, y_bracket], color="black",
                linewidth=1.2, clip_on=False, zorder=3)
        ax.plot([x0 - 0.4, x0 - 0.4], [y_bracket - tick, y_bracket], color="black",
                linewidth=1.2, clip_on=False, zorder=3)
        ax.plot([x1 + 0.4, x1 + 0.4], [y_bracket - tick, y_bracket], color="black",
                linewidth=1.2, clip_on=False, zorder=3)
        ax.text(mid, y_label, f"{clade} (n={n_genomes})", ha="center", va="bottom",
                fontsize=12, fontweight="bold", clip_on=False)

    ax.set_xlim(-1, x_positions[-1] + 1)
    ax.set_ylim(0, y_max)
    ax.set_ylabel("Percentage of genome length (%)", fontsize=11)
    ax.set_xticks(x_positions)
    ax.set_xticklabels(genome_order, rotation=90, fontsize=5.0)
    ax.tick_params(axis="y", labelsize=10)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)

    ax.set_title(
        "Genomic footprint of mobile genetic elements per genome, grouped by "
        "X. translucens pv. translucens clade (K0/K1/K2)",
        fontsize=13, pad=40,
    )

    legend_elements = [mpatches.Patch(color=CATEGORY_COLOR[c], label=CATEGORY_NAME[c])
                        for c in CATEGORY_ORDER]
    ax.legend(handles=legend_elements, loc="upper center",
              bbox_to_anchor=(0.5, -0.28), ncol=4, fontsize=10, frameon=False)

    fig.text(0.5, 0.005,
              "Bar height = % of that genome's total length occupied by MGEs from the SSN "
              "network node set (ssn_nodes_final.tsv, post redundancy-exclusion); NOT "
              "normalized to 100% -- unlike element-count composition, this reflects actual "
              "genomic real estate per category.",
              ha="center", fontsize=7.5, color="#666666", style="italic")

    plt.tight_layout(rect=[0, 0.04, 1, 1])
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_FILE, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"[DONE] Saved: {OUT_FILE} ({OUT_FILE.stat().st_size} bytes)")


def main():
    nodes, clade_meta = load_data()
    merged = build_fractions(nodes, clade_meta)
    genome_order = order_genomes(merged)
    plot(merged, genome_order)


if __name__ == "__main__":
    main()
