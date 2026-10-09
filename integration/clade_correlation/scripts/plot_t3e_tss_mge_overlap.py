#!/usr/bin/env python3
"""
plot_t3e_tss_mge_overlap.py -- heatmap(s) of the percentage of T3E/TSS
virulence-gene hits physically located inside an MGE region, from
t3e_tss_mge_physical_overlap_long.tsv (build_t3e_tss_mge_overlap.py).

Rows = MGE category (PPH/GI/IS/PLS/Any_MGE), columns = virulence category
(T3E, T2SS, T3SS, T4SS). Cell = % of that virulence category's gene hits
whose genomic position (resolved to the chromosome contig, see build
script's docstring) overlaps that MGE category's regions. Four panels:
pooled + K0/K1/K2, same layout conventions as
ds_mge_physical_overlap_heatmap.png (dedicated colorbar axis, thick block
borders, wrapped footnote, zero-signal columns dropped dynamically).

Output:
    integration/clade_correlation/plots/t3e_tss_mge_overlap_heatmap.png

Usage
-----
python plot_t3e_tss_mge_overlap.py
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# project root = first parent folder that contains metadata/ and integration/
# (works wherever the script is placed inside the repo)
PROJECT_ROOT = next(p for p in Path(__file__).resolve().parents
                    if (p / "metadata").is_dir() and (p / "integration").is_dir())

IN_FILE = PROJECT_ROOT / "integration" / "clade_correlation" / "matrices" / "t3e_tss_mge_physical_overlap_long.tsv"
OUT_DIR = PROJECT_ROOT / "integration" / "clade_correlation" / "plots"
OUT_FILE = OUT_DIR / "t3e_tss_mge_overlap_heatmap.png"

MGE_CATEGORY_ORDER = ["PPH", "GI", "IS", "PLS", "Any_MGE"]
MGE_LABEL = {"PPH": "PPH", "GI": "GI", "IS": "IS", "PLS": "PLS", "Any_MGE": "Any MGE"}

VIR_CATEGORY_ORDER = ["T3E", "T2SS", "T3SS", "T4SS", "TSS_other"]
VIR_LABEL = {"T3E": "T3E (effectors)", "T2SS": "T2SS", "T3SS": "T3SS",
             "T4SS": "T4SS", "TSS_other": "TSS (other)"}

GROUP_ORDER = ["pooled", "K0", "K1", "K2"]
GROUP_TITLE = {"pooled": "All genomes (pooled)", "K0": "Clade K0", "K1": "Clade K1", "K2": "Clade K2"}


def pivot_group(long_df: pd.DataFrame, group: str, vir_order: list):
    sub = long_df[long_df["Group"] == group]
    if sub.empty:
        return None, None
    pct = (sub.pivot(index="MGE_category", columns="Virulence_category", values="pct_contained")
           .reindex(index=MGE_CATEGORY_ORDER, columns=vir_order))
    n_genes = (sub[sub["MGE_category"] == "Any_MGE"]
               .set_index("Virulence_category")["n_genes"].reindex(vir_order))
    return pct, n_genes


def plot_panel(ax, pct, n_genes, title, vir_order: list):
    # GnBu (teal) on purpose: the DS x MGE overlap heatmap uses Blues and the
    # two figures sit side by side in the supplement -- keep them distinguishable.
    im = ax.imshow(pct.to_numpy(), cmap="GnBu", vmin=0, vmax=100, aspect="auto")

    vir_labels = [f"{VIR_LABEL[c]}\n(n={int(n_genes[c])})" if pd.notna(n_genes[c]) else
                  f"{VIR_LABEL[c]}\n(n=0)" for c in vir_order]
    ax.set_xticks(range(len(vir_order)))
    ax.set_xticklabels(vir_labels, rotation=45, ha="right", fontsize=7)
    ax.set_yticks(range(len(MGE_CATEGORY_ORDER)))
    ax.set_yticklabels([MGE_LABEL[c] for c in MGE_CATEGORY_ORDER], fontsize=8)
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("black")
        spine.set_linewidth(1.8)
    ax.set_xticks(np.arange(-0.5, len(vir_order), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(MGE_CATEGORY_ORDER), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=1)
    ax.tick_params(which="minor", length=0)

    # separator line above "Any_MGE" row -- summary row (union across
    # PPH/GI/IS/PLS), not a fifth independent category
    any_row_idx = MGE_CATEGORY_ORDER.index("Any_MGE")
    ax.axhline(any_row_idx - 0.5, color="black", linewidth=1.4)

    for i, mc in enumerate(MGE_CATEGORY_ORDER):
        for j, vc in enumerate(vir_order):
            v = pct.loc[mc, vc]
            if pd.isna(v):
                continue
            text_color = "white" if v > 55 else "black"
            ax.text(j, i, f"{v:.1f}%", ha="center", va="center",
                    fontsize=7, color=text_color)

    ax.set_title(title, fontsize=10)
    return im


def main():
    long_df = pd.read_csv(IN_FILE, sep="\t")
    print(f"[INFO] Loaded {len(long_df)} overlap rows from {IN_FILE}")

    has_signal = long_df.groupby("Virulence_category")["n_genes"].apply(lambda s: (s > 0).any())
    vir_order = [c for c in VIR_CATEGORY_ORDER if has_signal.get(c, False)]
    dropped = [c for c in VIR_CATEGORY_ORDER if c not in vir_order]
    if dropped:
        print(f"[INFO] Dropping virulence categories with 0 genes in every group: {dropped}")

    fig, axes = plt.subplots(2, 2, figsize=(14, 11), facecolor="white")
    axes = axes.flatten()
    im = None
    for ax, group in zip(axes, GROUP_ORDER):
        pct, n_genes = pivot_group(long_df, group, vir_order)
        if pct is None:
            ax.axis("off")
            ax.set_title(f"{GROUP_TITLE[group]} (no data)", fontsize=10)
            continue
        im = plot_panel(ax, pct, n_genes, GROUP_TITLE[group], vir_order)

    fig.subplots_adjust(left=0.08, right=0.87, top=0.90, bottom=0.20,
                         wspace=0.40, hspace=0.65)

    if im is not None:
        cax = fig.add_axes([0.90, 0.20, 0.02, 0.62])
        cbar = fig.colorbar(im, cax=cax)
        cbar.set_label("% of gene hits contained in MGE region", fontsize=9)

    fig.suptitle(
        "Physical overlap: T3E / TSS virulence genes located inside MGE regions "
        "(BLAST hits on the chromosome contig vs ssn_nodes_final.tsv)",
        fontsize=13, y=0.97,
    )
    dropped_note = f" Category(ies) with 0 genes dropped: {', '.join(dropped)}." if dropped else ""
    footnote = (
        "Virulence-gene coordinates: Xtt_2021_2023_T3E_TSS (positions).xlsx, resolved to "
        "each strain's chromosome contig (the one whose length equals the recorded "
        "query_length)." + dropped_note + "\n"
        "Cell = % of that column's gene hits whose contig+coordinates overlap that row's "
        "MGE regions (ssn_nodes_final.tsv, all-tools-agree consensus).\n"
        "'Any MGE' = union across PPH/GI/IS/PLS (not a sum -- avoids double-counting genes "
        "that overlap more than one MGE category). 2 strains (LG 2, LG 5) excluded: no "
        "matching genome FASTA found to resolve the chromosome contig."
    )
    fig.text(0.5, 0.01, footnote, ha="center", va="bottom", fontsize=7.5,
              color="#666666", style="italic", multialignment="center")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_FILE, dpi=300, facecolor="white")
    plt.close(fig)
    print(f"[DONE] Saved: {OUT_FILE} ({OUT_FILE.stat().st_size} bytes)")


if __name__ == "__main__":
    main()
