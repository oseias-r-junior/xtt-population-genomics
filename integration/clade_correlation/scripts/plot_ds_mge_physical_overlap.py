#!/usr/bin/env python3
"""
plot_ds_mge_physical_overlap.py (v2, SYSTEM-LEVEL) -- heatmap(s) of the
percentage of defense-system INSTANCES (loci) located inside an MGE region,
from ds_mge_physical_overlap_long_system.tsv (build_ds_mge_physical_overlap.py).

Rows = MGE category (PPH/GI/IS/PLS/Any_MGE), columns = DS mechanism
category. Cell = % of that DS category's system instances (PADLOC
system.number) having at least one gene overlapping that MGE category's
regions. Four panels: pooled + K0/K1/K2.

Output:
    integration/clade_correlation/plots/ds_mge_physical_overlap_heatmap_v2_system.png

Usage
-----
python plot_ds_mge_physical_overlap.py
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[3]
IN_FILE = PROJECT_ROOT / "integration" / "clade_correlation" / "matrices" / "ds_mge_physical_overlap_long_system.tsv"
OUT_DIR = PROJECT_ROOT / "integration" / "clade_correlation" / "plots"
OUT_FILE = OUT_DIR / "ds_mge_physical_overlap_heatmap_v2_system.png"

MGE_CATEGORY_ORDER = ["PPH", "GI", "IS", "PLS", "Any_MGE"]
MGE_LABEL = {"PPH": "PPH", "GI": "GI", "IS": "IS", "PLS": "PLS", "Any_MGE": "Any MGE"}
DS_CATEGORY_ORDER = [
    "Restriction-Modification (RM)",
    "CRISPR-Cas (adaptive immunity)",
    "Toxin-antitoxin / Abortive infection",
    "Retrons",
    "Anti-plasmid",
    "Cyclic-nucleotide signaling (CBASS)",
    "Other/recently discovered anti-phage",
]
DS_LABEL_SHORT = {
    "Restriction-Modification (RM)": "RM",
    "CRISPR-Cas (adaptive immunity)": "CRISPR-Cas",
    "Toxin-antitoxin / Abortive infection": "TA / Abortive",
    "Retrons": "Retrons",
    "Anti-plasmid": "Anti-plasmid",
    "Cyclic-nucleotide signaling (CBASS)": "CBASS",
    "Other/recently discovered anti-phage": "Other anti-phage",
}

GROUP_ORDER = ["pooled", "K0", "K1", "K2"]
GROUP_TITLE = {"pooled": "All genomes (pooled)", "K0": "Clade K0", "K1": "Clade K1", "K2": "Clade K2"}


def pivot_group(long_df: pd.DataFrame, group: str, ds_order: list):
    sub = long_df[long_df["Group"] == group]
    if sub.empty:
        return None, None
    pct = (sub.pivot(index="MGE_category", columns="DS_category", values="pct_contained")
           .reindex(index=MGE_CATEGORY_ORDER, columns=ds_order))
    n_ds = (sub[sub["MGE_category"] == "Any_MGE"]
            .set_index("DS_category")["n_ds_instances"].reindex(ds_order))
    return pct, n_ds


def plot_panel(ax, pct, n_ds, title, ds_order: list):
    im = ax.imshow(pct.to_numpy(), cmap="Blues", vmin=0, vmax=100, aspect="auto")

    ds_labels = [f"{DS_LABEL_SHORT[c]}\n(n={int(n_ds[c])})" if pd.notna(n_ds[c]) else
                 f"{DS_LABEL_SHORT[c]}\n(n=0)" for c in ds_order]
    ax.set_xticks(range(len(ds_order)))
    ax.set_xticklabels(ds_labels, rotation=45, ha="right", fontsize=6.5)
    ax.set_yticks(range(len(MGE_CATEGORY_ORDER)))
    ax.set_yticklabels([MGE_LABEL[c] for c in MGE_CATEGORY_ORDER], fontsize=8)
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("black")
        spine.set_linewidth(1.8)
    ax.set_xticks(np.arange(-0.5, len(ds_order), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(MGE_CATEGORY_ORDER), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=1)
    ax.tick_params(which="minor", length=0)

    any_row_idx = MGE_CATEGORY_ORDER.index("Any_MGE")
    ax.axhline(any_row_idx - 0.5, color="black", linewidth=1.4)

    for i, mc in enumerate(MGE_CATEGORY_ORDER):
        for j, dc in enumerate(ds_order):
            v = pct.loc[mc, dc]
            if pd.isna(v):
                continue
            text_color = "white" if v > 55 else "black"
            ax.text(j, i, f"{v:.0f}%", ha="center", va="center",
                    fontsize=7, color=text_color)

    ax.set_title(title, fontsize=10)
    return im


def main():
    long_df = pd.read_csv(IN_FILE, sep="\t")
    print(f"[INFO] Loaded {len(long_df)} overlap rows from {IN_FILE}")

    has_signal = long_df.groupby("DS_category")["n_ds_instances"].apply(lambda s: (s > 0).any())
    ds_order = [c for c in DS_CATEGORY_ORDER if has_signal.get(c, False)]
    dropped = [c for c in DS_CATEGORY_ORDER if c not in ds_order]
    if dropped:
        print(f"[INFO] Dropping DS categories with 0 instances in every group: {dropped}")

    fig, axes = plt.subplots(2, 2, figsize=(15, 11), facecolor="white")
    axes = axes.flatten()
    im = None
    for ax, group in zip(axes, GROUP_ORDER):
        pct, n_ds = pivot_group(long_df, group, ds_order)
        if pct is None:
            ax.axis("off")
            ax.set_title(f"{GROUP_TITLE[group]} (no data)", fontsize=10)
            continue
        im = plot_panel(ax, pct, n_ds, GROUP_TITLE[group], ds_order)

    fig.subplots_adjust(left=0.07, right=0.88, top=0.90, bottom=0.20,
                         wspace=0.45, hspace=0.65)

    if im is not None:
        cax = fig.add_axes([0.91, 0.20, 0.018, 0.62])
        cbar = fig.colorbar(im, cax=cax)
        cbar.set_label("% of defense-system instances contained in MGE region", fontsize=9)

    fig.suptitle(
        "Physical overlap: defense-system loci located inside MGE regions "
        "(PADLOC coordinates vs ssn_nodes_final.tsv)",
        fontsize=13, y=0.97,
    )
    dropped_note = f" Category(ies) with 0 instances dropped: {', '.join(dropped)}." if dropped else ""
    footnote = (
        "Unit = system instance (one PADLOC system.number = one locus); an instance counts as "
        "inside an MGE if any of its genes overlaps the region (same contig)." + dropped_note + "\n"
        "Only PADLOC provides coordinates, so systems detected only by DefenseFinder are not placed; "
        "restricted to the 25 both-tool-confirmed systems.\n"
        "MGE = ssn_nodes_final.tsv (all-tools-agree consensus). 'Any MGE' = union across "
        "PPH/GI/IS/PLS (not a sum -- avoids double-counting)."
    )
    fig.text(0.5, 0.01, footnote, ha="center", va="bottom", fontsize=7.5,
              color="#666666", style="italic", multialignment="center")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_FILE, dpi=300, facecolor="white")
    plt.close(fig)
    print(f"[DONE] Saved: {OUT_FILE} ({OUT_FILE.stat().st_size} bytes)")


if __name__ == "__main__":
    main()