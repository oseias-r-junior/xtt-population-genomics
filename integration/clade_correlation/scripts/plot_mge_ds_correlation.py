#!/usr/bin/env python3
"""
plot_mge_ds_correlation.py -- heatmap(s) of Spearman correlation between MGE
category burden and DS mechanism-category burden per genome, from
mge_ds_correlation_long.tsv (build_mge_ds_correlation.py).

Four panels: pooled (all genomes with clade assignment), and K0/K1/K2
stratified as a sensitivity check (per user's explicit request -- to see
whether any pooled MGE-DS association holds within clades or is confounded
by population structure). Cell = Spearman rho (diverging colormap, -1 to
+1); asterisks = FDR-corrected significance (Benjamini-Hochberg, within
each panel's own 5x8=40 tests): * q<0.05, ** q<0.01, *** q<0.001.

Output:
    integration/clade_correlation/plots/mge_ds_correlation_heatmap.png

Usage
-----
python plot_mge_ds_correlation.py
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[3]
IN_FILE = PROJECT_ROOT / "integration" / "clade_correlation" / "matrices" / "mge_ds_correlation_long_presence.tsv"
OUT_DIR = PROJECT_ROOT / "integration" / "clade_correlation" / "plots"
OUT_FILE = OUT_DIR / "mge_ds_correlation_heatmap_v5_system_presence.png"

MGE_CATEGORY_ORDER = ["PPH", "GI", "IS", "PLS", "Total"]
DS_CATEGORY_ORDER = [
    "Restriction-Modification (RM)",
    "CRISPR-Cas (adaptive immunity)",
    "Toxin-antitoxin / Abortive infection",
    "Retrons",
    "Anti-plasmid",
    "Cyclic-nucleotide signaling (CBASS)",
    "Other/recently discovered anti-phage",
    "Total",
]
DS_LABEL_SHORT = {
    "Restriction-Modification (RM)": "RM",
    "CRISPR-Cas (adaptive immunity)": "CRISPR-Cas",
    "Toxin-antitoxin / Abortive infection": "TA / Abortive",
    "Retrons": "Retrons",
    "Anti-plasmid": "Anti-plasmid",
    "Cyclic-nucleotide signaling (CBASS)": "CBASS",
    "Other/recently discovered anti-phage": "Other anti-phage",
    "Total": "Total",
}

GROUP_ORDER = ["pooled", "K0", "K1", "K2"]
GROUP_TITLE = {"pooled": "All genomes (pooled)", "K0": "Clade K0", "K1": "Clade K1", "K2": "Clade K2"}


def sig_stars(q):
    if pd.isna(q):
        return ""
    if q < 0.001:
        return "***"
    if q < 0.01:
        return "**"
    if q < 0.05:
        return "*"
    return ""


def pivot_group(long_df: pd.DataFrame, group: str, ds_order: list):
    sub = long_df[long_df["Group"] == group]
    if sub.empty:
        return None, None, None
    rho = (sub.pivot(index="MGE_category", columns="DS_category", values="rho")
           .reindex(index=MGE_CATEGORY_ORDER, columns=ds_order))
    qval = (sub.pivot(index="MGE_category", columns="DS_category", values="qval_fdr_bh")
            .reindex(index=MGE_CATEGORY_ORDER, columns=ds_order))
    n = sub["n"].iloc[0] if len(sub) else np.nan
    return rho, qval, n


def plot_panel(ax, rho, qval, n, title, ds_order: list):
    im = ax.imshow(rho.to_numpy(), cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")

    ax.set_xticks(range(len(ds_order)))
    ax.set_xticklabels([DS_LABEL_SHORT[c] for c in ds_order], rotation=45,
                        ha="right", fontsize=7)
    ax.set_yticks(range(len(MGE_CATEGORY_ORDER)))
    ax.set_yticklabels(MGE_CATEGORY_ORDER, fontsize=8)
    ax.tick_params(length=0)
    # thick border around the whole panel block, so the 4 blocks (pooled/
    # K0/K1/K2) read as clearly separate sections rather than one jumble
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("black")
        spine.set_linewidth(1.8)
    ax.set_xticks(np.arange(-0.5, len(ds_order), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(MGE_CATEGORY_ORDER), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=1)
    ax.tick_params(which="minor", length=0)

    for i, mc in enumerate(MGE_CATEGORY_ORDER):
        for j, dc in enumerate(ds_order):
            r = rho.loc[mc, dc]
            if pd.isna(r):
                continue
            q = qval.loc[mc, dc]
            stars = sig_stars(q)
            text_color = "white" if abs(r) > 0.6 else "black"
            ax.text(j, i, f"{r:.2f}{stars}", ha="center", va="center",
                    fontsize=6.5, color=text_color)

    n_note = f" (n={n})" if pd.notna(n) else ""
    ax.set_title(f"{title}{n_note}", fontsize=10)
    return im


def main():
    long_df = pd.read_csv(IN_FILE, sep="\t")
    print(f"[INFO] Loaded {len(long_df)} correlation rows from {IN_FILE}")

    # Drop DS categories with zero variance across ALL genomes in every group
    # (rho undefined -> NaN everywhere -> empty column in every panel). This
    # happens when a DS system was only ever detected in a genome that gets
    # excluded from the clade-restricted analysis -- e.g. CBASS here was only
    # found in LW16, whose Clade is "Xtu" (X. translucens pv. undulosa, an
    # outgroup reference genome, not a K0/K1/K2 Xtt genome), so once the
    # analysis is restricted to K0/K1/K2 the CBASS column is 0 for all 101
    # genomes -- correctly excluded, not a data bug.
    has_signal = long_df.groupby("DS_category")["rho"].apply(lambda s: s.notna().any())
    ds_order = [c for c in DS_CATEGORY_ORDER if has_signal.get(c, False)]
    dropped = [c for c in DS_CATEGORY_ORDER if c not in ds_order]
    if dropped:
        print(f"[INFO] Dropping DS categories with zero variance in every group (empty "
              f"columns, no genomes carrying them within K0/K1/K2): {dropped}")

    fig, axes = plt.subplots(2, 2, figsize=(16, 12), facecolor="white")
    axes = axes.flatten()
    im = None
    for ax, group in zip(axes, GROUP_ORDER):
        rho, qval, n = pivot_group(long_df, group, ds_order)
        if rho is None:
            ax.axis("off")
            ax.set_title(f"{GROUP_TITLE[group]} (no data)", fontsize=10)
            continue
        im = plot_panel(ax, rho, qval, n, GROUP_TITLE[group], ds_order)

    # explicit layout (no tight_layout -- it doesn't reliably account for a
    # shared colorbar spanning multiple axes, which was overlapping the
    # right-hand columns of the right-column panels). Reserve a real margin
    # on the right for a dedicated colorbar axis instead.
    fig.subplots_adjust(left=0.06, right=0.88, top=0.90, bottom=0.18,
                         wspace=0.45, hspace=0.55)

    if im is not None:
        cax = fig.add_axes([0.91, 0.19, 0.018, 0.64])  # [left, bottom, width, height]
        cbar = fig.colorbar(im, cax=cax)
        cbar.set_label("Spearman rho", fontsize=10)

    fig.suptitle(
        "MGE burden vs defense-system burden per genome (Spearman correlation, "
        "conservative consensus sets)",
        fontsize=14, y=0.97,
    )
    # NOTE: fig.text does not auto-wrap long strings -- a single long line
    # here overflowed the figure width and got clipped at both edges on
    # save (v2 bug). Broken into short explicit lines instead.
    dropped_note = (
        f" Category(ies) dropped (zero variance within K0/K1/K2): {', '.join(dropped)}."
        if dropped else ""
    )
    footnote = (
        "MGE = ssn_nodes_final.tsv (all-tools-agree consensus). DS = number of distinct "
        "both-tool-confirmed systems present per genome (each system counted once), per "
        "mechanism category." + dropped_note + "\n"
        "* q<0.05, ** q<0.01, *** q<0.001 (Benjamini-Hochberg FDR, corrected within each "
        "panel's 40 tests).\n"
        "K0/K1/K2 panels are a sensitivity check for confounding by population structure "
        "-- interpret low-n panels with caution."
    )
    fig.text(0.5, 0.01, footnote, ha="center", va="bottom", fontsize=8,
              color="#666666", style="italic", multialignment="center")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_FILE, dpi=300, facecolor="white")
    plt.close(fig)
    print(f"[DONE] Saved: {OUT_FILE} ({OUT_FILE.stat().st_size} bytes)")


if __name__ == "__main__":
    main()