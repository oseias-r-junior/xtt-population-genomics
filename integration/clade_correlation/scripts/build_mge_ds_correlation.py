#!/usr/bin/env python3
"""
build_mge_ds_correlation.py -- Spearman correlation between MGE burden and
defense-system (DS) burden per genome, using the conservative set on both
sides (per user's explicit instruction: "Conjunto conservador em ambos"):

  - MGE: integration/ssn/matrices/ssn_nodes_final.tsv (the "all tools agree"
    consensus node set, post redundancy-exclusion -- the same set used
    throughout this project's clade-correlation figures).
  - DS: defense_systems/matrices/ds_consensus_gene_counts.tsv, restricted to
    the 25 systems confirmed by BOTH DefenseFinder and PADLOC (build_ds_
    consensus_matrix.py's CERTAIN_SYSTEMS filter).

Granularity (per user decision): MGE by category (PPH/GI/IS/PLS/Total) x DS
by mechanism category (RM, CRISPR-Cas, TA/Abortive infection, Retrons,
Anti-plasmid, CBASS, Other/recently discovered anti-phage, Total) -- a 5x8
matrix of Spearman correlations, not per-individual-system (fewer tests,
more direct for the manuscript narrative).

Statistic (per user decision): Spearman rank correlation per (MGE category,
DS category) pair across genomes, with Benjamini-Hochberg FDR correction
applied across all pairs tested within each analysis block.

Stratification (per user decision): the primary analysis pools all genomes
with both MGE and DS data and a valid clade assignment; as a sensitivity
check, the same 5x8 correlation matrix is recomputed separately within each
clade (K0/K1/K2) to see whether any pooled association holds up within
clades or is confounded by population structure. Clade blocks with few
genomes have correspondingly low power -- reported explicitly, not hidden.

Output:
  integration/clade_correlation/matrices/mge_ds_correlation_long.tsv
  columns: Group (pooled/K0/K1/K2), MGE_category, DS_category, n, rho,
  pval, qval_fdr_bh

Usage
-----
python build_mge_ds_correlation.py
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

PROJECT_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(PROJECT_ROOT))
from genome_ids import canon  # noqa: E402  (unifies AB22010_2 / AB22010-2 etc.)
MGE_NODES_FILE = PROJECT_ROOT / "integration" / "ssn" / "matrices" / "ssn_nodes_final.tsv"
# v2: system-level presence (0/1 per system per genome, build_ds_presence_matrix.py);
# DS burden per category = number of DISTINCT systems present, not summed gene copies.
DS_MATRIX_FILE = PROJECT_ROOT / "defense_systems" / "matrices" / "ds_consensus_presence.tsv"
CLADE_FILE = PROJECT_ROOT / "metadata" / "genome_clades_master.tsv"
OUT_DIR = PROJECT_ROOT / "integration" / "clade_correlation" / "matrices"
OUT_FILE = OUT_DIR / "mge_ds_correlation_long_presence.tsv"

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

# Same mechanism-category assignment used in plot_ds_heatmap_by_clade.py,
# restricted to the 25 both-tool-confirmed systems.
SYSTEM_CATEGORY = {
    "RM_Type_I": "Restriction-Modification (RM)",
    "RM_Type_II": "Restriction-Modification (RM)",
    "RM_Type_IIG": "Restriction-Modification (RM)",
    "RM_Type_III": "Restriction-Modification (RM)",
    "Cas_Type_I-C": "CRISPR-Cas (adaptive immunity)",
    "AbiD": "Toxin-antitoxin / Abortive infection",
    "DarTG": "Toxin-antitoxin / Abortive infection",
    "PrrC": "Toxin-antitoxin / Abortive infection",
    "SoFic": "Toxin-antitoxin / Abortive infection",
    "Retron_XI": "Retrons",
    "Retron_XII": "Retrons",
    "Wadjet_I": "Anti-plasmid",
    "CBASS": "Cyclic-nucleotide signaling (CBASS)",
    "Gabija": "Other/recently discovered anti-phage",
    "Gao_Ppl": "Other/recently discovered anti-phage",
    "Gao_TerY": "Other/recently discovered anti-phage",
    "HEC-01": "Other/recently discovered anti-phage",
    "HEC-03": "Other/recently discovered anti-phage",
    "HEC-09": "Other/recently discovered anti-phage",
    "Lamassu": "Other/recently discovered anti-phage",
    "Mokosh_TypeII": "Other/recently discovered anti-phage",
    "PD-Lambda-5": "Other/recently discovered anti-phage",
    "Septu": "Other/recently discovered anti-phage",
    "Shango": "Other/recently discovered anti-phage",
    "Zorya_TypeIII": "Other/recently discovered anti-phage",
}

CLADE_ORDER = ["K0", "K1", "K2"]


def bh_fdr(pvals: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg FDR correction, implemented directly (no
    statsmodels dependency, to avoid an extra env requirement on CCAST)."""
    pvals = np.asarray(pvals, dtype=float)
    n = len(pvals)
    order = np.argsort(pvals)
    ranked = pvals[order]
    q = ranked * n / (np.arange(1, n + 1))
    # enforce monotonicity (BH step-up procedure)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)
    out = np.empty(n)
    out[order] = q
    return out


def load_mge_counts():
    nodes = pd.read_csv(MGE_NODES_FILE, sep="\t")
    nodes["Genome"] = nodes["Genome"].map(canon)
    counts = (nodes.groupby(["Genome", "Category"]).size()
              .unstack(fill_value=0)
              .reindex(columns=["PPH", "GI", "IS", "PLS"], fill_value=0))
    # A genome with no MGE in the consensus node set has burden 0 -- it is real
    # data, not missing data -- so every genome of the master table is kept.
    master = pd.read_csv(CLADE_FILE, sep="\t", usecols=["Genome"])["Genome"].astype(str).str.strip()
    zero_mge = sorted(set(master) - set(counts.index))
    counts = counts.reindex(sorted(set(counts.index) | set(master)), fill_value=0)
    if zero_mge:
        print(f"[INFO] {len(zero_mge)} genome(s) with no MGE node set to burden 0: {zero_mge}")
    counts["Total"] = counts.sum(axis=1)
    counts = counts.add_prefix("MGE_")
    print(f"[INFO] MGE (ssn_nodes_final.tsv): {len(nodes)} nodes, "
          f"{counts.shape[0]} distinct genomes")
    return counts


def load_ds_counts():
    matrix = pd.read_csv(DS_MATRIX_FILE, sep="\t", index_col=0)
    uncategorized = [c for c in matrix.columns if c not in SYSTEM_CATEGORY]
    if uncategorized:
        print(f"[WARN] {len(uncategorized)} DS system(s) in the matrix have no category "
              f"assignment and will be dropped from the correlation analysis: {uncategorized}")
        matrix = matrix.drop(columns=uncategorized)
    cat_of = pd.Series({c: SYSTEM_CATEGORY[c] for c in matrix.columns})
    by_cat = matrix.T.groupby(cat_of).sum().T
    by_cat = by_cat.reindex(columns=[c for c in DS_CATEGORY_ORDER if c != "Total"], fill_value=0)
    by_cat["Total"] = by_cat.sum(axis=1)
    by_cat = by_cat.add_prefix("DS_")
    print(f"[INFO] DS (ds_consensus_gene_counts.tsv, both-tool-confirmed systems): "
          f"{matrix.shape[0]} genomes x {matrix.shape[1]} systems -> "
          f"{by_cat.shape[1] - 1} mechanism categories")
    return by_cat


def build_merged(mge_counts, ds_counts):
    clades = pd.read_csv(CLADE_FILE, sep="\t", usecols=["Genome", "Clade"]).set_index("Genome")

    merged = mge_counts.join(ds_counts, how="inner")
    only_mge = sorted(set(mge_counts.index) - set(ds_counts.index))
    only_ds = sorted(set(ds_counts.index) - set(mge_counts.index))
    if only_mge:
        print(f"[WARN] {len(only_mge)} genome(s) have MGE data but no DS data, excluded: {only_mge}")
    if only_ds:
        print(f"[WARN] {len(only_ds)} genome(s) have DS data but no MGE data, excluded: {only_ds}")

    merged = merged.join(clades, how="left")
    no_clade = merged.index[merged["Clade"].isna()].tolist()
    if no_clade:
        print(f"[WARN] {len(no_clade)} genome(s) have MGE+DS data but no clade assignment, "
              f"excluded from analysis: {no_clade}")
    merged = merged[merged["Clade"].isin(CLADE_ORDER)]

    print(f"[INFO] Final merged genome count (MGE + DS + clade, all present): {len(merged)}")
    print(merged["Clade"].value_counts().reindex(CLADE_ORDER).to_string())
    return merged


def correlate_block(df: pd.DataFrame, group_label: str):
    mge_cols = [f"MGE_{c}" for c in MGE_CATEGORY_ORDER]
    ds_cols = [f"DS_{c}" for c in DS_CATEGORY_ORDER]
    n = len(df)

    rows = []
    for mc in mge_cols:
        for dc in ds_cols:
            x, y = df[mc].to_numpy(), df[dc].to_numpy()
            if n < 3 or np.all(x == x[0]) or np.all(y == y[0]):
                rho, pval = np.nan, np.nan
            else:
                rho, pval = spearmanr(x, y)
            rows.append({
                "Group": group_label,
                "MGE_category": mc.replace("MGE_", ""),
                "DS_category": dc.replace("DS_", ""),
                "n": n, "rho": rho, "pval": pval,
            })
    block = pd.DataFrame(rows)
    valid = block["pval"].notna()
    block.loc[valid, "qval_fdr_bh"] = bh_fdr(block.loc[valid, "pval"].to_numpy())
    return block


def main():
    mge_counts = load_mge_counts()
    ds_counts = load_ds_counts()
    merged = build_merged(mge_counts, ds_counts)

    if len(merged) < 3:
        print("[FATAL] Fewer than 3 genomes with complete MGE+DS+clade data -- cannot correlate.")
        return

    blocks = [correlate_block(merged, "pooled")]
    for clade in CLADE_ORDER:
        sub = merged[merged["Clade"] == clade]
        print(f"[INFO] Clade {clade}: n={len(sub)} genomes for stratified sensitivity check")
        if len(sub) < 5:
            print(f"[WARN] Clade {clade} has only {len(sub)} genomes -- correlation estimates "
                  f"here will be very low-power / unstable; interpret with caution.")
        blocks.append(correlate_block(sub, clade))

    result = pd.concat(blocks, ignore_index=True)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    result.to_csv(OUT_FILE, sep="\t", index=False)
    print(f"[DONE] Saved: {OUT_FILE} ({len(result)} rows)")

    sig = result[(result["qval_fdr_bh"] < 0.05)].sort_values("qval_fdr_bh")
    if len(sig):
        print(f"[INFO] {len(sig)} (Group, MGE_category, DS_category) pairs significant at FDR<0.05:")
        print(sig[["Group", "MGE_category", "DS_category", "n", "rho", "pval", "qval_fdr_bh"]]
              .to_string(index=False))
    else:
        print("[INFO] No pairs significant at FDR<0.05.")


if __name__ == "__main__":
    main()
