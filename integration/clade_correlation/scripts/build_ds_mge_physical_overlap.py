#!/usr/bin/env python3
"""
build_ds_mge_physical_overlap.py (v2, SYSTEM-LEVEL) -- checks whether
defense-system (DS) loci are physically located INSIDE an MGE region in the
same genome, as a follow-up to the genome-level Spearman correlation
(mge_ds_correlation_long_presence.tsv): a positive MGE x DS correlation could
reflect (a) real co-selection of two independent features, or (b) the DS
systems being cargo carried BY the MGE itself ("defense islands"; Bernheim &
Sorek 2020, Rousset et al. 2022).

v2 change: the unit is now the defense-system INSTANCE (one locus, i.e. one
PADLOC system.number), not the individual gene hit. An instance counts as
"in an MGE" if ANY of its genes overlaps an MGE region on the same contig.
This keeps a 3-gene system from weighing 3 times and avoids counting a
system's genes separately. A second, coarser view is also saved: the
(genome, system) PAIR is "in an MGE" if at least one of its instances is --
the same "count once per genome" convention as the correlation.

Coordinates: only PADLOC has real nucleotide positions (seqid, start, end);
DefenseFinder's outputs here do not. So systems detected ONLY by
DefenseFinder cannot be placed and are not in this analysis (they ARE in the
presence matrix used for the correlation). Only the 25 both-tool-confirmed
systems are used.

Output:
  integration/clade_correlation/matrices/ds_mge_physical_overlap_long_system.tsv
  columns: Group, MGE_category, DS_category, n_ds_instances, n_contained,
           pct_contained, n_genome_system_pairs, n_pairs_contained,
           pct_pairs_contained

Usage
-----
python build_ds_mge_physical_overlap.py
"""
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[3]
MGE_NODES_FILE = PROJECT_ROOT / "integration" / "ssn" / "matrices" / "ssn_nodes_final.tsv"
PADLOC_DIR = PROJECT_ROOT / "defense_systems" / "padloc" / "results"
CLADE_FILE = PROJECT_ROOT / "metadata" / "genome_clades_master.tsv"
OUT_DIR = PROJECT_ROOT / "integration" / "clade_correlation" / "matrices"
OUT_FILE = OUT_DIR / "ds_mge_physical_overlap_long_system.tsv"

MGE_CATEGORY_ORDER = ["PPH", "GI", "IS", "PLS"]
CLADE_ORDER = ["K0", "K1", "K2"]

PADLOC_TO_CONSENSUS = {
    "AbiD": "AbiD",
    "cbass_type_I": "CBASS",
    "cas_type_I-C": "Cas_Type_I-C",
    "darTG": "DarTG",
    "gabija": "Gabija",
    "ppl": "Gao_Ppl",
    "TerY-P": "Gao_TerY",
    "HEC-01": "HEC-01",
    "HEC-03": "HEC-03",
    "HEC-09": "HEC-09",
    "Lamassu_Family": "Lamassu",
    "Mokosh_TypeII": "Mokosh_TypeII",
    "PD-Lambda-5": "PD-Lambda-5",
    "PrrC": "PrrC",
    "RM_type_I": "RM_Type_I",
    "RM_type_II": "RM_Type_II",
    "RM_type_IIG": "RM_Type_IIG",
    "RM_type_III": "RM_Type_III",
    "retron_XI": "Retron_XI",
    "retron_XII": "Retron_XII",
    "septu_type_I": "Septu",
    "Shango": "Shango",
    "SoFic": "SoFic",
    "wadjet_type_I": "Wadjet_I",
    "zorya_type_III": "Zorya_TypeIII",
}

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
DS_CATEGORY_ORDER = [
    "Restriction-Modification (RM)",
    "CRISPR-Cas (adaptive immunity)",
    "Toxin-antitoxin / Abortive infection",
    "Retrons",
    "Anti-plasmid",
    "Cyclic-nucleotide signaling (CBASS)",
    "Other/recently discovered anti-phage",
]


def load_mge_regions():
    nodes = pd.read_csv(MGE_NODES_FILE, sep="\t")
    nodes["Contig"] = nodes["Contig"].astype(str)
    by_genome_contig = {}
    for _, row in nodes.iterrows():
        key = (row["Genome"], row["Contig"])
        by_genome_contig.setdefault(key, []).append(
            (row["Start"], row["End"], row["Category"])
        )
    print(f"[INFO] MGE regions loaded: {len(nodes)} nodes across "
          f"{nodes['Genome'].nunique()} genomes")
    return by_genome_contig


def load_ds_genes():
    files = sorted(PADLOC_DIR.glob("*/*_padloc.csv"))
    print(f"[INFO] Parsing {len(files)} PADLOC *_padloc.csv files for gene coordinates")
    rows = []
    unmapped = set()
    for f in files:
        genome = f.name.replace("_padloc.csv", "")
        df = pd.read_csv(f)
        needed = {"system", "system.number", "seqid", "start", "end"}
        if not needed.issubset(df.columns):
            print(f"[WARN] {f} missing one of {needed}, skipping")
            continue
        for _, row in df.iterrows():
            cs = PADLOC_TO_CONSENSUS.get(str(row["system"]))
            if cs is None:
                unmapped.add(str(row["system"]))
                continue
            rows.append({
                "Genome": genome, "Contig": str(row["seqid"]),
                "SystemNumber": str(row["system.number"]),
                "Start": row["start"], "End": row["end"],
                "consensus_system": cs, "DS_category": SYSTEM_CATEGORY[cs],
            })
    if unmapped:
        print(f"[INFO] {len(unmapped)} PADLOC system label(s) not in the both-tool-confirmed "
              f"set, excluded: {sorted(unmapped)}")
    genes = pd.DataFrame(rows)
    print(f"[INFO] DS gene hits (both-tool-confirmed, with coordinates): {len(genes)}")
    return genes


def annotate_containment(genes: pd.DataFrame, mge_by_genome_contig: dict):
    for cat in MGE_CATEGORY_ORDER:
        genes[f"in_{cat}"] = False
    genes["in_Any_MGE"] = False

    for idx, row in genes.iterrows():
        regions = mge_by_genome_contig.get((row["Genome"], row["Contig"]), [])
        g_start, g_end = row["Start"], row["End"]
        hit_any = False
        for m_start, m_end, m_cat in regions:
            if min(g_end, m_end) - max(g_start, m_start) > 0:
                genes.at[idx, f"in_{m_cat}"] = True
                hit_any = True
        genes.at[idx, "in_Any_MGE"] = hit_any
    return genes


def collapse_to_instances(genes: pd.DataFrame):
    """One row per system instance (genome + PADLOC system.number + system):
    in_<cat> = True if ANY of its genes overlaps that MGE category."""
    flag_cols = [f"in_{c}" for c in MGE_CATEGORY_ORDER] + ["in_Any_MGE"]
    keys = ["Genome", "Clade", "SystemNumber", "consensus_system", "DS_category"]
    inst = genes.groupby(keys, as_index=False)[flag_cols].any()
    return inst


def summarize(inst: pd.DataFrame, group_label: str):
    flag = {c: f"in_{c}" for c in MGE_CATEGORY_ORDER + ["Any_MGE"]}
    rows = []
    for ds_cat in DS_CATEGORY_ORDER:
        sub = inst[inst["DS_category"] == ds_cat]
        pairs = sub.groupby(["Genome", "consensus_system"], as_index=False)[list(flag.values())].any()
        n, n_pairs = len(sub), len(pairs)
        for mge_cat, col in flag.items():
            if n == 0:
                pct, nc, pctp, ncp = float("nan"), 0, float("nan"), 0
            else:
                nc = int(sub[col].sum())
                pct = 100.0 * nc / n
                ncp = int(pairs[col].sum())
                pctp = 100.0 * ncp / n_pairs
            rows.append({
                "Group": group_label, "MGE_category": mge_cat, "DS_category": ds_cat,
                "n_ds_instances": n, "n_contained": nc, "pct_contained": pct,
                "n_genome_system_pairs": n_pairs, "n_pairs_contained": ncp,
                "pct_pairs_contained": pctp,
            })
    return pd.DataFrame(rows)


def main():
    mge_by_genome_contig = load_mge_regions()
    genes = load_ds_genes()
    if genes.empty:
        print("[FATAL] No DS gene hits with coordinates found -- nothing to check.")
        return
    genes = annotate_containment(genes, mge_by_genome_contig)

    clades = pd.read_csv(CLADE_FILE, sep="\t", usecols=["Genome", "Clade"])
    genes["Clade"] = genes["Genome"].map(dict(zip(clades["Genome"], clades["Clade"])))

    no_clade = genes.loc[genes["Clade"].isna(), "Genome"].unique().tolist()
    if no_clade:
        print(f"[WARN] {len(no_clade)} genome(s) with DS genes have no clade assignment "
              f"and will be excluded: {no_clade}")
    non_k_mask = genes["Clade"].notna() & ~genes["Clade"].isin(CLADE_ORDER)
    if non_k_mask.any():
        print(f"[WARN] non-K0/K1/K2 genomes excluded: "
              f"{sorted(genes.loc[non_k_mask, 'Genome'].unique().tolist())} "
              f"(clades {sorted(genes.loc[non_k_mask, 'Clade'].unique().tolist())})")
    genes = genes[genes["Clade"].isin(CLADE_ORDER)]

    inst = collapse_to_instances(genes)
    print(f"[INFO] {len(genes)} gene hits collapsed into {len(inst)} system instances "
          f"(K0/K1/K2 only)")

    blocks = [summarize(inst, "pooled")]
    print(f"[INFO] pooled: n={len(inst)} system instances")
    for clade in CLADE_ORDER:
        sub = inst[inst["Clade"] == clade]
        blocks.append(summarize(sub, clade))
        print(f"[INFO] Clade {clade}: n={len(sub)} system instances")

    result = pd.concat(blocks, ignore_index=True)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    result.to_csv(OUT_FILE, sep="\t", index=False)
    print(f"[DONE] Saved: {OUT_FILE} ({len(result)} rows)")

    print("[INFO] Pooled Any_MGE containment by DS category (instance level | genome-system pair level):")
    print(result[(result["Group"] == "pooled") & (result["MGE_category"] == "Any_MGE")]
          [["DS_category", "n_ds_instances", "n_contained", "pct_contained",
            "n_genome_system_pairs", "n_pairs_contained", "pct_pairs_contained"]]
          .to_string(index=False))


if __name__ == "__main__":
    main()