#!/usr/bin/env python3
"""
build_ds_consensus_matrix.py -- builds the genome x defense-system gene-count
consensus matrix from DefenseFinder + PADLOC raw outputs, using the manually
curated name crosswalk (ds_system_crosswalk.tsv, no automated fuzzy matching
-- every "certain" row was confirmed by matching identical/near-identical
names and matching n_genomes/genes_count profiles across the two tools;
"tentative" candidate matches were explicitly NOT unified, per user decision,
and are kept as separate DF-only / PADLOC-only systems; PADLOC-only novel/
uncharacterized systems (PDC-*, DMS_other, etc.) are kept as their own
columns, contributing 0 from DefenseFinder).

Per-genome, per-consensus-system gene count:
  - DefenseFinder: sum of `genes_count` across ALL system instances in that
    genome whose (type, subtype) maps to this consensus_system (multiple
    instances/subunits count independently of genomic proximity, per the
    user's explicit definition, mirroring the Uceda-Campos Fig. 4 style
    heatmap: darkness ~ number of genes present, regardless of clustering).
  - PADLOC: total row count (= gene/protein hits) across ALL system.number
    instances in that genome for the mapped padloc_system.
  - Consensus value = max(DefenseFinder value, PADLOC value), per user's
    explicit choice (rather than requiring agreement or summing).

Defense systems never detected in ANY genome by either tool are already
absent from the crosswalk (it was built from an inventory of what WAS
detected), so there is nothing to omit here -- unlike the Uceda-Campos PICI
case, which was a system actively searched for and never found.

Data sources:
  - defense_systems/defensefinder/raw_outputs/<genome>/<genome>_defense_finder_systems.tsv
  - defense_systems/padloc/results/<genome>/<genome>_padloc.csv
  - CROSSWALK below (manually curated inline -- was originally a separate TSV,
    moved inline after a heredoc paste corrupted its tab characters)

Output:
  defense_systems/matrices/ds_consensus_gene_counts.tsv
  (rows = genomes, columns = consensus_system, values = consensus gene count;
   long-format companion also saved: ds_consensus_gene_counts_long.tsv with
   per-genome/system DF value, PADLOC value, and the max used)

Usage
-----
python build_ds_consensus_matrix.py
"""
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))
from genome_ids import canon  # noqa: E402  (unifies AB22010_2 / AB22010-2 etc.)
DF_DIR = PROJECT_ROOT / "defense_systems" / "defensefinder" / "raw_outputs"
PADLOC_DIR = PROJECT_ROOT / "defense_systems" / "padloc" / "results"
OUT_DIR = PROJECT_ROOT / "defense_systems" / "matrices"
OUT_WIDE = OUT_DIR / "ds_consensus_gene_counts.tsv"
OUT_LONG = OUT_DIR / "ds_consensus_gene_counts_long.tsv"

# ---------------------------------------------------------------------------
# Manually curated DefenseFinder <-> PADLOC name crosswalk.
# Each tuple: (consensus_system, df_type_or_None, df_subtype_or_None, padloc_system_or_None)
# "certain" matches: identical/near-identical name, or identical n_genomes +
# genes_count profile between the two tools (see chat discussion for the
# per-system reasoning). "tentative" candidate matches (CapRel/PDC-S33,
# Kongming/ietAS, Prometheus/PDC-S19, RloC/PDC-S64) were explicitly NOT
# unified per user decision -- kept as separate DF-only/PADLOC-only rows.
# PADLOC-only novel/uncharacterized systems (PDC-*, DMS_other, etc.) are
# included as their own columns per user decision, contributing 0 from
# DefenseFinder.
# ---------------------------------------------------------------------------
CROSSWALK = [
    # certain matches
    ("AbiD", "AbiD", "AbiD", None),
    ("AbiD", None, None, "AbiD"),
    ("CBASS", "CBASS", "CBASS_I_TM", None),
    ("CBASS", None, None, "cbass_type_I"),
    ("Cas_Type_I-C", "Cas", "CAS_Class1-Subtype-I-C", None),
    ("Cas_Type_I-C", None, None, "cas_type_I-C"),
    ("DarTG", "DarTG", "DarTG", None),
    ("DarTG", None, None, "darTG"),
    ("Gabija", "Gabija", "Gabija", None),
    ("Gabija", None, None, "gabija"),
    ("Gao_Ppl", "Gao_Ppl", "Gao_Ppl", None),
    ("Gao_Ppl", None, None, "ppl"),
    ("Gao_TerY", "Gao_TerY", "Gao_TerY", None),
    ("Gao_TerY", None, None, "TerY-P"),
    ("HEC-01", "HEC-01", "HEC-01", None),
    ("HEC-01", None, None, "HEC-01"),
    ("HEC-03", "HEC-03", "HEC-03", None),
    ("HEC-03", None, None, "HEC-03"),
    ("HEC-09", "HEC-09", "HEC-09", None),
    ("HEC-09", None, None, "HEC-09"),
    ("Lamassu", "Lamassu-Fam", "Lamassu_Cap4", None),
    ("Lamassu", "Lamassu-Fam", "Lamassu_Mrr", None),
    ("Lamassu", None, None, "Lamassu_Family"),
    ("Mokosh_TypeII", "Mokosh", "Mokosh_TypeII", None),
    ("Mokosh_TypeII", None, None, "Mokosh_TypeII"),
    ("PD-Lambda-5", "PD-Lambda-5", "PD-Lambda-5", None),
    ("PD-Lambda-5", None, None, "PD-Lambda-5"),
    ("PrrC", "PrrC", "PrrC", None),
    ("PrrC", None, None, "PrrC"),
    ("RM_Type_I", "RM_Type_I", "N6_RM_Type_I", None),
    ("RM_Type_I", None, None, "RM_type_I"),
    ("RM_Type_II", "RM_Type_II", "DNA_met_RM_Type_II_PDDEXK", None),
    ("RM_Type_II", "RM_Type_II", "N6_RM_Type_II_PDDEXK", None),
    ("RM_Type_II", None, None, "RM_type_II"),
    ("RM_Type_IIG", "RM_Type_IIG", "N6_RM_Type_IIG", None),
    ("RM_Type_IIG", None, None, "RM_type_IIG"),
    ("RM_Type_III", "RM_Type_III", "N6_N4_RM_Type_III", None),
    ("RM_Type_III", None, None, "RM_type_III"),
    ("Retron_XI", "Retron", "Retron_XI", None),
    ("Retron_XI", None, None, "retron_XI"),
    ("Retron_XII", "Retron", "Retron_XII", None),
    ("Retron_XII", None, None, "retron_XII"),
    ("Septu", "Septu", "Septu", None),
    ("Septu", None, None, "septu_type_I"),
    ("Shango", "Shango", "Shango", None),
    ("Shango", None, None, "Shango"),
    ("SoFic", "SoFIC", "SoFic", None),
    ("SoFic", None, None, "SoFic"),
    ("Wadjet_I", "Wadjet", "Wadjet_I", None),
    ("Wadjet_I", None, None, "wadjet_type_I"),
    ("Zorya_TypeIII", "Zorya", "Zorya_TypeIII", None),
    ("Zorya_TypeIII", None, None, "zorya_type_III"),
    # DefenseFinder-only (no PADLOC candidate)
    ("AbiH", "AbiH", "AbiH", None),
    ("Avs8", "Avs", "Avs8", None),
    ("DS-1", "DS-1", "DS-1", None),
    ("DS-27", "DS-27", "DS-27", None),
    ("Gao_Ape", "Gao_Ape", "Gao_Ape", None),
    ("Gao_RL", "Gao_RL", "Gao_RL", None),
    ("Hna", "Hna", "Hna", None),
    ("McrBC", "McrBC", "McrBC", None),
    ("dCTPdeaminase", "dCTPdeaminase", "dCTPdeaminase", None),
    # tentative matches -- kept SEPARATE per user decision
    ("CapRel", "CapRel", "CapRel", None),
    ("PDC-S33", None, None, "PDC-S33"),
    ("Kongming", "Kongming", "Kongming", None),
    ("ietAS", None, None, "ietAS"),
    ("Prometheus", "Prometheus", "Prometheus", None),
    ("PDC-S19", None, None, "PDC-S19"),
    ("RloC", "RloC", "RloC", None),
    ("PDC-S64", None, None, "PDC-S64"),
    # PADLOC-only (PADLOC-DB catalog systems, no published DF equivalent)
    ("DMS_other", None, None, "DMS_other"),
    ("PD-T4-6", None, None, "PD-T4-6"),
    ("PDC-M01", None, None, "PDC-M01"),
    ("PDC-M05", None, None, "PDC-M05"),
    ("PDC-S04", None, None, "PDC-S04"),
    ("PDC-S05", None, None, "PDC-S05"),
    ("PDC-S06", None, None, "PDC-S06"),
    ("PDC-S08", None, None, "PDC-S08"),
    ("PDC-S13", None, None, "PDC-S13"),
    ("PDC-S18", None, None, "PDC-S18"),
    ("PDC-S28", None, None, "PDC-S28"),
    ("PDC-S35", None, None, "PDC-S35"),
    ("PDC-S55", None, None, "PDC-S55"),
    ("PDC-S62", None, None, "PDC-S62"),
    ("PDC-S67", None, None, "PDC-S67"),
    ("SEFIR", None, None, "SEFIR"),
    ("cas_type_other", None, None, "cas_type_other"),
    ("wadjet_other", None, None, "wadjet_other"),
    ("zorya_other", None, None, "zorya_other"),
    ("retron_I-C", None, None, "retron_I-C"),
]

# ---------------------------------------------------------------------------
# Per user decision: only keep systems CONFIRMED BY BOTH TOOLS (i.e. every
# system that has at least one DefenseFinder row AND at least one PADLOC row
# in CROSSWALK above -- the "certain" matches). DF-only, PADLOC-only, and the
# 4 tentative unconfirmed pairs (CapRel/PDC-S33, Kongming/ietAS,
# Prometheus/PDC-S19, RloC/PDC-S64) are excluded from the output matrix.
# ---------------------------------------------------------------------------
def _both_tool_systems():
    has_df, has_padloc = set(), set()
    for cs, df_type, df_subtype, padloc_system in CROSSWALK:
        if df_type is not None:
            has_df.add(cs)
        if padloc_system is not None:
            has_padloc.add(cs)
    return has_df & has_padloc


CERTAIN_SYSTEMS = _both_tool_systems()


def load_crosswalk():
    df_map = {}   # (type, subtype) -> consensus_system
    padloc_map = {}  # padloc_system -> consensus_system
    consensus_systems = []
    for cs, df_type, df_subtype, padloc_system in CROSSWALK:
        if cs not in CERTAIN_SYSTEMS:
            continue  # DF-only / PADLOC-only / tentative -- excluded per user decision
        if cs not in consensus_systems:
            consensus_systems.append(cs)
        if df_type is not None:
            df_map[(df_type, df_subtype)] = cs
        if padloc_system is not None:
            padloc_map[padloc_system] = cs
    print(f"[INFO] Crosswalk loaded (both-tool-confirmed systems only): "
          f"{len(consensus_systems)} consensus systems, "
          f"{len(df_map)} DefenseFinder (type,subtype) mappings, "
          f"{len(padloc_map)} PADLOC system mappings")
    return df_map, padloc_map, consensus_systems


def parse_defensefinder(df_map):
    files = sorted(DF_DIR.glob("*/*_defense_finder_systems.tsv"))
    print(f"[INFO] Parsing {len(files)} DefenseFinder systems.tsv files")
    values = defaultdict(int)  # (genome, consensus_system) -> summed genes_count
    unmapped = defaultdict(int)
    for f in files:
        genome = canon(f.name.replace("_defense_finder_systems.tsv", ""))
        df = pd.read_csv(f, sep="\t")
        for _, row in df.iterrows():
            key = (str(row.get("type", "")), str(row.get("subtype", "")))
            cs = df_map.get(key)
            if cs is None:
                unmapped[key] += 1
                continue
            try:
                gc = int(row.get("genes_count", 0))
            except (ValueError, TypeError):
                gc = 0
            values[(genome, cs)] += gc
    if unmapped:
        print(f"[INFO] {len(unmapped)} DefenseFinder (type,subtype) combos are DF-only/"
              f"tentative and excluded per the both-tool-confirmed filter: {dict(unmapped)}")
    print(f"[INFO] DefenseFinder: {len(values)} (genome, consensus_system) cells with genes")
    return values


def parse_padloc(padloc_map):
    files = sorted(PADLOC_DIR.glob("*/*_padloc.csv"))
    print(f"[INFO] Parsing {len(files)} PADLOC *_padloc.csv files")
    values = defaultdict(int)  # (genome, consensus_system) -> total gene-hit rows
    unmapped = defaultdict(int)
    for f in files:
        genome = canon(f.name.replace("_padloc.csv", ""))
        df = pd.read_csv(f)
        if "system" not in df.columns:
            print(f"[WARN] {f} missing 'system' column, skipping")
            continue
        for sys_name, count in df["system"].astype(str).value_counts().items():
            cs = padloc_map.get(sys_name)
            if cs is None:
                unmapped[sys_name] += 1
                continue
            values[(genome, cs)] += int(count)
    if unmapped:
        print(f"[INFO] {len(unmapped)} PADLOC system labels are PADLOC-only/tentative and "
              f"excluded per the both-tool-confirmed filter: {dict(unmapped)}")
    print(f"[INFO] PADLOC: {len(values)} (genome, consensus_system) cells with genes")
    return values


def main():
    df_map, padloc_map, consensus_systems = load_crosswalk()
    df_values = parse_defensefinder(df_map)
    padloc_values = parse_padloc(padloc_map)

    all_genomes = sorted({g for (g, _) in df_values} | {g for (g, _) in padloc_values})
    print(f"[INFO] Total distinct genomes with >=1 DS detection: {len(all_genomes)}")

    long_rows = []
    for genome in all_genomes:
        for cs in consensus_systems:
            dfv = df_values.get((genome, cs), 0)
            plv = padloc_values.get((genome, cs), 0)
            if dfv == 0 and plv == 0:
                continue
            consensus_val = max(dfv, plv)
            long_rows.append({
                "Genome": genome, "consensus_system": cs,
                "defensefinder_genes": dfv, "padloc_genes": plv,
                "consensus_genes": consensus_val,
            })
    if not long_rows:
        print("[FATAL] No (genome, consensus_system) cells with genes found -- check that the "
              "crosswalk keys match the raw file values exactly (case/whitespace).")
        return
    long_df = pd.DataFrame(long_rows)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    long_df.to_csv(OUT_LONG, sep="\t", index=False)
    print(f"[DONE] Saved long-format: {OUT_LONG} ({len(long_df)} rows)")

    wide = (long_df.pivot(index="Genome", columns="consensus_system", values="consensus_genes")
            .reindex(columns=consensus_systems, fill_value=0)
            .fillna(0).astype(int))
    wide = wide.reindex(all_genomes)
    wide.to_csv(OUT_WIDE, sep="\t")
    print(f"[DONE] Saved wide matrix: {OUT_WIDE} ({wide.shape[0]} genomes x {wide.shape[1]} systems)")

    empty_cols = wide.columns[(wide == 0).all(axis=0)].tolist()
    if empty_cols:
        print(f"[WARN] {len(empty_cols)} consensus systems have 0 genes across ALL genomes "
              f"(unexpected, check crosswalk): {empty_cols}")
    else:
        print("[INFO] Every consensus system has >=1 genome with >=1 gene (as expected).")


if __name__ == "__main__":
    main()
