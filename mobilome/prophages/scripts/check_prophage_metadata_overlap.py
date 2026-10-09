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

# ======================================
# INPUT
# ======================================

PROPHAGE_FILE = (
    PROPHAGE_MATRICES
    / "prophage_summary.tsv"
)

CLADE_FILE = GENOME_CLADES

# ======================================
# OUTPUT
# ======================================

PROPHAGE_STATS.mkdir(
    parents=True,
    exist_ok=True
)

MISSING_IN_CLADES = (
    PROPHAGE_STATS
    / "genomes_missing_in_clades.tsv"
)

MISSING_IN_PROPHAGES = (
    PROPHAGE_STATS
    / "genomes_missing_in_prophages.tsv"
)

OVERLAP_SUMMARY = (
    PROPHAGE_STATS
    / "prophage_metadata_overlap_summary.tsv"
)

# ======================================
# LOAD
# ======================================

prophages = pd.read_csv(
    PROPHAGE_FILE,
    sep="\t"
)

clades = pd.read_csv(
    CLADE_FILE,
    sep="\t"
)

# ======================================
# NORMALIZE GENOME NAMES
# ======================================

prophages["Genome"] = (
    prophages["Genome"]
    .astype(str)
    .str.replace("_", "-", regex=False)
)

clades["Genome"] = (
    clades["Genome"]
    .astype(str)
    .str.replace("_", "-", regex=False)
)

# ======================================
# SETS
# ======================================

prophage_genomes = set(
    prophages["Genome"]
)

clade_genomes = set(
    clades["Genome"]
)

# ======================================
# COMPARISONS
# ======================================

missing_in_clades = sorted(
    prophage_genomes - clade_genomes
)

missing_in_prophages = sorted(
    clade_genomes - prophage_genomes
)

overlap = sorted(
    prophage_genomes & clade_genomes
)

# ======================================
# SAVE DETAILS
# ======================================

pd.DataFrame({
    "Genome": missing_in_clades
}).to_csv(
    MISSING_IN_CLADES,
    sep="\t",
    index=False
)

pd.DataFrame({
    "Genome": missing_in_prophages
}).to_csv(
    MISSING_IN_PROPHAGES,
    sep="\t",
    index=False
)

# ======================================
# SUMMARY
# ======================================

summary = pd.DataFrame([{

    "Prophage_Genomes":
        len(prophage_genomes),

    "Clade_Genomes":
        len(clade_genomes),

    "Overlap":
        len(overlap),

    "Missing_In_Clades":
        len(missing_in_clades),

    "Missing_In_Prophages":
        len(missing_in_prophages)

}])

summary.to_csv(
    OVERLAP_SUMMARY,
    sep="\t",
    index=False
)

# ======================================
# REPORT
# ======================================

print("\n========== OVERLAP REPORT ==========")

print(
    f"Prophage genomes: "
    f"{len(prophage_genomes)}"
)

print(
    f"Clade genomes: "
    f"{len(clade_genomes)}"
)

print(
    f"Overlap: "
    f"{len(overlap)}"
)

print(
    f"Missing in clades: "
    f"{len(missing_in_clades)}"
)

print(
    f"Missing in prophages: "
    f"{len(missing_in_prophages)}"
)

print("\nSaved:")

print(MISSING_IN_CLADES)
print(MISSING_IN_PROPHAGES)
print(OVERLAP_SUMMARY)

print("\nFinished.")