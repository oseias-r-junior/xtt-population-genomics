#!/usr/bin/env python3

from pathlib import Path
import sys

SCRIPT_DIR = Path(__file__).resolve().parent

PROJECT_ROOT = SCRIPT_DIR.parents[2]

sys.path.append(
    str(PROJECT_ROOT / "scripts")
)

from config import *
from normalize_genome_names import normalize_genome_name

import pandas as pd

# ======================================
# INPUT
# ======================================

REGIONS_FILE = (
    PROPHAGE_MATRICES
    / "prophage_regions.tsv"
)

CLADE_FILE = GENOME_CLADES

# ======================================
# OUTPUT
# ======================================

PROPHAGE_STATS.mkdir(
    parents=True,
    exist_ok=True
)

COUNTS_OUTFILE = (
    PROPHAGE_STATS
    / "prophage_family_counts.tsv"
)

CLADE_OUTFILE = (
    PROPHAGE_STATS
    / "prophage_family_by_clade.tsv"
)

PREVALENCE_OUTFILE = (
    PROPHAGE_STATS
    / "prophage_family_prevalence.tsv"
)

# ======================================
# LOAD
# ======================================

regions = pd.read_csv(
    REGIONS_FILE,
    sep="\t"
)

clades = pd.read_csv(
    CLADE_FILE,
    sep="\t"
)

# ======================================
# NORMALIZE
# ======================================

regions["Genome"] = (
    regions["Genome"]
    .astype(str)
    .apply(normalize_genome_name)
)

clades["Genome"] = (
    clades["Genome"]
    .astype(str)
    .apply(normalize_genome_name)
)

# ======================================
# MERGE
# ======================================

df = regions.merge(
    clades,
    on="Genome",
    how="inner"
)

print(
    f"Merged prophage regions: {len(df)}"
)

# ======================================
# FAMILY COUNTS
# ======================================

family_counts = (
    df["Most_Common_Phage"]
    .value_counts()
    .reset_index()
)

family_counts.columns = [
    "Family",
    "Count"
]

family_counts.to_csv(
    COUNTS_OUTFILE,
    sep="\t",
    index=False
)

# ======================================
# FAMILY x CLADE
# ======================================

family_clade = pd.crosstab(
    df["Most_Common_Phage"],
    df["Clade"]
)

family_clade = (
    family_clade
    .reset_index()
)

family_clade.to_csv(
    CLADE_OUTFILE,
    sep="\t",
    index=False
)

# ======================================
# PREVALENCE
# ======================================

n_per_clade = (
    clades["Clade"]
    .value_counts()
    .to_dict()
)

presence = (
    df[
        [
            "Genome",
            "Clade",
            "Most_Common_Phage"
        ]
    ]
    .drop_duplicates()
)

prev = pd.crosstab(
    presence["Most_Common_Phage"],
    presence["Clade"]
)

for clade in prev.columns:

    prev[clade] = (
        prev[clade]
        / n_per_clade[clade]
        * 100
    )

prev = (
    prev
    .round(2)
    .reset_index()
)

prev.to_csv(
    PREVALENCE_OUTFILE,
    sep="\t",
    index=False
)

# ======================================
# REPORT
# ======================================

print(
    f"Saved: {COUNTS_OUTFILE}"
)

print(
    f"Saved: {CLADE_OUTFILE}"
)

print(
    f"Saved: {PREVALENCE_OUTFILE}"
)

print("\nFinished.")