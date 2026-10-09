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

MATRIX_FILE = (
    DEFENSEFINDER_MATRICES
    / "defensefinder_matrix.tsv"
)

# ======================================
# OUTPUT
# ======================================

DEFENSEFINDER_STATS.mkdir(
    parents=True,
    exist_ok=True
)

PREVALENCE_OUT = (
    DEFENSEFINDER_STATS
    / "defense_prevalence.tsv"
)

GENOME_OUT = (
    DEFENSEFINDER_STATS
    / "defense_systems_per_genome.tsv"
)

# ======================================
# LOAD MATRIX
# ======================================

df = pd.read_csv(
    MATRIX_FILE,
    sep="\t"
)

feature_cols = [
    c for c in df.columns
    if c != "Genome"
]

n_genomes = len(df)

# ======================================
# PREVALENCE
# ======================================

prevalence = pd.DataFrame({

    "Subtype": feature_cols,

    "Count": [
        int(df[c].sum())
        for c in feature_cols
    ]

})

prevalence["Frequency_percent"] = (
    prevalence["Count"]
    / n_genomes
    * 100
).round(2)

prevalence = prevalence.sort_values(
    "Count",
    ascending=False
)

prevalence.to_csv(
    PREVALENCE_OUT,
    sep="\t",
    index=False
)

# ======================================
# SYSTEMS PER GENOME
# ======================================

genome_stats = pd.DataFrame({

    "Genome": df["Genome"],

    "Total_defense_systems":
        df[feature_cols].sum(axis=1)

})

genome_stats = genome_stats.sort_values(
    "Total_defense_systems",
    ascending=False
)

genome_stats.to_csv(
    GENOME_OUT,
    sep="\t",
    index=False
)

# ======================================
# SUMMARY
# ======================================

print("\n===== DefenseFinder summary =====")

print(f"Genomes: {n_genomes}")

print(
    f"Defense systems: {len(feature_cols)}"
)

print(
    f"Mean systems/genome: "
    f"{genome_stats['Total_defense_systems'].mean():.2f}"
)

print(
    f"Min systems/genome: "
    f"{genome_stats['Total_defense_systems'].min()}"
)

print(
    f"Max systems/genome: "
    f"{genome_stats['Total_defense_systems'].max()}"
)

print("\nSaved:")

print(PREVALENCE_OUT)

print(GENOME_OUT)