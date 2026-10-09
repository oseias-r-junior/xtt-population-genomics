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

# ======================================
# OUTPUT
# ======================================

OUTFILE = (
    PROPHAGE_MATRICES
    / "prophage_family_presence_absence.tsv"
)

# ======================================
# LOAD
# ======================================

df = pd.read_csv(
    REGIONS_FILE,
    sep="\t"
)

# ======================================
# NORMALIZE
# ======================================

df["Genome"] = (
    df["Genome"]
    .astype(str)
    .apply(normalize_genome_name)
)

# ======================================
# BUILD MATRIX
# ======================================

matrix = pd.crosstab(

    df["Genome"],
    df["Most_Common_Phage"]

)

# convert counts -> presence/absence

matrix = (
    matrix > 0
).astype(int)

# ======================================
# ADD GENOMES WITHOUT PROPHAGES
# ======================================

clades = pd.read_csv(
    GENOME_CLADES,
    sep="\t"
)

clades["Genome"] = (
    clades["Genome"]
    .astype(str)
    .apply(normalize_genome_name)
)

all_genomes = (
    clades["Genome"]
    .drop_duplicates()
)

matrix = matrix.reindex(
    all_genomes,
    fill_value=0
)

matrix = matrix.reset_index()

matrix.columns.values[0] = "Genome"

# ======================================
# SAVE
# ======================================

matrix.to_csv(

    OUTFILE,
    sep="\t",
    index=False

)

print(
    f"Genomes: {len(matrix)}"
)

print(
    f"Phage families: {len(matrix.columns)-1}"
)

print(
    f"Saved: {OUTFILE}"
)

print("\nFinished.")