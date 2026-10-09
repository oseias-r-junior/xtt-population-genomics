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

# =====================================
# PATHS
# =====================================

RAW_DIR = DEFENSEFINDER_RAW

OUT_DIR = DEFENSEFINDER_MATRICES

OUT_DIR.mkdir(
    parents=True,
    exist_ok=True
)

OUTFILE = (
    DEFENSEFINDER_MATRICES
    / "defensefinder_matrix.tsv"
)

# =====================================
# FIND FILES
# =====================================

files = sorted(
    RAW_DIR.glob(
        "*/*_defense_finder_systems.tsv"
    )
)

print(f"Found {len(files)} system files")

# =====================================
# COLLECT SYSTEMS
# =====================================

all_subtypes = set()

genome_to_subtypes = {}

for f in files:

    genome = normalize_genome_name(
        f.parent.name
    )

    df = pd.read_csv(
        f,
        sep="\t"
    )

    if df.empty:

        genome_to_subtypes[genome] = set()
        continue

    subtypes = set(
        df["subtype"]
        .dropna()
        .astype(str)
    )

    genome_to_subtypes[genome] = subtypes

    all_subtypes.update(
        subtypes
    )

all_subtypes = sorted(
    all_subtypes
)

print(
    f"Detected {len(all_subtypes)} unique defense systems"
)

# =====================================
# BUILD MATRIX
# =====================================

rows = []

for genome in sorted(
    genome_to_subtypes
):

    row = {
        "Genome": genome
    }

    present = genome_to_subtypes[
        genome
    ]

    for subtype in all_subtypes:

        row[subtype] = int(
            subtype in present
        )

    rows.append(
        row
    )

matrix = pd.DataFrame(
    rows
)

# =====================================
# SAVE
# =====================================

matrix.to_csv(
    OUTFILE,
    sep="\t",
    index=False
)

print(
    f"Saved: {OUTFILE}"
)

print(
    f"Genomes: {matrix.shape[0]}"
)

print(
    f"Defense systems: {matrix.shape[1] - 1}"
)