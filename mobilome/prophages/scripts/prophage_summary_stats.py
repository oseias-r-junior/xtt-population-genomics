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

REGIONS_FILE = (
    PROPHAGE_MATRICES
    / "prophage_regions.tsv"
)

CLADE_FILE = GENOME_CLADES

PROPHAGE_MATRICES.mkdir(
    parents=True,
    exist_ok=True
)

OUTFILE = (
    PROPHAGE_MATRICES
    / "prophage_summary.tsv"
)

regions = pd.read_csv(
    REGIONS_FILE,
    sep="\t"
)

clades = pd.read_csv(
    CLADE_FILE,
    sep="\t"
)

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

counts = (
    regions
    .groupby("Genome")
    .size()
    .rename("Total_Prophages")
    .reset_index()
)

all_genomes = (
    clades["Genome"]
    .drop_duplicates()
    .to_frame()
)

summary = all_genomes.merge(
    counts,
    on="Genome",
    how="left"
)

summary["Total_Prophages"] = (
    summary["Total_Prophages"]
    .fillna(0)
    .astype(int)
)

summary.to_csv(
    OUTFILE,
    sep="\t",
    index=False
)

print(f"Genomes: {len(summary)}")
print(f"Genomes with prophages: {(summary['Total_Prophages'] > 0).sum()}")
print(summary.head(10))
print(f"Saved summary to {OUTFILE}")
