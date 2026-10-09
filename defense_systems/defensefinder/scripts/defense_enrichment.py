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
from scipy.stats import fisher_exact

# =========================
# INPUT
# =========================

MATRIX_FILE = (
    DEFENSEFINDER_MATRICES
    / "defensefinder_matrix.tsv"
)

CLADE_FILE = GENOME_CLADES

OUTDIR = (
    DEFENSEFINDER_STATS
    / "enrichment"
)

OUTDIR.mkdir(
    parents=True,
    exist_ok=True
)

OUTFILE = (
    OUTDIR
    / "defense_enrichment.tsv"
)

# =========================
# LOAD
# =========================

df = pd.read_csv(
    MATRIX_FILE,
    sep="\t"
)

clades = pd.read_csv(
    CLADE_FILE,
    sep="\t"
)

df["Genome"] = (
    df["Genome"]
    .astype(str)
    .apply(normalize_genome_name)
)

clades["Genome"] = (
    clades["Genome"]
    .astype(str)
    .apply(normalize_genome_name)
)

merged = pd.merge(
    df,
    clades,
    on="Genome",
    how="inner"
)

print(
    f"Genomes after merge: {len(merged)}"
)

# =========================
# TESTS
# =========================

subtypes = [
    c
    for c in merged.columns
    if c not in ["Genome", "Clade"]
]

results = []

for subtype in subtypes:

    for clade in sorted(
        merged["Clade"].unique()
    ):

        in_clade = (
            merged["Clade"] == clade
        )

        present_in = (
            merged.loc[
                in_clade,
                subtype
            ].sum()
        )

        absent_in = (
            in_clade.sum()
            - present_in
        )

        present_out = (
            merged.loc[
                ~in_clade,
                subtype
            ].sum()
        )

        absent_out = (
            (~in_clade).sum()
            - present_out
        )

        table = [
            [present_in, absent_in],
            [present_out, absent_out]
        ]

        odds_ratio, p = fisher_exact(
            table
        )

        results.append({
            "Subtype": subtype,
            "Clade": clade,
            "Present_in_clade": present_in,
            "Present_outside": present_out,
            "Odds_ratio": odds_ratio,
            "P_value": p
        })

results_df = pd.DataFrame(
    results
)

results_df = results_df.sort_values(
    "P_value"
)

results_df.to_csv(
    OUTFILE,
    sep="\t",
    index=False
)

print("\nSaved:")
print(OUTFILE)

print(
    f"Rows: {len(results_df)}"
)