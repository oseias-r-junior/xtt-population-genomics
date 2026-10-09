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
import numpy as np

from scipy.stats import fisher_exact

from statsmodels.stats.multitest import (
    multipletests
)

# =====================================
# INPUT
# =====================================

DEFENSE_FILE = (
    INTEGRATION_MATRICES
    / "defensefinder_matrix.tsv"
)

IS_FILE = (
    IS_MATRICES
    / "is_family_matrix.tsv"
)

# =====================================
# OUTPUT
# =====================================

OUTDIR = (
    INTEGRATION_STATS
    / "defense_is_association"
)

OUTDIR.mkdir(
    parents=True,
    exist_ok=True
)

OUTFILE = (
    OUTDIR
    / "defense_vs_is_family.tsv"
)

# =====================================
# LOAD
# =====================================

defense = pd.read_csv(
    DEFENSE_FILE,
    sep="\t"
)

is_df = pd.read_csv(
    IS_FILE,
    sep="\t"
)

defense["Genome"] = (
    defense["Genome"]
    .astype(str)
    .apply(normalize_genome_name)
)

is_df["Genome"] = (
    is_df["Genome"]
    .astype(str)
    .apply(normalize_genome_name)
)

# =====================================
# MERGE
# =====================================

df = defense.merge(
    is_df,
    on="Genome"
)

print(
    f"Genomes after merge: {len(df)}"
)

# =====================================
# FEATURES
# =====================================

defense_cols = [
    c
    for c in defense.columns
    if c != "Genome"
]

is_cols = [
    c
    for c in is_df.columns
    if c != "Genome"
]

print(
    f"Defense systems: {len(defense_cols)}"
)

print(
    f"IS families: {len(is_cols)}"
)

# =====================================
# BINARIZE IS
# =====================================

for col in is_cols:

    df[col] = (
        df[col] > 0
    ).astype(int)

# =====================================
# TESTS
# =====================================

results = []

for defense_system in defense_cols:

    defense_present = df[
        defense_system
    ].sum()

    if defense_present < 5:
        continue

    for is_family in is_cols:

        is_present = df[
            is_family
        ].sum()

        if is_present < 5:
            continue

        a = (
            (df[defense_system] == 1)
            &
            (df[is_family] == 1)
        ).sum()

        b = (
            (df[defense_system] == 1)
            &
            (df[is_family] == 0)
        ).sum()

        c = (
            (df[defense_system] == 0)
            &
            (df[is_family] == 1)
        ).sum()

        d = (
            (df[defense_system] == 0)
            &
            (df[is_family] == 0)
        ).sum()

        table = [
            [a, b],
            [c, d]
        ]

        try:

            OR, p = fisher_exact(
                table
            )

        except Exception:
            continue

        results.append({
            "Defense_system": defense_system,
            "IS_family": is_family,
            "Cooccur": a,
            "Odds_ratio": OR,
            "P_value": p
        })

# =====================================
# DATAFRAME
# =====================================

results_df = pd.DataFrame(
    results
)

# =====================================
# FDR
# =====================================

results_df["FDR"] = multipletests(
    results_df["P_value"],
    method="fdr_bh"
)[1]

# =====================================
# SORT
# =====================================

results_df = results_df.sort_values(
    [
        "FDR",
        "Odds_ratio"
    ],
    ascending=[
        True,
        False
    ]
)

# =====================================
# SAVE
# =====================================

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