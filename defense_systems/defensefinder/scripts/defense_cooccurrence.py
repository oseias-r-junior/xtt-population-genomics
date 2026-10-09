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
from statsmodels.stats.multitest import multipletests

# =====================================
# INPUT
# =====================================

MATRIX_FILE = (
    DEFENSEFINDER_MATRICES
    / "defensefinder_matrix.tsv"
)

# =====================================
# OUTPUT
# =====================================

OUTDIR = (
    DEFENSEFINDER_STATS
    / "cooccurrence"
)

OUTDIR.mkdir(
    parents=True,
    exist_ok=True
)

OUT_ALL = (
    OUTDIR
    / "defense_cooccurrence.tsv"
)

OUT_SIG = (
    OUTDIR
    / "defense_cooccurrence_significant.tsv"
)

# =====================================
# LOAD
# =====================================

df = pd.read_csv(
    MATRIX_FILE,
    sep="\t"
)

systems = [
    c
    for c in df.columns
    if c != "Genome"
]

# =====================================
# FILTER RARE SYSTEMS
# =====================================

keep = []

for s in systems:

    prevalence = df[s].sum()

    if prevalence >= 5:
        keep.append(s)

systems = sorted(keep)

print(
    f"Systems retained: {len(systems)}"
)

# =====================================
# PAIRWISE TESTS
# =====================================

results = []

for i in range(len(systems)):

    s1 = systems[i]

    for j in range(i + 1, len(systems)):

        s2 = systems[j]

        a = (
            (df[s1] == 1)
            &
            (df[s2] == 1)
        ).sum()

        b = (
            (df[s1] == 1)
            &
            (df[s2] == 0)
        ).sum()

        c = (
            (df[s1] == 0)
            &
            (df[s2] == 1)
        ).sum()

        d = (
            (df[s1] == 0)
            &
            (df[s2] == 0)
        ).sum()

        table = [
            [a, b],
            [c, d]
        ]

        try:

            odds_ratio, p = fisher_exact(
                table,
                alternative="two-sided"
            )

        except:

            odds_ratio = float("nan")
            p = 1.0

        results.append({
            "System_A": s1,
            "System_B": s2,
            "Cooccur": a,
            "Odds_ratio": odds_ratio,
            "P_value": p
        })

# =====================================
# FDR
# =====================================

results_df = pd.DataFrame(
    results
)

results_df["FDR"] = multipletests(
    results_df["P_value"],
    method="fdr_bh"
)[1]

results_df = (
    results_df
    .sort_values(
        "FDR"
    )
)

# =====================================
# SAVE
# =====================================

results_df.to_csv(
    OUT_ALL,
    sep="\t",
    index=False
)

sig_df = results_df.loc[
    results_df["FDR"] < 0.05
]

sig_df.to_csv(
    OUT_SIG,
    sep="\t",
    index=False
)

print(
    f"Total pairs: {len(results_df)}"
)

print(
    f"Significant pairs: {len(sig_df)}"
)

print(
    f"Saved:\n{OUT_ALL}\n{OUT_SIG}"
)