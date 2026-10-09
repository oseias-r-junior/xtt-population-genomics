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
# from statsmodels.stats.multitest import multipletests

# ======================================
# INPUT
# ======================================

MATRIX_FILE = (
    PROPHAGE_MATRICES
    / "prophage_family_presence_absence.tsv"
)

CLADE_FILE = GENOME_CLADES

# ======================================
# OUTPUT
# ======================================

OUTFILE = (
    PROPHAGE_STATS
    / "prophage_family_enrichment.tsv"
)

# ======================================
# LOAD
# ======================================

matrix = pd.read_csv(
    MATRIX_FILE,
    sep="\t"
)

clades = pd.read_csv(
    CLADE_FILE,
    sep="\t"
)

# ======================================
# NORMALIZE
# ======================================

matrix["Genome"] = (
    matrix["Genome"]
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

df = matrix.merge(
    clades,
    on="Genome",
    how="inner"
)

# ======================================
# FISHER'S EXACT TEST
# ======================================

results = []

families = [
    c
    for c in df.columns
    if c not in [
        "Genome",
        "Clade"
    ]
]

for family in families:

    for target_clade in [
        "K0",
        "K1",
        "K2"
    ]:

        present_target = len(
            df[
                (df["Clade"] == target_clade)
                &
                (df[family] == 1)
            ]
        )

        absent_target = len(
            df[
                (df["Clade"] == target_clade)
                &
                (df[family] == 0)
            ]
        )

        present_other = len(
            df[
                (df["Clade"] != target_clade)
                &
                (df[family] == 1)
            ]
        )

        absent_other = len(
            df[
                (df["Clade"] != target_clade)
                &
                (df[family] == 0)
            ]
        )

        table = [

            [
                present_target,
                absent_target
            ],

            [
                present_other,
                absent_other
            ]

        ]

        odds_ratio, pvalue = fisher_exact(
            table
        )

        results.append({

            "Family":
                family,

            "Target_Clade":
                target_clade,

            "Present_Target":
                present_target,

            "Present_Other":
                present_other,

            "Odds_Ratio":
                odds_ratio,

            "P_value":
                pvalue

        })

# ======================================
# SAVE
# ======================================

out = pd.DataFrame(results)

out = out.sort_values(
    "P_value"
)

m = len(out)

out["Rank"] = range(
    1,
    m + 1
)

out["FDR"] = (
    out["P_value"]
    * m
    / out["Rank"]
)

out["FDR"] = (
    out["FDR"]
    .clip(upper=1)
)

out = out.drop(
    columns="Rank"
)

out.to_csv(
    OUTFILE,
    sep="\t",
    index=False
)

print(
    f"Saved: {OUTFILE}"
)

print("\nFinished.")