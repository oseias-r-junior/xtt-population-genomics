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
import statsmodels.formula.api as smf
from patsy import *

# =====================================
# INPUTS
# =====================================

DEFENSE_MATRIX = (
    INTEGRATION_MATRICES
    / "defensefinder_matrix.tsv"
)

PROPHAGE_FILE = (
    PROPHAGE_MATRICES
    / "prophage_family_presence_absence.tsv"
)

CLADE_FILE = GENOME_CLADES

# =====================================
# OUTPUT
# =====================================

OUTDIR = (
    INTEGRATION_STATS
    / "clade_corrected"
)

OUTDIR.mkdir(
    parents=True,
    exist_ok=True
)

OUTFILE = (
    OUTDIR
    / "defense_vs_prophage_clade_corrected.tsv"
)

# =====================================
# LOAD
# =====================================

defense = pd.read_csv(
    DEFENSE_MATRIX,
    sep="\t"
)

prophage = pd.read_csv(
    PROPHAGE_FILE,
    sep="\t"
)

clades = pd.read_csv(
    CLADE_FILE,
    sep="\t"
)

defense["Genome"] = (
    defense["Genome"]
    .astype(str)
    .apply(normalize_genome_name)
)

prophage["Genome"] = (
    prophage["Genome"]
    .astype(str)
    .apply(normalize_genome_name)
)

clades["Genome"] = (
    clades["Genome"]
    .astype(str)
    .apply(normalize_genome_name)
)

# =====================================
# TOTAL PROPHAGES
# =====================================

prophage["Total_prophages"] = (
    prophage.drop(columns=["Genome"])
    .sum(axis=1)
)

# =====================================
# MERGE
# =====================================

df = (
    defense
    .merge(
        prophage[
            ["Genome", "Total_prophages"]
        ],
        on="Genome"
    )
    .merge(
        clades,
        on="Genome"
    )
)

print(
    f"Genomes after merge: {len(df)}"
)

# =====================================
# SYSTEMS
# =====================================

systems = [
    c
    for c in defense.columns
    if c != "Genome"
]

results = []

# =====================================
# MODEL
# =====================================

safe_systems = [
    c
    for c in defense.columns
    if c != "Genome"
]

for system in systems:

    present = df[system].sum()

    if present < 5:
        continue

    if df[system].nunique() < 2:
        continue

    try:

        model = smf.ols(
            formula=f"Total_prophages ~ Q('{system}') + C(Clade)",
            data=df
        ).fit()

        term = f"Q('{system}')"

        coef = model.params.get(
            term,
            float("nan")
        )

        pval = model.pvalues.get(
            term,
            float("nan")
        )

        results.append({
            "Subtype": system,
            "N_present": int(present),
            "Effect_size": coef,
            "P_value": pval
        })

    except Exception as e:

        print(f"Skipping {system}: {e}")

# =====================================
# SAVE
# =====================================

results_df = (
    pd.DataFrame(results)
      .sort_values("P_value")
)

results_df.to_csv(
    OUTFILE,
    sep="\t",
    index=False
)

print("\nSaved:")
print(OUTFILE)
print(f"Rows: {len(results_df)}")