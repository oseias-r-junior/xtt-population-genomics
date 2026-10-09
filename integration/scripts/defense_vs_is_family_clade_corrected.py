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

import statsmodels.formula.api as smf

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
    / "defense_vs_is_family_clade_corrected.tsv"
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

clades = pd.read_csv(
    CLADE_FILE,
    sep="\t"
)

# =====================================
# NORMALIZE
# =====================================

for df in [defense, is_df, clades]:

    df["Genome"] = (
        df["Genome"]
        .astype(str)
        .apply(normalize_genome_name)
    )

# =====================================
# MERGE
# =====================================

df = (
    defense
    .merge(
        is_df,
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

N = len(df)

MIN_PRESENT = 5
MAX_PRESENT = N - 5

# =====================================
# VARIABLES
# =====================================

defense_systems = []

for system in defense.columns:

    if system == "Genome":
        continue

    present = (df[system] > 0).sum()

    if (
        present < MIN_PRESENT
        or
        present > MAX_PRESENT
    ):
        continue

    defense_systems.append(system)

is_families = []

for fam in is_df.columns:

    if fam == "Genome":
        continue

    present = (df[fam] > 0).sum()

    if (
        present < MIN_PRESENT
        or
        present > MAX_PRESENT
    ):
        continue

    is_families.append(fam)

results = []

print(
    f"Defense systems retained: {len(defense_systems)}"
)

print(
    f"IS families retained: {len(is_families)}"
)

print(defense_systems)
print(is_families)


# =====================================
# MODELS
# =====================================

for system in defense_systems:

    n_present = int(
        df[system].sum()
    )

    if n_present < 5:
        continue

    for family in is_families:

        try:

            model = smf.ols(
                formula=(
                    f"Q('{family}') ~ "
                    f"Q('{system}') + "
                    f"C(Clade)"
                ),
                data=df
            ).fit()

            coef_name = (
                f"Q('{system}')"
            )

            coef = model.params.get(
                coef_name,
                np.nan
            )

            pval = model.pvalues.get(
                coef_name,
                np.nan
            )

            results.append({
                "Defense_system": system,
                "IS_family": family,
                "N_present": n_present,
                "Effect_size": coef,
                "P_value": pval
            })

        except Exception as e:

            print(
                f"Skipping {system} vs {family}: {e}"
            )

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
print(
    f"Rows: {len(results_df)}"
)