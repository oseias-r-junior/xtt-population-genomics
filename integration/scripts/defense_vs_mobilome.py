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

import matplotlib.pyplot as plt
import seaborn as sns

# =====================================
# INPUT
# =====================================

DEFENSE_FILE = (
    INTEGRATION_MATRICES
    / "defensefinder_matrix.tsv"
)

PROPHAGE_FILE = (
    PROPHAGE_MATRICES
    / "prophage_family_presence_absence.tsv"
)

IS_FILE = (
    IS_MATRICES
    / "is_family_matrix.tsv"
)

# CLADE_FILE = (
#     CLADE_METADATA
#     / "genome_clades.tsv"
# )

CLADE_FILE = GENOME_CLADES

# =====================================
# OUTPUT
# =====================================

OUTFILE = (
    INTEGRATION_STATS
    / "clade_corrected"
    / "defense_vs_mobilome.tsv"
)

OUTPLOT = (
    INTEGRATION_PLOTS
    / "defense_vs_mobilome_effects.png"
)

OUTFILE.parent.mkdir(
    parents=True,
    exist_ok=True
)

# =====================================
# LOAD
# =====================================

defense = pd.read_csv(
    DEFENSE_FILE,
    sep="\t"
)

prophage = pd.read_csv(
    PROPHAGE_FILE,
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
# NORMALIZE NAMES
# =====================================

for df_ in [defense, prophage, is_df, clades]:

    df_["Genome"] = (
        df_["Genome"]
        .astype(str)
        .apply(normalize_genome_name)
    )

# =====================================
# TOTAL PROPHAGES
# =====================================

prophage["Total_prophages"] = (
    prophage
    .drop(columns=["Genome"])
    .sum(axis=1)
)

# =====================================
# KEEP ONLY TOTAL IS
# =====================================

is_df = is_df[
    ["Genome", "Total_IS"]
]

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

# =====================================
# MOBILOME INDEX
# =====================================

df["Total_Mobilome"] = (
    df["Total_prophages"]
    +
    df["Total_IS"]
)

# =====================================
# FILTER SYSTEMS
# =====================================

N = len(df)

MIN_PRESENT = 5
MAX_PRESENT = N - 5

systems = []

for system in defense.columns:

    if system == "Genome":
        continue

    present = (
        df[system] > 0
    ).sum()

    if present < MIN_PRESENT:
        continue

    if present > MAX_PRESENT:
        continue

    systems.append(system)

print(
    f"Systems retained: {len(systems)}"
)

# =====================================
# MODELS
# =====================================

results = []

for system in systems:

    present = (
        df[system] > 0
    ).sum()

    try:

        model = smf.ols(
            formula=(
                f"Total_Mobilome ~ Q('{system}') + C(Clade)"
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
            "N_present": int(present),
            "Mobilome_effect": coef,
            "P_value": pval
        })

    except Exception as e:

        print(
            f"Skipping {system}: {e}"
        )

# =====================================
# SAVE TABLE
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

# =====================================
# PLOT
# =====================================

plot_df = (
    results_df
    .dropna()
    .sort_values(
        "Mobilome_effect"
    )
)

plt.figure(
    figsize=(8, 6)
)

sns.barplot(
    data=plot_df,
    x="Mobilome_effect",
    y="Defense_system"
)

plt.axvline(
    0,
    color="black",
    linestyle="--",
    linewidth=1
)

plt.xlabel(
    "Effect on mobilome load"
)

plt.ylabel(
    "Defense system"
)

plt.title(
    "Defense systems associated with mobilome abundance"
)

plt.tight_layout()

plt.savefig(
    OUTPLOT,
    dpi=300
)

plt.close()

# =====================================
# REPORT
# =====================================

print("\nSaved:")
print(OUTFILE)
print(OUTPLOT)
print(f"Rows: {len(results_df)}")