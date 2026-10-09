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

# ======================================
# INPUT
# ======================================

PHASTEST_FILE = (
    PROPHAGE_MATRICES
    / "prophage_summary.tsv"
)

VIRSORTER_FILE = (
    VIRSORTER2_MATRICES
    / "virsorter2_summary_matrix.tsv"
)

# ======================================
# OUTPUT
# ======================================

OUTFILE = (
    PROPHAGE_STATS
    / "virsorter2_phastest_disagreement.tsv"
)

OUTFILE.parent.mkdir(
    parents=True,
    exist_ok=True
)

# ======================================
# LOAD
# ======================================

phastest = pd.read_csv(
    PHASTEST_FILE,
    sep="\t"
)

virsorter = pd.read_csv(
    VIRSORTER_FILE,
    sep="\t"
)

# ======================================
# NORMALIZE
# ======================================

for df in [phastest, virsorter]:

    df["Genome"] = (
        df["Genome"]
        .astype(str)
        .apply(normalize_genome_name)
    )

# # ======================================
# # KEEP COLUMNS
# # ======================================

# phastest = phastest[
#     [
#         "Genome",
#         "Total_Prophages",
#         "Intact",
#         "Questionable",
#         "Incomplete"
#     ]
# ]

# virsorter = virsorter[
#     [
#         "Genome",
#         "Total_sequences",
#         "Full",
#         "Partial",
#         "Short",
#         "Total_length",
#         "Mean_length",
#         "Mean_score",
#         "Max_score",
#         "Mean_hallmark",
#         "Mean_viral_pct",
#         "Mean_cellular_pct"
#     ]
# ]

# ======================================
# MERGE
# ======================================

df = (
    phastest
    .merge(
        virsorter,
        on="Genome",
        how="outer"
    )
)

# ======================================
# FILL MISSING
# ======================================

numeric_cols = (
    df
    .select_dtypes(include=np.number)
    .columns
)

df[numeric_cols] = (
    df[numeric_cols]
    .fillna(0)
)


# ======================================
# DIFFERENCE
# ======================================

df["Difference"] = (
    df["Total_sequences"]
    -
    df["Total_Prophages"]
)

df["Agreement"] = (
    df["Difference"] == 0
)

# ======================================
# DETECTOR BIAS
# ======================================

df["Detector_bias"] = np.where(

    df["Difference"] > 0,
    "VirSorter2",

    np.where(
        df["Difference"] < 0,
        "PHASTEST",
        "Equal"
    )

)

# ======================================
# DISCORDANT ONLY
# ======================================

discordant = (
    df[
        df["Agreement"] == False
    ]
    .copy()
)

discordant["Abs_difference"] = (
    discordant["Difference"]
    .abs()
)

discordant = (
    discordant
    .sort_values(
        [
            "Abs_difference",
            "Difference"
        ],
        ascending=False
    )
)

discordant["Reason"] = np.where(
    discordant["Intact"] > 0,
    "Contains intact prophage",
    np.where(
        discordant["Questionable"] > 0,
        "Only questionable prophages",
        np.where(
            discordant["Incomplete"] > 0,
            "Only incomplete prophages",
            "No PHASTEST prediction"
        )
    )
)

# ======================================
# SAVE
# ======================================

discordant.to_csv(
    OUTFILE,
    sep="\t",
    index=False
)

# ======================================
# REPORT
# ======================================

print()

print(
    f"Genomes compared: {len(df)}"
)

print(
    f"Discordant genomes: {len(discordant)}"
)

print(
    f"Agreement: {len(df)-len(discordant)} / {len(df)}"
)

print()

print(
    "VirSorter2 > PHASTEST:",
    (discordant["Difference"] > 0).sum()
)

print(
    "PHASTEST > VirSorter2:",
    (discordant["Difference"] < 0).sum()
)

print()

print(OUTFILE)