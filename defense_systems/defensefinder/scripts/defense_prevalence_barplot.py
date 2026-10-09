#!/usr/bin/env python3

from pathlib import Path
import sys

SCRIPT_DIR = Path(__file__).resolve().parent

PROJECT_ROOT = SCRIPT_DIR.parents[2]

sys.path.append(
    str(PROJECT_ROOT / "scripts")
)

from config import *

import pandas as pd
import matplotlib.pyplot as plt

# ======================================
# INPUT
# ======================================

INPUT_FILE = (
    DEFENSEFINDER_STATS
    / "defense_prevalence.tsv"
)

# ======================================
# OUTPUT
# ======================================

DEFENSEFINDER_PLOTS.mkdir(
    parents=True,
    exist_ok=True
)

OUTFILE = (
    DEFENSEFINDER_PLOTS
    / "defense_prevalence.png"
)

# ======================================
# LOAD DATA
# ======================================

df = pd.read_csv(
    PREVALENCE_FILE,
    sep="\t"
)

df = df.sort_values(
    "Frequency_percent",
    ascending=True
)

# ======================================
# PLOT
# ======================================

plt.figure(
    figsize=(8, 10)
)

plt.barh(
    df["Subtype"],
    df["Frequency_percent"]
)

plt.xlabel(
    "Frequency (%)"
)

plt.ylabel(
    "Defense system"
)

plt.title(
    "Defense system prevalence across Xtt genomes"
)

plt.tight_layout()

plt.savefig(
    OUTFILE,
    dpi=300,
    bbox_inches="tight"
)

plt.close()

print("\nSaved:")
print(OUTFILE)