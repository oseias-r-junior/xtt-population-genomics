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

# from xtt_population_genomics.scripts.config import PROPHAGE_PLOTS

import pandas as pd
import matplotlib.pyplot as plt

# ======================================
# INPUT
# ======================================

PREVALENCE_FILE = (
    PROPHAGE_STATS
    / "prophage_family_prevalence.tsv"
)

ENRICHMENT_FILE = (
    PROPHAGE_STATS
    / "prophage_family_enrichment.tsv"
)

# ======================================
# OUTPUT
# ======================================

OUTFILE = (
    PROPHAGE_PLOTS
    / "prophage_family_heatmap.png"
)

# ======================================
# LOAD
# ======================================

heat = pd.read_csv(
    PREVALENCE_FILE,
    sep="\t"
)

enrich = pd.read_csv(
    ENRICHMENT_FILE,
    sep="\t"
)
# ======================================
# FDR CORRECTION
# ======================================

order = (

    enrich
    .groupby("Family")["FDR"]
    .min()
    .sort_values()
    .index

)

heat = heat.set_index(
    "Most_Common_Phage"
)

valid_families = [

    fam
    for fam in order
    if fam in heat.index

]

heat = heat.loc[
    valid_families
]

# ======================================
# PLOT
# ======================================
fig, ax = plt.subplots(
    figsize=(8, 6)
)

im = ax.imshow(
    heat,
    aspect="auto"
)

ax.set_xticks(
    range(len(heat.columns))
)

ax.set_xticklabels(
    heat.columns
)

ax.set_yticks(
    range(len(heat.index))
)

ax.set_yticklabels(
    heat.index
)

plt.colorbar(
    im,
    ax=ax,
    label="Prevalence (%)"
)

plt.tight_layout()

plt.savefig(
    OUTFILE,
    dpi=300,
    bbox_inches="tight"
)

plt.close()

# ======================================
# SAVE
# ======================================

print(
    f"Saved: {OUTFILE}"
)

print("\nFinished.")
