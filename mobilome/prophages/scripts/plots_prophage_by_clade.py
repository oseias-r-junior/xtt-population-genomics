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

SUMMARY_FILE = (
    PROPHAGE_STATS
    / "prophage_clade_summary.tsv"
)

PRESENCE_FILE = (
    PROPHAGE_STATS
    / "prophage_clade_presence.tsv"
)

COMPLETENESS_FILE = (
    PROPHAGE_STATS
    / "prophage_completeness_by_clade.tsv"
)

LENGTH_FILE = (
    PROPHAGE_STATS
    / "prophage_length_by_clade.tsv"
)

# ======================================
# OUTPUT
# ======================================

PROPHAGE_PLOTS.mkdir(
    parents=True,
    exist_ok=True
)

# ======================================
# LOAD
# ======================================

summary = pd.read_csv(
    SUMMARY_FILE,
    sep="\t"
)

presence = pd.read_csv(
    PRESENCE_FILE,
    sep="\t"
)

# ======================================
# PLOT 1
# MEAN PROPHAGES / GENOME
# ======================================

plt.figure(figsize=(6,4))

plt.bar(
    summary["Clade"],
    summary["Mean_Prophages"],
    yerr=summary["SD_Prophages"],
    capsize=5
)

plt.ylabel(
    "Mean prophages per genome"
)

plt.xlabel(
    "Clade"
)

plt.tight_layout()

outfile = (
    PROPHAGE_PLOTS
    / "prophage_mean_count_by_clade.png"
)

plt.savefig(
    outfile,
    dpi=300
)

plt.close()

print(
    f"Saved: {outfile}"
)

# ======================================
# PLOT 2
# PRESENCE RATE
# ======================================

plt.figure(figsize=(6,4))

plt.bar(
    presence["Clade"],
    presence["Presence_Rate"]
)

plt.ylabel(
    "Prophage presence rate"
)

plt.xlabel(
    "Clade"
)

plt.ylim(0, 1)

plt.tight_layout()

outfile = (
    PROPHAGE_PLOTS
    / "prophage_presence_by_clade.png"
)

plt.savefig(
    outfile,
    dpi=300
)

plt.close()

print(
    f"Saved: {outfile}"
)

# ======================================
# OPTIONAL PLOTS
# ======================================

if LENGTH_FILE.exists():

    length_df = pd.read_csv(
        LENGTH_FILE,
        sep="\t"
    )

    plt.figure(figsize=(6,4))

    plt.bar(
        length_df["Clade"],
        length_df["Mean_Length_bp"],
        yerr=length_df["SD_Length_bp"],
        capsize=5
    )

    plt.ylabel(
        "Mean prophage length (bp)"
    )

    plt.xlabel(
        "Clade"
    )

    plt.tight_layout()

    outfile = (
        PROPHAGE_PLOTS
        / "prophage_length_by_clade.png"
    )

    plt.savefig(
        outfile,
        dpi=300
    )

    plt.close()

    print(
        f"Saved: {outfile}"
    )

# ======================================
# COMPLETENESS
# ======================================

if COMPLETENESS_FILE.exists():

    comp = pd.read_csv(
        COMPLETENESS_FILE,
        sep="\t"
    )

    categories = [

        "Intact",
        "Questionable",
        "Incomplete"

    ]

    comp = comp.set_index(
        "Clade"
    )

    comp[categories].plot(
        kind="bar",
        figsize=(7,4)
    )

    plt.ylabel(
        "Number of prophages"
    )

    plt.xlabel(
        "Clade"
    )

    plt.tight_layout()

    outfile = (
        PROPHAGE_PLOTS
        / "prophage_completeness_by_clade.png"
    )

    plt.savefig(
        outfile,
        dpi=300
    )

    plt.close()

    print(
        f"Saved: {outfile}"
    )

print("\nFinished.")