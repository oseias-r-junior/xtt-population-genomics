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

FAMILY_FILE = (
    PROPHAGE_STATS
    / "prophage_family_by_clade.tsv"
)

# ======================================
# OUTPUT
# ======================================

OUTFILE = (
    PROPHAGE_STATS
    / "prophage_key_findings.tsv"
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

families = pd.read_csv(
    FAMILY_FILE,
    sep="\t"
)

# ======================================
# DOMINANT FAMILY PER CLADE
# ======================================

dominant = {}

for clade in families.columns[1:]:

    tmp = families[
        [
            "Most_Common_Phage",
            clade
        ]
    ].copy()

    tmp = tmp.sort_values(
        clade,
        ascending=False
    )

    dominant[clade] = (
        tmp.iloc[0]["Most_Common_Phage"]
    )

# ======================================
# BUILD TABLE
# ======================================

out = summary.merge(
    presence,
    on="Clade"
)

out["Dominant_Family"] = (
    out["Clade"]
    .map(dominant)
)

# ======================================
# SAVE
# ======================================

out.to_csv(
    OUTFILE,
    sep="\t",
    index=False
)

print(
    f"Saved: {OUTFILE}"
)

print("\nFinished.")