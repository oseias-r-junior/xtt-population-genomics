#!/usr/bin/env python3

from pathlib import Path
import sys
import numpy as np

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

ENRICHMENT_FILE = (
    PROPHAGE_STATS
    / "prophage_family_enrichment.tsv"
)

# ======================================
# OUTPUT
# ======================================

OUTFILE = (
    PROPHAGE_PLOTS
    / "prophage_family_enrichment.png"
)

# ======================================
# LOAD
# ======================================

df = pd.read_csv(
    ENRICHMENT_FILE,
    sep="\t"
)

# ======================================
# FILTER
# ======================================

df = df[
    (df["FDR"] < 0.05)
    &
    (df["Odds_Ratio"] > 1)
]

if len(df) == 0:

    print(
        "No enriched families found."
    )

    sys.exit()

# ======================================
# LABELS
# ======================================

df["Label"] = (

    df["Family"]
    + "\n("
    + df["Target_Clade"]
    + ")"

)

# ======================================
# EFFECT SIZE
# ======================================

mapping = {
    "PHAGE_Stenot_Smp131_NC_023588": "Smp131",
    "PHAGE_Pseudo_phiPSA1_NC_024365": "phiPSA1",
    "PHAGE_Stenot_PSH1_NC_010429": "PSH1"
}

plot_df = df.sort_values(
    "FDR",
    ascending=True
)

plot_df["Label"] = (
    plot_df["Family"]
    .map(lambda x: mapping.get(x, x))
    + " ("
    + plot_df["Target_Clade"]
    + ")"
)

plot_df["Score"] = (
    -np.log10(plot_df["FDR"])
)

# ======================================
# PLOT
# ======================================

fig, ax = plt.subplots(
    figsize=(8, 5)
)

ax.barh(
    plot_df["Label"],
    plot_df["Score"]
)

ax.set_xlabel(
    "-log10(FDR)"
)

ax.set_ylabel(
    ""
)

ax.set_title(
    "Enriched prophage families by clade"
)

plt.tight_layout()

plt.savefig(
    OUTFILE,
    dpi=300,
    bbox_inches="tight"
)

plt.close()

print(
    f"Saved: {OUTFILE}"
)

print("\nFinished.")