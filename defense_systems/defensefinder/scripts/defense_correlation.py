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
import seaborn as sns
import matplotlib.pyplot as plt

# =========================
# INPUT
# =========================

MATRIX_FILE = (
    DEFENSEFINDER_MATRICES
    / "defensefinder_matrix.tsv"
)

# =========================
# OUTPUT
# =========================

CORR_DIR = (
    DEFENSEFINDER_STATS
    / "correlation"
)

CORR_DIR.mkdir(
    parents=True,
    exist_ok=True
)

CORR_OUTFILE = (
    CORR_DIR
    / "defense_correlation.tsv"
)

PLOT_OUTFILE = (
    DEFENSEFINDER_PLOTS
    / "defense_correlation_heatmap.png"
)

# =========================
# LOAD DATA
# =========================

df = pd.read_csv(
    MATRIX_FILE,
    sep="\t"
)

X = df.drop(
    columns=["Genome"]
)

# =========================
# REMOVE INVARIANT SYSTEMS
# =========================

variable_cols = [
    c
    for c in X.columns
    if X[c].nunique() > 1
]

X = X[variable_cols]

print(
    f"After invariant filtering: {X.shape[1]} systems"
)

# =========================
# REMOVE VERY RARE SYSTEMS
# =========================

MIN_PREVALENCE = 5

keep_cols = [
    c
    for c in X.columns
    if X[c].sum() >= MIN_PREVALENCE
]

removed_cols = sorted(
    set(X.columns) - set(keep_cols)
)

print(
    f"Removed rare systems (< {MIN_PREVALENCE} genomes): {len(removed_cols)}"
)

for c in removed_cols:
    print(c)

X = X[keep_cols]

print(
    f"Systems retained: {X.shape[1]}"
)

# =========================
# CORRELATION MATRIX
# =========================

corr = X.corr(
    method="pearson"
)

corr = corr.fillna(0)

print(
    f"Variable systems retained: {X.shape[1]}"
)

corr.to_csv(
    CORR_OUTFILE,
    sep="\t"
)

print(
    f"Saved: {CORR_OUTFILE}"
)

# =========================
# CLUSTERED HEATMAP
# =========================

g = sns.clustermap(
    corr,
    cmap="vlag",
    center=0,
    figsize=(14, 14),
    linewidths=0,
    vmin=-1,
    vmax=1
)

g.fig.suptitle(
    "Defense system co-occurrence correlations",
    y=1.02
)

plt.savefig(
    PLOT_OUTFILE,
    dpi=300,
    bbox_inches="tight"
)

plt.close()

print(
    f"Saved: {PLOT_OUTFILE}"
)