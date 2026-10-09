import pandas as pd

from pathlib import Path
import sys

SCRIPT_DIR = Path(__file__).resolve().parent

PROJECT_ROOT = SCRIPT_DIR.parents[2]

sys.path.append(
    str(PROJECT_ROOT / "scripts")
)

from config import *
from normalize_genome_names import normalize_genome_name


# ======================================
# INPUT
# ======================================

IS_FILE = (
    IS_MATRICES
    / "is_family_matrix.tsv"
)

META_FILE = GENOME_CLADES

OUTDIR = (
    IS_ITOL
)

OUTDIR.mkdir(
    parents=True,
    exist_ok=True
)

OUTDIR.mkdir(
    parents=True,
    exist_ok=True
)

# ======================================
# LOAD
# ======================================

is_df = pd.read_csv(IS_FILE, sep="\t")

is_df["Genome"] = (
    is_df["Genome"]
    .astype(str)
    .apply(normalize_genome_name)
)

meta_df = pd.read_csv(META_FILE, sep="\t")

meta_df["Genome"] = (
    meta_df["Genome"]
    .astype(str)
    .apply(normalize_genome_name)
)

df = pd.merge(
    is_df,
    meta_df,
    on="Genome",
    how="inner"
)

# ======================================
# COLORS
# ======================================

clade_colors = {
    "K0": "#1f77b4",
    "K1": "#ff7f0e",
    "K2": "#2ca02c"
}

# ======================================
# 1. CLADES
# ======================================

clade_file = OUTDIR / "itol_clades.txt"

with open(clade_file, "w") as f:

    f.write("DATASET_COLORSTRIP\n")
    f.write("SEPARATOR TAB\n")
    f.write("DATASET_LABEL\tANI_clades\n")
    f.write("COLOR\t#000000\n")
    f.write("DATA\n")

    for _, row in df.iterrows():

        genome = row["Genome"]
        clade = row["Clade"]
        color = clade_colors[clade]

        f.write(
            f"{genome}\t{color}\t{clade}\n"
        )

print(f"Saved: {clade_file}")

# ======================================
# 2. TOTAL IS BARS
# ======================================

bar_file = OUTDIR / "itol_total_is.txt"

with open(bar_file, "w") as f:

    f.write("DATASET_SIMPLEBAR\n")
    f.write("SEPARATOR TAB\n")
    f.write("DATASET_LABEL\tTotal_IS\n")
    f.write("COLOR\t#444444\n")
    f.write("DATA\n")

    for _, row in df.iterrows():

        genome = row["Genome"]
        total_is = row["Total_IS"]

        f.write(
            f"{genome}\t{total_is}\n"
        )

print(f"Saved: {bar_file}")

# ======================================
# 3. IS FAMILY HEATMAP
# ======================================

heatmap_file = OUTDIR / "itol_is_heatmap.txt"

is_families = [
    c for c in df.columns
    if c.startswith("IS")
]

with open(heatmap_file, "w") as f:

    f.write("DATASET_HEATMAP\n")
    f.write("SEPARATOR TAB\n")
    f.write("DATASET_LABEL\tIS_families\n")
    f.write("COLOR\t#ff0000\n")

    f.write(
        "FIELD_LABELS\t"
        + "\t".join(is_families)
        + "\n"
    )

    f.write("DATA\n")

    for _, row in df.iterrows():

        genome = row["Genome"]

        values = [
            str(row[x])
            for x in is_families
        ]

        f.write(
            genome
            + "\t"
            + "\t".join(values)
            + "\n"
        )

print(f"Saved: {heatmap_file}")