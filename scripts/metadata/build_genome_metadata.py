import pandas as pd
from pathlib import Path

# ======================================
# INPUT
# ======================================

EXCEL = Path(
    "~/diel_paper_colab/input_data_mge_ds/"
    "Supplemental_Table_1.xlsx"
).expanduser()

OUTFILE = Path(
    "~/diel_paper_colab/results/genome_clades_master.tsv"
).expanduser()

# ======================================
# LOAD EXCEL
# ======================================

df = pd.read_excel(EXCEL, header=1)

# ======================================
# CLEAN
# ======================================

df["Genome"] = df["Name"].astype(str) + "_out"

# normalize genome size
df["Genome_Length"] = (
    df["Genome Length"]
    .astype(str)
    .str.replace(".", "", regex=False)
    .astype(int)
)

# keep relevant columns
meta = df[[
    "Genome",
    "Genome_Length",
    "ANI-clade"
]].copy()

meta.columns = [
    "Genome",
    "Genome_Length",
    "Clade"
]

# ======================================
# SAVE
# ======================================

meta.to_csv(
    OUTFILE,
    sep="\t",
    index=False
)

print(meta.head())

print(f"\nSaved: {OUTFILE}")