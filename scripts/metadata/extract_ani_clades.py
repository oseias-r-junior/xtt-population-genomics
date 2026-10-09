#!/usr/bin/env python3

import pandas as pd
from xtt_population_genomics.scripts.metadata.normalize_genome_names import normalize_genome_name

# =========================
# INPUT
# =========================

EXCEL = "Supplemental_Table_1.xlsx"

# =========================
# LOAD EXCEL
# =========================

df = pd.read_excel(EXCEL, header=1)

print(df.head())

# =========================
# FORMAT GENOME NAMES
# =========================

df["Genome"] = (
    df["Name"]
    .astype(str)
    .apply(normalize_genome_name)
)

# =========================
# SELECT COLUMNS
# =========================

out = df[["Genome", "ANI-clade"]]

out.columns = ["Genome", "Clade"]

# =========================
# SAVE
# =========================

out.to_csv(
    "results/genome_clades.tsv",
    sep="\t",
    index=False
)

print("\nSaved:")
print("genome_clades.tsv")
