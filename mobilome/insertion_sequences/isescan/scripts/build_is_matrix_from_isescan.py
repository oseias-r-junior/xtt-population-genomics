#!/usr/bin/env python3
"""
mobilome/insertion_sequences/isescan/scripts/build_is_matrix_from_isescan.py

Build is_family_matrix.tsv from ISEScan's OWN native per-instance predictions
(isescan/results/{genome}/fna/{genome}.fna.tsv), replacing the deprecated
build_is_matrix.py which counted "transposase" keyword hits in Prokka-annotated
GFF files (annotations/gff/) -- a crude text-search proxy with no concept of
element boundaries, terminal inverted repeats, or ISEScan's cluster/score
model, and which completely missed the IS21 family (96 real instances)
because it wasn't in that script's hardcoded 7-family list.

Real family diversity confirmed across all genomes with ISEScan results:
IS110, IS5, IS3, IS1595, IS4, IS21, IS481 (main families, >10 instances each)
plus 4 singleton instances (new, ISL3, IS30, IS200/IS605) grouped into "Other".

Iterates over the OFFICIAL genome list (metadata/genome_clades_master.tsv via
workflow/config/genomes.py), NOT the isescan/results/ directory listing --
that listing contains one bogus entry ("fna/", debris from a malformed past
invocation, not a real genome) and, more importantly, ISEScan does not write
ANY output files (.tsv/.gff/.csv/.sum) for a genome where zero IS elements
were found, so a directory-listing-only approach would silently DROP those
genomes from the matrix instead of correctly recording them as Total_IS=0.
Dropping them would bias any downstream per-genome statistic (means, PCA,
clade comparisons) that assumes complete genome coverage.

Output schema is UNCHANGED from the deprecated script (Genome, <families...>,
Total_IS) so all downstream scripts (defense_vs_is_*, PCA, heatmaps, Mantel
test, iTol tracks) work without modification against the corrected data.
"""

from pathlib import Path
from collections import defaultdict
import sys

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parents[3]

sys.path.append(str(PROJECT_ROOT / "scripts"))
sys.path.append(str(PROJECT_ROOT / "workflow" / "config"))

from config import *
from genomes import GENOMES

import pandas as pd

# ======================================================
# INPUT
# ======================================================

RESULTS_DIR = ISESCAN / "results"

MAIN_FAMILIES = [
    "IS110",
    "IS5",
    "IS3",
    "IS1595",
    "IS4",
    "IS21",
    "IS481",
]

# ======================================================
# OUTPUT
# ======================================================

IS_MATRICES.mkdir(parents=True, exist_ok=True)

OUT_FILE = IS_MATRICES / "is_family_matrix.tsv"

# ======================================================
# PARSE
# ======================================================

results = []
rare_family_log = defaultdict(int)
zero_is_genomes = []
missing_results_dir = []

for genome in sorted(GENOMES):
    genome_dir = RESULTS_DIR / genome

    row = {"Genome": genome}
    for fam in MAIN_FAMILIES:
        row[fam] = 0
    row["Other"] = 0

    if not genome_dir.is_dir():
        missing_results_dir.append(genome)
        row["Total_IS"] = 0
        results.append(row)
        continue

    tsv_file = genome_dir / "fna" / f"{genome}.fna.tsv"

    if not tsv_file.exists():
        # ISEScan completed but wrote no summary files -- genuinely zero IS
        # elements found (confirmed manually for a sample: log shows normal
        # completion, only intermediate ORF/HMM files present, no .tsv/.gff).
        zero_is_genomes.append(genome)
        row["Total_IS"] = 0
        results.append(row)
        continue

    df = pd.read_csv(tsv_file, sep="\t")

    counts = defaultdict(int)
    for family in df["family"]:
        if family in MAIN_FAMILIES:
            counts[family] += 1
        else:
            counts["Other"] += 1
            rare_family_log[family] += 1

    for fam in MAIN_FAMILIES:
        row[fam] = counts.get(fam, 0)
    row["Other"] = counts.get("Other", 0)
    row["Total_IS"] = sum(row[fam] for fam in MAIN_FAMILIES) + row["Other"]

    results.append(row)

df_out = pd.DataFrame(results)
df_out = df_out.sort_values("Total_IS", ascending=False)
df_out.to_csv(OUT_FILE, sep="\t", index=False)

print(f"Genomes in official list: {len(GENOMES)}")
print(f"Genomes in matrix: {len(df_out)}")
print()
print("Family totals:")
for fam in MAIN_FAMILIES:
    print(f"  {fam}: {df_out[fam].sum()}")
print(f"  Other: {df_out['Other'].sum()}")
if rare_family_log:
    print()
    print("Rare families grouped into 'Other':")
    for fam, n in sorted(rare_family_log.items(), key=lambda x: -x[1]):
        print(f"  {fam}: {n}")
if zero_is_genomes:
    print()
    print(f"Genomes with 0 IS elements (ISEScan completed, no elements found): {len(zero_is_genomes)}")
    print(f"  {', '.join(zero_is_genomes)}")
if missing_results_dir:
    print()
    print(f"WARNING -- genomes with NO ISEScan results directory at all (never ran / failed?): {len(missing_results_dir)}")
    print(f"  {', '.join(missing_results_dir)}")
print()
print(f"Saved: {OUT_FILE}")
