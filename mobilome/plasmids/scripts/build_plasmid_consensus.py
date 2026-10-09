#!/usr/bin/env python3
"""
mobilome/plasmids/scripts/build_plasmid_consensus.py

2-way plasmid consensus: MOB-suite x geNomad, joined by exact contig ID
(unlike prophages, plasmids are whole assembled contigs, not sub-regions
within a larger contig -- both tools report calls against the SAME contig
IDs from the same annotations/fna assembly, so an exact-ID join is the
correct consensus operation here, not coordinate overlap).

Each plasmid contig called by either tool becomes one consensus record,
annotated with:
    - how many / which tools support it (1-2)
    - mobility classification (Conjugative / Mobilizable / Non-mobilizable /
      Unknown), derived from MOB-suite's relaxase_type + mpf_type
    - replicon type (Rep_type, MOB-suite) and conjugation/AMR gene content
      (geNomad) as secondary typing evidence
"""

from pathlib import Path
import sys

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parents[2]

sys.path.append(str(PROJECT_ROOT / "scripts"))

from config import *

import pandas as pd

# ======================================================
# INPUT
# ======================================================

MOBSUITE_FILE = MOB_SUITE_MATRICES / "mobsuite_plasmids.tsv"
GENOMAD_FILE = GENOMAD_PLASMIDS_MATRICES / "genomad_plasmids.tsv"

# ======================================================
# OUTPUT
# ======================================================

PLASMID_MATRICES.mkdir(parents=True, exist_ok=True)
PLASMID_STATS.mkdir(parents=True, exist_ok=True)

OUT_CONSENSUS = PLASMID_MATRICES / "plasmid_consensus.tsv"
OUT_GENOME_SUMMARY = PLASMID_STATS / "plasmid_consensus_genome_summary.tsv"

# ======================================================
# LOAD
# ======================================================

mobsuite = pd.read_csv(MOBSUITE_FILE, sep="\t")
genomad = pd.read_csv(GENOMAD_FILE, sep="\t")

mobsuite_key = set(zip(mobsuite["Genome"], mobsuite["Contig_id"]))
genomad_key = set(zip(genomad["Genome"], genomad["Contig_id"]))

all_keys = sorted(mobsuite_key | genomad_key)

mobsuite_idx = mobsuite.set_index(["Genome", "Contig_id"])
genomad_idx = genomad.set_index(["Genome", "Contig_id"])

records = []

for genome, contig_id in all_keys:
    in_mobsuite = (genome, contig_id) in mobsuite_key
    in_genomad = (genome, contig_id) in genomad_key

    tools = []
    if in_mobsuite:
        tools.append("MOB-suite")
    if in_genomad:
        tools.append("geNomad")

    m = mobsuite_idx.loc[(genome, contig_id)] if in_mobsuite else None
    g = genomad_idx.loc[(genome, contig_id)] if in_genomad else None

    size_bp = None
    mobility = "Unknown"
    rep_type = ""
    relaxase_type = ""
    mpf_type = ""
    if m is not None:
        size_bp = m["Size_bp"]
        mobility = m["Predicted_mobility_derived"]
        rep_type = m["Rep_type"]
        relaxase_type = m["Relaxase_type"]
        mpf_type = m["MPF_type"]

    plasmid_score = None
    conjugation_genes = ""
    amr_genes = ""
    if g is not None:
        if size_bp is None:
            size_bp = g["Length_bp"]
        plasmid_score = g["Plasmid_score"]
        conjugation_genes = g["Conjugation_genes"]
        amr_genes = g["AMR_genes"]

    records.append({
        "Genome": genome,
        "Contig_id": contig_id,
        "Size_bp": size_bp,
        "N_Tools_Support": len(tools),
        "Tools_Supporting": ",".join(tools),
        "Mobility_Class": mobility,
        "Rep_type": rep_type,
        "Relaxase_type": relaxase_type,
        "MPF_type": mpf_type,
        "GeNomad_Plasmid_Score": plasmid_score,
        "Conjugation_genes": conjugation_genes,
        "AMR_genes": amr_genes,
    })

result = pd.DataFrame(records)
result.to_csv(OUT_CONSENSUS, sep="\t", index=False)

genome_summary = (
    result
    .groupby("Genome")
    .agg(
        Total_plasmid_contigs=("Contig_id", "count"),
        N_supported_by_2_tools=("N_Tools_Support", lambda s: (s == 2).sum()),
        N_supported_by_1_tool=("N_Tools_Support", lambda s: (s == 1).sum()),
    )
    .reset_index()
)
genome_summary.to_csv(OUT_GENOME_SUMMARY, sep="\t", index=False)

print(f"Genomes with >=1 plasmid contig: {genome_summary['Genome'].nunique()}")
print(f"Total consensus plasmid contigs: {len(result)}")
print()
print("Support level:")
print(result["N_Tools_Support"].value_counts().sort_index().to_string())
print()
print("Mobility classification:")
print(result["Mobility_Class"].value_counts().to_string())
print()
print(f"Saved: {OUT_CONSENSUS}")
print(f"Saved: {OUT_GENOME_SUMMARY}")
