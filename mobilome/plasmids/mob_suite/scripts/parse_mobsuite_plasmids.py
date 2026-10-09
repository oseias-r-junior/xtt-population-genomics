#!/usr/bin/env python3
"""
mobilome/plasmids/mob_suite/scripts/parse_mobsuite_plasmids.py

Parse MOB-suite's contig_report.txt (one per genome, under raw_outputs/)
into a single plasmid-contig-level matrix, keeping only rows where
molecule_type == "plasmid" (dropping the chromosome row).

MOB-suite's own predicted_mobility column is sometimes empty even when
relaxase_type and mpf_type are both present (observed on real data), so a
mobility class is ALSO derived independently here, using the standard
MOB-suite logic:
    relaxase present + MPF present  -> Conjugative
    relaxase present, no MPF        -> Mobilizable
    no relaxase                     -> Non-mobilizable
"""

from pathlib import Path
import sys

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parents[3]

sys.path.append(str(PROJECT_ROOT / "scripts"))
sys.path.append(str(PROJECT_ROOT / "scripts" / "metadata"))

from config import *
from normalize_genome_names import normalize_genome_name

import pandas as pd

# ======================================================
# INPUT
# ======================================================

RAW_DIR = MOB_SUITE_RAW

# ======================================================
# OUTPUT
# ======================================================

MOB_SUITE_MATRICES.mkdir(parents=True, exist_ok=True)

OUT_FILE = MOB_SUITE_MATRICES / "mobsuite_plasmids.tsv"


def derive_mobility(relaxase, mpf):
    has_relaxase = isinstance(relaxase, str) and relaxase.strip() not in ("", "-")
    has_mpf = isinstance(mpf, str) and mpf.strip() not in ("", "-")

    if has_relaxase and has_mpf:
        return "Conjugative"
    if has_relaxase:
        return "Mobilizable"
    return "Non-mobilizable"


records = []

for genome_dir in sorted(RAW_DIR.iterdir()):
    if not genome_dir.is_dir():
        continue

    report_file = genome_dir / "contig_report.txt"
    if not report_file.exists():
        continue

    df = pd.read_csv(report_file, sep="\t")
    if df.empty:
        continue

    plasmid_rows = df[df["molecule_type"] == "plasmid"]

    for _, row in plasmid_rows.iterrows():
        genome = normalize_genome_name(genome_dir.name)

        records.append({
            "Genome": genome,
            "Contig_id": row.get("contig_id", ""),
            "Size_bp": row.get("size", None),
            "GC": row.get("gc", None),
            "Rep_type": row.get("rep_type(s)", ""),
            "Relaxase_type": row.get("relaxase_type(s)", ""),
            "MPF_type": row.get("mpf_type", ""),
            "Predicted_mobility_MOBsuite": row.get("predicted_mobility", ""),
            "Predicted_mobility_derived": derive_mobility(
                row.get("relaxase_type(s)", ""),
                row.get("mpf_type", "")
            ),
            "Primary_cluster_id": row.get("primary_cluster_id", ""),
            "Mash_nearest_neighbor": row.get("mash_nearest_neighbor", ""),
            "Mash_identification": row.get("mash_neighbor_identification", ""),
        })

result = pd.DataFrame(records)
result.to_csv(OUT_FILE, sep="\t", index=False)

print(f"Genomes scanned: {sum(1 for d in RAW_DIR.iterdir() if d.is_dir())}")
print(f"Plasmid contigs extracted: {len(result)}")
if not result.empty:
    print(f"Genomes with >=1 plasmid: {result['Genome'].nunique()}")
    print(result["Predicted_mobility_derived"].value_counts().to_string())
print(f"Saved: {OUT_FILE}")
