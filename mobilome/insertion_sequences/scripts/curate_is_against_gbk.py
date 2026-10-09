#!/usr/bin/env python3
"""
mobilome/insertion_sequences/scripts/curate_is_against_gbk.py

Manual-curation-style QC for ISEScan IS calls, mirroring the curation
Guillermo performed on ISEScan-derived IS calls in the Xylella paper: for
every IS instance ISEScan detected, check whether there is an
independently annotated IS/transposase-related gene in the GBK at the
same coordinates (same contig, overlapping interval).

ISEScan is a purpose-built, more sensitive/complete detector than a
keyword search over Prokka/NCBI annotations, so most instances are
expected to be confirmed. Instances with NO corresponding annotated gene
are flagged as GBK_Unconfirmed and excluded from the SSN network (but
kept on record here, not silently dropped).

Match criterion for a GBK gene "counting" as IS-related: its /product
qualifier matches (case-insensitive) any of: "transposase",
"insertion sequence", or "ISxxx family" (xxx = digits).

Coordinates: ISEScan's isBegin/isEnd are on the same coordinate system as
the GBK contig they were called on (seqID == GBK record id/name for that
contig); overlap is any interval intersection > 0.
"""

from pathlib import Path
import re
import sys

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parents[2]

sys.path.append(str(PROJECT_ROOT / "scripts"))
sys.path.append(str(PROJECT_ROOT / "workflow" / "config"))

from config import *
from genomes import GENOMES

from Bio import SeqIO
import pandas as pd

# ======================================================
# CONFIG
# ======================================================

ISESCAN_RESULTS = ISESCAN / "results"
GBK_DIR = PROJECT_ROOT / "annotations" / "gbk"

IS_PRODUCT_RE = re.compile(
    r"transposase|insertion sequence|IS\d+\s+family",
    re.IGNORECASE,
)

OUT_LOG = IS_MATRICES / "is_gbk_curation_log.tsv"
OUT_CONFIRMED = IS_MATRICES / "is_gbk_confirmed_instances.tsv"
OUT_SUMMARY = IS_STATS / "is_gbk_curation_summary.tsv"

IS_STATS.mkdir(parents=True, exist_ok=True)
IS_MATRICES.mkdir(parents=True, exist_ok=True)


def load_gbk_is_genes(gbk_path: Path) -> dict:
    """Return {contig_id: [(start, end), ...]} for IS/transposase-related
    CDS features, across every contig/record in the (possibly
    multi-record) GBK file."""
    genes_by_contig = {}

    if not gbk_path.exists():
        return genes_by_contig

    for record in SeqIO.parse(str(gbk_path), "genbank"):
        contig_id = record.name or record.id
        intervals = []

        for feature in record.features:
            if feature.type != "CDS":
                continue

            product = feature.qualifiers.get("product", [""])[0]

            if IS_PRODUCT_RE.search(product):
                start = int(feature.location.start) + 1  # GBK 1-based
                end = int(feature.location.end)
                intervals.append((start, end))

        genes_by_contig[contig_id] = intervals

    return genes_by_contig


def overlaps(a_start, a_end, intervals) -> bool:
    for b_start, b_end in intervals:
        if min(a_end, b_end) - max(a_start, b_start) > 0:
            return True
    return False


# ======================================================
# MAIN
# ======================================================

log_rows = []

for genome in sorted(GENOMES):
    tsv_file = ISESCAN_RESULTS / genome / "fna" / f"{genome}.fna.tsv"

    if not tsv_file.exists():
        continue  # genome with 0 IS instances, nothing to curate

    is_df = pd.read_csv(tsv_file, sep="\t")

    if is_df.empty:
        continue

    gbk_path = GBK_DIR / f"{genome}.gbk"
    genes_by_contig = load_gbk_is_genes(gbk_path)

    for _, row in is_df.iterrows():
        seqid = str(row["seqID"])
        is_begin = int(row["isBegin"])
        is_end = int(row["isEnd"])
        family = row["family"]

        contig_genes = genes_by_contig.get(seqid, [])
        confirmed = overlaps(is_begin, is_end, contig_genes)

        log_rows.append({
            "Genome": genome,
            "Contig": seqid,
            "Family": family,
            "isBegin": is_begin,
            "isEnd": is_end,
            "isLen": int(row["isLen"]),
            "GBK_Confirmed": confirmed,
            "GBK_file_found": gbk_path.exists(),
        })

log_df = pd.DataFrame(log_rows)
log_df.to_csv(OUT_LOG, sep="\t", index=False)

confirmed_df = log_df[log_df["GBK_Confirmed"]].copy()
confirmed_df.to_csv(OUT_CONFIRMED, sep="\t", index=False)

# ======================================================
# SUMMARY
# ======================================================

total = len(log_df)
n_confirmed = int(log_df["GBK_Confirmed"].sum())
n_unconfirmed = total - n_confirmed

summary = (
    log_df
    .groupby("Family")
    .agg(
        Total_Instances=("GBK_Confirmed", "count"),
        Confirmed=("GBK_Confirmed", "sum"),
    )
    .assign(Unconfirmed=lambda d: d["Total_Instances"] - d["Confirmed"])
    .assign(Pct_Confirmed=lambda d: 100 * d["Confirmed"] / d["Total_Instances"])
    .reset_index()
    .sort_values("Total_Instances", ascending=False)
)
summary.to_csv(OUT_SUMMARY, sep="\t", index=False)

print(f"Total ISEScan instances scanned: {total}")
print(f"GBK-confirmed: {n_confirmed} ({100*n_confirmed/total:.1f}%)")
print(f"GBK-UNconfirmed (flagged, excluded from SSN): {n_unconfirmed} ({100*n_unconfirmed/total:.1f}%)")
print()
print("Missing GBK files (genome scanned but no annotation available):",
      int((~log_df["GBK_file_found"]).sum()))
print()
print("Breakdown by family:")
print(summary.to_string(index=False))
print()
print(f"Saved: {OUT_LOG}")
print(f"Saved: {OUT_CONFIRMED}")
print(f"Saved: {OUT_SUMMARY}")
