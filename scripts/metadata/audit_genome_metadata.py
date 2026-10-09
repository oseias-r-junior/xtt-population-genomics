#!/usr/bin/env python3
"""
audit_genome_metadata.py

Audit genome metadata used throughout the Xtt Population Genomics Pipeline.

This script compares genome lengths recovered from all available annotation
formats (GBK, FNA, FSA and GFF) and generates a comprehensive audit report.

Unlike update_genome_clades_master.py, this script NEVER modifies the official
metadata.

Outputs
-------
genome_length_audit.tsv

Author
------
Xtt Population Genomics Pipeline
"""

from __future__ import annotations

import logging
import re
from pathlib import Path

import pandas as pd
from Bio import SeqIO

from normalize_genome_names import normalize_genome_name

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))


from workflow.config.paths import (
    GBK,
    FNA,
    GFF,
    METADATA,
    XTT_PROJECT,
)

# ==========================================================
# Optional directories
# ==========================================================

FSA = XTT_PROJECT / "annotations" / "fsa"

MASTER = METADATA / "genome_clades_master.tsv"

AUDIT = METADATA / "genome_length_audit.tsv"

# ==========================================================
# Logging
# ==========================================================

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
)

LOGGER = logging.getLogger(__name__)

# ==========================================================
# GBK
# ==========================================================


def gbk_seqio_length(path: Path):

    if not path.exists():
        return None

    try:

        records = list(
            SeqIO.parse(path, "genbank")
        )

        if not records:
            return None

        return sum(len(r.seq) for r in records)

    except Exception:

        return None


def gbk_locus_length(path: Path):

    if not path.exists():
        return None

    with open(path) as f:

        for line in f:

            if line.startswith("LOCUS"):

                m = re.search(r"(\d+)\s+bp", line)

                if m:

                    return int(m.group(1))

    return None


def gbk_origin_length(path: Path):

    if not path.exists():
        return None

    origin = False

    total = 0

    with open(path) as f:

        for line in f:

            if line.startswith("ORIGIN"):

                origin = True

                continue

            if origin:

                if line.startswith("//"):

                    break

                seq = re.sub("[^ATGCNatgcn]", "", line)

                total += len(seq)

    if total == 0:

        return None

    return total


# ==========================================================
# FASTA
# ==========================================================


def fasta_length(path: Path):

    if not path.exists():
        return None

    length = 0

    with open(path) as f:

        for line in f:

            if line.startswith(">"):

                continue

            length += len(line.strip())

    return length if length else None


# ==========================================================
# GFF
# ==========================================================


def gff_length(path: Path):

    if not path.exists():
        return None

    total = 0

    found = False

    with open(path) as f:

        for line in f:

            if line.startswith("##sequence-region"):

                found = True

                parts = line.split()

                total += int(parts[3]) - int(parts[2]) + 1

    if found:

        return total

    return None


# ==========================================================
# Confidence
# ==========================================================


def choose_best(row):

    candidates = [
        ("GBK_SEQIO", row["GBK_SEQIO"], "HIGH"),
        ("GBK_LOCUS", row["GBK_LOCUS"], "HIGH"),
        ("GBK_ORIGIN", row["GBK_ORIGIN"], "HIGH"),
        ("FNA", row["FNA"], "HIGH"),
        ("FSA", row["FSA"], "MEDIUM"),
        ("GFF", row["GFF"], "LOW"),
    ]

    for source, value, confidence in candidates:

        if pd.notna(value):

            return source, int(value), confidence

    return "FAILED", pd.NA, "FAILED"





# ==========================================================
# Agreement
# ==========================================================


def agreement(values):

    vals = {
        int(v)
        for v in values
        if pd.notna(v)
    }

    return len(vals) <= 1


def qc_status(row):

    if pd.isna(row["Chosen_Length"]):
        return "FAILED"

    if not row["Primary_Agreement"]:
        return "PRIMARY_CONFLICT"

    # SeqIO indisponível, mas demais fontes concordam
    if (
        pd.isna(row["GBK_SEQIO"])
        and row["Primary_Agreement"]
    ):
        return "SEQIO_FAILED"

    # LOCUS representa apenas um contig
    if (
        pd.notna(row["GBK_SEQIO"])
        and pd.notna(row["GBK_LOCUS"])
        and row["GBK_SEQIO"] != row["GBK_LOCUS"]
    ):
        return "PASS_PARTIAL_LOCUS"

    return "PASS"


# ==========================================================
# Main
# ==========================================================


def main():

    metadata = pd.read_csv(
        MASTER,
        sep="\t",
    )

    rows = []

    LOGGER.info(
        "Auditing %d genomes.",
        len(metadata),
    )

    for genome in metadata["Genome"]:

        genome = normalize_genome_name(genome)

        gbk = GBK / f"{genome}.gbk"

        fna = FNA / f"{genome}.fna"

        fsa = FSA / f"{genome}.fsa"

        gff = GFF / f"{genome}.gff"

        row = dict(

            Genome=genome,

            GBK_SEQIO=gbk_seqio_length(gbk),

            GBK_LOCUS=gbk_locus_length(gbk),

            GBK_ORIGIN=gbk_origin_length(gbk),

            FNA=fasta_length(fna),

            FSA=fasta_length(fsa),

            GFF=gff_length(gff),

        )

        # ----------------------------------------------------------
        # Agreement among primary data sources
        # ----------------------------------------------------------

        row["Primary_Agreement"] = agreement(

            [

                row["GBK_SEQIO"],

                row["FNA"],

                row["FSA"],

                row["GFF"],

            ]

        )

        # ----------------------------------------------------------
        # Agreement among diagnostic fields
        # ----------------------------------------------------------

        row["Diagnostic_Agreement"] = agreement(

            [

                row["GBK_LOCUS"],

                row["GBK_ORIGIN"],

            ]

        )

        source, length, confidence = choose_best(row)

        row["Chosen_Length"] = length

        row["Chosen_Source"] = source

        row["Confidence"] = confidence

        row["QC_Status"] = qc_status(row)

        rows.append(row)

    audit = pd.DataFrame(rows)

    audit.to_csv(

        AUDIT,

        sep="\t",

        index=False,

    )

    LOGGER.info(

        "Audit written to %s",

        AUDIT,

    )

    LOGGER.info(

        "Primary agreement: %.1f%%",

        100 * audit["Primary_Agreement"].mean(),

    )

    LOGGER.info(

        "Diagnostic agreement: %.1f%%",

        100 * audit["Diagnostic_Agreement"].mean(),

    )

    LOGGER.info(

        audit["QC_Status"].value_counts().to_string()

    )

    LOGGER.info("")

    LOGGER.info("QC summary")

    for status, n in audit["QC_Status"].value_counts().items():

        LOGGER.info(
            "%-22s %3d",
            status,
            n,
        )

    LOGGER.info(

        "Finished.",

    )


if __name__ == "__main__":

    main()