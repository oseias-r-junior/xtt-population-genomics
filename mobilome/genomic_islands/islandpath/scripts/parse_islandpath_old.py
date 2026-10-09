#!/usr/bin/env python3
"""
Parse IslandPath-DIMOB GFF3 predictions.

Converts the raw GFF3 output produced by IslandPath into a standardized
tab-separated table used throughout the Xtt Population Genomics pipeline.

Input
-----
IslandPath GFF3

Output
------
TSV with one genomic island per row.

Author
------
Xtt Population Genomics Pipeline
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import pandas as pd


# ============================================================
# Logging
# ============================================================

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
)

LOGGER = logging.getLogger(__name__)


# ============================================================
# Helpers
# ============================================================

def parse_attributes(attributes: str) -> dict:
    """
    Parse GFF3 attribute column.
    """

    result = {}

    for item in attributes.split(";"):

        if "=" in item:

            key, value = item.split("=", 1)

            result[key] = value

    return result


# ============================================================
# Parser
# ============================================================

def parse_gff(gff_file: Path) -> pd.DataFrame:

    genome = gff_file.parent.name

    rows = []

    with open(gff_file) as handle:

        for line in handle:

            if line.startswith("#"):
                continue

            fields = line.rstrip().split("\t")

            if len(fields) != 9:
                continue

            (
                seqid,
                source,
                feature,
                start,
                end,
                score,
                strand,
                phase,
                attributes,
            ) = fields

            attr = parse_attributes(attributes)

            start = int(start)
            end = int(end)

            rows.append(
                {
                    "Genome": genome,
                    "Island_ID": attr.get("ID", ""),
                    "SeqID": seqid,
                    "Source": source,
                    "Feature": feature,
                    "Start": start,
                    "End": end,
                    "Length": end - start + 1,
                    "Score": score,
                    "Strand": strand,
                    "Phase": phase,
                    "Detector": "IslandPath",
                }
            )

    return pd.DataFrame(rows)


# ============================================================
# Main
# ============================================================

def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--input",
        required=True,
        type=Path,
        help="IslandPath GFF3 file.",
    )

    parser.add_argument(
        "--output",
        required=True,
        type=Path,
        help="Output TSV.",
    )

    args = parser.parse_args()

    LOGGER.info("Reading %s", args.input.name)

    df = parse_gff(args.input)

    args.output.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    df.to_csv(
        args.output,
        sep="\t",
        index=False,
    )

    LOGGER.info(
        "Parsed %d genomic islands.",
        len(df),
    )

    LOGGER.info("Finished.")


if __name__ == "__main__":

    main()