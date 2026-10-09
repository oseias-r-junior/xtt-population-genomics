#!/usr/bin/env python3
"""
Build IslandPath summary matrices.

This script parses all IslandPath GFF3 files generated for the Xtt population
and builds standardized matrices used throughout the Genomic Islands module.

Outputs
-------
islandpath_summary_matrix.tsv
islandpath_coordinates.tsv
islandpath_presence_absence.tsv

Author
------
Xtt Population Genomics Pipeline
"""

from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd

from workflow.config.paths import (
    GENOME_CLADES,
    ISLANDPATH_RAW,
    ISLANDPATH_MATRICES,
)

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

def parse_attributes(attribute_string: str) -> dict:

    attributes = {}

    for item in attribute_string.split(";"):

        if "=" not in item:
            continue

        key, value = item.split("=", 1)

        attributes[key] = value

    return attributes


# ============================================================
# Read all IslandPath GFF files
# ============================================================

def load_predictions(official_genomes: set) -> pd.DataFrame:

    rows = []

    gff_files = sorted(
        f for f in ISLANDPATH_RAW.glob("*/islandpath.gff")
        if f.parent.name in official_genomes
        # Skip anything not in the official genome list (e.g. a stray
        # "gbk" directory left over from a malformed past invocation --
        # not a real genome).
    )

    LOGGER.info(
        "Found %d IslandPath prediction files.",
        len(gff_files),
    )

    for gff in gff_files:

        genome = gff.parent.name

        with open(gff) as handle:

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
                    dict(
                        Genome=genome,
                        Region_ID=attr.get("ID", ""),
                        Start=start,
                        End=end,
                        Length=end - start + 1,
                        Strand=strand,
                        Detector="IslandPath",
                    )
                )

    return pd.DataFrame(rows)


# ============================================================
# Main
# ============================================================

def main():

    ISLANDPATH_MATRICES.mkdir(
        parents=True,
        exist_ok=True,
    )

    metadata = pd.read_csv(
        GENOME_CLADES,
        sep="\t",
    )

    official_genomes = set(metadata["Genome"])

    predictions = load_predictions(official_genomes)

    # --------------------------------------------------------
    # Summary matrix
    # --------------------------------------------------------

    summary = (
        predictions
        .groupby("Genome")
        .agg(
            Number_of_Islands=("Region_ID", "count"),
            Total_GI_Length=("Length", "sum"),
            Mean_GI_Length=("Length", "mean"),
            Median_GI_Length=("Length", "median"),
            Longest_GI=("Length", "max"),
        )
        .reset_index()
    )

    summary = metadata.merge(
        summary,
        on="Genome",
        how="left",
    )

    summary.fillna(
        {
            "Number_of_Islands": 0,
            "Total_GI_Length": 0,
        },
        inplace=True,
    )

    summary.to_csv(
        ISLANDPATH_MATRICES /
        "islandpath_summary_matrix.tsv",
        sep="\t",
        index=False,
    )

    LOGGER.info(
        "Summary matrix written."
    )

    # --------------------------------------------------------
    # Coordinates
    # --------------------------------------------------------

    predictions.to_csv(
        ISLANDPATH_MATRICES /
        "islandpath_coordinates.tsv",
        sep="\t",
        index=False,
    )

    LOGGER.info(
        "Coordinate matrix written."
    )

    # --------------------------------------------------------
    # Presence/absence
    # --------------------------------------------------------

    presence = (
        predictions.assign(Present=1)
        .pivot_table(
            index="Genome",
            columns="Region_ID",
            values="Present",
            fill_value=0,
        )
        .reset_index()
    )

    presence.to_csv(
        ISLANDPATH_MATRICES /
        "islandpath_presence_absence.tsv",
        sep="\t",
        index=False,
    )

    LOGGER.info(
        "Presence/absence matrix written."
    )

    LOGGER.info(
        "Finished."
    )


if __name__ == "__main__":

    main()