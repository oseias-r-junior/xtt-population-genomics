#!/usr/bin/env python3

"""
Summarize IslandPath predictions.

Produces

- summary statistics
- descriptive plots

Author
------
Xtt Population Genomics Pipeline
"""

from pathlib import Path
import logging

import matplotlib.pyplot as plt
import pandas as pd

from workflow.config.paths import (
    ISLANDPATH_MATRICES,
    ISLANDPATH_STATS,
    ISLANDPATH_PLOTS,
)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
)

LOGGER = logging.getLogger(__name__)


def main():

    ISLANDPATH_STATS.mkdir(
        parents=True,
        exist_ok=True,
    )

    ISLANDPATH_PLOTS.mkdir(
        parents=True,
        exist_ok=True,
    )

    summary = pd.read_csv(
        ISLANDPATH_MATRICES /
        "islandpath_summary_matrix.tsv",
        sep="\t",
    )

    coordinates = pd.read_csv(
        ISLANDPATH_MATRICES /
        "islandpath_coordinates.tsv",
        sep="\t",
    )

    stats = pd.DataFrame(
        [
            {
                "Genomes": len(summary),
                "Genomes_with_GIs":
                    (summary["Number_of_Islands"] > 0).sum(),
                "Total_GIs":
                    len(coordinates),
                "Mean_GIs_per_Genome":
                    summary["Number_of_Islands"].mean(),
                "Median_GIs_per_Genome":
                    summary["Number_of_Islands"].median(),
                "Mean_GI_Length":
                    coordinates["Length"].mean(),
                "Median_GI_Length":
                    coordinates["Length"].median(),
                "Longest_GI":
                    coordinates["Length"].max(),
                "Total_GI_bp":
                    coordinates["Length"].sum(),
            }
        ]
    )

    stats.to_csv(
        ISLANDPATH_STATS /
        "islandpath_summary.tsv",
        sep="\t",
        index=False,
    )

    LOGGER.info("Summary table written.")

    # --------------------------------------------------
    # Number of GIs per genome
    # --------------------------------------------------

    plt.figure(figsize=(6,4))

    plt.hist(
        summary["Number_of_Islands"],
        bins=20,
    )

    plt.xlabel("Genomic islands")

    plt.ylabel("Genomes")

    plt.tight_layout()

    plt.savefig(
        ISLANDPATH_PLOTS /
        "islandpath_prevalence.png",
        dpi=300,
    )

    plt.close()

    # --------------------------------------------------
    # GI lengths
    # --------------------------------------------------

    plt.figure(figsize=(6,4))

    plt.hist(
        coordinates["Length"],
        bins=30,
    )

    plt.xlabel("Island length (bp)")

    plt.ylabel("Count")

    plt.tight_layout()

    plt.savefig(
        ISLANDPATH_PLOTS /
        "islandpath_length_distribution.png",
        dpi=300,
    )

    plt.close()

    # --------------------------------------------------
    # Total bp per genome
    # --------------------------------------------------

    plt.figure(figsize=(6,4))

    plt.hist(
        summary["Total_GI_Length"],
        bins=20,
    )

    plt.xlabel("Total genomic island length")

    plt.ylabel("Genomes")

    plt.tight_layout()

    plt.savefig(
        ISLANDPATH_PLOTS /
        "islandpath_total_length_distribution.png",
        dpi=300,
    )

    plt.close()

    # --------------------------------------------------
    # Longest island
    # --------------------------------------------------

    plt.figure(figsize=(6,4))

    plt.hist(
        summary["Longest_GI"],
        bins=20,
    )

    plt.xlabel("Longest genomic island")

    plt.ylabel("Genomes")

    plt.tight_layout()

    plt.savefig(
        ISLANDPATH_PLOTS /
        "islandpath_longest_distribution.png",
        dpi=300,
    )

    plt.close()

    # --------------------------------------------------
    # Density
    # --------------------------------------------------

    summary["GI_density"] = (
        summary["Number_of_Islands"]
        /
        (summary["Genome_Length"] / 1e6)
    )

    plt.figure(figsize=(6,4))

    plt.hist(
        summary["GI_density"],
        bins=20,
    )

    plt.xlabel("Genomic islands / Mb")

    plt.ylabel("Genomes")

    plt.tight_layout()

    plt.savefig(
        ISLANDPATH_PLOTS /
        "islandpath_density.png",
        dpi=300,
    )

    plt.close()

    LOGGER.info("Finished.")


if __name__ == "__main__":

    main()