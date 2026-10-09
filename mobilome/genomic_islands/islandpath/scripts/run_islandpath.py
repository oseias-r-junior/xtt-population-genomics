#!/usr/bin/env python3
"""
Run IslandPath-DIMOB on a single genome.

This wrapper isolates IslandPath execution inside a temporary directory,
and automatically organizes outputs by genome, following the same
architecture used by VirSorter2, DefenseFinder and other modules.

Input
-----
GenBank (.gbk)

Output
------
raw_outputs/<GENOME_NAME>/islandpath.gff

Author
------
Xtt Population Genomics Pipeline
"""

from __future__ import annotations

import argparse
import logging
import shutil
import subprocess
import tempfile
from pathlib import Path


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

def run_command(command: list[str]) -> None:
    """
    Execute external command.
    """

    LOGGER.info("Running:")
    LOGGER.info(" ".join(command))

    result = subprocess.run(
        command,
        capture_output=True,
        text=True,
    )

    if result.stdout:
        LOGGER.info(result.stdout)

    if result.returncode != 0:
        LOGGER.error(result.stderr)
        raise RuntimeError("IslandPath execution failed.")


# ============================================================
# Main
# ============================================================

def main():

    parser = argparse.ArgumentParser(
        description="Run IslandPath-DIMOB."
    )

    parser.add_argument(
        "--genome",
        required=True,
        type=Path,
        help="Input GenBank file."
    )

    parser.add_argument(
        "--outdir",
        required=True,
        type=Path,
        help="Directory where per-genome outputs will be created."
    )

    parser.add_argument(
        "--threads",
        type=int,
        default=1,
        help="Reserved for future compatibility."
    )

    args = parser.parse_args()

    genome = args.genome.resolve()
    outdir = args.outdir.resolve()

    if not genome.exists():
        raise FileNotFoundError(genome)

    # --------------------------------------------
    # Derive genome name and output directory
    # --------------------------------------------
    genome_name = genome.stem  # AB22007.gbk → AB22007

    genome_outdir = outdir / genome_name
    genome_outdir.mkdir(parents=True, exist_ok=True)

    output_gff = genome_outdir / "islandpath.gff"

    LOGGER.info("Genome : %s", genome.name)
    LOGGER.info("Output directory : %s", genome_outdir)

    # --------------------------------------------
    # Temporary execution directory
    # --------------------------------------------
    with tempfile.TemporaryDirectory(prefix="islandpath_") as tmpdir:

        tmpdir = Path(tmpdir)

        tmp_gbk = tmpdir / genome.name
        shutil.copy2(genome, tmp_gbk)

        tmp_output = tmpdir / "islandpath_output.gff"

        command = [
            "islandpath",
            str(tmp_gbk),
            str(tmp_output),
        ]

        run_command(command)

        if not tmp_output.exists():
            raise FileNotFoundError(
                "IslandPath did not produce the expected GFF output."
            )

        shutil.copy2(tmp_output, output_gff)

    LOGGER.info("Finished successfully.")
    LOGGER.info("IslandPath output written to: %s", output_gff)


if __name__ == "__main__":
    main()
