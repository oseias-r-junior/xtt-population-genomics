#!/usr/bin/env python3
"""
Run AlienHunter on a single annotated genome.

This wrapper executes AlienHunter for one genome and stores all generated
output files in a dedicated genome directory.

Example
-------
python run_alienhunter.py \
    --genome annotations/gbk/AB22007.gbk \
    --outdir mobilome/genomic_islands/alienhunter/raw_outputs

Author
------
Xtt Population Genomics Pipeline
"""

from __future__ import annotations

import argparse
import logging
import shutil
import subprocess
import time
from pathlib import Path

# ==========================================================
# Logging
# ==========================================================

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
)

LOGGER = logging.getLogger(__name__)

# ==========================================================
# Locate AlienHunter
# ==========================================================

PROJECT_ROOT = Path(__file__).resolve().parents[5]

DEFAULT_ALIENHUNTER = (
    PROJECT_ROOT
    / "third_party"
    / "alienhunter"
    / "alien_hunter-1.7"
    / "alien_hunter"
)

FNA_DIR = (
    PROJECT_ROOT
    / "annotations"
    / "fna"
)

# ==========================================================
# Arguments
# ==========================================================

parser = argparse.ArgumentParser()

parser.add_argument("--genome", required=True, type=Path)
parser.add_argument("--outdir", required=True, type=Path)
parser.add_argument("--alienhunter", default=DEFAULT_ALIENHUNTER, type=Path)

parser.add_argument(
    "--optimize",
    action="store_true",
    help="Enable AlienHunter boundary optimization (-c).",
)

args = parser.parse_args()

# ==========================================================
# Validate input
# ==========================================================

if not args.genome.exists():
    raise FileNotFoundError(args.genome)

if not args.alienhunter.exists():
    raise FileNotFoundError(args.alienhunter)

genome_name = args.genome.stem

# FNA correspondente ao GBK informado
fna_source = FNA_DIR / f"{genome_name}.fna"

if not fna_source.exists():
    raise FileNotFoundError(fna_source)

genome_dir = args.outdir / genome_name
genome_dir.mkdir(parents=True, exist_ok=True)

# ==========================================================
# Temporary input copy (.fna)
# ==========================================================

working_fna = genome_dir / f"{genome_name}.fna"

if not working_fna.exists():
    shutil.copy2(fna_source, working_fna)

output_prefix = genome_dir / "alienhunter"

# ==========================================================
# Build command
# ==========================================================

command = [
    str(args.alienhunter),
    str(working_fna),
    str(output_prefix),
]

if args.optimize:
    command.append("-c")

LOGGER.info("Genome: %s", genome_name)
LOGGER.info("Running AlienHunter...")
LOGGER.info("Command:\n%s", " ".join(command))

# ==========================================================
# Execute
# ==========================================================

start = time.perf_counter()

result = subprocess.run(
    command,
    capture_output=True,
    text=True,
)

(genome_dir / "stdout.log").write_text(result.stdout)
(genome_dir / "stderr.log").write_text(result.stderr)

# ==========================================================
# Check return code
# ==========================================================

if result.returncode != 0:
    LOGGER.error(result.stderr)
    raise RuntimeError(
        f"AlienHunter exited with code {result.returncode}"
    )

elapsed = time.perf_counter() - start
LOGGER.info("Elapsed time: %.2f min", elapsed / 60)

LOGGER.info("Finished successfully.")
LOGGER.info("Output directory: %s", genome_dir)

# ==========================================================
# Normalize main EMBL output name
# ==========================================================

candidates = [
    genome_dir / "alienhunter",
    genome_dir / "alienhunter.",
]

main_output = None
for c in candidates:
    if c.exists():
        main_output = c
        break

if main_output is not None:
    target = genome_dir / "alienhunter.embl"

    if target.exists():
        LOGGER.warning("Target file already exists, overwriting: %s", target)
        target.unlink()

    LOGGER.info("Renaming AlienHunter EMBL output to: %s", target.name)
    main_output.rename(target)
else:
    LOGGER.warning(
        "Expected main AlienHunter output file not found (checked: %s)",
        ", ".join(str(c) for c in candidates),
    )

# ==========================================================
# Verify expected output files
# ==========================================================

expected_files = [
    genome_dir / "alienhunter.embl",
    genome_dir / "alienhunter.plot",
    genome_dir / "alienhunter.sco",
]

missing = [f.name for f in expected_files if not f.exists()]

if missing:
    LOGGER.warning("Missing AlienHunter outputs: %s", ", ".join(missing))
else:
    LOGGER.info("All AlienHunter outputs successfully generated.")

# ==========================================================
# Count predicted genomic islands
# ==========================================================

embl_file = genome_dir / "alienhunter.embl"

n_regions = 0

if embl_file.exists():
    with embl_file.open() as f:
        for line in f:
            if line.startswith("FT   misc_feature"):
                n_regions += 1

    LOGGER.info("Predicted genomic islands: %d", n_regions)
else:
    LOGGER.warning("Cannot count regions — EMBL file missing.")