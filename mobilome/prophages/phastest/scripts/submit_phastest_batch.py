#!/usr/bin/env python3

from pathlib import Path
import sys

import pandas as pd
import requests
import time
import json

SCRIPT_DIR = Path(__file__).resolve().parent

PROJECT_ROOT = SCRIPT_DIR.parents[3]

sys.path.append(
    str(PROJECT_ROOT / "scripts")
)

from config import *
from normalize_genome_names import normalize_genome_name

# ======================================
# INPUT
# ======================================

GENOME_DIR = RAW_GENOMES

OUTFILE = (
    PHASTEST_JOBS
    / "phastest_jobs.tsv"
)

submitted = set()

if OUTFILE.exists():

    old_df = pd.read_csv(
        OUTFILE,
        sep="\t"
    )

    submitted = set(
        old_df["Genome"]
        .astype(str)
    )

    print(
        f"{len(submitted)} genomes already submitted"
    )

API_URL = "https://phastest.ca/phastest_api"

# ======================================
# PARAMETERS
# ======================================

BATCH_SIZE = 5
SLEEP_BETWEEN = 180  # seconds

# ======================================
# GENOMES
# ======================================

all_genomes = sorted(
    GENOME_DIR.glob("*.fasta")
)

genomes = []

for genome_file in all_genomes:

    genome = normalize_genome_name(
        genome_file.stem
    )

    if genome in submitted:

        continue

    genomes.append(
        genome_file
    )

print(
    f"{len(genomes)} genomes pending submission"
)

if OUTFILE.exists():

    old_df = pd.read_csv(
        OUTFILE,
        sep="\t"
    )

    results = old_df.to_dict(
        "records"
    )

else:

    results = []

# ======================================
# SUBMIT
# ======================================

for i, genome_file in enumerate(genomes):

    genome_name = normalize_genome_name(
        genome_file.stem
    )

    if genome_name in submitted:

        print(
            f"Skipping {genome_name}"
        )

        continue

    print(
        f"\nSubmitting: {genome_name}"
    )

    try:
        with open(genome_file, "rb") as f:

            response = requests.post(
                API_URL,
                data=f
            )

        data = response.json()

        print(data)

        results.append({
            "Genome": genome_name,
            "Job_ID": data.get("job_id"),
            "Status": data.get("status")
        })
        submitted.add(genome_name)

    except Exception as e:

        print(f"ERROR: {genome_name}")
        print(e)

        results.append({
            "Genome": genome_name,
            "Job_ID": "ERROR",
            "Status": str(e)
        })

    # ==================================
    # SAVE PROGRESS
    # ==================================

    df = pd.DataFrame(results)

    OUTFILE.parent.mkdir(
        parents=True,
        exist_ok=True
    )

    df.to_csv(
        OUTFILE,
        sep="\t",
        index=False
    )

    # ==================================
    # RATE LIMIT CONTROL
    # ==================================

    if (i + 1) % BATCH_SIZE == 0:

        print(
            f"\nSleeping {SLEEP_BETWEEN} seconds..."
        )

        time.sleep(SLEEP_BETWEEN)

print("\nFinished.")
print(f"Saved: {OUTFILE}")