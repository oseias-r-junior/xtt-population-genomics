#!/usr/bin/env python3

import requests
import pandas as pd
import time
from pathlib import Path

from xtt_population_genomics.scripts.metadata.normalize_genome_names import normalize_genome_name

# ======================================
# INPUT
# ======================================

GENOME_DIR = Path(
    "~/diel_paper_colab/input_data_mge_ds/Xtt_project/raw_genomes"
).expanduser()

OUTFILE = Path(
    "~/diel_paper_colab/results/phastest/phastest_jobs.tsv"
).expanduser()

API_URL = "https://phastest.ca/phastest_api"

# ======================================
# PARAMETERS
# ======================================

BATCH_SIZE = 5
SLEEP_BETWEEN = 180  # seconds

# ======================================
# GENOMES
# ======================================

genomes = sorted(
    GENOME_DIR.glob("*.fasta")
)

results = []

# ======================================
# SUBMIT
# ======================================

for i, genome_file in enumerate(genomes):

    genome_name = normalize_genome_name(
        genome_file.stem
    )

    print(f"\nSubmitting: {genome_name}")

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