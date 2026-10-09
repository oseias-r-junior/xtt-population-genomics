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
# PATHS
# ======================================

JOBS_FILE = (
    PHASTEST_JOBS
    / "phastest_jobs.tsv"
)

OUTFILE = (
    PHASTEST_JOBS
    / "phastest_status.tsv"
)

API_URL = "https://phastest.ca/phastest_api?acc="

# ======================================
# LOAD JOBS
# ======================================

jobs = pd.read_csv(
    JOBS_FILE,
    sep="\t"
)

print(f"Loaded {len(jobs)} jobs")

# ======================================
# LOAD PREVIOUS STATUS
# ======================================

completed = set()

if OUTFILE.exists():

    old_status = pd.read_csv(
        OUTFILE,
        sep="\t"
    )

    completed = set(
        old_status.loc[
            old_status["Status"].isin([
                "COMPLETE_WITH_PHAGE",
                "COMPLETE_NO_PHAGE"
            ]),
            "Genome"
        ]
    )

    print(
        f"{len(completed)} completed genomes"
    )

# ======================================
# CHECK STATUS
# ======================================

results = []

for idx, row in jobs.iterrows():

    genome = row["Genome"]
    if genome in completed:

        print(
            f"Skipping completed: {genome}"
        )

        old_row = old_status[
            old_status["Genome"] == genome
        ].iloc[0]

        results.append(
            old_row.to_dict()
        )

        continue

    job_id = row["Job_ID"]
    print(f"\nChecking: {genome}")

    # Skip failed jobs
    if job_id == "ERROR":

        results.append({
            "Genome": genome,
            "Job_ID": job_id,
            "Status": "ERROR"
        })

        continue

    try:

        response = requests.get(
            API_URL + str(job_id),
            timeout=60
        )

        data = response.json()

        print(data)

        raw_status = str(
            data.get("status", "")
        )

        zip_url = data.get("zip")

        if raw_status == "Complete":

            if zip_url and zip_url != "N/A":

                status = "COMPLETE_WITH_PHAGE"

            else:

                status = "COMPLETE_NO_PHAGE"

        else:

            status = raw_status

        results.append({

            "Genome": genome,

            "Job_ID": job_id,

            "Status": status,

            "URL": data.get("url"),

            "ZIP": zip_url

        })

    except Exception as e:

        print(f"ERROR: {genome}")
        print(e)

        results.append({

            "Genome": genome,

            "Job_ID": job_id,

            "Status": "QUERY_ERROR",

            "URL": "",

            "ZIP": "",

            "Summary": str(e)

        })

    # polite delay
    time.sleep(2)

# ======================================
# SAVE
# ======================================

status_df = pd.DataFrame(results)

status_df.to_csv(
    OUTFILE,
    sep="\t",
    index=False
)

print("\nSaved:")
print(OUTFILE)