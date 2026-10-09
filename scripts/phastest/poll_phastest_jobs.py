#!/usr/bin/env python3

import requests
import pandas as pd
import time
import sys

from pathlib import Path

# ======================================
# PATHS
# ======================================

JOBS_FILE = Path(
    "~/diel_paper_colab/results/phastest/phastest_jobs.tsv"
).expanduser()

OUTFILE = Path(
    "~/diel_paper_colab/results/phastest/phastest_status.tsv"
).expanduser()

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
# CHECK STATUS
# ======================================

results = []

for idx, row in jobs.iterrows():

    genome = row["Genome"]
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
            API_URL + str(job_id)
        )

        data = response.json()

        print(data)

        results.append({

            "Genome": genome,

            "Job_ID": job_id,

            "Status": data.get("status"),

            "URL": data.get("url"),

            "ZIP": data.get("zip"),

            "Summary": data.get("summary")

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