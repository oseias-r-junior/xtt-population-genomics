import requests
import pandas as pd
from pathlib import Path
import time

# =====================================
# INPUT
# =====================================

JOB_FILE = Path(
    "~/diel_paper_colab/results/phastest/phastest_jobs.tsv"
).expanduser()

OUTDIR = Path(
    "~/diel_paper_colab/results/phastest"
).expanduser()

# =====================================
# LOAD JOBS
# =====================================

jobs = pd.read_csv(JOB_FILE, sep="\t")

results = []

# =====================================
# QUERY JOBS
# =====================================

for _, row in jobs.iterrows():

    genome = row["Genome"]
    job_id = row["job_id"]

    print(f"\nChecking: {genome} ({job_id})")

    url = f"https://phastest.ca/phastest_api?acc={job_id}"

    response = requests.get(url)

    try:
        data = response.json()

    except Exception:

        print("ERROR:")
        print(response.text)

        continue

    print(data)

    results.append({
        "Genome": genome,
        "job_id": job_id,
        "status": data.get("status"),
        "url": data.get("url"),
        "zip": data.get("zip"),
        "summary": data.get("summary")
    })

    time.sleep(2)

# =====================================
# SAVE
# =====================================

outdf = pd.DataFrame(results)

outfile = OUTDIR / "phastest_status.tsv"

outdf.to_csv(
    outfile,
    sep="\t",
    index=False
)

print(f"\nSaved: {outfile}")