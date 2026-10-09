import requests
import json
from pathlib import Path
import pandas as pd
import time

# ====================================
# INPUT DIRECTORY
# ====================================

GENOME_DIR = Path(
    "~/diel_paper_colab/input_data_mge_ds/Xtt_project/raw_genomes"
).expanduser()

OUTDIR = Path(
    "~/diel_paper_colab/results/phastest"
).expanduser()

OUTDIR.mkdir(exist_ok=True)

# ====================================
# SELECT TEST GENOMES
# ====================================

genomes = [
    "AB22007.fasta",
    "BW23007.fasta",
    "AH23019.fasta"
]

# ====================================
# SUBMIT
# ====================================

results = []

for genome in genomes:

    genome_path = GENOME_DIR / genome

    print(f"\nSubmitting: {genome}")

    with open(genome_path, "rb") as f:

        response = requests.post(
            "https://phastest.ca/phastest_api",
            data=f
        )

    try:
        data = response.json()

    except Exception:

        print("ERROR:")
        print(response.text)

        continue

    print(data)

    results.append({
        "Genome": genome,
        "job_id": data.get("job_id"),
        "status": data.get("status")
    })

    time.sleep(5)

# ====================================
# SAVE
# ====================================

df = pd.DataFrame(results)

outfile = OUTDIR / "phastest_jobs.tsv"

df.to_csv(
    outfile,
    sep="\t",
    index=False
)

print(f"\nSaved: {outfile}")