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

STATUS_FILE = (
    PHASTEST_JOBS
    / "phastest_status.tsv"
)

OUTDIR = PHASTEST_DOWNLOADS

OUTDIR.mkdir(
    parents=True,
    exist_ok=True
)

# ======================================
# LOAD STATUS
# ======================================

df = pd.read_csv(
    STATUS_FILE,
    sep="\t"
)

print(f"Loaded {len(df)} records")

# ======================================
# DOWNLOAD
# ======================================

for idx, row in df.iterrows():

    genome = row["Genome"]

    status = str(row["Status"])

    zip_url = row["ZIP"]

    # ==================================
    # ONLY COMPLETE JOBS WITH PHAGES
    # ==================================

    if status != "COMPLETE_WITH_PHAGE":

        print(
            f"Skipping {genome} ({status})"
        )

        continue

    # ==================================
    # CHECK ZIP URL
    # ==================================

    if pd.isna(zip_url):

        print(f"No ZIP for {genome}")

        continue

    zip_url = str(zip_url)

    if not zip_url.startswith("http"):

        zip_url = "https://" + zip_url.lstrip("/")

    outfile = OUTDIR / f"{genome}.zip"

    # ==================================
    # SKIP EXISTING
    # ==================================

    if outfile.exists():

        print(f"Already downloaded: {genome}")

        continue

    print(f"\nDownloading: {genome}")

    print(zip_url)

    try:

        response = requests.get(
            zip_url,
            stream=True
        )

        response.raise_for_status()

        with open(outfile, "wb") as f:

            for chunk in response.iter_content(
                chunk_size=8192
            ):

                f.write(chunk)

        print(f"Saved: {outfile}")

    except Exception as e:

        print(f"ERROR downloading {genome}")

        print(e)

    # polite delay
    time.sleep(5)

print("\nFinished downloads.")