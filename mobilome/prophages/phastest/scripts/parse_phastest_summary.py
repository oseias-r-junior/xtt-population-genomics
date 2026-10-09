#!/usr/bin/env python3

from pathlib import Path
import sys

SCRIPT_DIR = Path(__file__).resolve().parent

PROJECT_ROOT = SCRIPT_DIR.parents[3]

sys.path.append(
    str(PROJECT_ROOT / "scripts")
)

from config import *

import json
import zipfile
import pandas as pd

# ======================================
# INPUT
# ======================================

ZIP_DIR = PHASTEST_DOWNLOADS

# ======================================
# OUTPUT
# ======================================

PROPHAGE_MATRICES.mkdir(
    parents=True,
    exist_ok=True
)

REGIONS_OUTFILE = (
    PROPHAGE_MATRICES
    / "prophage_regions.tsv"
)

SUMMARY_OUTFILE = (
    PROPHAGE_MATRICES
    / "prophage_summary.tsv"
)

# ======================================
# PARSE ZIP FILES
# ======================================

regions = []

for zip_file in sorted(
    ZIP_DIR.glob("*.zip")
):

    genome = zip_file.stem

    print(f"Processing {genome}")

    try:

        with zipfile.ZipFile(zip_file) as z:

            if (
                "predicted_phage_regions.json"
                not in z.namelist()
            ):

                print(
                    f"Missing JSON: {genome}"
                )

                continue

            data = json.loads(
                z.read(
                    "predicted_phage_regions.json"
                )
            )

            for region in data:

                regions.append({

                    "Genome":
                        genome,

                    "Region":
                        region.get(
                            "region"
                        ),

                    "Start":
                        region.get(
                            "start"
                        ),

                    "Stop":
                        region.get(
                            "stop"
                        ),

                    "Length_bp":
                        (
                            region.get(
                                "stop"
                            )
                            -
                            region.get(
                                "start"
                            )
                            + 1
                        ),

                    "Completeness":
                        region.get(
                            "completeness"
                        ),

                    "GC":
                        region.get(
                            "GC"
                        ),

                    "Most_Common_Phage":
                        region.get(
                            "most_common_phage"
                        )

                })

    except Exception as e:

        print(
            f"ERROR {genome}: {e}"
        )

# ======================================
# REGIONS TABLE
# ======================================

regions_df = pd.DataFrame(regions)

regions_df.to_csv(
    REGIONS_OUTFILE,
    sep="\t",
    index=False
)

print(
    f"Saved: {REGIONS_OUTFILE}"
)

STATUS_FILE = (
    PHASTEST_JOBS
    / "phastest_status.tsv"
)

status_df = pd.read_csv(
    STATUS_FILE,
    sep="\t"
)


# ======================================
# SUMMARY TABLE
# ======================================

summary_rows = []

for genome in status_df["Genome"]:

    subdf = regions_df[
        regions_df["Genome"] == genome
    ]

    if len(subdf) == 0:

        summary_rows.append({

            "Genome": genome,

            "Total_Prophages": 0,

            "Intact": 0,

            "Questionable": 0,

            "Incomplete": 0,

            "Mean_GC": None

        })

        continue

    row = {

        "Genome":
            genome,

        "Total_Prophages":
            len(subdf),

        "Intact":
            (
                subdf["Completeness"]
                .eq("intact")
                .sum()
            ),

        "Questionable":
            (
                subdf["Completeness"]
                .eq("questionable")
                .sum()
            ),

        "Incomplete":
            (
                subdf["Completeness"]
                .eq("incomplete")
                .sum()
            ),

        "Mean_GC":
            round(
                subdf["GC"].mean(),
                2
            )

    }

    summary_rows.append(row)

summary_df = pd.DataFrame(
    summary_rows
)

summary_df.to_csv(
    SUMMARY_OUTFILE,
    sep="\t",
    index=False
)

print(
    f"Saved: {SUMMARY_OUTFILE}"
)

print("\nFinished.")