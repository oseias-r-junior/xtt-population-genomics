#!/usr/bin/env python3

from pathlib import Path
from datetime import datetime
import subprocess
import sys
import time
import pandas as pd

cycle = 1

print(
    f"\n[{datetime.now()}]",
    flush=True
)

print(
    f"PHASTEST cycle {cycle}",
    flush=True
)

SCRIPT_DIR = Path(__file__).resolve().parent

# ======================================
# CONFIG
# ======================================

SLEEP_TIME = 1800

# ======================================
# SCRIPTS
# ======================================

SUBMIT_SCRIPT = (
    SCRIPT_DIR
    / "submit_phastest_batch.py"
)

POLL_SCRIPT = (
    SCRIPT_DIR
    / "poll_phastest_jobs.py"
)

DOWNLOAD_SCRIPT = (
    SCRIPT_DIR
    / "download_phastest_results.py"
)

PARSE_SCRIPT = (
    SCRIPT_DIR
    / "parse_phastest_summary.py"
)

STATUS_FILE = (
    SCRIPT_DIR.parent
    / "jobs"
    / "phastest_status.tsv"
)

# ======================================
# RUNNER
# ======================================

def run_script(script):

    print("\n" + "=" * 60)
    print(
        f"Running {script.name}",
        flush=True
    )
    print("=" * 60)

    result = subprocess.run(
        [sys.executable, str(script)]
    )

    if result.returncode == 0:

        print(
            f"SUCCESS: {script.name}",
            flush=True
        )

    else:

        print(
            f"FAILED: {script.name}",
            flush=True
        )

    return result.returncode == 0

# ======================================
# MAIN LOOP
# ======================================

while True:

    print("\n", flush=True)
    print("#" * 80, flush=True)
    print(f"PHASTEST cycle {cycle}", flush=True)
    print("#" * 80, flush=True)

    run_script(SUBMIT_SCRIPT)

    run_script(POLL_SCRIPT)

    run_script(DOWNLOAD_SCRIPT)

    run_script(PARSE_SCRIPT)

    # ==================================
    # REPORT
    # ==================================

    if STATUS_FILE.exists():

        df = pd.read_csv(
            STATUS_FILE,
            sep="\t"
        )

        with_phage = (
            df["Status"]
            .eq("COMPLETE_WITH_PHAGE")
            .sum()
        )

        no_phage = (
            df["Status"]
            .eq("COMPLETE_NO_PHAGE")
            .sum()
        )

        total = len(df)

        print("\nProgress report")

        print(
            f"Total jobs: {total}"
        )

        print(
            f"With prophage: {with_phage}"
        )

        print(
            f"No prophage: {no_phage}"
        )

    print(
        f"\nSleeping {SLEEP_TIME} seconds..."
    )

    time.sleep(
        SLEEP_TIME
    )

    cycle += 1