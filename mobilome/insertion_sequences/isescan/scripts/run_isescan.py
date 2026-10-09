#!/usr/bin/env python3
"""
run_isescan.py — Robust wrapper for ISEScan

States:
    done    — IS elements detected
    empty   — successful execution, no IS detected
    skipped — genome already processed
    failed  — real execution error

Author:
    Oseias Rodrigues Feitosa Junior
"""

from pathlib import Path
from typing import List
import argparse
import subprocess
import time
import shutil
import sys


# ------------------------------------------------------------
# utilities
# ------------------------------------------------------------

def timestamp():
    return time.strftime("%Y-%m-%d %H:%M:%S")


def log(msg):
    print(f"[{timestamp()}] {msg}", flush=True)


# ------------------------------------------------------------
# Detect whether ISEScan produced complete output (IS found)
# ------------------------------------------------------------

def completed(outdir: Path, genome: str) -> bool:
    if not outdir.exists():
        return False

    # Recursively search for the TSV file produced by ISEScan
    tsv_matches = list(outdir.rglob(f"{genome}.fna.tsv"))
    if not tsv_matches:
        return False

    # Check all possible result directories
    for tsv in tsv_matches:
        result_dir = tsv.parent
        expected = [
            result_dir / f"{genome}.fna.tsv",
            result_dir / f"{genome}.fna.csv",
            result_dir / f"{genome}.fna.sum",
            result_dir / f"{genome}.fna.gff",
        ]
        if all(x.exists() and x.stat().st_size > 0 for x in expected):
            return True

    return False


# ------------------------------------------------------------
# Execute ISEScan for a single genome
# ------------------------------------------------------------

def run_one(fna: Path,
            output: Path,
            threads: int,
            force: bool) -> str:

    genome = fna.stem
    outdir = output / genome
    logfile = outdir / f"{genome}.log"

    if force and outdir.exists():
        shutil.rmtree(outdir)

    outdir.mkdir(parents=True, exist_ok=True)

    # State: SKIPPED
    if completed(outdir, genome):
        log(f"[SKIP ] {genome}")
        return "skipped"

    cmd = [
        "isescan.py",
        "--seqfile", str(fna),
        "--output", str(outdir),
        "--no-FragGeneScan",
        "--nthread", str(threads),
    ]

    log(f"[RUN  ] {genome}")

    start = time.time()

    # Run with timeout protection
    try:
        with logfile.open("w") as fout:
            proc = subprocess.run(
                cmd,
                stdout=fout,
                stderr=subprocess.STDOUT,
                timeout=21600,  # 6 hours
            )
    except subprocess.TimeoutExpired:
        log(f"[FAIL ] {genome} (execution timed out)")
        return "failed"

    elapsed = time.time() - start

    # State: FAILED (real error)
    if proc.returncode != 0:
        log(f"[FAIL ] {genome} (return code {proc.returncode})")
        return "failed"

    # State: DONE (IS found)
    if completed(outdir, genome):
        log(f"[ DONE] {genome} ({elapsed/60:.1f} min)")
        return "done"

    # ------------------------------------------------------------
    # State: EMPTY (successful run, no IS detected)
    # Confirmed by log content — prevents false EMPTY
    # ------------------------------------------------------------

    try:
        text = logfile.read_text(errors="ignore")
    except Exception as e:
        log(f"[WARN ] Unable to read log file for {genome}: {e}")
        text = ""

    if "No IS element was found" in text:
        log(f"[EMPTY] {genome} ({elapsed/60:.1f} min; no IS detected)")
        return "empty"

    # Silent failure: no TSV, no "No IS element was found"
    log(f"[FAIL ] {genome} (no output and no 'No IS element was found')")
    return "failed"


# ------------------------------------------------------------
# Discover valid genomes (requires matching GBK file)
# ------------------------------------------------------------

def discover_genomes(input_dir: Path) -> List[Path]:
    genomes = []
    gbk_dir = input_dir.parent / "gbk"

    all_fna = sorted(input_dir.glob("*.fna"))

    for fna in all_fna:
        gbk = gbk_dir / f"{fna.stem}.gbk"

        if not gbk.exists():
            log(f"[WARN ] Missing GBK for {fna.stem}; genome skipped.")
            continue

        genomes.append(fna)

    log(f"[INFO ] Valid genomes   : {len(genomes)}")
    log(f"[INFO ] Ignored genomes : {len(all_fna) - len(genomes)}")

    return genomes


# ------------------------------------------------------------
# main
# ------------------------------------------------------------

def main():

    parser = argparse.ArgumentParser()

    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--threads", default=4, type=int)
    parser.add_argument("--force", action="store_true")

    args = parser.parse_args()

    valid_genomes = discover_genomes(args.input)

    if not valid_genomes:
        sys.exit("No valid genomes found (all inputs invalid or missing GBK).")

    log("=========================================")
    log("ISEScan Execution")
    log(f"Valid genomes : {len(valid_genomes)}")
    log(f"Threads       : {args.threads}")
    log("=========================================")

    done = 0
    empty = 0
    skipped = 0
    failed = 0

    overall = time.time()

    # ------------------------------------------------------------
    # Main loop
    # ------------------------------------------------------------

    for i, genome in enumerate(valid_genomes, start=1):

        log(f"[{i}/{len(valid_genomes)}] Processing {genome.stem}")

        result = run_one(
            genome,
            args.output,
            args.threads,
            args.force,
        )

        if result == "done":
            done += 1
        elif result == "empty":
            empty += 1
        elif result == "skipped":
            skipped += 1
        elif result == "failed":
            failed += 1
        else:
            failed += 1  # defensive fallback

    elapsed = time.time() - overall

    log("")
    log("=========================================")
    log("Summary Report")
    log("-----------------------------------------")
    log(f"IS detected        : {done}")
    log(f"No IS detected     : {empty}")
    log(f"Skipped            : {skipped}")
    log(f"Failed             : {failed}")
    log(f"Total runtime      : {elapsed/3600:.2f} hours")
    log("=========================================")


if __name__ == "__main__":
    main()