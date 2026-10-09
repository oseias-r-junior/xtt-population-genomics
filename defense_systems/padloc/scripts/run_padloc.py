#!/usr/bin/env python3
"""
run_padloc.py

Run PADLOC on all genomes from the XTT Population Genomics pipeline.

Author
------
Oseias Rodrigues Feitosa Junior

Pipeline
--------
XTT Population Genomics

Description
-----------
Executes PADLOC independently for every genome.

Expected input
--------------

annotations/
        fna/
            *.fna

Results
-------

defense_systems/
        padloc/
            results/
                <genome>/

Each genome receives its own directory containing:

    <genome>_padloc.csv
    <genome>_padloc.gff
    <genome>.domtblout
    <genome>_prodigal.faa
    <genome>_prodigal.gff

Features
--------

✓ Resume interrupted runs
✓ Skip completed genomes
✓ Automatic validation of input genomes
✓ Individual logs
✓ Runtime statistics
✓ Summary report
"""

from pathlib import Path
from typing import List

import argparse
import shutil
import subprocess
import sys
import time


# ------------------------------------------------------------
# utilities
# ------------------------------------------------------------

def timestamp():
    return time.strftime("%Y-%m-%d %H:%M:%S")


def log(msg):
    print(f"[{timestamp()}] {msg}", flush=True)


# ------------------------------------------------------------
# completed?
# ------------------------------------------------------------

def completed(outdir: Path, genome: str) -> bool:

    expected = [
        outdir / f"{genome}_padloc.csv",
        outdir / f"{genome}_padloc.gff",
        outdir / f"{genome}.domtblout",
        outdir / f"{genome}_prodigal.faa",
        outdir / f"{genome}_prodigal.gff",
    ]

    return all(
        f.exists() and f.stat().st_size > 0
        for f in expected
    )


# ------------------------------------------------------------
# discover genomes
# ------------------------------------------------------------

def discover_genomes(input_dir: Path) -> List[Path]:

    gbk_dir = input_dir.parent / "gbk"

    genomes = []

    all_fna = sorted(input_dir.glob("*.fna"))

    for fna in all_fna:

        gbk = gbk_dir / f"{fna.stem}.gbk"

        if not gbk.exists():
            log(f"[WARN ] Missing GBK for {fna.stem}; skipped.")
            continue

        genomes.append(fna)

    log(f"[INFO ] Valid genomes   : {len(genomes)}")
    log(f"[INFO ] Ignored genomes : {len(all_fna)-len(genomes)}")

    return genomes


# ------------------------------------------------------------
# run one genome
# ------------------------------------------------------------

def run_one(
        fna: Path,
        output: Path,
        threads: int,
        force: bool,
):

    genome = fna.stem

    outdir = output / genome

    logfile = outdir / f"{genome}.log"

    if force and outdir.exists():
        shutil.rmtree(outdir)

    outdir.mkdir(parents=True, exist_ok=True)

    if completed(outdir, genome):
        log(f"[SKIP ] {genome}")
        return "skipped"

    cmd = [
        "padloc",
        "--fna", str(fna),
        "--outdir", str(outdir),
        "--cpu", str(threads),
    ]

    log(f"[RUN  ] {genome}")

    start = time.time()

    try:

        with logfile.open("w") as fout:

            proc = subprocess.run(
                cmd,
                stdout=fout,
                stderr=subprocess.STDOUT,
                timeout=21600,
            )

    except subprocess.TimeoutExpired:

        log(f"[FAIL ] {genome} (timeout)")
        return "failed"

    elapsed = time.time() - start

    if proc.returncode != 0:

        log(f"[FAIL ] {genome} (return code {proc.returncode})")
        return "failed"

    if not completed(outdir, genome):

        log(f"[FAIL ] {genome} (missing output)")
        return "failed"

    csv = outdir / f"{genome}_padloc.csv"

    try:

        n_hits = max(
            0,
            sum(1 for _ in csv.open()) - 1
        )

    except Exception:

        n_hits = -1

    if n_hits == 0:

        log(f"[EMPTY] {genome} ({elapsed/60:.1f} min)")
        return "empty"

    log(
        f"[ DONE] {genome} "
        f"({n_hits} systems; {elapsed/60:.1f} min)"
    )

    return "done"


# ------------------------------------------------------------
# main
# ------------------------------------------------------------

def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--input",
        required=True,
        type=Path,
    )

    parser.add_argument(
        "--output",
        required=True,
        type=Path,
    )

    parser.add_argument(
        "--threads",
        type=int,
        default=4,
    )

    parser.add_argument(
        "--force",
        action="store_true",
    )

    parser.add_argument(
        "--max-genomes",
        type=int,
        default=None,
        help="Process only the first N genomes (for testing).",
    )

    args = parser.parse_args()

    genomes = discover_genomes(args.input)

    if args.max_genomes is not None:
        genomes = genomes[:args.max_genomes]

    if not genomes:
        sys.exit("No valid genomes found.")

    done = 0
    empty = 0
    skipped = 0
    failed = 0

    overall = time.time()

    log("=========================================")
    log("PADLOC")
    log(f"Genomes : {len(genomes)}")
    log(f"Threads : {args.threads}")
    log("=========================================")

    for i, genome in enumerate(genomes, start=1):

        log(f"[{i}/{len(genomes)}] {genome.stem}")

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

        else:
            failed += 1

    elapsed = time.time() - overall

    log("")
    log("=========================================")
    log("Summary")
    log("-----------------------------------------")
    log(f"Processed (systems found) : {done}")
    log(f"Processed (no systems)    : {empty}")
    log(f"Skipped                   : {skipped}")
    log(f"Failed                    : {failed}")
    log(f"Elapsed                   : {elapsed/3600:.2f} h")
    log("=========================================")


if __name__ == "__main__":
    main()