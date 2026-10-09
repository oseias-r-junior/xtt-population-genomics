#!/usr/bin/env python3
"""
mobilome/genomic_islands/alienhunter/scripts/parse_alienhunter_islands.py

Parse AlienHunter's native alienhunter.embl output (NOT alienhunter.sco) into
a per-region matrix, one row per predicted island.

IMPORTANT: alienhunter.sco is the RAW per-window IVOM/KL sliding-window score
track (all windows, dense, unthresholded -- used only to draw alienhunter.plot).
alienhunter.embl is AlienHunter's OWN native post-processed output: it has
already applied its own internal score threshold (reported per-genome in each
record's /note="threshold: X.XXX") and already merged adjacent above-threshold
windows into contiguous candidate regions. No additional threshold/merge logic
is needed or appropriate here -- building a second, redundant threshold+merge
step on top of the tool's own would double-apply stringency with no principled
basis. This mirrors why IslandPath's own island calls are used directly rather
than re-derived from a lower-level signal.

AlienHunter also flags some regions as "probably region overlapping rRNA
operon" -- a well-known AlienHunter false-positive mode, since rRNA operons
have atypical, highly conserved k-mer composition that mimics horizontally
acquired DNA without actually being foreign. These are NOT dropped (the
consensus step downstream can use tool agreement to filter them out
naturally), but are flagged so they can be treated with appropriate caution.
"""

from pathlib import Path
import re
import sys

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parents[3]

sys.path.append(str(PROJECT_ROOT / "scripts"))
sys.path.append(str(PROJECT_ROOT / "workflow" / "config"))

from config import *
from genomes import GENOMES

import pandas as pd

# ======================================================
# INPUT
# ======================================================

RAW_DIR = ALIENHUNTER_RAW

# ======================================================
# OUTPUT
# ======================================================

ALIENHUNTER_MATRICES = ALIENHUNTER / "matrices"
ALIENHUNTER_MATRICES.mkdir(parents=True, exist_ok=True)

OUT_FILE = ALIENHUNTER_MATRICES / "alienhunter_islands.tsv"

LOCUS_RE = re.compile(r"^FT\s+misc_feature\s+(\d+)\.\.(\d+)")
SCORE_RE = re.compile(r"/score=([\d.]+)")
THRESHOLD_RE = re.compile(r'/note="threshold:\s*([\d.]+)')


def parse_embl(path: Path):
    """Parse one alienhunter.embl file into a list of region dicts."""
    records = []
    current = None

    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")

            m = LOCUS_RE.match(line)
            if m:
                if current is not None:
                    records.append(current)
                current = {
                    "start": int(m.group(1)),
                    "end": int(m.group(2)),
                    "score": None,
                    "threshold": None,
                    "rrna_operon_flag": False,
                }
                continue

            if current is None:
                continue

            m = SCORE_RE.search(line)
            if m:
                current["score"] = float(m.group(1))
                continue

            m = THRESHOLD_RE.search(line)
            if m:
                current["threshold"] = float(m.group(1))
                if "rRNA operon" in line:
                    current["rrna_operon_flag"] = True
                continue

    if current is not None:
        records.append(current)

    return records


records = []

OFFICIAL_GENOMES = set(GENOMES)

for genome_dir in sorted(RAW_DIR.iterdir()):
    if not genome_dir.is_dir():
        continue

    genome = genome_dir.name  # AlienHunter raw_outputs dirs already use normalized names

    # Skip anything not in the official genome list (e.g. stray/debris
    # directories left over from malformed past invocations, such as an
    # "fna" directory -- not a real genome).
    if genome not in OFFICIAL_GENOMES:
        continue

    embl_file = genome_dir / "alienhunter.embl"
    if not embl_file.exists():
        continue

    islands = parse_embl(embl_file)

    for i, isl in enumerate(islands, start=1):
        records.append({
            "Genome": genome,
            "Region_ID": f"alienhunter_gi{i}",
            "Start": isl["start"],
            "End": isl["end"],
            "Length": isl["end"] - isl["start"] + 1,
            "Score": isl["score"],
            "Threshold": isl["threshold"],
            "RRNA_Operon_Flag": isl["rrna_operon_flag"],
        })

result = pd.DataFrame(records)
result.to_csv(OUT_FILE, sep="\t", index=False)

print(f"Genomes scanned: {sum(1 for d in RAW_DIR.iterdir() if d.is_dir())}")
print(f"Islands extracted: {len(result)}")
if not result.empty:
    print(f"Genomes with >=1 island: {result['Genome'].nunique()}")
    print(f"Regions flagged as rRNA-operon-overlap: {int(result['RRNA_Operon_Flag'].sum())}")
print(f"Saved: {OUT_FILE}")
