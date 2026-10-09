#!/usr/bin/env python3
"""
mobilome/scripts/split_consensus_confidence.py

Split each MGE consensus table (prophage / genomic island / plasmid) into
two derived tables based on how many of the tools used for that category
agree on a region:

    - {category}_high_confidence.tsv : N_Tools_Support == TOTAL tools used
      for that category (ALL tools agree). This is the conservative,
      "all-or-nothing" set used for the main discussion/analysis (feeds
      the SSN, any downstream biological interpretation).
    - {category}_supplementary.tsv   : everything else (partial support --
      1 of N, up to N-1 of N). Kept for transparency/completeness, listed
      in a supplementary table, but NOT used in the main analysis.

Rationale (per project PI decision, 2026-09-28): rather than presenting a
consensus that mixes strong and weak multi-tool agreement together, the
main discussion should only include MGEs validated by the maximum number
of tools available for that category. Partial-support regions are not
claimed to be false; they are set aside under a deliberately conservative
policy.

Tool counts per category:
    prophage -- 3 tools (PHASTEST x VirSorter2 x geNomad)
    GI       -- 3 tools (AlienHunter x IslandPath x PanISLE)
    plasmid  -- 2 tools (MOB-suite x geNomad)
"""

from pathlib import Path
import sys

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parents[1]

sys.path.append(str(PROJECT_ROOT / "scripts"))

from config import *

import pandas as pd

# ======================================================
# CONFIG: (input file, max tools, output basename, output dir)
# ======================================================

TARGETS = [
    {
        "name": "prophage",
        "input": PROPHAGE_MATRICES / "prophage_consensus_regions.tsv",
        "max_tools": 3,
        "out_dir": PROPHAGE_MATRICES,
    },
    {
        "name": "gi",
        "input": GI_MATRICES / "gi_consensus_regions.tsv",
        "max_tools": 3,
        "out_dir": GI_MATRICES,
    },
    {
        "name": "plasmid",
        "input": PLASMID_MATRICES / "plasmid_consensus.tsv",
        "max_tools": 2,
        "out_dir": PLASMID_MATRICES,
    },
]

for t in TARGETS:
    df = pd.read_csv(t["input"], sep="\t")

    high_conf = df[df["N_Tools_Support"] == t["max_tools"]].copy()
    supplementary = df[df["N_Tools_Support"] < t["max_tools"]].copy()

    out_high = t["out_dir"] / f"{t['name']}_high_confidence.tsv"
    out_supp = t["out_dir"] / f"{t['name']}_supplementary.tsv"

    high_conf.to_csv(out_high, sep="\t", index=False)
    supplementary.to_csv(out_supp, sep="\t", index=False)

    print(f"--- {t['name']} (max_tools={t['max_tools']}) ---")
    print(f"Total regions: {len(df)}")
    print(f"High-confidence (all {t['max_tools']} tools agree): {len(high_conf)} "
          f"({100*len(high_conf)/len(df):.1f}%)")
    print(f"Supplementary (partial support): {len(supplementary)} "
          f"({100*len(supplementary)/len(df):.1f}%)")
    print(f"Genomes represented in high-confidence set: {high_conf['Genome'].nunique()}")
    print(f"Saved: {out_high}")
    print(f"Saved: {out_supp}")
    print()
