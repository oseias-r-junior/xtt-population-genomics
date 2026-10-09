#!/usr/bin/env python3
from pathlib import Path
from normalize_genome_names import normalize_genome_name

PROJECT_ROOT = Path(__file__).resolve().parents[2]
RAW = PROJECT_ROOT / "mobilome" / "genomic_islands" / "islandpath" / "raw_outputs"

for d in RAW.iterdir():
    if d.is_dir():
        original = d.name
        normalized = normalize_genome_name(original)
        new_path = RAW / normalized

        if d != new_path:
            print(f"{d.name} → {new_path.name}")
            d.rename(new_path)