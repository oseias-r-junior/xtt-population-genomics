#!/usr/bin/env python3
from pathlib import Path
from normalize_genome_names import normalize_genome_name

PROJECT_ROOT = Path(__file__).resolve().parents[2]

FNA_DIR = PROJECT_ROOT / "annotations" / "fna"
FSA_DIR = PROJECT_ROOT / "annotations" / "fsa"

def rename_all(dirpath: Path, ext: str):
    for file in dirpath.glob(f"*.{ext}"):
        original = file.stem
        normalized = normalize_genome_name(original)
        new_path = dirpath / f"{normalized}.{ext}"

        if file != new_path:
            print(f"{file.name} → {new_path.name}")
            file.rename(new_path)

rename_all(FNA_DIR, "fna")
rename_all(FSA_DIR, "fsa")