#!/usr/bin/env python3
from pathlib import Path
from normalize_genome_names import normalize_genome_name

PROJECT_ROOT = Path(__file__).resolve().parents[2]

ANNOT_DIR = PROJECT_ROOT / "annotations"

EXTENSIONS = ["faa", "ffn", "fna", "fsa", "gff", "tbl", "sqn"]

def rename_all(dirpath: Path, ext: str):
    for file in dirpath.glob(f"*.{ext}"):
        original = file.stem
        normalized = normalize_genome_name(original)
        new_path = dirpath / f"{normalized}.{ext}"

        if file != new_path:
            print(f"{file.name} → {new_path.name}")
            file.rename(new_path)

for ext in EXTENSIONS:
    subdir = ANNOT_DIR / ext
    if subdir.exists():
        rename_all(subdir, ext)