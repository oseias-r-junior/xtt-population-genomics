#!/usr/bin/env python3
from pathlib import Path
from normalize_genome_names import normalize_genome_name

GBK_DIR = Path(__file__).resolve().parents[2] / "annotations" / "gbk"

for gbk in GBK_DIR.glob("*.gbk"):
    original = gbk.stem
    normalized = normalize_genome_name(original)

    new_path = GBK_DIR / f"{normalized}.gbk"

    if gbk != new_path:
        print(f"{gbk.name} → {new_path.name}")
        gbk.rename(new_path)
