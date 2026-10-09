#!/usr/bin/env python3
"""
list_ds_systems.py -- inventory of unique defense-system labels used by
DefenseFinder and PADLOC across all genomes in the project, as a first step
toward building a manual name crosswalk between the two tools (no
crosswalk exists yet, confirmed with the user).

DefenseFinder: reads every */*_defense_finder_systems.tsv under
defense_systems/defensefinder/raw_outputs/, collects unique (type, subtype)
pairs, counts how many genomes and how many system instances carry each,
and the range of genes_count observed for that (type, subtype).

PADLOC: reads every */*_padloc.csv under defense_systems/padloc/results/,
collects unique `system` values, counts genomes and instances (system.number
is the per-genome instance id), and range of gene-count per instance
(rows per system.number).

Output (printed, and also written to defense_systems/matrices/ for the
manual crosswalk step):
    defense_systems/matrices/defensefinder_systems_inventory.tsv
    defense_systems/matrices/padloc_systems_inventory.tsv

Usage
-----
python list_ds_systems.py
"""
from collections import defaultdict
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
DF_DIR = PROJECT_ROOT / "defense_systems" / "defensefinder" / "raw_outputs"
PADLOC_DIR = PROJECT_ROOT / "defense_systems" / "padloc" / "results"
OUT_DIR = PROJECT_ROOT / "defense_systems" / "matrices"


def inventory_defensefinder():
    files = sorted(DF_DIR.glob("*/*_defense_finder_systems.tsv"))
    print(f"[INFO] DefenseFinder systems.tsv files found: {len(files)}")
    if not files:
        print(f"[FATAL] No files matched under {DF_DIR}")
        return None

    stats = defaultdict(lambda: {"genomes": set(), "instances": 0, "genes_counts": []})
    for f in files:
        genome = f.name.replace("_defense_finder_systems.tsv", "")
        try:
            df = pd.read_csv(f, sep="\t")
        except Exception as e:
            print(f"[WARN] Could not read {f}: {e}")
            continue
        for _, row in df.iterrows():
            key = (str(row.get("type", "")), str(row.get("subtype", "")))
            stats[key]["genomes"].add(genome)
            stats[key]["instances"] += 1
            try:
                stats[key]["genes_counts"].append(int(row.get("genes_count", 0)))
            except (ValueError, TypeError):
                pass

    rows = []
    for (type_, subtype), d in stats.items():
        gc = d["genes_counts"]
        rows.append({
            "type": type_,
            "subtype": subtype,
            "n_genomes": len(d["genomes"]),
            "n_instances": d["instances"],
            "genes_count_min": min(gc) if gc else "",
            "genes_count_max": max(gc) if gc else "",
        })
    out = pd.DataFrame(rows).sort_values(["type", "subtype"]).reset_index(drop=True)
    print(f"[INFO] DefenseFinder: {len(out)} unique (type, subtype) combinations across "
          f"{len(files)} genomes")
    return out


def inventory_padloc():
    files = sorted(PADLOC_DIR.glob("*/*_padloc.csv"))
    print(f"[INFO] PADLOC *_padloc.csv files found: {len(files)}")
    if not files:
        print(f"[FATAL] No files matched under {PADLOC_DIR}")
        return None

    stats = defaultdict(lambda: {"genomes": set(), "instance_keys": set(), "gene_counts": defaultdict(int)})
    for f in files:
        genome = f.name.replace("_padloc.csv", "")
        try:
            df = pd.read_csv(f)
        except Exception as e:
            print(f"[WARN] Could not read {f}: {e}")
            continue
        if "system" not in df.columns or "system.number" not in df.columns:
            print(f"[WARN] {f} missing expected columns, skipping")
            continue
        for _, row in df.iterrows():
            sys_name = str(row["system"])
            inst_key = (genome, row["system.number"])
            stats[sys_name]["genomes"].add(genome)
            stats[sys_name]["instance_keys"].add(inst_key)
            stats[sys_name]["gene_counts"][inst_key] += 1

    rows = []
    for sys_name, d in stats.items():
        gc = list(d["gene_counts"].values())
        rows.append({
            "system": sys_name,
            "n_genomes": len(d["genomes"]),
            "n_instances": len(d["instance_keys"]),
            "genes_count_min": min(gc) if gc else "",
            "genes_count_max": max(gc) if gc else "",
        })
    out = pd.DataFrame(rows).sort_values("system").reset_index(drop=True)
    print(f"[INFO] PADLOC: {len(out)} unique system labels across {len(files)} genomes")
    return out


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    df_inv = inventory_defensefinder()
    if df_inv is not None:
        out_path = OUT_DIR / "defensefinder_systems_inventory.tsv"
        df_inv.to_csv(out_path, sep="\t", index=False)
        print(f"[DONE] Saved: {out_path}")
        print(df_inv.to_string(index=False))

    print()
    padloc_inv = inventory_padloc()
    if padloc_inv is not None:
        out_path = OUT_DIR / "padloc_systems_inventory.tsv"
        padloc_inv.to_csv(out_path, sep="\t", index=False)
        print(f"[DONE] Saved: {out_path}")
        print(padloc_inv.to_string(index=False))


if __name__ == "__main__":
    main()
