#!/usr/bin/env python3
"""
list_ds_gene_names.py -- diagnostic for the gene-by-gene defense-system
heatmap: for each of the 25 both-tool-confirmed systems, lists which GENE
names each tool reports (DefenseFinder: `name_of_profiles_in_sys` from
systems.tsv; PADLOC: `protein.name` from padloc.csv) and in how many genomes
each gene name appears, plus max copies per genome.

Reuses CROSSWALK / CERTAIN_SYSTEMS from build_ds_consensus_matrix.py.

Output:
  defense_systems/matrices/ds_gene_name_inventory.tsv
  columns: consensus_system, tool, gene_name, n_genomes, max_copies_in_a_genome

Usage (from anywhere)
---------------------
python defense_systems/scripts/list_ds_gene_names.py
"""
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from build_ds_consensus_matrix import CROSSWALK, CERTAIN_SYSTEMS  # noqa: E402

PROJECT_ROOT = HERE.parents[1]
DF_DIR = PROJECT_ROOT / "defense_systems" / "defensefinder" / "raw_outputs"
PADLOC_DIR = PROJECT_ROOT / "defense_systems" / "padloc" / "results"
OUT_FILE = PROJECT_ROOT / "defense_systems" / "matrices" / "ds_gene_name_inventory.tsv"


def norm(x):
    return None if pd.isna(x) else str(x)


def main():
    df_map, pl_map = {}, {}
    for cs, dt, ds, ps in CROSSWALK:
        if cs not in CERTAIN_SYSTEMS:
            continue
        if dt is not None:
            df_map[(dt, ds)] = cs
        if ps is not None:
            pl_map[ps] = cs

    # (tool, system, gene) -> {genome: copies}
    copies = defaultdict(lambda: defaultdict(int))

    for f in sorted(DF_DIR.glob("*/*_defense_finder_systems.tsv")):
        genome = f.name.replace("_defense_finder_systems.tsv", "")
        t = pd.read_csv(f, sep="\t")
        for _, r in t.iterrows():
            cs = df_map.get((norm(r.get("type")), norm(r.get("subtype"))))
            if cs is None:
                continue
            genes = str(r.get("name_of_profiles_in_sys", "")).split(",")
            for g in genes:
                g = g.strip()
                if g:
                    copies[("DefenseFinder", cs, g)][genome] += 1

    for f in sorted(PADLOC_DIR.glob("*/*_padloc.csv")):
        genome = f.name.replace("_padloc.csv", "")
        t = pd.read_csv(f)
        for _, r in t.iterrows():
            cs = pl_map.get(str(r.get("system")))
            if cs is None:
                continue
            g = norm(r.get("protein.name")) or norm(r.get("hmm.name")) or "NA"
            copies[("PADLOC", cs, g)][genome] += 1

    rows = [{"consensus_system": cs, "tool": tool, "gene_name": g,
             "n_genomes": len(gm), "max_copies_in_a_genome": max(gm.values())}
            for (tool, cs, g), gm in copies.items()]
    out = pd.DataFrame(rows).sort_values(["consensus_system", "tool", "n_genomes"],
                                         ascending=[True, True, False])
    OUT_FILE.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT_FILE, sep="\t", index=False)
    print(f"[DONE] Saved: {OUT_FILE} ({len(out)} rows)\n")
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()