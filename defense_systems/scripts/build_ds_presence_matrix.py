#!/usr/bin/env python3
"""
build_ds_presence_matrix.py -- genome x system PRESENCE matrix (0/1) for the
25 both-tool-confirmed defense systems, counting each system ONCE per genome.

Why: ds_consensus_gene_counts.tsv holds gene counts (DefenseFinder genes_count
summed, vs PADLOC hit rows, max of the two). Summed by mechanism category in
the MGE x DS correlation, that metric conflates number of systems with number
of gene copies and with the number of subunits each system has. Here a
system is "present" in a genome if DefenseFinder reported at least one
instance of it OR PADLOC reported at least one instance of it. Both tools only
report a system when it satisfies their model's gene-quorum rules and group
genes by genomic proximity, so a lone, distant copy of one gene does not
count as a system, and a genome with several copies/instances of the same
system still counts it once.

Also saves instance counts (DefenseFinder instances, PADLOC distinct
system.number) in a long file, for a possible sensitivity analysis.

Outputs:
  defense_systems/matrices/ds_consensus_presence.tsv        (Genome x system, 0/1)
  defense_systems/matrices/ds_consensus_instances_long.tsv  (per genome/system:
      df_instances, padloc_instances, present)

Usage
-----
python defense_systems/scripts/build_ds_presence_matrix.py
"""
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from build_ds_consensus_matrix import CROSSWALK, CERTAIN_SYSTEMS  # noqa: E402

PROJECT_ROOT = HERE.parents[1]
sys.path.insert(0, str(PROJECT_ROOT))
from genome_ids import canon  # noqa: E402  (unifies AB22010_2 / AB22010-2 etc.)

DF_DIR = PROJECT_ROOT / "defense_systems" / "defensefinder" / "raw_outputs"
PADLOC_DIR = PROJECT_ROOT / "defense_systems" / "padloc" / "results"
OUT_DIR = PROJECT_ROOT / "defense_systems" / "matrices"
OUT_WIDE = OUT_DIR / "ds_consensus_presence.tsv"
OUT_LONG = OUT_DIR / "ds_consensus_instances_long.tsv"


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
    systems = sorted(CERTAIN_SYSTEMS)

    df_inst = defaultdict(int)       # (genome, system) -> n DefenseFinder instances
    pl_inst = defaultdict(set)       # (genome, system) -> {PADLOC system.number}
    genomes = set()

    for f in sorted(DF_DIR.glob("*/*_defense_finder_systems.tsv")):
        genome = canon(f.name.replace("_defense_finder_systems.tsv", ""))
        genomes.add(genome)
        t = pd.read_csv(f, sep="\t")
        for _, r in t.iterrows():
            cs = df_map.get((norm(r.get("type")), norm(r.get("subtype"))))
            if cs is not None:
                df_inst[(genome, cs)] += 1

    for f in sorted(PADLOC_DIR.glob("*/*_padloc.csv")):
        genome = canon(f.name.replace("_padloc.csv", ""))
        genomes.add(genome)
        t = pd.read_csv(f)
        for _, r in t.iterrows():
            cs = pl_map.get(str(r.get("system")))
            if cs is not None:
                pl_inst[(genome, cs)].add(r.get("system.number"))

    genomes = sorted(genomes)
    print(f"[INFO] {len(genomes)} genomes with DefenseFinder and/or PADLOC output")

    rows = []
    for g in genomes:
        for s in systems:
            d = df_inst.get((g, s), 0)
            p = len(pl_inst.get((g, s), ()))
            rows.append({"Genome": g, "system": s, "df_instances": d,
                         "padloc_instances": p, "present": int(d > 0 or p > 0)})
    long_df = pd.DataFrame(rows)
    wide = long_df.pivot(index="Genome", columns="system", values="present").reindex(columns=systems)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    wide.to_csv(OUT_WIDE, sep="\t")
    long_df.to_csv(OUT_LONG, sep="\t", index=False)
    print(f"[DONE] Saved: {OUT_WIDE} ({wide.shape[0]} genomes x {wide.shape[1]} systems)")
    print(f"[DONE] Saved: {OUT_LONG}")

    multi = long_df[(long_df["df_instances"] > 1) | (long_df["padloc_instances"] > 1)]
    print(f"\n[INFO] (genome, system) pairs with >1 instance in either tool "
          f"(now counted once): {len(multi)}")
    print(multi.groupby("system").size().sort_values(ascending=False).to_string())
    print("\n[INFO] Genomes carrying each system (presence):")
    print(wide.sum().sort_values(ascending=False).to_string())


if __name__ == "__main__":
    main()
