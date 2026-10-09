import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import sys

SCRIPT_DIR = Path(__file__).resolve().parent

PROJECT_ROOT = SCRIPT_DIR.parents[2]

sys.path.append(
    str(PROJECT_ROOT / "scripts")
)

from config import *
from normalize_genome_names import normalize_genome_name



from scipy.stats import kruskal
from pathlib import Path

# ======================================
# INPUT
# ======================================

IS_FILE = (
    IS_MATRICES
    / "is_family_matrix.tsv"
)

META_FILE = GENOME_CLADES

OUTDIR = (
    IS_STATS
    / "enrichment"
)

OUTDIR.mkdir(
    parents=True,
    exist_ok=True
)

# ======================================
# LOAD
# ======================================

is_df = pd.read_csv(IS_FILE, sep="\t")
meta_df = pd.read_csv(META_FILE, sep="\t")

merged = pd.merge(
    is_df,
    meta_df,
    on="Genome",
    how="inner"
)
print(merged["Clade"].value_counts(dropna=False))

# ======================================
# IS families
# ======================================

is_families = [
    c for c in merged.columns
    if c.startswith("IS")
]

print(is_families)

# ======================================
# KRUSKAL TESTS
# ======================================

results = []

for fam in is_families:

    groups = []

    for clade in sorted(merged["Clade"].unique()):

        vals = merged.loc[
            merged["Clade"] == clade,
            fam
        ]

        groups.append(vals)

    # Remove empty groups

    groups = [
        g for g in groups
        if len(g) > 0
    ]

    # Need at least 2 groups

    if len(groups) < 2:

        print(f"Skipping {fam}: fewer than 2 groups")

        continue

    # Skip invariant values

    all_values = np.concatenate(groups)

    if len(set(all_values)) <= 1:

        print(f"Skipping {fam}: invariant values")

        continue

    stat, p = kruskal(*groups)

    means = (
        merged
        .groupby("Clade")[fam]
        .mean()
        .to_dict()
    )

    results.append({
        "IS_family": fam,
        "Kruskal_H": stat,
        "pvalue": p,
        **means
    })

results_df = pd.DataFrame(results)

results_df = results_df.sort_values("pvalue")

print(results_df)

# ======================================
# SAVE TABLE
# ======================================

outfile = OUTDIR / "is_family_enrichment.tsv"

results_df.to_csv(
    outfile,
    sep="\t",
    index=False
)

print(f"\nSaved: {outfile}")

# ======================================
# HEATMAP OF CLADE MEANS
# ======================================

heatmap_df = (
    results_df
    .set_index("IS_family")[["K0", "K1", "K2"]]
)

plt.figure(figsize=(6, 5))

sns.heatmap(
    heatmap_df,
    annot=True,
    cmap="viridis"
)

plt.title(
    "Mean IS family abundance by ANI clade"
)

plt.tight_layout()

heatmap_out = OUTDIR / "is_family_clade_means_heatmap.png"

plt.savefig(
    heatmap_out,
    dpi=300
)

print(f"Saved: {heatmap_out}")