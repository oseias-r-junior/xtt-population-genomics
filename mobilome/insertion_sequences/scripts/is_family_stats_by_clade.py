import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import kruskal
from pathlib import Path
import numpy as np
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[4]

sys.path.insert(0, str(PROJECT_ROOT / "scripts" / "metadata"))
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))


SCRIPT_DIR = Path(__file__).resolve().parent

# PROJECT_ROOT = SCRIPT_DIR.parents[2]

# sys.path.append(
#     str(PROJECT_ROOT / "scripts")
# )

from config import *
from normalize_genome_names import normalize_genome_name




# =========================
# INPUT FILES
# =========================

MATRIX_FILE = (
    IS_MATRICES
    / "is_family_matrix.tsv"
)

CLADE_FILE = GENOME_CLADES

OUTDIR = (
    IS_STATS
    / "is_clade_stats"
)

OUTDIR.mkdir(
    parents=True,
    exist_ok=True
)

# =========================
# LOAD DATA
# =========================

df = pd.read_csv(MATRIX_FILE, sep="\t")
clades = pd.read_csv(CLADE_FILE, sep="\t")

merged = df.merge(clades, on="Genome")

print(merged.head())

# =========================
# IS FAMILIES
# =========================

is_families = [
    "IS110",
    "IS1595",
    "IS3",
    "IS4",
    "IS481",
    "IS5",
    "IS605"
]

# =========================
# STATS
# =========================

results = []

for family in is_families:

    groups = []

    for clade in sorted(merged["Clade"].unique()):

        vals = merged.loc[
            merged["Clade"] == clade,
            family
        ]

        groups.append(vals)

    # Skip invariant features

    all_values = np.concatenate(groups)

    if len(set(all_values)) <= 1:

        print(f"Skipping {family}: invariant values")

        continue


    H, p = kruskal(*groups)

    results.append({
        "Family": family,
        "H_statistic": H,
        "p_value": p
    })

stats_df = pd.DataFrame(results)

stats_df.to_csv(
    OUTDIR / "kruskal_results.tsv",
    sep="\t",
    index=False
)

print("\n===== Kruskal-Wallis results =====")
print(stats_df)

# =========================
# LONG FORMAT
# =========================

long_df = merged.melt(
    id_vars=["Genome", "Clade"],
    value_vars=is_families,
    var_name="IS_family",
    value_name="Count"
)

# =========================
# BOXPLOTS
# =========================

plt.figure(figsize=(12, 6))

sns.boxplot(
    data=long_df,
    x="IS_family",
    y="Count",
    hue="Clade"
)

plt.title("IS family abundance across ANI clades")
plt.tight_layout()

plt.savefig(
    OUTDIR / "is_family_boxplots_by_clade.png",
    dpi=300
)

print("\nSaved boxplots.")

# =========================
# SUMMARY TABLE
# =========================

summary = merged.groupby("Clade")[is_families].mean()

summary.to_csv(
    OUTDIR / "mean_is_by_clade.tsv",
    sep="\t"
)

print("\nSaved summary table.")
