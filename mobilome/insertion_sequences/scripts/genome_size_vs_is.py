import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from pathlib import Path
import sys

SCRIPT_DIR = Path(__file__).resolve().parent

PROJECT_ROOT = SCRIPT_DIR.parents[2]

sys.path.append(
    str(PROJECT_ROOT / "scripts")
)

from config import *




from scipy.stats import pearsonr, spearmanr
from pathlib import Path

# ======================================
# INPUT
# ======================================

IS_FILE = (
    IS_MATRICES
    / "is_family_matrix.tsv"
)

META_FILE = Path(
    GENOME_METADATA
).expanduser()

OUTDIR = IS_PLOTS

OUTDIR.mkdir(
    parents=True,
    exist_ok=True
)
# ======================================
# LOAD
# ======================================

is_df = pd.read_csv(IS_FILE, sep="\t")

meta_df = pd.read_csv(META_FILE, sep="\t")

# ======================================
# MERGE
# ======================================

df = pd.merge(
    is_df,
    meta_df,
    on="Genome",
    how="inner"
)

print(df.head())

# ======================================
# CORRELATIONS
# ======================================

pearson_r, pearson_p = pearsonr(
    df["Genome_Length"],
    df["Total_IS"]
)

spearman_r, spearman_p = spearmanr(
    df["Genome_Length"],
    df["Total_IS"]
)

print("\n===== Correlations =====")

print(f"Pearson r = {pearson_r:.4f}")
print(f"Pearson p = {pearson_p:.6f}")

print(f"Spearman rho = {spearman_r:.4f}")
print(f"Spearman p = {spearman_p:.6f}")

# ======================================
# PLOT
# ======================================

plt.figure(figsize=(10, 8))

sns.scatterplot(
    data=df,
    x="Genome_Length",
    y="Total_IS",
    hue="Clade",
    s=120
)

sns.regplot(
    data=df,
    x="Genome_Length",
    y="Total_IS",
    scatter=False,
    color="black"
)

plt.xlabel("Genome size (bp)")
plt.ylabel("Total IS count")

plt.title(
    "Genome size vs insertion sequence abundance"
)

plt.text(
    0.05,
    0.95,
    (
        f"Pearson r = {pearson_r:.3f}\n"
        f"Spearman rho = {spearman_r:.3f}"
    ),
    transform=plt.gca().transAxes,
    verticalalignment="top",
    fontsize=12
)

plt.tight_layout()

outfile = OUTDIR / "genome_size_vs_is.png"

plt.savefig(
    outfile,
    dpi=300
)

print(f"\nSaved: {outfile}")