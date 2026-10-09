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
from normalize_genome_names import normalize_genome_name


# ======================================
# INPUT
# ======================================

INPUT_MATRIX = (
    IS_MATRICES
    / "is_family_matrix.tsv"
)

# ======================================
# OUTPUT
# ======================================

IS_PLOTS.mkdir(
    parents=True,
    exist_ok=True
)

OUTPUT_PLOT = (
    IS_PLOTS
    / "is_family_clustered_heatmap.png"
)

# ======================================
# LOAD MATRIX
# ======================================

df = pd.read_csv(
    INPUT_MATRIX,
    sep="\t"
)

# Set genome names as index
df = df.set_index("Genome")

# Remove Total_IS for clustering
features = df.drop(columns=["Total_IS"])

# Generate clustered heatmap
g = sns.clustermap(
    features,
    cmap="viridis",
    figsize=(10, 20),
    metric="euclidean",
    method="average",
    standard_scale=1
)

# Save
g.savefig(
    OUTPUT_PLOT,
    dpi=300
)

print("Clustered heatmap generated.")