from pathlib import Path
import pandas as pd

try:
    # When this file is loaded via Snakemake's `include:` directive,
    # __file__ is NOT reliable (Snakemake's include mechanism does not
    # guarantee it points to this file's own location -- it can resolve
    # to Snakemake's own internal package instead). Snakemake exposes
    # `workflow.basedir`, the absolute path of the directory containing
    # the top-level Snakefile (workflow/), which IS reliable here.
    PROJECT_ROOT = Path(workflow.basedir).parent
except NameError:
    # Plain Python import (e.g. `from genomes import GENOMES` in the
    # various mobilome/*/scripts/*.py analysis scripts, via
    # sys.path.append(.../workflow/config)) -- here __file__ IS
    # reliable, anchor on it instead.
    SCRIPT_DIR = Path(__file__).resolve().parent
    PROJECT_ROOT = SCRIPT_DIR.parents[1]

GENOME_CLADES = (
    PROJECT_ROOT
    / "metadata"
    / "genome_clades_master.tsv"
)

GENOMES = pd.read_csv(
    GENOME_CLADES,
    sep="\t"
)["Genome"].tolist()
