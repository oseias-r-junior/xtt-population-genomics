# *Xanthomonas translucens* pv. *translucens* population genomics pipeline

[![License: MIT](https://img.shields.io/badge/Code-MIT-blue.svg)](LICENSE-CODE)

Analysis workflow for the manuscript (in preparation; working title):
> Velasco, D., Belen, G., Feitosa-Junior, O., Schachterle, J., Friskop, A., Liu, Z., and Baldwin, T. *Exploring the population diversity of* Xanthomonas translucens *pv.* translucens*, causing bacterial leaf streak on barley in North Dakota and other barley-producing states.*

Snakemake workflow for mobilome, defense-system and virulence-gene analyses of *X. translucens* pv. *translucens* (Xtt) population genomes, run on an HPC cluster (SLURM).

> **Status: manuscript in preparation.** This repository shares the **analysis workflow** (code, environments, run order). **Data, results, figures and derived tables are not included**; they will be released together with the associated publication.

## Overview

Pipeline order:

**annotation → pangenome → recombination and phylogeny → MGE detection → MGE similarity network → defense systems → integration**

| Stage | Folder | Tools |
| --- | --- | --- |
| Annotation | run outside Snakemake | Prokka |
| Pangenome | `pangenome/` | Roary |
| Recombination and phylogeny | `phylogeny/` | Gubbins, ClonalFrameML, IQ-TREE |
| MGE detection | `mobilome/{prophages,plasmids,genomic_islands,insertion_sequences}/` | PHASTEST, VirSorter2, geNomad (prophages); MOB-suite, geNomad (plasmids); AlienHunter, IslandPath-DIMOB (genomic islands); ISEScan (IS) |
| MGE network | `integration/ssn/` | BLASTN all-vs-all, networkx |
| Defense systems | `defense_systems/` | DefenseFinder, PADLOC (consensus of both tools) |
| Integration | `integration/clade_correlation/` | MGE per clade, MGE × defense correlation, physical overlaps |

MGE calls use a consensus rule (all tools of a class must agree). The genomic-island step also uses an in-house pangenome-based classifier that is **not distributed** here; see [THIRD_PARTY_LICENSES.md](THIRD_PARTY_LICENSES.md).

## Repository structure

```
.
├── workflow/            # Snakefile, rules/*.smk, config/, envs/
├── pangenome/           # Roary helper scripts
├── phylogeny/           # recombination / tree scripts and SLURM jobs
├── mobilome/            # per-class MGE detection scripts and jobs
├── defense_systems/     # DefenseFinder / PADLOC consensus scripts
├── integration/         # MGE network and clade-level integration scripts
├── scripts/             # metadata and maintenance helpers
├── data/README.md       # expected inputs (data not included)
├── metadata/README.md   # expected metadata (not included)
├── palette.py           # shared figure colors
└── THIRD_PARTY_LICENSES.md
```

## Setup

```bash
git clone https://github.com/oseias-r-junior/xtt-population-genomics.git
cd xtt-population-genomics

mamba create -n xtt_snakemake -c conda-forge -c bioconda snakemake=9
conda activate xtt_snakemake

# dry-run of the manuscript targets (needs the input data, see data/README.md)
snakemake -s workflow/Snakefile -n paper --cores 1
```

Job scripts assume they are submitted from the project root (or set `PROJECT_DIR`).

Rule environments are defined in `workflow/envs/` and built with `--use-conda`.
Steps that are heavy or run on SLURM (annotation, pangenome, recombination, BLAST all-vs-all, defense-system tools) are launched outside Snakemake from the `jobs/` or `slurm/` folders and enter the workflow as plain inputs.

## Data

Input data are not distributed here. See [data/README.md](data/README.md) and [metadata/README.md](metadata/README.md).

## Third-party tools

See [THIRD_PARTY_LICENSES.md](THIRD_PARTY_LICENSES.md).

## Citation

Citation information will be added upon publication.

## License

Code: [MIT](LICENSE-CODE).
