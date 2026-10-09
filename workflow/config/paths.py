#!/usr/bin/env python3
"""
workflow/config/paths.py

Central definition of project directories used by the Snakemake workflow.

This file contains ONLY filesystem locations.
No metadata, no genome lists, no analysis-specific constants.
"""

from pathlib import Path

# ======================================
# PROJECT ROOT
# ======================================

# xtt_population_genomics/
PROJECT_ROOT = Path(__file__).resolve().parents[2]

# ======================================
# WORKFLOW_ROOT
# ======================================

WORKFLOW_ROOT = PROJECT_ROOT / "workflow"

WORKFLOW_CONFIG = WORKFLOW_ROOT / "config"
WORKFLOW_RULES = WORKFLOW_ROOT / "rules"
WORKFLOW_ENVS = WORKFLOW_ROOT / "envs"

# ======================================
# GENERAL
# ======================================

RESULTS = PROJECT_ROOT / "results"
LOGS = PROJECT_ROOT / "logs"
DATABASES = PROJECT_ROOT / "databases"
THIRD_PARTY = PROJECT_ROOT / "third_party"

# ======================================
# MAIN PROJECT
# ======================================

XTT_PROJECT = PROJECT_ROOT  # alias legado; camada xtt_project/ removida


# ======================================
# METADATA
# ======================================

METADATA = XTT_PROJECT / "metadata"
GENOME_CLADES = METADATA / "genome_clades_master.tsv"

# ======================================
# RAW DATA
# ======================================

RAW_GENOMES = XTT_PROJECT / "raw_genomes"
ANNOTATIONS = XTT_PROJECT / "annotations"

GFF = ANNOTATIONS / "gff"
FAA = ANNOTATIONS / "faa"
FFN = ANNOTATIONS / "ffn"
FNA = ANNOTATIONS / "fna"
GBK = ANNOTATIONS / "gbk"

# ======================================
# WORKFLOW MODULES
# ======================================

PHYLOGENY = XTT_PROJECT / "phylogeny"
MOBILOME = XTT_PROJECT / "mobilome"
DEFENSE_SYSTEMS = XTT_PROJECT / "defense_systems"
INTEGRATION = XTT_PROJECT / "integration"
GENOMIC_ISLANDS = MOBILOME / "genomic_islands"
QC = XTT_PROJECT / "qc"

# ======================================
# PHYLOGENY
# ======================================

PHYLO_RAW = PHYLOGENY / "raw_outputs"
PHYLO_MATRICES = PHYLOGENY / "matrices"
PHYLO_STATS = PHYLOGENY / "stats"
PHYLO_PLOTS = PHYLOGENY / "plots"
PHYLO_LOGS = PHYLOGENY / "logs"
PHYLO_SCRIPTS = PHYLOGENY / "scripts"
PHYLO_ITOL = PHYLOGENY / "itol"

FINAL_TREES = PHYLOGENY / "final_trees"
FINAL_ALIGNMENTS = PHYLOGENY / "final_alignments"
RECOMBINATION = PHYLOGENY / "recombination_intermediates"

# ======================================
# MOBILOME
# ======================================

# ------------------------------
# INSERTION SEQUENCES
# ------------------------------

INSERTION_SEQUENCES = MOBILOME / "insertion_sequences"

IS_RAW = INSERTION_SEQUENCES / "raw_outputs"
IS_MATRICES = INSERTION_SEQUENCES / "matrices"
IS_STATS = INSERTION_SEQUENCES / "stats"
IS_PLOTS = INSERTION_SEQUENCES / "plots"
IS_LOGS = INSERTION_SEQUENCES / "logs"
IS_ITOL = INSERTION_SEQUENCES / "itol"
IS_SCRIPTS = INSERTION_SEQUENCES / "scripts"

# ------------------------------
# ISESCAN
# ------------------------------

ISESCAN = INSERTION_SEQUENCES / "isescan"

ISESCAN_DATABASE = ISESCAN / "database"
ISESCAN_RAW = ISESCAN / "raw_outputs"
ISESCAN_MATRICES = ISESCAN / "matrices"
ISESCAN_STATS = ISESCAN / "stats"
ISESCAN_PLOTS = ISESCAN / "plots"
ISESCAN_LOGS = ISESCAN / "logs"
ISESCAN_ITOL = ISESCAN / "itol"
ISESCAN_SCRIPTS = ISESCAN / "scripts"

# ------------------------------
# PROPHAGES
# ------------------------------

PROPHAGES = MOBILOME / "prophages"

############################################################
# PLASMIDS
############################################################

PLASMIDS = MOBILOME / "plasmids"

PLASMID_MATRICES = PLASMIDS / "matrices"
PLASMID_STATS = PLASMIDS / "stats"
PLASMID_PLOTS = PLASMIDS / "plots"
PLASMID_LOGS = PLASMIDS / "logs"
PLASMID_SCRIPTS = PLASMIDS / "scripts"

MOB_SUITE = PLASMIDS / "mob_suite"
MOB_SUITE_RAW = MOB_SUITE / "raw_outputs"
MOB_SUITE_MATRICES = MOB_SUITE / "matrices"
MOB_SUITE_SCRIPTS = MOB_SUITE / "scripts"

GENOMAD_PLASMIDS = PLASMIDS / "genomad"
GENOMAD_PLASMIDS_MATRICES = GENOMAD_PLASMIDS / "matrices"
GENOMAD_PLASMIDS_SCRIPTS = GENOMAD_PLASMIDS / "scripts"

PROPHAGE_RAW = PROPHAGES / "raw_outputs"
PROPHAGE_MATRICES = PROPHAGES / "matrices"
PROPHAGE_STATS = PROPHAGES / "stats"
PROPHAGE_PLOTS = PROPHAGES / "plots"
PROPHAGE_LOGS = PROPHAGES / "logs"
PROPHAGE_ITOL = PROPHAGES / "itol"
PROPHAGE_SCRIPTS = PROPHAGES / "scripts"

# ------------------------------
# PHASTEST
# ------------------------------

PHASTEST = PROPHAGES / "phastest"

PHASTEST_RAW = PHASTEST / "raw_outputs"
PHASTEST_JOBS = PHASTEST / "jobs"
PHASTEST_DOWNLOADS = PHASTEST / "downloads"
PHASTEST_SUMMARIES = PHASTEST / "summaries"
PHASTEST_STATS = PHASTEST / "stats"
PHASTEST_PLOTS = PHASTEST / "plots"
PHASTEST_LOGS = PHASTEST / "logs"
PHASTEST_SCRIPTS = PHASTEST / "scripts"

# ------------------------------
# VIRSORTER2
# ------------------------------

VIRSORTER2 = PROPHAGES / "virsorter2"

VIRSORTER2_DB = VIRSORTER2 / "db"
VIRSORTER2_RAW = VIRSORTER2 / "raw_outputs"
VIRSORTER2_MATRICES = VIRSORTER2 / "matrices"
VIRSORTER2_STATS = VIRSORTER2 / "stats"
VIRSORTER2_PLOTS = VIRSORTER2 / "plots"
VIRSORTER2_LOGS = VIRSORTER2 / "logs"
VIRSORTER2_SCRIPTS = VIRSORTER2 / "scripts"

# ------------------------------
# CHECKV
# ------------------------------

CHECKV = PROPHAGES / "checkv"

CHECKV_RAW = CHECKV / "raw_outputs"
CHECKV_STATS = CHECKV / "stats"
CHECKV_PLOTS = CHECKV / "plots"
CHECKV_LOGS = CHECKV / "logs"
CHECKV_SCRIPTS = CHECKV / "scripts"

# ------------------------------
# GENOMAD
# ------------------------------

GENOMAD = PROPHAGES / "genomad"

GENOMAD_DATABASE = GENOMAD / "database"
GENOMAD_RAW = GENOMAD / "raw_outputs"
GENOMAD_MATRICES = GENOMAD / "matrices"
GENOMAD_STATS = GENOMAD / "stats"
GENOMAD_PLOTS = GENOMAD / "plots"
GENOMAD_LOGS = GENOMAD / "logs"
GENOMAD_SCRIPTS = GENOMAD / "scripts"

# ======================================
# DEFENSE SYSTEMS
# ======================================

# ------------------------------
# DEFENSEFINDER
# ------------------------------

DEFENSEFINDER = DEFENSE_SYSTEMS / "defensefinder"

DEFENSEFINDER_RAW = DEFENSEFINDER / "raw_outputs"
DEFENSEFINDER_MATRICES = DEFENSEFINDER / "matrices"
DEFENSEFINDER_STATS = DEFENSEFINDER / "stats"
DEFENSEFINDER_PLOTS = DEFENSEFINDER / "plots"
DEFENSEFINDER_LOGS = DEFENSEFINDER / "logs"
DEFENSEFINDER_ITOL = DEFENSEFINDER / "itol"
DEFENSEFINDER_SCRIPTS = DEFENSEFINDER / "scripts"

# ------------------------------
# PADLOC
# ------------------------------

PADLOC = DEFENSE_SYSTEMS / "padloc"

PADLOC_DATABASE = PADLOC / "database"
PADLOC_RAW = PADLOC / "raw_outputs"
PADLOC_MATRICES = PADLOC / "matrices"
PADLOC_STATS = PADLOC / "stats"
PADLOC_PLOTS = PADLOC / "plots"
PADLOC_LOGS = PADLOC / "logs"
PADLOC_ITOL = PADLOC / "itol"
PADLOC_SCRIPTS = PADLOC / "scripts"

# ======================================
# INTEGRATION
# ======================================

INTEGRATION_RAW = INTEGRATION / "raw_outputs"
INTEGRATION_MATRICES = INTEGRATION / "matrices"
INTEGRATION_STATS = INTEGRATION / "stats"
INTEGRATION_PLOTS = INTEGRATION / "plots"
INTEGRATION_LOGS = INTEGRATION / "logs"
INTEGRATION_SCRIPTS = INTEGRATION / "scripts"

# ------------------------------
# SSN (Sequence Similarity Network, MGE nucleotide-level, style of
# Uceda-Campos et al. 2022 Microorganisms Fig. -- prophages, genomic
# islands, insertion sequences, plasmids as nodes; BLASTN nucleotide
# similarity >=50% identity / >=80% coverage as edges)
# ------------------------------

SSN = INTEGRATION / "ssn"

SSN_FASTA = SSN / "fasta"
SSN_BLAST = SSN / "blast"
SSN_MATRICES = SSN / "matrices"
SSN_PLOTS = SSN / "plots"
SSN_LOGS = SSN / "logs"
SSN_SCRIPTS = SSN / "scripts"

# ======================================
# GENOMIC ISLANDS
# ======================================

GI_RAW = GENOMIC_ISLANDS / "raw_outputs"
GI_MATRICES = GENOMIC_ISLANDS / "matrices"
GI_STATS = GENOMIC_ISLANDS / "stats"
GI_PLOTS = GENOMIC_ISLANDS / "plots"
GI_LOGS = GENOMIC_ISLANDS / "logs"
GI_ITOL = GENOMIC_ISLANDS / "itol"
GI_SCRIPTS = GENOMIC_ISLANDS / "scripts"

# ------------------------------
# ALIENHUNTER
# ------------------------------

ALIENHUNTER = GENOMIC_ISLANDS / "alienhunter"

ALIENHUNTER_RAW = ALIENHUNTER / "raw_outputs"
ALIENHUNTER_STATS = ALIENHUNTER / "stats"
ALIENHUNTER_PLOTS = ALIENHUNTER / "plots"
ALIENHUNTER_LOGS = ALIENHUNTER / "logs"
ALIENHUNTER_SCRIPTS = ALIENHUNTER / "scripts"

# ------------------------------
# GIPSY2
# ------------------------------

GIPSY2 = GENOMIC_ISLANDS / "gipsy2"

GIPSY2_RAW = GIPSY2 / "raw_outputs"
GIPSY2_STATS = GIPSY2 / "stats"
GIPSY2_PLOTS = GIPSY2 / "plots"
GIPSY2_LOGS = GIPSY2 / "logs"
GIPSY2_SCRIPTS = GIPSY2 / "scripts"

# ------------------------------
# ISLANDPATH
# ------------------------------

ISLANDPATH = GENOMIC_ISLANDS / "islandpath"

ISLANDPATH_RAW = ISLANDPATH / "raw_outputs"
ISLANDPATH_STATS = ISLANDPATH / "stats"
ISLANDPATH_PLOTS = ISLANDPATH / "plots"
ISLANDPATH_LOGS = ISLANDPATH / "logs"
ISLANDPATH_SCRIPTS = ISLANDPATH / "scripts"
ISLANDPATH_MATRICES = ISLANDPATH / "matrices"

# ------------------------------
# PANISLE (isolated code, inspired by GIPSy2 -- see THIRD_PARTY_LICENSES.md;
# results/ produced manually via SLURM array job, not orchestrated here)
# ------------------------------

PANISLE = GENOMIC_ISLANDS / "panisle"

PANISLE_RAW = PANISLE / "results"
PANISLE_MATRICES = PANISLE / "matrices"
PANISLE_STATS = PANISLE / "stats"
PANISLE_LOGS = PANISLE / "logs"
PANISLE_SCRIPTS = PANISLE / "scripts"


