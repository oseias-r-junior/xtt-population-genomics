#!/usr/bin/env python3
"""
workflow.config

Central configuration package for the Xtt Population Genomics Pipeline.

This package provides:

    paths.py
        Filesystem locations.

    genomes.py
        Genome collections and sample lists.

    resources.py
        Snakemake resource definitions.
"""

from .paths import *
from .genomes import *
from .resources import *