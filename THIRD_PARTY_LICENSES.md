# Third-party tools and licensing notes

This repository documents an analysis pipeline. Some third-party tools it
depends on are **not distributed in this repository** because their
licenses do not permit redistribution, or because code related to them
carries the same restriction. This file lists those exclusions.

## GIPSy2 / PanISLE genomic island classifier

**Kept out of this repository.**

PanISLE, the genomic island classifier developed during this project, was
slightly inspired by GIPSy2. GIPSy2 is distributed under a license that
allows adaptation of the software but does not allow distribution of the
adapted version. In addition, the GIPSy2 distributions hosted on Zenodo are
not actively supported, so a different approach was needed to work around
that limitation and make the method usable in this pipeline.

PanISLE differs substantially from GIPSy2's original method: it classifies
islands from Roary pangenome output, rather than by comparing the genome of
a virulent strain against that of an avirulent one. Even so, because of its
inspiration in GIPSy2 and the license terms above, we keep the PanISLE
code, together with its development history, out of this pipeline and off
GitHub.

GIPSy2 reference and license:

- Rodrigues, D. L. N., Azevedo, V. A. de C., Soares, S. de C., & Aburjaile, F. F.
  (2024-2025). *GIPSy2 - Genomic Island Prediction Software 2* (Version 2.0.7)
  [Computer software]. Zenodo. https://doi.org/10.5281/zenodo.14969252
- License: Creative Commons Attribution-NonCommercial-NoDerivatives 4.0
  International (CC BY-NC-ND 4.0).

## Other third-party tools (installed, not vendored)

The following tools are installed locally on CCAST (see `third_party/`
and the conda environments under `workflow/envs/`) and are referenced by
this pipeline's Snakemake rules, but their source code is not included in
this repository. Install them from their own distribution channels under
each tool's own license:

- AlienHunter
- IslandPath-DIMOB
- ISEScan
- MMseqs2
- ClonalFrameML
- IQ-TREE