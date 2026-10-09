#!/usr/bin/env python3
"""
palette.py -- single source of truth for figure colors in this project.

Import from any plotting script (put this file in integration/scripts/ and
add that folder to sys.path, or copy it next to the script):

    from palette import MGE_COLOR, GI_CLASS_COLOR, CLADE_COLOR, CMAP

Design rules
------------
1. MGE category identity is the same color in every figure:
   PPH black, GI purple, IS orange, PLS grey.
2. GI sub-classes (functional classification) get ONE color each (5 classes,
   no pairing and "unknown" is its own pale tone), all outside the warm
   YlOrRd ramp used for IS families and outside black/grey, so a GI node is
   never confused with an IS or PPH/PLS node.
3. Heatmaps of different analyses use different colormaps so figures that sit
   next to each other in the supplement are not mistaken for one another.
4. Clade colors (K0/K1/K2) are used only for clade annotation.
"""

MGE_COLOR = {
    "PPH": "#1a1a1a",
    "GI": "#807dba",
    "IS": "#fd8d3c",
    "PLS": "#808080",
}

# One color per GI functional class (ssn_nodes_final.tsv, Category == GI).
# Order = legend order. "unknown" is the not-assigned class.
GI_CLASS_ORDER = ["virulence", "resistance", "metabolic", "symbiotic", "unknown"]
GI_CLASS_COLOR = {
    "virulence": "#4a1486",   # deep purple
    "resistance": "#c51b8a",  # magenta
    "metabolic": "#2c7fb8",   # blue
    "symbiotic": "#1b9e77",   # teal green
    "unknown": "#cbc9e2",     # pale lavender (not assigned)
}

# IS families keep their existing warm ramp (light -> dark).
IS_FAMILY_COLOR = {
    "IS110": "#ffeda0",
    "IS1595": "#feb24c",
    "IS3": "#fd8d3c",
    "IS4": "#f03b20",
    "IS481": "#bd0026",
    "IS5": "#800026",
}

CLADE_COLOR = {"K0": "#1f77b4", "K1": "#ff7f0e", "K2": "#2ca02c"}

# One colormap family per analysis.
CMAP = {
    "ds_gene_content": None,        # discrete reds, defined in plot_ds_gene_heatmap_by_clade.py
    "mge_ds_correlation": "RdBu_r",  # diverging, -1..+1
    "ds_mge_overlap": "Blues",       # % of DS system loci inside MGE
    "txss_mge_overlap": "GnBu",      # % of T3E/TSS genes inside MGE
}
