#!/usr/bin/env python3
"""
plot_ds_gene_heatmap_by_clade.py -- GENE-BY-GENE defense-system heatmap per
genome, grouped by Xtt clade (K0/K1/K2). Replaces the system-level
"total gene count" heatmap (plot_ds_heatmap_by_clade.py), which mixed three
different things in one number (distinct genes present, extra copies of the
same gene, and separate partial instances), and so could not distinguish a
complete multi-gene system from a lone fragment of one.

Columns = individual genes/subunits ("components") of each of the 25
both-tool-confirmed systems, grouped by system and then by mechanism
category. A component counts as present in a genome if DefenseFinder OR
PADLOC reported it there; copies = max(DefenseFinder copies, PADLOC copies)
(same max rule as the consensus matrix). Components are the union of gene
names observed anywhere in the dataset (the tools' outputs only list
detected genes, not each model's full mandatory gene list), so "complete"
below means "all components ever observed for that system", a proxy.

DefenseFinder gene names come from `name_of_profiles_in_sys` in
*_defense_finder_systems.tsv; PADLOC from `protein.name` in *_padloc.csv.
Gene-name mapping between tools is the COMPONENTS table below (built from
ds_gene_name_inventory.tsv).

ASSUMPTION TO CONFIRM: Gao_TerY -- DefenseFinder reports TerYA/TerYB/TerYC
while PADLOC reports TerY/PP2C/STK_OB; they are paired in that order here
(2 genomes only), but the pairing is not verified.

Outputs:
  defense_systems/plots/ds_gene_heatmap_by_clade.png
  defense_systems/matrices/ds_gene_level_matrix.tsv           (Genome x System|gene, copies)
  defense_systems/matrices/ds_gene_level_system_completeness.tsv

Usage
-----
python defense_systems/scripts/plot_ds_gene_heatmap_by_clade.py
"""
import re
import sys
import warnings
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.patches import Patch
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from build_ds_consensus_matrix import CROSSWALK, CERTAIN_SYSTEMS  # noqa: E402

PROJECT_ROOT = HERE.parents[1]
sys.path.insert(0, str(PROJECT_ROOT))
from genome_ids import canon  # noqa: E402  (unifies AB22010_2 / AB22010-2 etc.)

DF_DIR = PROJECT_ROOT / "defense_systems" / "defensefinder" / "raw_outputs"
PADLOC_DIR = PROJECT_ROOT / "defense_systems" / "padloc" / "results"
CLADE_FILE = PROJECT_ROOT / "metadata" / "genome_clades_master.tsv"
OUT_PLOT = PROJECT_ROOT / "defense_systems" / "plots" / "ds_gene_heatmap_by_clade.png"
OUT_MATRIX = PROJECT_ROOT / "defense_systems" / "matrices" / "ds_gene_level_matrix.tsv"
OUT_COMPLETE = PROJECT_ROOT / "defense_systems" / "matrices" / "ds_gene_level_system_completeness.tsv"

CLADE_ORDER = ["K0", "K1", "K2"]
CLADE_GAP_ROWS = 2

# --- gene/component mapping: system -> [(label, DefenseFinder regex, PADLOC regex)]
COMPONENTS = {
    "RM_Type_I": [("MTase", r"Type_I_MTases", r"^MTase_I$"),
                  ("REase", r"Type_I_REases", r"^REase_I$"),
                  ("S (specificity)", r"Type_I_S_", r"^Specificity_I$")],
    "RM_Type_II": [("MTase", r"Type_II_MTases", r"^MTase_II$"),
                   ("REase", r"Type_II_REase", r"^REase_II$")],
    "RM_Type_IIG": [("REase-MTase", r"Type_IIG", r"^REase_MTase_IIG$")],
    "RM_Type_III": [("MTase", r"Type_III_MTases", r"^MTase_III$"),
                    ("REase", r"Type_III_REases", r"^REase_III$")],
    "Cas_Type_I-C": [("cas1", r"^cas1_", r"^Cas1c$"),
                     ("cas2", r"^cas2_", r"^Cas2c$"),
                     ("cas3", r"^cas3(HD)?_", r"^Cas3c$"),
                     ("cas4", r"^cas4_", r"^Cas4c$"),
                     ("cas5", r"^cas5_", r"^Cas5c$"),
                     ("cas7", r"^cas7_", r"^Cas7c$"),
                     ("cas8c", r"^cas8c_", r"^Cas8c$")],
    "AbiD": [("AbiD", r"^AbiD__AbiD$", r"^AbiD$")],
    "DarTG": [("DarT", r"__DarT$", r"^DarT$"), ("DarG", r"__DarG$", r"^DarG$")],
    "PrrC": [("PrrC", r"^PrrC__PrrC$", r"^PrrC$")],
    "SoFic": [("SoFic", r"__SoFic$", r"^SoFic$")],
    "Retron_XI": [("RT", r"RT_11$", r"^RT-Protease_XI$")],
    "Retron_XII": [("RT", r"RT_12$", r"^RT-TIR_XII$")],
    "Wadjet_I": [("JetA", r"__JetA", r"^JetA\d*$"), ("JetB", r"__JetB", r"^JetB\d*$"),
                 ("JetC", r"__JetC", r"^JetC\d*$"), ("JetD", r"__JetD", r"^JetD\d*$")],
    "CBASS": [("Cyclase", r"Cyclase", r"^Cyclase$"), ("Effector", r"Effector", r"^Effector$")],
    "Gabija": [("GajA", r"__GajA$", r"^GajA$"), ("GajB", r"__GajB", r"^GajB$")],
    "Gao_Ppl": [("Ppl", r"__PplA$", r"^Ppl$")],
    "Gao_TerY": [("TerYA/TerY", r"__TerYA$", r"^TerY$"),
                 ("TerYB/PP2C", r"__TerYB$", r"^PP2C$"),
                 ("TerYC/STK_OB", r"__TerYC$", r"^STK_OB$")],
    "HEC-01": [("HEC-01A", r"HEC-01A$", r"^HEC-01A$"), ("HEC-01B", r"HEC-01B$", r"^HEC-01B$")],
    "HEC-03": [("HEC-03A", r"HEC-03A$", r"^HEC-03A$"), ("HEC-03B", r"HEC-03B$", r"^HEC-03B$")],
    "HEC-09": [("HEC-09", r"^HEC-09__HEC-09$", r"^HEC-09$")],
    "Lamassu": [("LmuA", r"__LmuA", r"^LmuA$"), ("LmuB", r"__LmuB", r"^LmuB$"),
                ("LmuC", r"__LmuC", r"^LmuC$")],
    "Mokosh_TypeII": [("MkoC", r"__MkoC$", r"^MkoC$")],
    "PD-Lambda-5": [("PD-Lambda-5_A", r"_A$", r"^PD-Lambda-5_A$"),
                    ("PD-Lambda-5_B", r"_B$", r"^PD-Lambda-5_B$")],
    "Septu": [("PtuA", r"__PtuA$", r"^PtuA\d*$"), ("PtuB", r"__PtuB$", r"^PtuB\d*$")],
    "Shango": [("SngA", r"__SngA$", r"^SngA$"), ("SngB", r"__SngB$", r"^SngB$"),
               ("SngC", r"__SngC$", r"^SngC$")],
    "Zorya_TypeIII": [("ZorA", r"^Zorya__ZorA$", r"^ZorA3$"),
                      ("ZorB", r"^Zorya__ZorB$", r"^ZorB3$"),
                      ("ZorF3", r"__ZorF3$", r"^ZorF3$"),
                      ("ZorG3", r"__ZorG3$", r"^ZorG3$")],
}

CATEGORY_ORDER = [
    "Restriction-Modification (RM)",
    "CRISPR-Cas (adaptive immunity)",
    "Toxin-antitoxin / Abortive infection",
    "Retrons",
    "Anti-plasmid",
    "Cyclic-nucleotide signaling (CBASS)",
    "Other/recently discovered anti-phage",
]
SYSTEM_CATEGORY = {
    "RM_Type_I": CATEGORY_ORDER[0], "RM_Type_II": CATEGORY_ORDER[0],
    "RM_Type_IIG": CATEGORY_ORDER[0], "RM_Type_III": CATEGORY_ORDER[0],
    "Cas_Type_I-C": CATEGORY_ORDER[1],
    "AbiD": CATEGORY_ORDER[2], "DarTG": CATEGORY_ORDER[2],
    "PrrC": CATEGORY_ORDER[2], "SoFic": CATEGORY_ORDER[2],
    "Retron_XI": CATEGORY_ORDER[3], "Retron_XII": CATEGORY_ORDER[3],
    "Wadjet_I": CATEGORY_ORDER[4],
    "CBASS": CATEGORY_ORDER[5],
    "Gabija": CATEGORY_ORDER[6], "Gao_Ppl": CATEGORY_ORDER[6], "Gao_TerY": CATEGORY_ORDER[6],
    "HEC-01": CATEGORY_ORDER[6], "HEC-03": CATEGORY_ORDER[6], "HEC-09": CATEGORY_ORDER[6],
    "Lamassu": CATEGORY_ORDER[6], "Mokosh_TypeII": CATEGORY_ORDER[6],
    "PD-Lambda-5": CATEGORY_ORDER[6], "Septu": CATEGORY_ORDER[6],
    "Shango": CATEGORY_ORDER[6], "Zorya_TypeIII": CATEGORY_ORDER[6],
}
CATEGORY_COLOR = {
    "Restriction-Modification (RM)": "#1b9e77",
    "CRISPR-Cas (adaptive immunity)": "#d95f02",
    "Toxin-antitoxin / Abortive infection": "#7570b3",
    "Retrons": "#e7298a",
    "Anti-plasmid": "#66a61e",
    "Cyclic-nucleotide signaling (CBASS)": "#e6ab02",
    "Other/recently discovered anti-phage": "#666666",
}

# Discrete red scale: 0 absent (white), 1 copy, 2 copies, 3+ copies. Discrete
# steps with strong luminance jumps make presence/absence and copy number
# readable by eye (a continuous ramp makes 1 vs 2 hard to tell apart).
LEVEL_COLORS = ["#ffffff", "#fc9272", "#de2d26", "#67000d"]
LEVEL_LABELS = ["0 (absent)", "1 copy", "2 copies", "3+ copies"]


def norm(x):
    return None if pd.isna(x) else str(x)


def classify(system, tool, gene, unmatched):
    idx = 1 if tool == "DF" else 2
    hits = [c[0] for c in COMPONENTS[system] if re.search(c[idx], gene)]
    if len(hits) == 1:
        return hits[0]
    unmatched[(system, tool, gene, len(hits))] += 1
    return None


def collect_copies(genomes):
    df_map, pl_map = {}, {}
    for cs, dt, ds, ps in CROSSWALK:
        if cs not in CERTAIN_SYSTEMS:
            continue
        if dt is not None:
            df_map[(dt, ds)] = cs
        if ps is not None:
            pl_map[ps] = cs

    copies = defaultdict(lambda: {"DF": 0, "PL": 0})
    unmatched = defaultdict(int)

    for f in sorted(DF_DIR.glob("*/*_defense_finder_systems.tsv")):
        genome = canon(f.name.replace("_defense_finder_systems.tsv", ""))
        if genome not in genomes:
            continue
        t = pd.read_csv(f, sep="\t")
        for _, r in t.iterrows():
            cs = df_map.get((norm(r.get("type")), norm(r.get("subtype"))))
            if cs is None:
                continue
            for g in str(r.get("name_of_profiles_in_sys", "")).split(","):
                g = g.strip()
                if not g:
                    continue
                label = classify(cs, "DF", g, unmatched)
                if label is not None:
                    copies[(genome, cs, label)]["DF"] += 1

    for f in sorted(PADLOC_DIR.glob("*/*_padloc.csv")):
        genome = canon(f.name.replace("_padloc.csv", ""))
        if genome not in genomes:
            continue
        t = pd.read_csv(f)
        for _, r in t.iterrows():
            cs = pl_map.get(str(r.get("system")))
            if cs is None:
                continue
            g = norm(r.get("protein.name")) or norm(r.get("hmm.name")) or "NA"
            label = classify(cs, "PL", g, unmatched)
            if label is not None:
                copies[(genome, cs, label)]["PL"] += 1

    if unmatched:
        print("[WARN] gene names not assigned to exactly one component "
              "(system, tool, gene, n_matches) -> count:")
        for k, v in sorted(unmatched.items()):
            print(f"    {k}: {v}")
    return copies


def build_matrix(genomes, copies):
    cols = [(s, c[0]) for cat in CATEGORY_ORDER
            for s in COMPONENTS if SYSTEM_CATEGORY[s] == cat for c in COMPONENTS[s]]
    data = {}
    for (s, label) in cols:
        data[f"{s}|{label}"] = [max(copies[(g, s, label)]["DF"], copies[(g, s, label)]["PL"])
                                if (g, s, label) in copies else 0 for g in genomes]
    mat = pd.DataFrame(data, index=genomes)
    return mat


def system_completeness(mat):
    rows = []
    systems = sorted({c.split("|")[0] for c in mat.columns})
    for s in systems:
        cols = [c for c in mat.columns if c.split("|")[0] == s]
        pres = (mat[cols] > 0)
        n_comp = pres.sum(axis=1)
        rows.append({
            "system": s, "category": SYSTEM_CATEGORY[s], "n_components": len(cols),
            "n_genomes_any": int((n_comp > 0).sum()),
            "n_genomes_complete": int((n_comp == len(cols)).sum()),
            "n_genomes_partial": int(((n_comp > 0) & (n_comp < len(cols))).sum()),
        })
    return pd.DataFrame(rows)


def order_rows(mat, clade_map):
    mat = mat.copy()
    mat["Clade"] = mat.index.map(clade_map)
    comp_cols = [c for c in mat.columns if c != "Clade"]
    mat["_npres"] = (mat[comp_cols] > 0).sum(axis=1)
    mat["_tot"] = mat[comp_cols].sum(axis=1)
    ordered, spans, row = [], {}, 0
    for i, clade in enumerate(CLADE_ORDER):
        sub = mat[mat["Clade"] == clade].sort_values(["_npres", "_tot"], ascending=False)
        if sub.empty:
            continue
        spans[clade] = (row, row + len(sub) - 1)
        ordered.extend(sub.index.tolist())
        row += len(sub)
        if i < len(CLADE_ORDER) - 1:
            row += CLADE_GAP_ROWS
    return ordered, spans


def order_columns(mat):
    """Category -> system (most prevalent first) -> components in table order."""
    comp_cols = list(mat.columns)
    sys_prev = {}
    for s in {c.split("|")[0] for c in comp_cols}:
        cols = [c for c in comp_cols if c.startswith(s + "|")]
        sys_prev[s] = int((mat[cols] > 0).any(axis=1).sum())
    ordered, sys_spans, cat_of_col = [], {}, {}
    col = 0
    for cat in CATEGORY_ORDER:
        systems = [s for s in sys_prev if SYSTEM_CATEGORY[s] == cat]
        for s in sorted(systems, key=lambda x: -sys_prev[x]):
            labels = [c for c in comp_cols if c.startswith(s + "|")]
            labels.sort(key=lambda c: [x[0] for x in COMPONENTS[s]].index(c.split("|", 1)[1]))
            sys_spans[s] = (col, col + len(labels) - 1)
            ordered.extend(labels)
            col += len(labels)
    return ordered, sys_spans


def plot(mat, ordered_genomes, clade_spans, ordered_cols, sys_spans):
    n_rows = clade_spans[CLADE_ORDER[-1]][1] + 1
    n_cols = len(ordered_cols)

    grid = np.full((n_rows, n_cols), np.nan)
    row_labels = [""] * n_rows
    gi = 0
    for clade in CLADE_ORDER:
        if clade not in clade_spans:
            continue
        s, e = clade_spans[clade]
        for r in range(s, e + 1):
            g = ordered_genomes[gi]
            grid[r, :] = np.minimum(mat.loc[g, ordered_cols].to_numpy(), 3)
            row_labels[r] = g
            gi += 1

    fig_w = max(14, n_cols * 0.30 + 4)
    fig_h = max(10, n_rows * 0.16 + 3)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h), facecolor="white")

    cmap = ListedColormap(LEVEL_COLORS)
    cmap.set_bad("white")
    norm_ = BoundaryNorm([-0.5, 0.5, 1.5, 2.5, 3.5], cmap.N)
    ax.imshow(np.ma.masked_invalid(grid), aspect="auto", cmap=cmap, norm=norm_,
              interpolation="none")

    # light grey cell grid so absent (white) cells still show the column structure
    ax.set_xticks(np.arange(-0.5, n_cols, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, n_rows, 1), minor=True)
    ax.grid(which="minor", color="#d9d9d9", linewidth=0.3)
    ax.tick_params(which="minor", bottom=False, left=False)

    gene_labels = [c.split("|", 1)[1] for c in ordered_cols]
    ax.set_xticks(range(n_cols))
    ax.set_xticklabels(gene_labels, rotation=90, fontsize=6)
    ax.set_yticks(range(n_rows))
    ax.set_yticklabels(row_labels, fontsize=4.0)
    ax.tick_params(axis="both", length=0)
    for i, c in enumerate(ordered_cols):
        ax.get_xticklabels()[i].set_color(CATEGORY_COLOR[SYSTEM_CATEGORY[c.split("|")[0]]])

    # frame around the heatmap body, thicker separators between systems
    # (and heavier between categories); drawn only over the data rows.
    for sp in ax.spines.values():
        sp.set_visible(True)
        sp.set_color("black")
        sp.set_linewidth(1.2)
    prev_cat = None
    for s, (c0, c1) in sys_spans.items():
        cat = SYSTEM_CATEGORY[s]
        if c0 > 0:
            heavy = cat != prev_cat
            ax.plot([c0 - 0.5, c0 - 0.5], [-0.5, n_rows - 0.5],
                    color="black" if heavy else "#555555",
                    linewidth=1.8 if heavy else 0.9, zorder=4)
        prev_cat = cat
    # horizontal separators between clade blocks (gap rows are blank)
    for clade in CLADE_ORDER[:-1]:
        if clade in clade_spans:
            y = clade_spans[clade][1] + 0.5 + CLADE_GAP_ROWS / 2
            ax.axhline(y, color="black", linewidth=1.0)

    # system names + category-colored brackets above
    y_br = -0.5 - n_rows * 0.012
    tk = n_rows * 0.004
    for s, (c0, c1) in sys_spans.items():
        color = CATEGORY_COLOR[SYSTEM_CATEGORY[s]]
        ax.plot([c0 - 0.4, c1 + 0.4], [y_br, y_br], color=color, linewidth=2.2,
                clip_on=False, zorder=3)
        ax.plot([c0 - 0.4, c0 - 0.4], [y_br, y_br + tk], color=color, linewidth=1.4, clip_on=False)
        ax.plot([c1 + 0.4, c1 + 0.4], [y_br, y_br + tk], color=color, linewidth=1.4, clip_on=False)
        ax.text((c0 + c1) / 2, y_br - tk * 1.5, s, ha="center", va="bottom", rotation=90,
                fontsize=7, fontweight="bold", color=color, clip_on=False)

    # clade brackets on the left
    xb = -0.5 - n_cols * 0.012
    t = n_cols * 0.004
    for clade in CLADE_ORDER:
        if clade not in clade_spans:
            continue
        y0, y1 = clade_spans[clade]
        ax.plot([xb, xb], [y0 - 0.4, y1 + 0.4], color="black", linewidth=1.2, clip_on=False)
        ax.plot([xb, xb - t], [y0 - 0.4, y0 - 0.4], color="black", linewidth=1.2, clip_on=False)
        ax.plot([xb, xb - t], [y1 + 0.4, y1 + 0.4], color="black", linewidth=1.2, clip_on=False)
        ax.text(xb - t * 2.5, (y0 + y1) / 2, f"{clade} (n={y1 - y0 + 1})", ha="right",
                va="center", fontsize=11, fontweight="bold", rotation=90, clip_on=False)

    # legends (outside, right)
    copy_handles = [Patch(facecolor=c, edgecolor="#888888", label=l)
                    for c, l in zip(LEVEL_COLORS, LEVEL_LABELS)]
    leg1 = ax.legend(handles=copy_handles, title="Gene copies", loc="upper left",
                     bbox_to_anchor=(1.01, 1.0), fontsize=8, title_fontsize=9, frameon=True)
    ax.add_artist(leg1)
    cat_handles = [Patch(facecolor=CATEGORY_COLOR[c], label=c) for c in CATEGORY_ORDER
                   if any(SYSTEM_CATEGORY[s] == c for s in sys_spans)]
    ax.legend(handles=cat_handles, title="Mechanism category", loc="upper left",
              bbox_to_anchor=(1.01, 0.82), fontsize=8, title_fontsize=9, frameon=True)

    ax.set_title("Defense-system content per genome, gene by gene, grouped by "
                 "X. translucens pv. translucens clade (K0/K1/K2)", fontsize=13, pad=95)

    footnote = (
        "Columns = individual genes/subunits of each system (union of genes observed in any genome); "
        "system names on top, colored by mechanism category.\n"
        "A gene is present if DefenseFinder or PADLOC detected it; shade = copies "
        "(max of the two tools, capped at 3+). Only systems confirmed by both tools are shown; "
        "genes never detected in K0/K1/K2 are omitted.\n"
        "Rows ordered within clade by number of distinct genes present. A system is complete "
        "only if every column of its block is red (proxy: genes observed in this dataset)."
    )
    fig.text(0.5, 0.0, footnote, ha="center", va="top", fontsize=7.5, color="#666666",
             style="italic", multialignment="center")

    OUT_PLOT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PLOT, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"[DONE] Saved: {OUT_PLOT} ({OUT_PLOT.stat().st_size} bytes)")


def main():
    clades = pd.read_csv(CLADE_FILE, sep="\t", usecols=["Genome", "Clade"])
    clades = clades[clades["Clade"].isin(CLADE_ORDER)]
    clade_map = clades.set_index("Genome")["Clade"]
    genomes = clades["Genome"].tolist()
    print(f"[INFO] {len(genomes)} genomes in K0/K1/K2 "
          f"({clades['Clade'].value_counts().reindex(CLADE_ORDER).to_dict()})")

    missing = [s for s in CERTAIN_SYSTEMS if s not in COMPONENTS]
    if missing:
        print(f"[FATAL] systems in CERTAIN_SYSTEMS without a COMPONENTS entry: {missing}")
        return

    copies = collect_copies(set(genomes))
    mat = build_matrix(genomes, copies)

    never = [c for c in mat.columns if not (mat[c] > 0).any()]
    if never:
        print(f"[INFO] Dropping genes never detected in K0/K1/K2 genomes: {never}")
        mat = mat.drop(columns=never)
    sys_gone = sorted(set(COMPONENTS) - {c.split('|')[0] for c in mat.columns})
    if sys_gone:
        print(f"[INFO] Systems with no gene detected in K0/K1/K2 (omitted): {sys_gone}")

    OUT_MATRIX.parent.mkdir(parents=True, exist_ok=True)
    mat.to_csv(OUT_MATRIX, sep="\t")
    comp = system_completeness(mat).sort_values(["category", "n_genomes_any"],
                                                ascending=[True, False])
    comp.to_csv(OUT_COMPLETE, sep="\t", index=False)
    print(f"[INFO] Saved {OUT_MATRIX} and {OUT_COMPLETE}\n")
    print("[INFO] Per-system completeness (K0/K1/K2 genomes):")
    print(comp.to_string(index=False))

    ordered_genomes, clade_spans = order_rows(mat, clade_map)
    ordered_cols, sys_spans = order_columns(mat)
    plot(mat, ordered_genomes, clade_spans, ordered_cols, sys_spans)


if __name__ == "__main__":
    main()
