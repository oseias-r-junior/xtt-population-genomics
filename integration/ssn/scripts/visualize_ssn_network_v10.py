#!/usr/bin/env python3
"""
visualize_ssn_network_v10.py -- SSN visualization v10.

Fix vs v8, per PI review:
  9. PLS SHELF: v8's fixed PLS_SHELF_Y=0.0 clipped slightly into the tail
     end of the GI dataset (the GI cluster's convex-hull halo extends a bit
     further down than the raw node positions). Per the PI's own framing
     ("bring it halfway back from where it was to where it is now"), PLS
     now sits at the midpoint between the v7 position (near the ground/
     legend) and the v8 position (0.0), computed from this run's actual
     y_ground rather than another guessed constant. x unchanged throughout.

Fix vs v7, per PI review (final polish before considering the size-legend
task closed):
  8. PLS SHELF REPOSITIONED: even after moving the legend to lower-right
     (v7), the PLS shelf was still low enough (just above the ground rows)
     to fall inside the legend's vertical footprint -- 2 of the 4 plasmid
     circles ended up hidden, and the "Plasmids" category title (positioned
     from the actual PLS node y) collided with the "plasmids (singletons)"
     annotation, rendering as garbled overlapping text. Fix: PLS singletons
     now sit at a fixed y (PLS_SHELF_Y = 0.0) in the empty gap below the GI
     cluster, with x unchanged (still centered on CATEGORY_CENTERS["PLS"][0]
     = 1.8, per the request to move only vertically) -- comfortably clear of
     both the GI cluster above and the legend/ground rows far below.

Fixes vs v6, per PI + lead-author review:
  6. LEGEND POSITION (overview only): moved from lower-left to lower-right,
     per the lead author's suggestion -- that quadrant is where PLS lives
     (only 4 singleton plasmids, no clusters), confirmed empty otherwise.
     The shared PPH/GI/IS singleton ground rows are now confined to a
     left-biased x-range (GROUND_X_MAX_FRAC) instead of spanning the full
     figure width, so they don't run underneath the relocated legend, and
     spread over 3 rows instead of 2 to keep density reasonable in the
     narrower range.
  7. NODE-SIZE LEGEND FIXED TO MATCH REALITY: in the overview, the old
     size legend showed a single GLOBAL bp min/max (across all 4
     categories combined) next to two dots -- but actual marker size is
     normalized INDEPENDENTLY per category (each category log-stretched
     across the same visual size range). That global number wasn't tied
     to any one category's real scale, which is exactly the disconnect
     the PI flagged. The two proxy dots are now generic ("smallest/
     largest node IN ITS CATEGORY"), and a per-category bp range is listed
     below them (one line per PPH/GI/IS/PLS, color-matched), which is what
     each dot actually corresponds to for that category. Standalone
     per-category plots are unaffected (already showed one unambiguous
     real bp range).

Fixes vs v5 (still present, from the previous round):
  1. LEGEND SIZE MARKERS: the two size-legend proxy dots (min/max bp) used
     sqrt-scaled markersize directly, which for the real MIN/MAX_MARKER_SIZE
     range produced a huge marker that overflowed the legend box and
     overlapped the smaller one. Legend proxy markers are now capped to a
     fixed, legend-safe size pair (visual only -- the real network markers
     are unaffected) with extra row spacing.
  2. SINGLETON VISIBILITY: singleton nodes were placed in a single ground
     row, and some ended up hidden behind the legend box (bottom-left).
     Mirrors the Uceda-Campos et al. 2022 approach: singleton "shared"
     nodes (PPH/GI/IS) are now split across TWO ground rows instead of one,
     roughly halving per-row density.
  3. CLUSTER LABEL SIZE: group-ID labels (e.g. IS-G27, GI-G13) increased
     ~2.5x (6.5pt -> 16pt) for readability at first glance, per feedback
     from a non-specialist reader.
  4. ALL CLUSTERS LABELED: previously only clusters at/above a per-category
     minimum size got a label (e.g. GI >=4, IS >=6), silently omitting
     smaller (2-3 member) clusters. Every cluster (>=2 members, i.e. every
     entry in cluster_groups) is now labeled.
  5. IS-FAMILY-COLORED CARRIER RINGS: host nodes (PPH/GI/PLS) that carry an
     internal IS no longer get a flat category-based ring color (gold/
     purple/dark-gray). The ring is now colored by the FAMILY of the IS
     found inside that host (same warm YlOrRd palette used for IS fill),
     giving an approximate visual correlation between host and the type of
     IS it carries. This uses a new IS_carrier_family column produced by
     add_is_carrier_family.py (run that script BEFORE this one). Hosts
     with an unresolved/unknown family fall back to a neutral dark-grey
     ring; hosts predating that column (older nodes files) fall back to
     the old flat per-category ring color.

Outputs (integration/ssn/plots/):
    ssn_prophages_v10.png
    ssn_genomic_islands_v10.png
    ssn_plasmids_v10.png
    ssn_insertion_sequences_v10.png
    ssn_overview_v10.png

Usage
-----
python add_is_carrier_family.py   # once, before this script (skip if already run)
python visualize_ssn_network_v10.py
"""

import math
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
from matplotlib.collections import LineCollection
from matplotlib.colors import to_rgba
from matplotlib.lines import Line2D
from scipy.spatial import ConvexHull

warnings.filterwarnings("ignore")

try:
    from adjustText import adjust_text
    HAS_ADJUST_TEXT = True
except ImportError:
    HAS_ADJUST_TEXT = False
    print("[WARN] adjustText not installed -- labels will not be auto-adjusted.")

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
PROJECT_ROOT = Path(__file__).resolve().parents[3]
SSN_MATRICES = PROJECT_ROOT / "integration" / "ssn" / "matrices"
PLOTS_DIR    = PROJECT_ROOT / "integration" / "ssn" / "plots"

NODES_FILE = SSN_MATRICES / "ssn_nodes_final.tsv"
EDGES_FILE = SSN_MATRICES / "ssn_edges_final.tsv"

# ---------------------------------------------------------------------------
# Visual constants
# ---------------------------------------------------------------------------
CATEGORY_COLOR = {
    "PPH": "#1a1a1a",   # near-black: no classification variance to encode
    "GI":  "#807dba",   # representative purple (real fill varies by function)
    "IS":  "#fd8d3c",   # representative orange (real fill varies by family)
    "PLS": "#808080",   # grey: no classification variance to encode
}
CATEGORY_NAME = {
    "PPH": "Prophages",
    "GI":  "Genomic islands",
    "IS":  "Insertion sequences",
    "PLS": "Plasmids",
}
# v10: ONE color per GI functional class (v9 paired resistance+virulence and
# metabolic+symbiotic). All five sit outside the warm IS-family ramp and
# outside black/grey (PPH/PLS), so a GI node is never mistaken for another
# MGE. Same values as palette.py.
GI_CLASS_COLOR = {
    "virulence":  "#4a1486",   # deep purple
    "resistance": "#c51b8a",   # magenta
    "metabolic":  "#2c7fb8",   # blue
    "symbiotic":  "#1b9e77",   # teal green
    "unknown":    "#cbc9e2",   # pale lavender (function not assigned)
}
# ColorBrewer YlOrRd-6 (yellow -> orange -> red -> wine); also reused for
# host-node carrier rings (fix #5), so the ring color and the IS fill color
# for that family always match.
IS_FAMILY_COLOR = {
    "IS110":  "#ffeda0",
    "IS1595": "#feb24c",
    "IS3":    "#fd8d3c",
    "IS4":    "#f03b20",
    "IS481":  "#bd0026",
    "IS5":    "#800026",
}
BORDER_DEFAULT = "#666666"
COL_EDGE       = "#999999"
# Fallback flat ring colors, used only when IS_carrier_family is missing
# (older nodes file) or unresolved (no coordinate match found).
BORDER_CARRIER_BY_HOST = {
    "PPH": "#FFD700",
    "GI":  "#9b59b6",
    "PLS": "#404040",
}
BORDER_CARRIER_UNKNOWN_FAMILY = "#555555"

CATEGORY_CENTERS = {
    "PPH": np.array([-1.3,  1.1]),
    "GI":  np.array([ 1.3,  1.1]),
    "PLS": np.array([ 1.8, -1.2]),
    "IS":  np.array([-0.8, -1.2]),
}
CATEGORY_AREA_SCALE = {"PPH": 0.55, "GI": 0.95, "PLS": 0.25, "IS": 1.3}

# Fix #4: every cluster (>=2 members) gets a label now -- no per-category
# minimum size threshold that could silently omit small clusters.
MIN_LABEL_SIZE = {"PPH": 2, "GI": 2, "IS": 2, "PLS": 2}
CLUSTER_LABEL_FONTSIZE = 16  # fix #3: was 6.5, ~2.5x increase

TOP_K_EDGES_STANDALONE = {"PPH": None, "GI": None, "PLS": None, "IS": 4}
TOP_K_EDGES_OVERVIEW   = {"PPH": None, "GI": 5,    "PLS": None, "IS": 2}

GOLDEN_ANGLE = 2.399963229728653  # radians

MIN_MARKER_SIZE = 25.0
MAX_MARKER_SIZE = 380.0

# fix #1: fixed, legend-safe proxy marker sizes (points), independent of the
# real MIN/MAX_MARKER_SIZE range used on the actual network nodes.
LEGEND_SIZE_MARKER_PT = (6.0, 13.0)

# fix #2: shared singleton ground now spans this many rows.
N_SINGLETON_ROWS = 2            # standalone per-category plots
N_SINGLETON_ROWS_OVERVIEW = 3   # overview: narrower x-range (see below) needs more rows
SINGLETON_ROW_GAP = 0.16

# v7 fix: reserve the bottom-right quadrant of the overview for the legend
# (that's where PLS already lives -- only 4 singleton plasmids, confirmed
# empty otherwise). Shared PPH/GI/IS ground-row x-range is capped at this
# fraction of x_span instead of spanning the full +/-x_span width.
GROUND_X_MAX_FRAC = 0.30
LEGEND_LOC_STANDALONE = ("lower left", (0.0, 0.0))
LEGEND_LOC_OVERVIEW = ("lower right", (1.0, 0.0))

# v8 fix: the PLS shelf (only 4 singleton plasmids) was still low enough to
# sit inside the relocated legend's vertical footprint, and the "Plasmids"
# category title (positioned from the PLS nodes' own y) ended up overlapping
# the "plasmids (singletons)" annotation as a result. Move PLS straight up
# (same x = CATEGORY_CENTERS["PLS"][0] = 1.8, unchanged) into the empty gap
# below the GI cluster.
#
# v9 fix: PLS_SHELF_Y=0.0 (v8) turned out to sit slightly UNDER the tail end
# of the GI dataset -- the GI cluster's convex-hull halo (drawn with padding
# beyond the raw node positions, see draw_group_hulls) extends further down
# than the raw "GI center 1.1 - area_scale 0.95 = 0.15" estimate suggested.
# Rather than guess another fixed number, the PI's own framing ("bring it
# halfway back from where it was before to where it is now") is applied
# literally at runtime: pls_shelf_y is now the midpoint between the v7
# formula's position (near the ground/legend) and the v8 fixed value (0.0),
# computed from THIS run's actual y_ground -- see PLS_SHELF_Y_FALLBACK below
# and its use in main().
PLS_SHELF_Y_V8 = 0.0


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load_data():
    nodes = pd.read_csv(NODES_FILE, sep="\t")
    edges = pd.read_csv(EDGES_FILE, sep="\t")
    print(f"[INFO] Loaded {len(nodes)} nodes, {len(edges)} edges")
    if "IS_carrier_family" not in nodes.columns:
        print("[WARN] IS_carrier_family column not found -- run add_is_carrier_family.py "
              "first for family-colored carrier rings. Falling back to flat per-category rings.")
    return nodes, edges


def top_k_per_node(edges: pd.DataFrame, k) -> pd.DataFrame:
    if k is None:
        return edges
    long = pd.concat([
        edges.rename(columns={"Source": "node", "Target": "other"}),
        edges.rename(columns={"Target": "node", "Source": "other"}),
    ], ignore_index=True)
    long = long.sort_values("bitscore", ascending=False)
    top = long.groupby("node").head(k)
    pairs, kept = set(), []
    for _, row in top.iterrows():
        key = tuple(sorted((row["node"], row["other"])))
        if key not in pairs:
            pairs.add(key)
            kept.append(row)
    if not kept:
        return edges.iloc[0:0]
    out = pd.DataFrame(kept).rename(columns={"node": "Source", "other": "Target"})
    return out[["Source", "Target", "pident", "bitscore", "qcovhsp", "alignment_length"]]


# ---------------------------------------------------------------------------
# Node size: per-category log-normalized min-max
# ---------------------------------------------------------------------------

def compute_size_stats(nodes_subset: pd.DataFrame) -> dict:
    stats = {}
    for cat, grp in nodes_subset.groupby("Category"):
        logs = np.log10(grp["Length_bp"].clip(lower=1))
        stats[cat] = (float(logs.min()), float(logs.max()),
                      int(grp["Length_bp"].min()), int(grp["Length_bp"].max()))
    return stats


def node_size(length_bp: float, cat: str, stats: dict) -> float:
    length_bp = max(length_bp, 1)
    lo, hi, _, _ = stats.get(cat, (0.0, 1.0, 1, 1))
    val = math.log10(length_bp)
    if hi - lo < 1e-9:
        norm = 0.5
    else:
        norm = (val - lo) / (hi - lo)
    norm = min(max(norm, 0.0), 1.0)
    return MIN_MARKER_SIZE + norm * (MAX_MARKER_SIZE - MIN_MARKER_SIZE)


def size_legend_handles(nodes_subset: pd.DataFrame, single_cat=None):
    """Fix #1: use fixed legend-safe marker sizes (LEGEND_SIZE_MARKER_PT)
    instead of a literal sqrt-area conversion of MIN/MAX_MARKER_SIZE -- the
    literal conversion made the "max" dot large enough to overflow the
    legend box and overlap the "min" dot.

    v7 fix: node size is normalized INDEPENDENTLY per category (each
    category is stretched log10 min-max across the same marker-size range
    -- see node_size()), so a single global bp min/max label (as in v5/v6)
    was misleading in the overview: it wasn't tied to any one category's
    actual scale. For a single-category plot (single_cat set) the real
    bp min/max for that one category is shown, which is unambiguous. For
    the overview (single_cat=None) these two dots are now generic
    ("smallest/largest node IN ITS OWN CATEGORY"), and the actual bp
    ranges per category are listed separately by category_size_range_handles()."""
    lo_ms, hi_ms = LEGEND_SIZE_MARKER_PT
    if single_cat is not None:
        sub = nodes_subset[nodes_subset["Category"] == single_cat]
        if sub.empty:
            return []
        lo_bp, hi_bp = int(sub["Length_bp"].min()), int(sub["Length_bp"].max())
        if lo_bp == hi_bp:
            return []
        return [
            Line2D([0], [0], marker='o', linestyle='none', markerfacecolor='#aaaaaa',
                   markeredgecolor='#555555', markersize=lo_ms,
                   label=f"{lo_bp:,} bp"),
            Line2D([0], [0], marker='o', linestyle='none', markerfacecolor='#aaaaaa',
                   markeredgecolor='#555555', markersize=hi_ms,
                   label=f"{hi_bp:,} bp"),
        ]
    if nodes_subset.empty:
        return []
    return [
        Line2D([0], [0], marker='o', linestyle='none', markerfacecolor='#aaaaaa',
               markeredgecolor='#555555', markersize=lo_ms,
               label="smallest node IN ITS CATEGORY"),
        Line2D([0], [0], marker='o', linestyle='none', markerfacecolor='#aaaaaa',
               markeredgecolor='#555555', markersize=hi_ms,
               label="largest node IN ITS CATEGORY"),
    ]


def category_size_range_handles(nodes_subset: pd.DataFrame):
    """v7: per-category bp min-max, shown only in the overview (where size
    scaling is relative/independent per category -- see size_legend_handles
    docstring). Answers directly: what does 'smallest/largest in its
    category' mean in actual bp, for EACH of the 4 categories."""
    handles = []
    for cat in ("PPH", "GI", "IS", "PLS"):
        sub = nodes_subset[nodes_subset["Category"] == cat]
        if sub.empty:
            continue
        lo_bp, hi_bp = int(sub["Length_bp"].min()), int(sub["Length_bp"].max())
        color = CATEGORY_COLOR.get(cat, "#888888")
        name = CATEGORY_NAME.get(cat, cat)
        if lo_bp == hi_bp:
            label = f"{name}: {lo_bp:,} bp (all equal)"
        else:
            label = f"{name}: {lo_bp:,}\u2013{hi_bp:,} bp"
        handles.append(mpatches.Patch(color=color, label=label))
    return handles


# ---------------------------------------------------------------------------
# Bubble packing of group centers
# ---------------------------------------------------------------------------

def pack_group_centers(sizes: dict, padding: float = 1.35):
    order = sorted(sizes.items(), key=lambda x: -x[1])
    radii = {g: math.sqrt(n) for g, n in sizes.items()}
    placed = {}
    for idx, (g, _) in enumerate(order):
        r = radii[g]
        if idx == 0:
            placed[g] = np.array([0.0, 0.0])
            continue
        theta0 = idx * GOLDEN_ANGLE
        step = 0.28 * (r + 1.0)
        found = False
        dist = step
        for k in range(1, 500):
            dist = step * k
            cand = np.array([dist * math.cos(theta0), dist * math.sin(theta0)])
            ok = True
            for og, opos in placed.items():
                if np.linalg.norm(cand - opos) < (radii[og] + r) * padding:
                    ok = False
                    break
            if ok:
                placed[g] = cand
                found = True
                break
        if not found:
            placed[g] = np.array([dist * math.cos(theta0), dist * math.sin(theta0)])
    return placed, radii


# ---------------------------------------------------------------------------
# Nested layout
# ---------------------------------------------------------------------------

def compute_nested_layout(cat_nodes, meta, full_G):
    groups = {}
    for n in cat_nodes:
        g = meta.loc[n, "Group"] or f"__singleton__{n}"
        groups.setdefault(g, []).append(n)

    cluster_groups = {g: ns for g, ns in groups.items() if len(ns) > 1}
    singleton_nodes = [ns[0] for ns in groups.values() if len(ns) == 1]

    sizes = {g: len(ns) for g, ns in cluster_groups.items()}
    centers, radii = pack_group_centers(sizes) if sizes else ({}, {})

    pos = {}
    for g, ns in cluster_groups.items():
        subG = full_G.subgraph(ns)
        if subG.number_of_edges() > 0:
            sub_pos = nx.spring_layout(subG, seed=abs(hash(g)) % (2**31),
                                        k=0.6, iterations=150)
        else:
            rng = np.random.default_rng(abs(hash(g)) % (2**31))
            sub_pos = {n: rng.standard_normal(2) * 0.3 for n in ns}
        coords = np.array(list(sub_pos.values()))
        if len(coords) > 1:
            maxabs = np.abs(coords - coords.mean(axis=0)).max()
            if maxabs < 1e-9:
                maxabs = 1.0
            coords = (coords - coords.mean(axis=0)) / maxabs
        target_r = radii[g] * 0.85
        coords = coords * target_r
        for node, c in zip(sub_pos.keys(), coords):
            pos[node] = c + centers[g]

    return pos, singleton_nodes, centers, radii, cluster_groups


def place_singleton_row(singleton_nodes, y_level, x_span, jitter_seed=0):
    """Single-row placement -- still used for the small PLS shelf."""
    n = len(singleton_nodes)
    if n == 0:
        return {}
    xs = np.linspace(-x_span, x_span, n) if n > 1 else np.array([0.0])
    rng = np.random.default_rng(jitter_seed)
    ys = y_level + rng.uniform(-0.03, 0.03, size=n)
    return {node: np.array([x, y]) for node, x, y in zip(singleton_nodes, xs, ys)}


def place_singleton_rows(singleton_nodes, y_base, x_min, x_max, n_rows=2,
                          row_gap=SINGLETON_ROW_GAP, jitter_seed=0):
    """Fix #2: distribute singletons across n_rows stacked ground rows
    instead of one, so a dense set of singletons doesn't end up
    hidden/omitted behind the legend box at the bottom of the figure.
    v7: takes an explicit (x_min, x_max) range instead of a symmetric
    x_span, so the overview can reserve the bottom-right quadrant for the
    legend (see GROUND_X_MAX_FRAC) while standalone plots stay symmetric."""
    n = len(singleton_nodes)
    if n == 0:
        return {}
    rng = np.random.default_rng(jitter_seed)
    order = list(singleton_nodes)
    rng.shuffle(order)  # interleave so rows aren't biased by input ordering
    rows = [order[i::n_rows] for i in range(n_rows)]
    pos = {}
    for r_idx, row_nodes in enumerate(rows):
        if not row_nodes:
            continue
        m = len(row_nodes)
        xs = np.linspace(x_min, x_max, m) if m > 1 else np.array([(x_min + x_max) / 2])
        y_level = y_base + r_idx * row_gap
        ys = y_level + rng.uniform(-0.025, 0.025, size=m)
        for node, x, y in zip(row_nodes, xs, ys):
            pos[node] = np.array([x, y])
    return pos


# ---------------------------------------------------------------------------
# Graph construction
# ---------------------------------------------------------------------------

def build_graph(nodes: pd.DataFrame, edges: pd.DataFrame, categories=None):
    if categories is not None:
        nodes = nodes[nodes["Category"].isin(categories)]
    node_ids = set(nodes["NodeID"])
    edges = edges[edges["Source"].isin(node_ids) & edges["Target"].isin(node_ids)]
    G = nx.Graph()
    for _, row in nodes.iterrows():
        G.add_node(row["NodeID"])
    for _, row in edges.iterrows():
        G.add_edge(row["Source"], row["Target"],
                    weight=float(row["bitscore"]), bitscore=float(row["bitscore"]))
    return G, nodes


# ---------------------------------------------------------------------------
# Draw helpers
# ---------------------------------------------------------------------------

def node_visuals(nodes_subset: pd.DataFrame, node_order):
    meta = nodes_subset.set_index("NodeID")
    stats = compute_size_stats(nodes_subset)
    has_family_col = "IS_carrier_family" in nodes_subset.columns
    fills, borders, bwidths, sizes = [], [], [], []
    for n in node_order:
        row = meta.loc[n]
        cat = row["Category"]
        classification = str(row.get("Classification", ""))
        if cat == "GI":
            fills.append(GI_CLASS_COLOR.get(classification, "#cbc9e2"))
        elif cat == "IS":
            fills.append(IS_FAMILY_COLOR.get(classification, CATEGORY_COLOR["IS"]))
        else:
            fills.append(CATEGORY_COLOR.get(cat, "#bbbbbb"))

        if bool(row.get("IS_carrier_flag", False)):
            # Fix #5: ring colored by the family of the IS the host carries,
            # not a flat per-category color. Falls back gracefully.
            fam = str(row.get("IS_carrier_family", "") or "").strip() if has_family_col else ""
            if fam and fam in IS_FAMILY_COLOR:
                borders.append(IS_FAMILY_COLOR[fam])
            elif fam:
                borders.append(BORDER_CARRIER_UNKNOWN_FAMILY)
            else:
                borders.append(BORDER_CARRIER_BY_HOST.get(cat, "#000000"))
            bwidths.append(2.6)
        else:
            borders.append(BORDER_DEFAULT)
            bwidths.append(0.4)
        sizes.append(node_size(row["Length_bp"], cat, stats))
    return fills, borders, bwidths, sizes


def draw_edges_with_shadow(ax, G, pos):
    edge_bitscores = [d["bitscore"] for _, _, d in G.edges(data=True)]
    max_bs = max(edge_bitscores) if edge_bitscores else 1.0
    segments, halo_colors, halo_widths = [], [], []
    sharp_colors, sharp_widths = [], []
    for u, v, d in G.edges(data=True):
        if u not in pos or v not in pos:
            continue
        frac = d["bitscore"] / max_bs
        seg = [pos[u], pos[v]]
        segments.append(seg)
        halo_colors.append(to_rgba(COL_EDGE, 0.035))
        halo_widths.append(2.0 + 3.0 * frac)
        sharp_colors.append(to_rgba(COL_EDGE, min(0.08 + 0.3 * frac, 0.55)))
        sharp_widths.append(0.25 + 1.0 * frac)
    if segments:
        ax.add_collection(LineCollection(segments, colors=halo_colors,
                                          linewidths=halo_widths, zorder=1))
        ax.add_collection(LineCollection(segments, colors=sharp_colors,
                                          linewidths=sharp_widths, zorder=2))


def draw_group_hulls(ax, cluster_groups, pos, category_color, meta,
                      min_label_size, label_texts):
    """Halo behind each non-trivial group. Padding is PURELY proportional to
    the group's own point spread (no fixed additive term) so it scales
    correctly whether this is drawn at native (standalone) scale or after
    CATEGORY_AREA_SCALE compression (overview)."""
    for g, ns in cluster_groups.items():
        pts = np.array([pos[n] for n in ns if n in pos])
        if len(pts) < 2:
            continue
        color = category_color
        if len(pts) == 2:
            d = pts[1] - pts[0]
            seglen = np.linalg.norm(d)
            if seglen < 1e-9:
                d = np.array([1.0, 0.0])
                seglen = 1.0
            perp = np.array([-d[1], d[0]]) / seglen
            pad_w = max(seglen * 0.35, seglen * 0.15 + 1e-4)
            p0, p1 = pts[0], pts[1]
            poly_pts = [p0 + perp * pad_w, p1 + perp * pad_w,
                        p1 - perp * pad_w, p0 - perp * pad_w]
            poly = plt.Polygon(poly_pts, closed=True, facecolor=color,
                                alpha=0.10, edgecolor=color, linewidth=0.6,
                                zorder=0)
            ax.add_patch(poly)
            centroid = pts.mean(axis=0)
        else:
            centroid = pts.mean(axis=0)
            dists = np.linalg.norm(pts - centroid, axis=1)
            mean_d = dists.mean()
            pad = mean_d * 0.55 + 1e-6
            expanded = centroid + (pts - centroid) * (1 + pad / (mean_d + 1e-6))
            try:
                hull = ConvexHull(expanded)
                hull_pts = expanded[hull.vertices]
                poly = plt.Polygon(hull_pts, closed=True, facecolor=color,
                                    alpha=0.10, edgecolor=color, linewidth=0.8,
                                    zorder=0)
                ax.add_patch(poly)
            except Exception:
                pass
        if len(ns) >= min_label_size:
            t = ax.text(centroid[0], centroid[1] + 0.02, g,
                        fontsize=CLUSTER_LABEL_FONTSIZE,
                        color="black", fontweight="bold", ha="center",
                        va="center", zorder=5,
                        bbox=dict(boxstyle="round,pad=0.18", facecolor="white",
                                  edgecolor="none", alpha=0.7))
            label_texts.append(t)


def draw_figure(G, pos, nodes_subset, cluster_groups_by_cat, singleton_nodes,
                 title, output: Path, category_titles=False, figsize=(16, 14),
                 y_ground=None, pls_shelf_y=None, single_cat=None,
                 legend_loc_bbox=None):
    meta = nodes_subset.set_index("NodeID")
    node_order = list(G.nodes())
    fills, borders, bwidths, sizes = node_visuals(nodes_subset, node_order)

    fig, ax = plt.subplots(figsize=figsize, facecolor="white")
    ax.set_facecolor("white")
    ax.axis("off")

    draw_edges_with_shadow(ax, G, pos)

    label_texts = []
    for cat, cluster_groups in cluster_groups_by_cat.items():
        draw_group_hulls(ax, cluster_groups, pos, CATEGORY_COLOR.get(cat, "#888888"),
                          meta, MIN_LABEL_SIZE.get(cat, 2), label_texts)

    nx.draw_networkx_nodes(G, pos, nodelist=node_order, ax=ax, node_color=fills,
                            node_size=sizes, edgecolors=borders,
                            linewidths=bwidths, alpha=0.9)
    ax.autoscale_view()

    if y_ground is not None:
        x0, x1 = ax.get_xlim()
        ax.axhline(y=y_ground - 0.05, color="#cccccc", linewidth=0.8,
                   linestyle="--", zorder=0)
        ax.text((x0 + x1) / 2, y_ground - 0.12, "singletons (no similarity partner)",
                fontsize=8, color="#888888", ha="center", style="italic")

    if pls_shelf_y is not None:
        ax.text(CATEGORY_CENTERS["PLS"][0], pls_shelf_y + 0.12,
                "plasmids (singletons)", fontsize=7.5, color="#888888",
                ha="center", style="italic")

    if HAS_ADJUST_TEXT and label_texts:
        adjust_text(label_texts, ax=ax, expand=(1.2, 1.3),
                    force_text=(0.3, 0.4), force_points=(0.1, 0.2))

    if category_titles:
        for cat, center in CATEGORY_CENTERS.items():
            cat_nodes = [n for n in node_order if meta.loc[n, "Category"] == cat]
            if not cat_nodes:
                continue
            top_y = max(pos[n][1] for n in cat_nodes) + 0.15
            ax.text(center[0], top_y, CATEGORY_NAME.get(cat, cat),
                    color=CATEGORY_COLOR.get(cat, "black"), fontsize=13,
                    fontweight="bold", ha="center",
                    bbox=dict(boxstyle="round,pad=0.25", facecolor="white",
                              edgecolor="none", alpha=0.85))

    present_cats = set(meta["Category"].unique())
    legend_elements = [
        mpatches.Patch(color="none", label="-- Category / classification --"),
    ]
    if "PPH" in present_cats:
        legend_elements.append(mpatches.Patch(color=CATEGORY_COLOR["PPH"],
                                label="Prophage (100% Caudoviricetes, intact)"))
    if "PLS" in present_cats:
        legend_elements.append(mpatches.Patch(color=CATEGORY_COLOR["PLS"],
                                label="Plasmid (100% conjugative)"))
    if "GI" in present_cats:
        legend_elements += [
            mpatches.Patch(color=GI_CLASS_COLOR["virulence"], label="GI: virulence"),
            mpatches.Patch(color=GI_CLASS_COLOR["resistance"], label="GI: resistance"),
            mpatches.Patch(color=GI_CLASS_COLOR["metabolic"], label="GI: metabolic"),
            mpatches.Patch(color=GI_CLASS_COLOR["symbiotic"], label="GI: symbiotic"),
            mpatches.Patch(color=GI_CLASS_COLOR["unknown"], label="GI: unassigned function"),
        ]
    if "IS" in present_cats:
        for fam, col in IS_FAMILY_COLOR.items():
            legend_elements.append(mpatches.Patch(color=col, label=f"IS family: {fam}"))

    legend_elements += [
        mpatches.Patch(color="none", label=" "),
        mpatches.Patch(color="none", label="-- Node border (carries internal IS) --"),
    ]
    if "IS" in present_cats and any(cat in present_cats for cat in ("PPH", "GI", "PLS")):
        # Fix #5: single explanatory line instead of 3 flat host-category
        # swatches -- the ring color now reuses the IS family palette above.
        legend_elements.append(
            mpatches.Patch(color="none",
                            label="Ring color = family of the internal IS\n(same palette as IS fill above)"))
    else:
        legend_elements += [
            mpatches.Patch(facecolor="none", edgecolor=BORDER_CARRIER_BY_HOST["PPH"],
                            linewidth=2.4, label="IS inside a prophage"),
            mpatches.Patch(facecolor="none", edgecolor=BORDER_CARRIER_BY_HOST["GI"],
                            linewidth=2.4, label="IS inside a genomic island"),
            mpatches.Patch(facecolor="none", edgecolor=BORDER_CARRIER_BY_HOST["PLS"],
                            linewidth=2.4, label="IS inside a plasmid"),
        ]
    if single_cat is not None:
        legend_elements.append(mpatches.Patch(color="none", label="-- Node size (length, bp) --"))
    else:
        legend_elements.append(mpatches.Patch(color="none",
                                label="-- Node size (log-scaled, RELATIVE per category) --"))
    legend_elements += size_legend_handles(nodes_subset, single_cat=single_cat)
    if single_cat is None:
        legend_elements.append(mpatches.Patch(color="none", label="   actual bp range per category:"))
        legend_elements += category_size_range_handles(nodes_subset)

    loc, bbox = legend_loc_bbox if legend_loc_bbox is not None else LEGEND_LOC_STANDALONE
    leg = ax.legend(handles=legend_elements, loc=loc,
                     bbox_to_anchor=bbox, fontsize=8.0, framealpha=1.0,
                     facecolor="white", edgecolor="#aaaaaa", labelcolor="black",
                     title="Legend", title_fontsize=9, labelspacing=1.05,
                     handletextpad=0.8)
    for handle, text in zip(leg.legend_handles, leg.get_texts()):
        if text.get_text().startswith("--") or text.get_text().strip() == "":
            handle.set_visible(False)

    ax.set_title(title, color="black", fontsize=14, pad=15)
    plt.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"[INFO] Saved: {output} ({output.stat().st_size} bytes) | "
          f"{G.number_of_nodes()} nodes, {G.number_of_edges()} edges, "
          f"{len(singleton_nodes)} singletons")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    nodes, edges = load_data()
    meta_all = nodes.set_index("NodeID")

    standalone_specs = [
        ("PPH", "ssn_prophages_v10.png",           "Prophage SSN (high-confidence)"),
        ("GI",  "ssn_genomic_islands_v10.png",      "Genomic island SSN (high-confidence)"),
        ("PLS", "ssn_plasmids_v10.png",             "Plasmid SSN (high-confidence)"),
        ("IS",  "ssn_insertion_sequences_v10.png",  "Insertion sequence SSN (GBK-curated, host-overlap excluded)"),
    ]

    for cat, fname, title in standalone_specs:
        cat_nodes_df = nodes[nodes["Category"] == cat]
        if cat_nodes_df.empty:
            print(f"[WARN] No nodes for {cat}, skipping.")
            continue
        cat_edges_full = edges[
            edges["Source"].str.startswith(cat) & edges["Target"].str.startswith(cat)
        ]
        full_G, _ = build_graph(nodes, cat_edges_full, categories=[cat])

        cat_nodes = list(cat_nodes_df["NodeID"])
        pos, singleton_nodes, centers, radii, cluster_groups = compute_nested_layout(
            cat_nodes, meta_all, full_G)

        if pos:
            all_y = [p[1] for p in pos.values()]
            y_ground = min(all_y) - 1.2
        else:
            y_ground = -1.2
        all_x_vals = [p[0] for p in pos.values()] if pos else [0.0]
        x_span = max(1.0, (max(all_x_vals) - min(all_x_vals)) / 2 if len(all_x_vals) > 1 else 1.0)
        pos.update(place_singleton_rows(singleton_nodes, y_ground, -x_span, x_span,
                                         n_rows=N_SINGLETON_ROWS,
                                         jitter_seed=hash(cat) % 1000))

        draw_edges_G, sub_nodes = build_graph(
            nodes, top_k_per_node(cat_edges_full, TOP_K_EDGES_STANDALONE.get(cat)),
            categories=[cat])
        draw_edges_G.add_nodes_from(full_G.nodes())

        draw_figure(draw_edges_G, pos, cat_nodes_df,
                    {cat: cluster_groups}, singleton_nodes, title,
                    PLOTS_DIR / fname, y_ground=y_ground if singleton_nodes else None,
                    single_cat=cat)

    # ---------------- Overview ----------------
    all_pos = {}
    all_cluster_groups = {}
    singletons_by_cat = {}
    node_frames = []

    for cat in CATEGORY_COLOR:
        cat_nodes_df = nodes[nodes["Category"] == cat]
        if cat_nodes_df.empty:
            continue
        node_frames.append(cat_nodes_df)
        cat_edges_full = edges[
            edges["Source"].str.startswith(cat) & edges["Target"].str.startswith(cat)
        ]
        full_G, _ = build_graph(nodes, cat_edges_full, categories=[cat])
        cat_nodes = list(cat_nodes_df["NodeID"])
        pos, singleton_nodes, centers, radii, cluster_groups = compute_nested_layout(
            cat_nodes, meta_all, full_G)

        if pos:
            coords = np.array(list(pos.values()))
            maxabs = np.abs(coords).max() or 1.0
            scale = CATEGORY_AREA_SCALE.get(cat, 0.5) / maxabs
            pos = {n: c * scale + CATEGORY_CENTERS[cat] for n, c in pos.items()}
        all_pos.update(pos)
        all_cluster_groups[cat] = cluster_groups
        singletons_by_cat[cat] = singleton_nodes

    all_nodes_df = pd.concat(node_frames, ignore_index=True)

    # Shared ground rows for PPH/GI/IS singletons; PLS keeps its own single
    # elevated shelf (only ~4 plasmids -- a second row would be overkill).
    shared_singletons = []
    for cat in ("PPH", "GI", "IS"):
        shared_singletons.extend(singletons_by_cat.get(cat, []))
    pls_singletons = singletons_by_cat.get("PLS", [])

    if all_pos:
        y_ground = min(p[1] for p in all_pos.values()) - 0.8
    else:
        y_ground = -2.5
    # v9: pls_shelf_y is the midpoint between the v7 position (near the
    # ground/legend) and the v8 fixed value (0.0, which turned out to
    # slightly clip the GI hull's lower edge) -- computed from THIS run's
    # real y_ground, not a re-guessed constant.
    pls_shelf_y_v7 = y_ground + (N_SINGLETON_ROWS - 1) * SINGLETON_ROW_GAP + 0.55
    pls_shelf_y = (pls_shelf_y_v7 + PLS_SHELF_Y_V8) / 2.0
    print(f"[INFO] PLS shelf y: v7={pls_shelf_y_v7:.3f}, v8={PLS_SHELF_Y_V8:.3f}, "
          f"v9 (midpoint)={pls_shelf_y:.3f}")

    x_vals = [p[0] for p in all_pos.values()] if all_pos else [0.0]
    x_span = max(1.5, (max(x_vals) - min(x_vals)) / 2 if len(x_vals) > 1 else 1.5)
    # v7: reserve the bottom-right quadrant (x > GROUND_X_MAX_FRAC*x_span) for
    # the legend -- that's where PLS already lives with only 4 singletons and
    # no clusters, confirmed empty otherwise. Shared PPH/GI/IS singletons are
    # confined to a left-biased range and spread over more rows to compensate.
    all_pos.update(place_singleton_rows(shared_singletons, y_ground,
                                         -x_span, x_span * GROUND_X_MAX_FRAC,
                                         n_rows=N_SINGLETON_ROWS_OVERVIEW, jitter_seed=99))
    # PLS shelf: narrow strip centered on the PLS macro-position, not spanning
    # the whole figure width, so it reads as "belonging" to the plasmid area
    pls_x_span = 0.4
    pls_pos = place_singleton_row(pls_singletons, pls_shelf_y, pls_x_span, jitter_seed=7)
    pls_pos = {n: p + np.array([CATEGORY_CENTERS["PLS"][0], 0.0]) for n, p in pls_pos.items()}
    all_pos.update(pls_pos)

    all_singletons = shared_singletons + pls_singletons

    overview_edge_parts = []
    for cat in CATEGORY_COLOR:
        cat_edges = edges[
            edges["Source"].str.startswith(cat) & edges["Target"].str.startswith(cat)
        ]
        overview_edge_parts.append(top_k_per_node(cat_edges, TOP_K_EDGES_OVERVIEW.get(cat)))
    overview_edges = pd.concat(overview_edge_parts, ignore_index=True)

    draw_edges_G, _ = build_graph(nodes, overview_edges, categories=None)
    draw_edges_G.add_nodes_from(all_nodes_df["NodeID"].tolist())

    draw_figure(draw_edges_G, all_pos, all_nodes_df, all_cluster_groups,
                all_singletons, "MGE similarity network overview (all categories)",
                PLOTS_DIR / "ssn_overview_v10.png", category_titles=True,
                figsize=(22, 20), y_ground=y_ground if shared_singletons else None,
                pls_shelf_y=pls_shelf_y if pls_singletons else None,
                legend_loc_bbox=LEGEND_LOC_OVERVIEW)

    print("[DONE] All v10 visualizations saved to:", PLOTS_DIR)


if __name__ == "__main__":
    main()