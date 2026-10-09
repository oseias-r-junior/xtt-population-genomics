"""genome_ids.py -- single source of truth for genome identifiers.

Several tools wrote genome names with different punctuation (AB22010-2 vs
AB22010_2, LG_2 vs L_G_2 vs "LG 2", Km8 vs XtKm8, Roary's "<id>_out").
`canon(x)` maps any of those spellings to the identifier used in
metadata/genome_clades_master.tsv, so tables from different tools join on the
same key. Unknown names are returned unchanged (stripped).

Place this file in the project root (next to palette.py). Scripts do:
    sys.path.insert(0, str(PROJECT_ROOT))
    from genome_ids import canon
"""
import re
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent
MASTER = ROOT / "metadata" / "genome_clades_master.tsv"

# file-name stems that differ from the master by more than punctuation
ALIASES = {"KM8": "XTKM8", "KM9": "XTKM9"}


def key(x):
    s = re.sub(r"_out$", "", str(x).strip())
    k = re.sub(r"[^A-Z0-9]", "", s.upper())
    return ALIASES.get(k, k)


_CANON = None


def _load():
    global _CANON
    if _CANON is None:
        ids = pd.read_csv(MASTER, sep="\t", usecols=["Genome"])["Genome"].astype(str).str.strip()
        _CANON = {}
        for g in ids:
            k = key(g)
            if k in _CANON and _CANON[k] != g:
                raise ValueError(f"Two master IDs collapse to the same key {k}: {_CANON[k]} / {g}")
            _CANON[k] = g
    return _CANON


def canon(x):
    """Return the master-table spelling of genome id x (or x itself if unknown)."""
    return _load().get(key(x), str(x).strip())
