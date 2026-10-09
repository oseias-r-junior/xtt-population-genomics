#!/usr/bin/env python3
__version__ = "1.0.0"
__author__ = "Oseias Rodrigues Feitosa Junior"

"""
migrate_defense_naming.py

Safely migrate Defense Systems naming conventions.

This script ONLY updates:

1) filenames written/read by scripts
2) figure names
3) output table names
4) user-facing text

It NEVER modifies:

- df["subtype"]
- df['subtype']
- dataframe column names
- variable names
- Python identifiers

Default = dry-run

Usage
-----

Preview

    python migrate_defense_naming.py .

Apply

    python migrate_defense_naming.py . --apply

Apply + backups

    python migrate_defense_naming.py . --apply --backup


Verify that no deprecated names remain

    python migrate_defense_naming.py . --verify
"""

from pathlib import Path
import argparse
import shutil
import re
import sys

###############################################################################
# SAFE REPLACEMENTS
###############################################################################

SAFE_REPLACEMENTS = {

    #
    # Matrix
    #
    "defensefinder_subtype_matrix.tsv":
        "defensefinder_matrix.tsv",

    #
    # Summary
    #
    "defense_subtype_prevalence.tsv":
        "defense_prevalence.tsv",

    "defense_subtype_prevalence.png":
        "defense_prevalence.png",

    #
    # Correlation
    #
    "defense_subtype_correlation.tsv":
        "defense_correlation.tsv",

    "defense_subtype_correlation_heatmap.png":
        "defense_correlation_heatmap.png",

    #
    # Cooccurrence
    #
    "defense_subtype_cooccurrence.tsv":
        "defense_cooccurrence.tsv",

    "defense_subtype_cooccurrence_significant.tsv":
        "defense_cooccurrence_significant.tsv",

    "defense_subtype_network.png":
        "defense_cooccurrence_network.png",

    #
    # Network outputs
    #
    "defense_network_edges.tsv":
        "defense_cooccurrence_edges.tsv",

    "defense_network_centrality.tsv":
        "defense_cooccurrence_centrality.tsv",

    #
    # IS
    #
    "defense_ISfamily_clade_corrected.tsv":
        "defense_is_family_clade_corrected.tsv",

    "defense_is_family_association.tsv":
        "defense_vs_is_family.tsv",

    #
    # Prophages
    #
    "prophage_defense_clade_corrected.tsv":
        "defense_vs_prophage_clade_corrected.tsv",

    #
    # Plot labels
    #
    "Defense subtype":
        "Defense system",

    "Defense subtypes":
        "Defense systems",

    "unique defense subtypes":
        "unique defense systems",

    "Defense subtype co-occurrence network":
        "Defense system co-occurrence network",
}

###############################################################################
# DEPRECATED TERMS
###############################################################################

DEPRECATED_TERMS = [
    "defensefinder_subtype_matrix.tsv",
    "defense_subtype_prevalence",
    "defense_subtype_correlation",
    "defense_subtype_cooccurrence",
    "defense_subtype_cooccurrence_significant",
    "defense_subtype_network",
    "defense_network_edges",
    "defense_network_centrality",
    "defense_ISfamily",
    "prophage_defense_clade_corrected",
]

###############################################################################
# PROTECTED PATTERNS
###############################################################################

PROTECTED = [

    r'df\["subtype"\]',
    r"df\['subtype'\]",

    r'\["subtype"\]',
    r"\['subtype'\]",

    r"\.loc\[.*subtype",

]

###############################################################################

TOKEN = "__PROTECTED_SUBTYPE__"

###############################################################################


def protect(text):

    matches = []

    for pattern in PROTECTED:

        while True:

            m = re.search(pattern, text)

            if m is None:
                break

            matches.append(m.group(0))

            text = (
                text[:m.start()]
                + TOKEN
                + text[m.end():]
            )

    return text, matches


###############################################################################

def restore(text, matches):

    for m in matches:

        text = text.replace(
            TOKEN,
            m,
            1
        )

    return text


###############################################################################

def migrate(text):

    protected_text, matches = protect(text)

    changes = []

    for old, new in SAFE_REPLACEMENTS.items():

        n = protected_text.count(old)

        if n:

            protected_text = protected_text.replace(
                old,
                new
            )

            changes.append((old, new, n))

    protected_text = restore(
        protected_text,
        matches
    )

    return protected_text, changes


###############################################################################

def process(path, apply=False, backup=False):

    original = path.read_text(
        encoding="utf-8"
    )

    updated, changes = migrate(original)

    if not changes:
        return False, 0

    print("=" * 80)
    print(path)

    total = 0

    for old, new, n in changes:

        total += n

        print(f"  {old}")
        print(f"      -> {new}")
        print(f"      ({n} replacement{'s' if n > 1 else ''})")

    if apply:

        if backup:

            shutil.copy2(
                path,
                path.with_suffix(
                    path.suffix + ".bak"
                )
            )

        path.write_text(
            updated,
            encoding="utf-8"
        )

    return True, total


###############################################################################

def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "directory",
        nargs="?",
        default=".",
        type=Path
    )

    parser.add_argument(
        "--apply",
        action="store_true"
    )

    parser.add_argument(
        "--backup",
        action="store_true"
    )

    parser.add_argument(
        "--verify",
        action="store_true",
        help="Verify that no deprecated names remain."
    )

    args = parser.parse_args()

    # --------------------------------------------------
    # Verify directory
    # --------------------------------------------------

    if not args.directory.exists():
        print("Directory not found.")
        sys.exit(1)

    # --------------------------------------------------
    # Verify mode
    # --------------------------------------------------

    if args.verify:
        ok = verify(args.directory)
        sys.exit(0 if ok else 1)

    # --------------------------------------------------
    # Normal migration mode
    # --------------------------------------------------

    files = sorted(args.directory.rglob("*.py"))

    modified = 0
    replacements = 0

###############################################################################
# VERIFY MODE
###############################################################################

def verify(directory):

    print("=" * 80)
    print("VERIFYING DEFENSE NAMING")
    print("=" * 80)

    failed = False

    files = sorted(directory.rglob("*.py"))

    for path in files:

        text = path.read_text(encoding="utf-8")

        for term in DEPRECATED_TERMS:

            if term in text:

                failed = True

                print(f"\n{path}")

                for i, line in enumerate(text.splitlines(), start=1):

                    if term in line:

                        print(f"  line {i}: {term}")

    print("\n" + "=" * 80)

    if failed:
        print("FAILED")
        print("Deprecated names still exist.")
        return False

    print("PASSED")
    print("No deprecated names found.")
    return True
###############################################################################

if __name__ == "__main__":
    main()