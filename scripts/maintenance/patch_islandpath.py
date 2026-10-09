#!/usr/bin/env python3
"""
patch_islandpath.py

Apply compatibility patches to IslandPath installations.

Currently implemented patches
-----------------------------

1. Remove the obsolete Bio::Perl import from GenomeUtils.pm.

Recent BioPerl releases distributed through Bioconda no longer include
Bio::Perl.pm, although IslandPath still imports it.

The remaining Bio::* modules (Bio::SeqIO, Bio::Seq, etc.) are sufficient
for IslandPath execution.

The patch is idempotent:
- creates a backup (*.bak)
- does nothing if already patched
- prints informative messages
"""

from pathlib import Path
import argparse
import shutil
import sys

PATCH_LINE = "use Bio::Perl;"


def locate_genomeutils(env_prefix: Path) -> Path:
    """
    Locate GenomeUtils.pm inside a Micromamba/Conda environment.
    """

    candidates = list(env_prefix.rglob("GenomeUtils.pm"))

    if not candidates:
        raise FileNotFoundError(
            "GenomeUtils.pm could not be found inside the supplied environment."
        )

    if len(candidates) > 1:
        print("Warning: multiple GenomeUtils.pm files found.")
        print("Using:", candidates[0])

    return candidates[0]


def apply_patch(file_path: Path) -> None:
    """
    Remove obsolete 'use Bio::Perl;' import.
    """

    text = file_path.read_text()

    if PATCH_LINE not in text:
        print("✓ Patch already applied.")
        return

    backup = file_path.with_suffix(file_path.suffix + ".bak")

    if not backup.exists():
        shutil.copy2(file_path, backup)
        print(f"Backup created: {backup}")

    new_text = text.replace(PATCH_LINE + "\n", "")

    file_path.write_text(new_text)

    print(f"Patched: {file_path}")


def main():

    parser = argparse.ArgumentParser(
        description="Patch IslandPath for compatibility with modern BioPerl."
    )

    parser.add_argument(
        "--env",
        required=True,
        type=Path,
        help="Path to Micromamba/Conda environment.",
    )

    args = parser.parse_args()

    if not args.env.exists():
        sys.exit(f"Environment not found:\n{args.env}")

    genomeutils = locate_genomeutils(args.env)

    print(f"GenomeUtils.pm: {genomeutils}")

    apply_patch(genomeutils)

    print("\nDone.")


if __name__ == "__main__":
    main()