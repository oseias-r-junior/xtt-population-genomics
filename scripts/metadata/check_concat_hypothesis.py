#!/usr/bin/env python3
"""
check_concat_hypothesis.py -- validates the hypothesis that the T3E/TSS BLAST
query sequences were built by simple concatenation of each genome's contigs,
in FASTA record order, with NO spacer/gap inserted between them.

Logic: if that hypothesis is true, then for any genome:
    query_length (as recorded in the T3E/TSS spreadsheet)
      == sum(len(contig) for contig in genome_fasta, in file order)

This script:
  1. Locates each genome's assembly FASTA file (tries a few common
     locations/patterns used elsewhere in this project -- edit
     FASTA_CANDIDATES below if none match).
  2. For a sample of strains (a mix of single-contig and multi-contig
     assemblies, picked automatically), compares the FASTA's total length
     and per-contig lengths/order against the spreadsheet's query_length.
  3. If they match exactly for multi-contig genomes too, the concatenation
     hypothesis is confirmed and we additionally print the per-contig
     cumulative offsets for one multi-contig example, so you can sanity
     check the reconstruction logic before it's used for the full overlap
     script.

Usage
-----
python check_concat_hypothesis.py
(run from anywhere; PROJECT_ROOT is auto-detected relative to this script's
location under <project_root>/metadata/scripts/ or wherever you placed it --
edit PROJECT_ROOT below if it doesn't resolve correctly)
"""
from pathlib import Path
import sys

import pandas as pd

# ---- EDIT THIS if the script isn't placed 2 levels under project root ----
PROJECT_ROOT = Path(__file__).resolve().parents[2]
if not (PROJECT_ROOT / "metadata").exists():
    # fallback: try parents[2] in case this is under project_root/x/scripts/
    alt = Path(__file__).resolve().parents[3]
    if (alt / "metadata").exists():
        PROJECT_ROOT = alt

XLSX_FILE = PROJECT_ROOT / "metadata" / "Xtt_2021_2023_T3E_TSS (positions).xlsx"

# Common places genome assemblies might live in this project -- add more
# patterns here if none of these hit.
FASTA_GLOB_PATTERNS = [
    "genomes/*.fasta", "genomes/*.fa", "genomes/*.fna",
    "genomes/assemblies/*.fasta", "genomes/assemblies/*.fa", "genomes/assemblies/*.fna",
    "assemblies/*.fasta", "assemblies/*.fa", "assemblies/*.fna",
    "raw_genomes/*.fasta", "raw_genomes/*.fa", "raw_genomes/*.fna",
    "**/*.fasta", "**/*.fa", "**/*.fna",
]


def find_fasta_files():
    found = {}
    for pattern in FASTA_GLOB_PATTERNS:
        for f in PROJECT_ROOT.glob(pattern):
            stem = f.stem
            if stem not in found:
                found[stem] = f
        if found:
            # stop widening the search once the narrower patterns hit something
            if pattern != "**/*.fasta" and pattern != "**/*.fa" and pattern != "**/*.fna":
                continue
    return found


def parse_fasta_contigs(path: Path):
    """Returns list of (contig_id, length) in file order, no Biopython dependency."""
    contigs = []
    name = None
    length = 0
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    contigs.append((name, length))
                name = line[1:].split()[0]
                length = 0
            else:
                length += len(line.strip())
    if name is not None:
        contigs.append((name, length))
    return contigs


def main():
    if not XLSX_FILE.exists():
        print(f"[FATAL] Cannot find {XLSX_FILE}")
        sys.exit(1)

    t3e = pd.read_excel(XLSX_FILE, sheet_name="T3E")
    t3e.columns = [c.strip() for c in t3e.columns]
    qlen_by_strain = t3e.groupby("Strain")["query length"].first()

    print(f"[INFO] {len(qlen_by_strain)} strains with a recorded query_length in T3E sheet")

    fasta_files = find_fasta_files()
    print(f"[INFO] Found {len(fasta_files)} candidate FASTA files under {PROJECT_ROOT}")
    if not fasta_files:
        print("[FATAL] No FASTA files found with the patterns tried. Edit "
              "FASTA_GLOB_PATTERNS in this script to point at wherever the "
              "genome assemblies actually live, then re-run.")
        sys.exit(1)

    matched_strains = [s for s in qlen_by_strain.index if s in fasta_files]
    print(f"[INFO] {len(matched_strains)} strains have both a query_length and a matching FASTA file")
    if not matched_strains:
        print("[FATAL] No strain names matched between the spreadsheet and the FASTA "
              "filenames found. FASTA filenames may use a different naming convention "
              "than the 'Strain' column -- check manually.")
        sys.exit(1)

    # classify by contig count so we test at least one multi-contig genome
    n_contigs = {}
    for s in matched_strains:
        contigs = parse_fasta_contigs(fasta_files[s])
        n_contigs[s] = len(contigs)

    multi = [s for s in matched_strains if n_contigs[s] > 1]
    single = [s for s in matched_strains if n_contigs[s] == 1]
    print(f"[INFO] Of matched strains: {len(single)} single-contig, {len(multi)} multi-contig")

    sample = single[:3] + multi[:5]
    print(f"\n[INFO] Testing concatenation hypothesis on {len(sample)} sample strain(s):\n")

    all_ok = True
    example_multi_offsets = None
    for s in sample:
        contigs = parse_fasta_contigs(fasta_files[s])
        total_len = sum(length for _, length in contigs)
        qlen = int(qlen_by_strain[s])
        ok = (total_len == qlen)
        all_ok &= ok
        print(f"  {s}: n_contigs={len(contigs)}  sum(FASTA contig lengths)={total_len}  "
              f"query_length(spreadsheet)={qlen}  {'MATCH' if ok else 'MISMATCH'}")
        if not ok:
            print(f"    -> per-contig lengths (file order): {[l for _, l in contigs]}")
        if len(contigs) > 1 and example_multi_offsets is None and ok:
            example_multi_offsets = (s, contigs)

    print()
    if all_ok:
        print("[RESULT] All sampled strains MATCH exactly: sum of FASTA contig lengths "
              "(in file order) == query_length recorded in the spreadsheet.")
        print("[RESULT] This confirms simple concatenation, FASTA record order, no spacer.")
        if example_multi_offsets:
            s, contigs = example_multi_offsets
            print(f"\n[INFO] Example reconstructed offsets for multi-contig strain '{s}':")
            cum = 0
            for cid, length in contigs:
                print(f"    contig {cid}: local range [1, {length}] -> concatenated range "
                      f"[{cum + 1}, {cum + length}]")
                cum += length
    else:
        print("[RESULT] MISMATCH found -- the simple-concatenation-no-spacer hypothesis is "
              "WRONG for at least one sampled genome. Do not proceed with coordinate "
              "conversion until this is resolved (check with the author: was a spacer/gap "
              "inserted between contigs? Was a different contig order used for the BLAST "
              "query than the FASTA file order?).")


if __name__ == "__main__":
    main()