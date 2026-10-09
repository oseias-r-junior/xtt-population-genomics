#!/usr/bin/env python3
"""
check_chromosome_hypothesis.py -- revised hypothesis after check_concat_hypothesis.py
ruled out simple concatenation: T3E/TSS query_length appears to equal the
length of the SINGLE LARGEST contig (the main chromosome) in the genome
assembly, not the sum of all contigs. I.e. the BLAST query for each strain
was the chromosome contig alone (plasmids / small contigs excluded).

If true: q.start/q.end in the T3E/TSS sheets are ALREADY coordinates
relative to that one chromosome contig -- no offset reconstruction needed
at all. We just need, per genome, which contig IS "the chromosome" (the
one whose FASTA length matches query_length), so we know which contig any
candidate overlapping MGE region must also be on.

This script checks, for ALL matched strains (not just a sample):
  1. query_length == max(contig lengths) ?
  2. is that max unique (no tie between two contigs of the same length)?
  3. reports any strain where either check fails, so we know exactly how
     general the hypothesis is before building the real overlap script.

Usage
-----
python check_chromosome_hypothesis.py
"""
from pathlib import Path
import sys

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if not (PROJECT_ROOT / "metadata").exists():
    alt = Path(__file__).resolve().parents[3]
    if (alt / "metadata").exists():
        PROJECT_ROOT = alt

XLSX_FILE = PROJECT_ROOT / "metadata" / "Xtt_2021_2023_T3E_TSS (positions).xlsx"

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
    return found


def parse_fasta_contigs(path: Path):
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

    fasta_files = find_fasta_files()
    matched_strains = sorted(s for s in qlen_by_strain.index if s in fasta_files)
    print(f"[INFO] Testing {len(matched_strains)} strains: "
          f"query_length == max(contig_length), and max is unique\n")

    n_ok, n_fail_value, n_fail_unique = 0, 0, 0
    chrom_contig_id = {}  # strain -> contig id identified as the chromosome

    for s in matched_strains:
        contigs = parse_fasta_contigs(fasta_files[s])
        lengths = [l for _, l in contigs]
        max_len = max(lengths)
        n_at_max = lengths.count(max_len)
        qlen = int(qlen_by_strain[s])

        if max_len != qlen:
            n_fail_value += 1
            print(f"  [FAIL:value] {s}: max(contig_len)={max_len}  query_length={qlen}  "
                  f"n_contigs={len(contigs)}  lengths={lengths}")
            continue
        if n_at_max > 1:
            n_fail_unique += 1
            print(f"  [FAIL:unique] {s}: {n_at_max} contigs tied at max length {max_len} "
                  f"(ambiguous which is 'the chromosome') lengths={lengths}")
            continue
        n_ok += 1
        chrom_id = [cid for cid, l in contigs if l == max_len][0]
        chrom_contig_id[s] = chrom_id

    print(f"\n[SUMMARY] OK: {n_ok}/{len(matched_strains)}   "
          f"value-mismatch: {n_fail_value}   ambiguous-tie: {n_fail_unique}")

    if n_ok == len(matched_strains):
        print("\n[RESULT] CONFIRMED for all matched strains: query_length == length of the "
              "single largest contig, and that contig is unique. The T3E/TSS q.start/q.end "
              "coordinates are directly usable against that one contig per genome -- no "
              "offset math needed. Any MGE node on a DIFFERENT contig in the same genome "
              "can never overlap a T3E/TSS hit (correctly, since T3E/TSS effectors here were "
              "only searched for on the chromosome).")
        print("\n[INFO] Example chromosome-contig IDs identified (first 5):")
        for s in list(chrom_contig_id)[:5]:
            print(f"    {s} -> {chrom_contig_id[s]}")
    else:
        print("\n[RESULT] Hypothesis holds for most but not all strains -- inspect the "
              "FAIL lines above before proceeding; those genomes may need a different rule "
              "or manual exclusion from the T3E/TSS x MGE overlap check.")


if __name__ == "__main__":
    main()
