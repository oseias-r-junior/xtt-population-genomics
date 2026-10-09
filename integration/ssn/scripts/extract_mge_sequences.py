#!/usr/bin/env python3
"""
integration/ssn/scripts/extract_mge_sequences.py

Extract nucleotide sequences for every consensus MGE instance (prophages,
genomic islands, plasmids, insertion sequences) into a single combined
multi-FASTA, for the all-vs-all BLASTN Sequence Similarity Network (SSN)
mirroring Uceda-Campos et al. 2022 (Microorganisms) Fig. 5 methodology.

Sources:
    PPH  -- mobilome/prophages/matrices/prophage_consensus_regions.tsv
             (Genome, Start, Stop) + annotations/fna/{genome}.fna
    GI   -- mobilome/genomic_islands/matrices/gi_consensus_regions.tsv
             (Genome, Start, End) + annotations/fna/{genome}.fna
    PLS  -- mobilome/plasmids/matrices/plasmid_consensus.tsv
             (Genome, Contig_id) -- whole contig, looked up directly by
             header in annotations/fna/{genome}.fna (no coordinate ambiguity)
    IS   -- mobilome/insertion_sequences/isescan/results/{genome}/fna/
             {genome}.fna.is.fna -- ISEScan's OWN per-instance sequences,
             already correctly extracted by ISEScan itself; just relabeled.

CONTIG RESOLUTION (PPH/GI only): the consensus scripts for prophages and
genomic islands track (Genome, Start, End) but NOT which contig the region
came from (a known limitation of those scripts, harmless for their own
clustering logic but NOT harmless here, where the actual sequence must be
sliced correctly). ~33% of genomes (37/113) have >1 contig. Resolution
heuristic: a region's End coordinate must fit within its source contig, so
any contig shorter than End cannot be the right one. In near-complete
bacterial assemblies contigs differ substantially in size (one large
chromosome + several small leftover contigs), so this usually identifies
the contig uniquely. When 0 or >1 contigs qualify, the region is AMBIGUOUS
and is skipped (logged), rather than risk extracting the wrong sequence.
"""

from pathlib import Path
import sys

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parents[2]

sys.path.append(str(PROJECT_ROOT / "scripts"))

from config import *

from Bio import SeqIO
import pandas as pd

# ======================================================
# INPUT
# ======================================================

# 2026-09-28: switched to the conservative "all-tools-agree" high-confidence
# subsets (see mobilome/scripts/split_consensus_confidence.py) and the
# GBK-curated IS instance list (see
# mobilome/insertion_sequences/scripts/curate_is_against_gbk.py) rather
# than the full consensus tables -- per PI decision, only MGEs supported
# by every tool used for that category (and, for IS, independently
# confirmed by GBK annotation) feed the main discussion/network. Partial
# and unconfirmed regions remain available in the *_supplementary.tsv /
# is_gbk_curation_log.tsv files but are intentionally excluded here.
PPH_FILE = PROPHAGE_MATRICES / "prophage_high_confidence.tsv"
GI_FILE = GI_MATRICES / "gi_high_confidence.tsv"
PLS_FILE = PLASMID_MATRICES / "plasmid_high_confidence.tsv"
IS_CONFIRMED_FILE = IS_MATRICES / "is_gbk_confirmed_instances.tsv"

# ======================================================
# OUTPUT
# ======================================================

SSN_FASTA.mkdir(parents=True, exist_ok=True)
SSN_MATRICES.mkdir(parents=True, exist_ok=True)

OUT_FASTA = SSN_FASTA / "all_mge_sequences.fna"
OUT_NODES = SSN_MATRICES / "ssn_nodes.tsv"

# ======================================================
# GENOME FASTA CACHE
# ======================================================

_genome_records_cache = {}


def load_genome_contigs(genome: str):
    """Return dict of {contig_id: SeqRecord} for a genome, cached."""
    if genome not in _genome_records_cache:
        fna_file = FNA / f"{genome}.fna"
        if not fna_file.exists():
            _genome_records_cache[genome] = {}
        else:
            _genome_records_cache[genome] = {
                rec.id: rec for rec in SeqIO.parse(fna_file, "fasta")
            }
    return _genome_records_cache[genome]


def resolve_contig(genome: str, end_coord: int):
    """Find the unique contig whose length >= end_coord. Returns contig_id
    or None if zero or multiple candidates qualify (ambiguous)."""
    contigs = load_genome_contigs(genome)
    if len(contigs) == 1:
        return next(iter(contigs))

    candidates = [cid for cid, rec in contigs.items() if len(rec.seq) >= end_coord]
    if len(candidates) == 1:
        return candidates[0]
    return None


nodes = []
fasta_records = []
ambiguous_log = []


def add_region_node(category, genome, region_id, start, end, classification,
                     n_tools, tools_supporting, extra):
    contig_id = resolve_contig(genome, end)
    if contig_id is None:
        ambiguous_log.append((category, genome, region_id, start, end))
        return

    contigs = load_genome_contigs(genome)
    seq = contigs[contig_id].seq[start - 1:end]
    if len(seq) == 0:
        ambiguous_log.append((category, genome, region_id, start, end))
        return

    node_id = f"{category}|{genome}|{region_id}|{len(seq)}"

    fasta_records.append((node_id, str(seq)))
    nodes.append({
        "NodeID": node_id,
        "Category": category,
        "Genome": genome,
        "Contig": contig_id,
        "Start": start,
        "End": end,
        "Length_bp": len(seq),
        "Classification": classification,
        "N_Tools_Support": n_tools,
        "Tools_Supporting": tools_supporting,
        "Extra": extra,
    })


# ======================================================
# PROPHAGES
# ======================================================

pph = pd.read_csv(PPH_FILE, sep="\t")
for _, row in pph.iterrows():
    add_region_node(
        category="PPH",
        genome=row["Genome"],
        region_id=row["Region"],
        start=int(row["Start"]),
        end=int(row["Stop"]),
        classification=row.get("Classification", ""),
        n_tools=row.get("N_Tools_Support", ""),
        tools_supporting=row.get("Tools_Supporting", ""),
        extra=row.get("PHASTEST_Most_Common_Phage", ""),
    )

print(f"Prophages processed: {len(pph)}")

# ======================================================
# GENOMIC ISLANDS
# ======================================================

gi = pd.read_csv(GI_FILE, sep="\t")
for _, row in gi.iterrows():
    add_region_node(
        category="GI",
        genome=row["Genome"],
        region_id=row["Region"],
        start=int(row["Start"]),
        end=int(row["End"]),
        classification=row.get("Classification", ""),
        n_tools=row.get("N_Tools_Support", ""),
        tools_supporting=row.get("Tools_Supporting", ""),
        extra="",
    )

print(f"Genomic islands processed: {len(gi)}")

# ======================================================
# PLASMIDS (whole contig, looked up directly by header -- no ambiguity)
# ======================================================

pls = pd.read_csv(PLS_FILE, sep="\t")
pls_skipped = 0
for _, row in pls.iterrows():
    genome = row["Genome"]
    contig_id = row["Contig_id"]
    contigs = load_genome_contigs(genome)

    if contig_id not in contigs:
        pls_skipped += 1
        continue

    seq = contigs[contig_id].seq
    node_id = f"PLS|{genome}|{contig_id}|{len(seq)}"

    fasta_records.append((node_id, str(seq)))
    nodes.append({
        "NodeID": node_id,
        "Category": "PLS",
        "Genome": genome,
        "Contig": contig_id,
        "Start": 1,
        "End": len(seq),
        "Length_bp": len(seq),
        "Classification": row.get("Mobility_Class", ""),
        "N_Tools_Support": row.get("N_Tools_Support", ""),
        "Tools_Supporting": row.get("Tools_Supporting", ""),
        "Extra": row.get("Rep_type", ""),
    })

print(f"Plasmids processed: {len(pls)} ({pls_skipped} skipped -- contig header not found)")

# ======================================================
# INSERTION SEQUENCES (already extracted by ISEScan -- just relabel)
# ======================================================

is_confirmed = pd.read_csv(IS_CONFIRMED_FILE, sep="\t")

is_count = 0
is_skipped = 0
for idx, row in is_confirmed.iterrows():
    genome = row["Genome"]
    contig_id = row["Contig"]
    start = int(row["isBegin"])
    end = int(row["isEnd"])
    family = row["Family"]

    contigs = load_genome_contigs(genome)
    if contig_id not in contigs:
        is_skipped += 1
        continue

    seq = contigs[contig_id].seq[start - 1:end]
    if len(seq) == 0:
        is_skipped += 1
        continue

    node_id = f"IS|{genome}|{idx}|{len(seq)}"

    fasta_records.append((node_id, str(seq)))
    nodes.append({
        "NodeID": node_id,
        "Category": "IS",
        "Genome": genome,
        "Contig": contig_id,
        "Start": start,
        "End": end,
        "Length_bp": len(seq),
        "Classification": family,
        "N_Tools_Support": 1,
        "Tools_Supporting": "ISEScan+GBK_curated",
        "Extra": "",
    })
    is_count += 1

print(f"IS instances processed: {is_count} ({is_skipped} skipped -- contig not found in genome fna)")

# ======================================================
# WRITE OUTPUTS
# ======================================================

with open(OUT_FASTA, "w") as fh:
    for node_id, seq in fasta_records:
        fh.write(f">{node_id}\n{seq}\n")

nodes_df = pd.DataFrame(nodes)
nodes_df.to_csv(OUT_NODES, sep="\t", index=False)

print()
print(f"Total nodes/sequences written: {len(fasta_records)}")
print(f"Ambiguous regions skipped (contig could not be resolved): {len(ambiguous_log)}")
if ambiguous_log:
    for cat, genome, region, start, end in ambiguous_log[:20]:
        print(f"  SKIPPED: {cat} {genome} {region} ({start}-{end})")
    if len(ambiguous_log) > 20:
        print(f"  ... and {len(ambiguous_log) - 20} more")
print()
print("Nodes by category:")
print(nodes_df["Category"].value_counts().to_string())
print()
print(f"Saved: {OUT_FASTA}")
print(f"Saved: {OUT_NODES}")
