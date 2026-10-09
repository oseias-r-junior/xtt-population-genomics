#!/usr/bin/env python3
"""
build_t3e_tss_mge_overlap.py -- checks whether T3E (type III effector) and
TSS (secretion system) gene hits, from Xtt_2021_2023_T3E_TSS (positions).xlsx,
are physically located inside an MGE region (task #24 / #55).

Coordinate system (validated by check_concat_hypothesis.py then
check_chromosome_hypothesis.py against the real assemblies): the BLAST
query for each strain was NOT the whole genome concatenated -- it was a
SINGLE CONTIG, specifically the one whose length exactly equals the
'query length' column recorded in the spreadsheet (confirmed for 98/98
matched strains; for 97 of them that contig also happens to be the largest
contig, but strain GR23016 is a case where it's the second-largest -- the
exact-length-match rule below handles that generally, with no special
case needed). This means q.start/q.end are ALREADY coordinates on that one
contig; no offset math is needed, but any MGE node on a DIFFERENT contig in
the same genome can never overlap (and correctly so).

Categories on the "virulence" side:
  - T3E sheet: every row is a single type-III effector gene hit -> one
    category, "T3E".
  - TSS sheet: 'subject id' values look like '<SYSTEM>_<gene>' (e.g.
    'T2SS_xcsE') -- the leading token before the first underscore is taken
    as the secretion-system subtype (T2SS/T3SS/T4SS/T6SS/...). Any row that
    doesn't match this pattern is bucketed as 'TSS_other' and logged.

Output:
  integration/clade_correlation/matrices/t3e_tss_mge_physical_overlap_long.tsv
  columns: Group, MGE_category, Virulence_category, n_genes, n_contained, pct_contained

Also writes a per-genome contig-resolution log
  integration/clade_correlation/matrices/t3e_tss_chromosome_contig_map.tsv
so the contig-matching step is auditable.

Usage
-----
python build_t3e_tss_mge_overlap.py
"""
from pathlib import Path
import re
import sys

import pandas as pd

# project root = first parent folder that contains metadata/ and integration/
# (works wherever the script is placed inside the repo)
PROJECT_ROOT = next(p for p in Path(__file__).resolve().parents
                    if (p / "metadata").is_dir() and (p / "integration").is_dir())
sys.path.insert(0, str(PROJECT_ROOT))
from genome_ids import canon  # noqa: E402  (unifies "LG 2" / LG_2 / L_G_2, MO22001-1 / MO22001_1 ...)

XLSX_FILE = PROJECT_ROOT / "metadata" / "Xtt_2021_2023_T3E_TSS (positions).xlsx"
MGE_NODES_FILE = PROJECT_ROOT / "integration" / "ssn" / "matrices" / "ssn_nodes_final.tsv"
CLADE_FILE = PROJECT_ROOT / "metadata" / "genome_clades_master.tsv"
OUT_DIR = PROJECT_ROOT / "integration" / "clade_correlation" / "matrices"
OUT_FILE = OUT_DIR / "t3e_tss_mge_physical_overlap_long.tsv"
CONTIG_MAP_FILE = OUT_DIR / "t3e_tss_chromosome_contig_map.tsv"

FASTA_GLOB_PATTERNS = [
    "genomes/*.fasta", "genomes/*.fa", "genomes/*.fna",
    "genomes/assemblies/*.fasta", "genomes/assemblies/*.fa", "genomes/assemblies/*.fna",
    "assemblies/*.fasta", "assemblies/*.fa", "assemblies/*.fna",
    "raw_genomes/*.fasta", "raw_genomes/*.fa", "raw_genomes/*.fna",
    "**/*.fasta", "**/*.fa", "**/*.fna",
]

MGE_CATEGORY_ORDER = ["PPH", "GI", "IS", "PLS"]
CLADE_ORDER = ["K0", "K1", "K2"]

TSS_PREFIX_RE = re.compile(r"^(T\d?SS)[_\-]", re.IGNORECASE)


def find_fasta_files():
    found = {}
    for pattern in FASTA_GLOB_PATTERNS:
        for f in PROJECT_ROOT.glob(pattern):
            stem = canon(f.stem)
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


def resolve_chromosome_contigs(strains, qlen_by_strain, fasta_files):
    """For each strain, find the contig whose FASTA length exactly equals
    the recorded query_length. Returns dict strain -> contig_id, and logs
    any strain that fails (no match, or >1 match)."""
    resolved = {}
    rows = []
    n_fail = 0
    for s in strains:
        if s not in fasta_files:
            print(f"[WARN] No FASTA found for strain '{s}' -- excluded.")
            rows.append({"Strain": s, "query_length": qlen_by_strain.get(s), "contig_id": None, "status": "no_fasta"})
            n_fail += 1
            continue
        contigs = parse_fasta_contigs(fasta_files[s])
        qlen = int(qlen_by_strain[s])
        matches = [cid for cid, l in contigs if l == qlen]
        if len(matches) == 1:
            resolved[s] = matches[0]
            rows.append({"Strain": s, "query_length": qlen, "contig_id": matches[0], "status": "ok"})
        elif len(matches) == 0:
            print(f"[FAIL] {s}: no contig with length == query_length ({qlen}); "
                  f"available lengths={[l for _, l in contigs]}. Excluded.")
            rows.append({"Strain": s, "query_length": qlen, "contig_id": None, "status": "no_length_match"})
            n_fail += 1
        else:
            print(f"[FAIL] {s}: {len(matches)} contigs tied at length {qlen} "
                  f"({matches}) -- ambiguous. Excluded.")
            rows.append({"Strain": s, "query_length": qlen, "contig_id": None, "status": "ambiguous"})
            n_fail += 1

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(CONTIG_MAP_FILE, sep="\t", index=False)
    print(f"[INFO] Chromosome-contig resolution: {len(resolved)}/{len(strains)} strains OK "
          f"({n_fail} excluded). Full map saved: {CONTIG_MAP_FILE}")
    return resolved


def load_virulence_genes():
    xl = pd.ExcelFile(XLSX_FILE)
    frames = []

    t3e = pd.read_excel(xl, sheet_name="T3E")
    t3e.columns = [c.strip() for c in t3e.columns]
    t3e["Virulence_category"] = "T3E"
    frames.append(t3e[["Strain", "q. start", "q. end", "Virulence_category"]]
                  .rename(columns={"q. start": "q_start", "q. end": "q_end"}))

    tss = pd.read_excel(xl, sheet_name="TSS")
    tss.columns = [c.strip() for c in tss.columns]

    def tss_category(subject_id):
        m = TSS_PREFIX_RE.match(str(subject_id).strip())
        return m.group(1).upper() if m else "TSS_other"

    tss["Virulence_category"] = tss["subject id"].apply(tss_category)
    other_labels = sorted(tss.loc[tss["Virulence_category"] == "TSS_other", "subject id"].unique().tolist())
    if other_labels:
        print(f"[INFO] TSS 'subject id' values not matching the T#SS_ prefix pattern "
              f"(bucketed as TSS_other, {len(other_labels)} unique): {other_labels[:20]}"
              f"{' ...' if len(other_labels) > 20 else ''}")
    cat_counts = tss["Virulence_category"].value_counts()
    print(f"[INFO] TSS category breakdown (row counts):\n{cat_counts.to_string()}")

    frames.append(tss[["Strain", "q. start", "q. end", "Virulence_category"]]
                  .rename(columns={"q. start": "q_start", "q. end": "q_end"}))

    genes = pd.concat(frames, ignore_index=True)
    genes["Strain"] = genes["Strain"].astype(str).map(canon)
    print(f"[INFO] Total virulence gene hits loaded (T3E + TSS): {len(genes)}")
    return genes


def load_mge_regions():
    nodes = pd.read_csv(MGE_NODES_FILE, sep="\t")
    nodes["Contig"] = nodes["Contig"].astype(str)
    by_genome_contig = {}
    for _, row in nodes.iterrows():
        key = (row["Genome"], row["Contig"])
        by_genome_contig.setdefault(key, []).append(
            (row["Start"], row["End"], row["Category"])
        )
    print(f"[INFO] MGE regions loaded: {len(nodes)} nodes across "
          f"{nodes['Genome'].nunique()} genomes")
    return by_genome_contig


def annotate_containment(genes: pd.DataFrame, mge_by_genome_contig: dict):
    for cat in MGE_CATEGORY_ORDER:
        genes[f"in_{cat}"] = False
    genes["in_Any_MGE"] = False

    for idx, row in genes.iterrows():
        key = (row["Genome"], row["Contig"])
        regions = mge_by_genome_contig.get(key, [])
        g_start, g_end = row["Start"], row["End"]
        hit_any = False
        for m_start, m_end, m_cat in regions:
            overlap = min(g_end, m_end) - max(g_start, m_start)
            if overlap > 0:
                genes.at[idx, f"in_{m_cat}"] = True
                hit_any = True
        genes.at[idx, "in_Any_MGE"] = hit_any
    return genes


def summarize(genes: pd.DataFrame, group_label: str, vir_categories: list):
    rows = []
    for vc in vir_categories:
        sub = genes[genes["Virulence_category"] == vc]
        n = len(sub)
        for mge_cat in MGE_CATEGORY_ORDER + ["Any_MGE"]:
            if n == 0:
                pct, n_contained = float("nan"), 0
            else:
                n_contained = int(sub[f"in_{mge_cat}"].sum())
                pct = 100.0 * n_contained / n
            rows.append({
                "Group": group_label, "MGE_category": mge_cat, "Virulence_category": vc,
                "n_genes": n, "n_contained": n_contained, "pct_contained": pct,
            })
    return pd.DataFrame(rows)


def main():
    if not XLSX_FILE.exists():
        print(f"[FATAL] Cannot find {XLSX_FILE}")
        sys.exit(1)

    genes = load_virulence_genes()
    fasta_files = find_fasta_files()

    qlen_file = pd.read_excel(XLSX_FILE, sheet_name="T3E")
    qlen_file.columns = [c.strip() for c in qlen_file.columns]
    qlen_file["Strain"] = qlen_file["Strain"].astype(str).map(canon)
    qlen_by_strain = qlen_file.groupby("Strain")["query length"].first()
    # TSS sheet may have extra strains T3E doesn't -- fold those in too
    tss_qlen = pd.read_excel(XLSX_FILE, sheet_name="TSS")
    tss_qlen.columns = [c.strip() for c in tss_qlen.columns]
    tss_qlen["Strain"] = tss_qlen["Strain"].astype(str).map(canon)
    tss_qlen_by_strain = tss_qlen.groupby("Strain")["query length"].first()
    qlen_by_strain = qlen_by_strain.combine_first(tss_qlen_by_strain)

    strains = sorted(genes["Strain"].unique())
    chrom_contig = resolve_chromosome_contigs(strains, qlen_by_strain, fasta_files)

    genes = genes[genes["Strain"].isin(chrom_contig)].copy()
    genes["Genome"] = genes["Strain"]
    genes["Contig"] = genes["Strain"].map(chrom_contig)
    genes["Start"] = genes[["q_start", "q_end"]].min(axis=1)
    genes["End"] = genes[["q_start", "q_end"]].max(axis=1)

    mge_by_genome_contig = load_mge_regions()
    genes = annotate_containment(genes, mge_by_genome_contig)

    clades = pd.read_csv(CLADE_FILE, sep="\t", usecols=["Genome", "Clade"])
    clade_map = dict(zip(clades["Genome"], clades["Clade"]))
    genes["Clade"] = genes["Genome"].map(clade_map)

    no_clade = genes.loc[genes["Clade"].isna(), "Genome"].unique().tolist()
    if no_clade:
        print(f"[WARN] {len(no_clade)} genome(s) with virulence-gene hits have no clade "
              f"assignment, excluded: {no_clade}")
    non_k_mask = genes["Clade"].notna() & ~genes["Clade"].isin(CLADE_ORDER)
    if non_k_mask.any():
        excluded_genomes = sorted(genes.loc[non_k_mask, "Genome"].unique().tolist())
        excluded_clades = sorted(genes.loc[non_k_mask, "Clade"].unique().tolist())
        print(f"[WARN] {len(excluded_genomes)} genome(s) have a non-K0/K1/K2 clade "
              f"({excluded_clades}) and are excluded: {excluded_genomes}")
    genes = genes[genes["Clade"].isin(CLADE_ORDER)]
    print(f"[INFO] Virulence gene hits after clade restriction (K0/K1/K2 only): {len(genes)}")

    vir_categories = sorted(genes["Virulence_category"].unique().tolist())
    print(f"[INFO] Virulence categories found: {vir_categories}")

    blocks = [summarize(genes, "pooled", vir_categories)]
    print(f"[INFO] pooled: n={len(genes)} virulence gene hits")
    for clade in CLADE_ORDER:
        sub = genes[genes["Clade"] == clade]
        blocks.append(summarize(sub, clade, vir_categories))
        print(f"[INFO] Clade {clade}: n={len(sub)} virulence gene hits")

    result = pd.concat(blocks, ignore_index=True)
    result.to_csv(OUT_FILE, sep="\t", index=False)
    print(f"[DONE] Saved: {OUT_FILE} ({len(result)} rows)")

    print("\n[INFO] Pooled Any_MGE containment by virulence category:")
    print(result[(result["Group"] == "pooled") & (result["MGE_category"] == "Any_MGE")]
          [["Virulence_category", "n_genes", "n_contained", "pct_contained"]]
          .to_string(index=False))


if __name__ == "__main__":
    main()
