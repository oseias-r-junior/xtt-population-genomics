import pathlib
ROOT = str(pathlib.Path(__file__).resolve().parents[2]) + "/"
from pathlib import Path
import pandas as pd

# ==========================================
# SETTINGS
# ==========================================

panseq_dir = ROOT + "input_data_hrr_new_phylo_tree/pan_genome_sequences"

presence_absence = ROOT + "input_data_hrr_new_phylo_tree/gene_presence_absence.csv"

output_xmfa = "core_genes_corrected.xmfa"

expected_nseq = 113

# ==========================================
# LOAD GENE PRESENCE/ABSENCE TABLE
# ==========================================

df = pd.read_csv(
    presence_absence,
    low_memory=False
)

# ==========================================
# IDENTIFY STRAIN COLUMNS
# ==========================================

metadata_cols = 14

strain_columns = list(df.columns[metadata_cols:])

print(f"Detected {len(strain_columns)} strains")

# ==========================================
# BUILD LOCUS -> STRAIN MAP
# ==========================================

locus_to_strain = {}

for _, row in df.iterrows():

    for strain in strain_columns:

        value = row[strain]

        if pd.isna(value):
            continue

        # Sometimes multiple loci separated by tabs/semicolons
        loci = str(value).replace("\t", ";").split(";")

        for locus in loci:

            locus = locus.strip()

            if locus:
                locus_to_strain[locus] = strain

print(f"Mapped {len(locus_to_strain)} locus tags")

# ==========================================
# FIND ALIGNMENTS
# ==========================================

alignment_files = sorted(
    Path(panseq_dir).glob("*.fa.aln")
)

# ==========================================
# BUILD XMFA
# ==========================================

kept = 0
removed = 0

with open(output_xmfa, "w") as out:

    for aln_file in alignment_files:

        with open(aln_file) as f:

            lines = f.readlines()

        # Count sequences
        nseq = sum(
            1 for line in lines
            if line.startswith(">")
        )

        # Keep only complete/core genes
        if nseq != expected_nseq:

            removed += 1
            continue

        kept += 1

        for line in lines:

            if line.startswith(">"):

                locus = line.strip()[1:]

                strain = locus_to_strain.get(
                    locus,
                    locus
                )

                out.write(f">{strain}\n")

            else:
                out.write(line)

        out.write("=\n")

# ==========================================
# REPORT
# ==========================================

print(f"Core genes kept: {kept}")

print(f"Genes removed: {removed}")

print(f"XMFA written to: {output_xmfa}")
