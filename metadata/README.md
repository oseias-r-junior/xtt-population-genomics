# Metadata

Metadata files are not included. The workflow expects:

- `genome_clades_master.tsv`: tab-separated table with one row per genome, giving the genome identifier and its clade assignment. See `workflow/config/genomes.py` for how it is read.
- A workbook with T3E / secretion-system gene positions per strain (provided by a collaborator; not distributed).

Exact column names are defined in `workflow/config/genomes.py` and in the scripts that read these files.
