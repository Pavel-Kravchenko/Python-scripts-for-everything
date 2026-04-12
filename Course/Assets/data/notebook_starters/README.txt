Notebook starter data pack
==========================

These files are intentionally tiny and synthetic so every notebook can run quickly.

Files:
- dna_sequences.csv: toy DNA fragments with labels for binary tasks.
- protein_sequences.fasta: short protein sequences for embedding/folding demos.
- variants.csv: simple SNV table with toy effect metadata.
- expression_matrix.csv: tiny gene-by-cell expression matrix.
- sample_metadata.json: metadata and example mappings for joins.
- small_molecules.csv: miniature molecule table with SMILES and toy activity.

Usage (Python):
    import pandas as pd
    base = "Course/Assets/data/notebook_starters"
    df = pd.read_csv(f"{base}/dna_sequences.csv")
    print(df.head())

Usage (R):
    base <- "Course/Assets/data/notebook_starters"
    df <- read.csv(file.path(base, "dna_sequences.csv"))
    head(df)
