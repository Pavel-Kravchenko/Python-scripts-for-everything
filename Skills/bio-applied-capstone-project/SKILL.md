---
name: bio-applied-capstone-project
description: "Integrative bioinformatics pipeline — sequence QC, BLAST identification, MSA, phylogenetics, structure analysis, GO enrichment, and publication figures"
tool_type: python
primary_tool: Matplotlib
---

# From Sequence to Discovery: Integrative Bioinformatics Pipeline

## Pipeline Overview

```
Unknown DNA → QC/clean → BLAST → MSA → Phylogeny → Structure → Domains → GO/Pathways → Figures → Report
```

## Setup

```python
import numpy as np, pandas as pd, matplotlib.pyplot as plt, seaborn as sns
from Bio import SeqIO, AlignIO, Entrez, Phylo
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import gc_fraction
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor

Entrez.email = "your.email@example.com"
```

## Step 1: QC and Cleaning

```python
cleaned_records = []
for rec in seq_records:
    raw_seq = str(rec.seq)
    stripped = raw_seq.strip('N')
    leading_n = len(raw_seq) - len(raw_seq.lstrip('N'))
    trailing_n = len(raw_seq) - len(raw_seq.rstrip('N'))
    clean_seq = Seq(stripped)
    gc = gc_fraction(clean_seq) * 100
    starts_atg = str(clean_seq)[:3] == 'ATG'
    clean_rec = SeqRecord(clean_seq, id=rec.id, description=f"Cleaned {rec.id}")
    cleaned_records.append(clean_rec)
```

**Validation checklist**: trim leading/trailing N's, verify ATG start, confirm length % 3 == 0, calculate GC%.

## Step 2: BLAST Identification

NCBI BLAST queries take 1-5 min. BLAST 2-3 representatives, infer the rest from alignment similarity.

## Step 3-4: MSA and Phylogenetics

Use ClustalOmega or MUSCLE for alignment, then `DistanceCalculator('identity')` + `DistanceTreeConstructor().upgma()` for tree.

## Step 5-6: Structure and Domain Analysis

Use PDB/UniProt for structure, Pfam/InterPro for domain annotation.

## Step 7: GO/Pathway Enrichment

Run enrichment on identified gene set using `goatools` or Enrichr API.

## Pitfalls

- **Coordinate systems**: BED uses 0-based half-open; VCF/GFF use 1-based inclusive — mixing them causes off-by-one errors
- **Batch effects**: Always check for batch confounding before interpreting biological signal
- **Multiple testing**: Apply FDR correction (Benjamini-Hochberg) when testing thousands of features simultaneously
