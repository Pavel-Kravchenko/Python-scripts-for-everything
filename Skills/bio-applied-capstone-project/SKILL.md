---
name: bio-applied-capstone-project
description: "**Estimated time:** 8-12 hours (can be split across multiple sessions)"
tool_type: python
source_notebook: "Tier_3_Applied_Bioinformatics/08_Capstone_Project/01_capstone_project.ipynb"
---

# From Sequence to Discovery: An Integrative Bioinformatics Project

*Source: Course notebook `Tier_3_Applied_Bioinformatics/08_Capstone_Project/01_capstone_project.ipynb`*

# From Sequence to Discovery: An Integrative Bioinformatics Project

---

### Capstone Project -- Tier 3: Applied Bioinformatics

**Estimated time:** 8-12 hours (can be split across multiple sessions)

---

## Overview

In this capstone project, you will analyze a set of **unknown DNA sequences** end-to-end,
applying skills from every tier of the course. You will identify the sequences, align them,
build phylogenetic trees, analyze protein structure, find functional domains, and produce
publication-quality figures -- exactly the workflow a bioinformatician follows when
confronted with new sequence data.

### What you will do

```
Unknown DNA sequences
       |
       v
  [Step 1] Quality control & cleaning
       |
       v
  [Step 2] Sequence identification (BLAST)
       |
       v
  [Step 3] Multiple sequence alignment
       |
       v
  [Step 4] Phylogenetic tree construction
       |
       v
  [Step 5] Protein structure analysis
       |
       v
  [Step 6] Functional domains & motifs
       |
       v
  [Step 7] GO / Pathway enrichment
       |
       v
  [Step 8] Publication-quality figures
       |
       v
  [Step 9] Scientific report
```

### How this notebook is organized

Each step has:
1. **Task description** -- what you need to accomplish
2. **Empty cells** -- where you write your code and analysis
3. **Hints** (collapsible) -- guidance if you get stuck
4. **Solution walkthrough** (collapsible) -- a complete reference implementation

Try each step yourself before opening the hints or solution.

---

### Learning Objectives

By completing this project, you will demonstrate your ability to:

- Assess raw sequence data quality and apply appropriate cleaning strategies
- Identify unknown sequences using BLAST and interpret the results critically
- Perform and visualize multiple sequence alignments
- Construct phylogenetic trees and evaluate their biological meaning
- Analyze protein structure and map functional residues
- Identify conserved domains and motifs using pattern matching
- Place genes in functional context using Gene Ontology and pathway databases
- Create publication-quality, multi-panel scientific figures
- Communicate findings in a clear, structured scientific report

### Prerequisites

This project assumes you have completed (or have equivalent knowledge from):
- **Tier 1**: Python fundamentals, data structures, regex, pandas, matplotlib
- **Tier 2**: BioPython, BLAST, sequence alignment, phylogenetics, protein structure, motifs, GO
- **Tier 3**: NGS fundamentals (for QC concepts)

### Estimated Time

| Step | Topic | Time |
|------|-------|------|
| 1 | QC & Cleaning | 30-45 min |
| 2 | BLAST Identification | 45-60 min |
| 3 | Multiple Alignment | 45-60 min |
| 4 | Phylogenetic Tree | 45-60 min |
| 5 | Protein Structure | 60-90 min |
| 6 | Domains & Motifs | 45-60 min |
| 7 | GO & Pathways | 30-45 min |
| 8 | Publication Figures | 60-90 min |
| 9 | Scientific Report | 60-90 min |
| **Total** | | **8-12 hours** |

## Setup and Dependencies

Run the cell below to install and import all required packages.

```python
# Install dependencies (uncomment if needed)
# !pip install biopython matplotlib seaborn pandas numpy requests

import warnings
warnings.filterwarnings('ignore')

# Core
import os
import re
import json
import time
from io import StringIO
from collections import Counter, defaultdict

# Data & plotting
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns

# BioPython
from Bio import SeqIO, AlignIO, Entrez, SearchIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import gc_fraction
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Align.Applications import ClustalOmegaCommandline, MuscleCommandline
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor

# Plotting defaults
sns.set_style('whitegrid')
plt.rcParams['figure.figsize'] = (10, 6)
plt.rcParams['figure.dpi'] = 120
plt.rcParams['font.size'] = 11

# Set your email for NCBI Entrez queries
Entrez.email = "your.email@example.com"  # <-- CHANGE THIS

print("All imports successful.")
```

---

## Input Data: Unknown DNA Sequences

Below are 7 DNA sequences isolated from environmental samples. Your mission is to determine
what they are, where they come from, and what they can tell us about the evolutionary
relationships among the source organisms.

These are real **cytochrome c** coding sequences from different species -- but you are not
supposed to know that yet. Pretend you received them from a sequencing facility with only
sample IDs.

Run the cell below to load the sequences into memory.

```python
# ============================================================
# Input sequences -- cytochrome c from 7 species
# Realistic CDS sequences with minor "sequencing noise" added
# ============================================================

raw_sequences = {
    "SAMPLE_001": (
        "NNNatgggtgatgttgagaaaggtaccaagcacggtaaagagttgactgaattgaagaat"
        "atcaagaaagctggtttcaaagctggtttgggtgatactttgaagaagcacactggtgaa"
        "gctccattcaagtacattacaaagaagatcccagatgctcctgctaagaagattttggct"
        "gcacctactaaatcttctaaattgtctactattgatgcaaagggttataaggatgtttca"
        "tgggctaaagaattgaatgaatacatcNNN"
    ),
    "SAMPLE_002": (
        "NNatgggagatgttgaaaaagggacaaagcatggaaagGaattgactgaactgaagaac"
        "atcaagaaggccggtttcaaagcgggtctcggtgatacattgaagaaacacacaggtgaa"
        "gcaccattcaaatacattactaaaaagatcccagatgcaccagcaaagaagattctcgca"
        "gcaccaaccaaatcatcaaaattgtccaccattgatgcgaagggctataaagacgtttcc"
        "tgggccaaagaattgaatgaatacattNN"
    ),
    "SAMPLE_003": (
        "atgggtgatgttgagaaaggtactaagcatggcaaagaactaactgaattaaagaat"
        "ataaaaagagctggctttaaagcaggcttgggagatactttgaagaaacatacaggagaa"
        "gcaccgttcaaatatattacaaaaaagataccagatgcgcctgctaagaagatcttggca"
        "gcaccaaccaaatcttcaaaattatccactattgatgcaaagggttacaaagatgtttca"
        "tgggctaaagaactaaatgaatacata"
    ),
    "SAMPLE_004": (
        "NNNNatgggcgacgtggagaaggggaccaagcacggcaaggagctgaccgagctgaagaac"
        "atcaagaaggccggcttcaaggccggcctgggcgacaccctgaagaagcacaccggcgag"
        "gccccgttcaagtacatcaccaagaagatcccggacgcgccggcgaagaagatcctggcg"
        "gcgccgaccaagtcgtcgaagctgtcgaccatcgacgcgaagggctacaaggacgtgtcg"
        "tgggcgaaggagctgaacgagtacatcNNNN"
    ),
    "SAMPLE_005": (
        "Natgggtgatgttgaaaaaggtactaaacatggtaaagaattgactgatttaaagaat"
        "ataaagaaagcaggttttaaagctggtttgggagatactttaaagaaacatacaggtgaa"
        "gctccatttaaatatataacaaaaaagataccagatgctcctgcaaagaagattttagct"
        "gcaccaactaaatcttcaaaattatctactattgatgcaaaaggatataaggatgtttca"
        "tgggctaaagaattaaatgaatatataN"
    ),
    "SAMPLE_006": (
        "NNNatgggtgatgtcgaaaagggcaccaagcatggcaaagagctcaccgaactaaagaat"
        "atcaagaaggctggttttaaggccggcctgggtgatactctaaagaaacatactggcgaa"
        "gccccattcaagtatatcactaagaagatcccggatgcccctgctaagaagattctggcg"
        "gcaccgactaagtcgtcgaagctctcgaccatcgatgcgaagggttataaggacgtctcg"
        "tgggccaaggagctgaacgaatacatcNNN"
    ),
    "SAMPLE_007": (
        "atgggtgatgttgaaaagggtactaaacatggtaaagaattaactgaattgaagaat"
        "atcaaaaaagctggtttcaaagctggcctgggagatactttgaaaaaacacaccggtgaa"
        "gcaccattcaaatatattactaaaaagataccagatgcaccagcaaaaaagattttggca"
        "gcaccaactaaatcttctaaattatctaccattgatgcaaagggctataaagatgtttca"
        "tgggctaaagagttaaacgagtacata"
    ),
}

# Convert to SeqRecord objects
seq_records = []
for name, seq_str in raw_sequences.items():
    record = SeqRecord(Seq(seq_str.upper()), id=name, description=f"Unknown sample {name}")
    seq_records.append(record)

print(f"Loaded {len(seq_records)} sequences:")
for rec in seq_records:
    print(f"  {rec.id}: {len(rec.seq)} bp")
```

```python
# Quick exploration of the raw input data
print("Raw Sequence Summary")
print("=" * 50)
for name, seq_str in raw_sequences.items():
    seq_upper = seq_str.upper()
    n_count = seq_upper.count('N')
    gc = (seq_upper.count('G') + seq_upper.count('C')) / len(seq_upper) * 100
    print(f"{name}: {len(seq_upper):>4} bp, "
          f"N's: {n_count:>2}, "
          f"GC: {gc:.1f}%, "
          f"Starts: {seq_upper[:6]}... "
          f"Ends: ...{seq_upper[-6:]}")
```

**Initial observations** (fill in after running the cell above):

- How many sequences contain ambiguous bases?
- Are the sequences roughly similar in length?
- Do they appear to be coding sequences (look at the first 3 bases after trimming N's)?

Write your observations here before proceeding to Step 1:

*Your observations:*

---

## Step 1: Quality Control and Sequence Cleaning

### Task

Before any analysis, you must inspect and clean the raw sequences:

1. Check each sequence for **ambiguous bases** (N characters) and report their positions
2. **Trim** leading and trailing N's from each sequence
3. Calculate **basic statistics** for each cleaned sequence: length, GC content
4. Verify that the cleaned sequences look like valid **coding sequences** (start with ATG, length divisible by 3)
5. Store the cleaned sequences for downstream analysis

Write your QC code in the cells below.

```python
# YOUR CODE: Step 1 -- Quality Control
# Inspect the raw sequences for ambiguous bases, trim N's, calculate stats
```

```python
# YOUR CODE: Step 1 continued -- Verify coding sequence properties
```

<details>
<summary><b>Hint: Step 1</b> (click to expand)</summary>

- Use `str.strip('N')` to remove leading/trailing N characters
- Count internal N's with `seq.count('N')` after stripping
- `gc_fraction()` from `Bio.SeqUtils` gives you GC content
- A valid CDS starts with ATG, has length divisible by 3, and ideally ends with a stop codon (TAA, TAG, TGA)
- Create a new list of cleaned `SeqRecord` objects for downstream use

</details>

<details>
<summary><b>Solution: Step 1</b> (click to expand)</summary>

```python
# ----- Step 1 Solution: Quality Control -----

cleaned_records = []

print(f"{'Sample':<14} {'Raw len':>8} {'N (lead)':>9} {'N (trail)':>10} "
      f"{'Clean len':>10} {'GC%':>6} {'Starts ATG':>11} {'Len%3':>6}")
print('-' * 85)

for rec in seq_records:
    raw_seq = str(rec.seq)
    raw_len = len(raw_seq)
    
    # Count leading/trailing N's
    stripped = raw_seq.strip('N')
    leading_n = len(raw_seq) - len(raw_seq.lstrip('N'))
    trailing_n = len(raw_seq) - len(raw_seq.rstrip('N'))
    
    # Clean sequence
    clean_seq = Seq(stripped)
    clean_len = len(clean_seq)
    gc = gc_fraction(clean_seq) * 100
    starts_atg = str(clean_seq)[:3] == 'ATG'
    divisible = clean_len % 3
    
    # Store cleaned record
    clean_rec = SeqRecord(clean_seq, id=rec.id, description=f"Cleaned {rec.id}")
    cleaned_records.append(clean_rec)
    
    print(f"{rec.id:<14} {raw_len:>8} {leading_n:>9} {trailing_n:>10} "
          f"{clean_len:>10} {gc:>5.1f}% {'Yes' if starts_atg else 'NO':>11} {divisible:>6}")

print(f"\nAll {len(cleaned_records)} sequences cleaned and validated.")
```

</details>

### Step 1 Checkpoint

Before moving on, verify you have:
- [ ] Identified and counted all ambiguous bases in each sequence
- [ ] Trimmed leading and trailing N's from all sequences
- [ ] Calculated length and GC content for each cleaned sequence
- [ ] Verified that all cleaned sequences start with ATG
- [ ] Stored cleaned sequences in a list of `SeqRecord` objects

---

---

## Step 2: Sequence Identification with BLAST

### Task

Now that you have clean sequences, identify what they are:

1. Run a **BLAST search** for at least one of your sequences against NCBI's `nt` or `nr` database
2. Parse the BLAST results to find the **top hits** (organism, gene name, E-value, identity%)
3. Determine the **gene family** and the **species** each sequence comes from
4. Summarize your identification results in a table

**Note:** NCBI BLAST queries can take 1-5 minutes. You may BLAST just 2-3 representative
sequences and infer the rest from alignment similarity.

```python
# YOUR CODE: Step 2 -- BLAST search
# Run BLAST for one or more sequences, parse and display results
```

```python
# YOUR CODE: Step 2 continued -- Summarize identifications
```
