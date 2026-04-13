---
name: bio-core-blast-searching
description: "1. Explain how BLAST works (the seed-and-extend heuristic) and why it is needed 2. Choose the correct BLAST variant for your data (blastn, blastp, blastx, tblastn, tblastx) 3. Run BLAST searches via t"
tool_type: python
source_notebook: "Tier_2_Core_Bioinformatics/04_BLAST_Searching/01_blast_searching.ipynb"
primary_tool: BLAST+
---

## Version Compatibility

Reference examples tested with: biopython 1.83+, matplotlib 3.8+, numpy 1.26+, pandas 2.1+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.


# BLAST: Sequence Similarity Searching

*Source: Course notebook `Tier_2_Core_Bioinformatics/04_BLAST_Searching/01_blast_searching.ipynb`*


**Tier 2 -- Core Bioinformatics**

---

## Learning Objectives

By the end of this notebook you will be able to:

1. Explain how BLAST works (the seed-and-extend heuristic) and why it is needed
2. Choose the correct BLAST variant for your data (blastn, blastp, blastx, tblastn, tblastx)
3. Run BLAST searches via the NCBI web interface and programmatically with BioPython
4. Parse BLAST results in both XML and tabular formats
5. Interpret E-values, bit scores, percent identity, and query coverage
6. Use PSI-BLAST for detecting remote homologs
7. Apply best practices for database selection, filtering, and parameter tuning

---

## How to use this notebook

1. Run cells top-to-bottom in order — later cells depend on earlier ones
2. Run the environment check cell first to confirm all imports work
3. Read the explanatory text before each code cell — the context matters
4. The exercises at the end are designed to be done after reading each section
5. If a code cell requires internet access (Entrez, PDB download), it is marked — these can be skipped if offline

## Complicated moments explained

- **Choosing the right BLAST program**: The five programs differ by query type and database type. blastn: nucleotide vs nucleotide. blastp: protein vs protein. blastx: nucleotide query (translated) vs protein database — use when you have a novel nucleotide sequence and want to find protein homologs. tblastn: protein query vs nucleotide database (translated) — use when you want to find unannotated protein-coding genes. tblastx: nucleotide query (translated) vs nucleotide database (translated) — slowest, rarely needed.
- **E-value depends on database size**: The same alignment gives a much lower (better) E-value against a small database than against nr. A hit with E=1e-5 in SwissProt is not the same as E=1e-5 in nr. Always report which database you searched.
- **E-value is not a p-value**: E-value is the *expected number of hits* with that score by chance. E=0.01 means you expect 0.01 random hits of that quality — not that there is a 1% chance of a false positive. For database searches, use E < 1e-5 as a typical significance threshold, but always consider alignment length and percent identity.
- **Low-complexity filtering (SEG/DUST)**: Sequences with low complexity (poly-A stretches, coiled-coil repeats, transmembrane helices) will generate thousands of spurious hits. Always leave low-complexity filtering enabled (the default) unless you have a specific reason to disable it.
- **PSI-BLAST convergence**: PSI-BLAST iterates until the set of hits stops changing. Watch for 'contamination' — if a false positive enters the PSSM in an early iteration, subsequent iterations become biased. Always verify hits biologically.

## Environment check (run this first)

```python
# Environment check
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Align import substitution_matrices
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

blosum62 = substitution_matrices.load("BLOSUM62")
print("BioPython BLAST modules loaded.")

# Show the 5 BLAST programs at a glance
programs = {
    "blastn":  ("Nucleotide", "Nucleotide", "Standard DNA/RNA search"),
    "blastp":  ("Protein",    "Protein",    "Standard protein search"),
    "blastx":  ("Nucleotide (translated)", "Protein", "Annotate novel nucleotide seqs"),
    "tblastn": ("Protein",    "Nucleotide (translated)", "Find unannotated protein genes"),
    "tblastx": ("Nucleotide (translated)", "Nucleotide (translated)", "Cross-species DNA comparison"),
}
print(f"\n{'Program':<10} {'Query':<30} {'Database':<30} {'Use case'}")
print("-" * 100)
for prog, (query, db, use) in programs.items():
    print(f"{prog:<10} {query:<30} {db:<30} {use}")
print("\nProceed to Section 1.")
```python

---

## 2. How BLAST Works: The Seed-and-Extend Heuristic

BLAST achieves its speed by **not** performing full Smith-Waterman on every query-database pair. Instead, it uses a three-phase heuristic:

### Phase 1: Find Seeds (Word Matching)

BLAST breaks the query into overlapping **words** of length $W$ (word size):
- Protein BLAST (blastp): $W = 3$ (default) — every 3-mer in the query
- Nucleotide BLAST (blastn): $W = 11$ (default) — every 11-mer

For protein BLAST, BLAST does not just look for exact word matches. It looks for **word neighborhoods** — all database words whose alignment score with the query word exceeds a threshold $T$. This is what makes blastp sensitive to divergent sequences even at the seeding stage.

```python
Query word: LEW
Neighborhood (score >= T=11 with BLOSUM62):
  LEW (exact, score=14), LEF (score=12), LEY (score=11), LKW (score=11), ...

This allows BLAST to find seeds even when the database sequence differs slightly.
```python

### Phase 2: Extend Seeds (Ungapped Extension)

When a seed is found, BLAST extends it in both directions **without allowing gaps**. Extension continues as long as the score does not drop more than $X$ below the best score seen so far (the **X-drop** heuristic). This produces a **high-scoring pair (HSP)**.

```python
Query:    ...MVLSPADKTNVKAAWGKVGAHAG...
               |||||||||||||||||||||||
Database: ...MVLSGEDKSNIKAAWGKIGGHGAE...
                    ^ seed ^
              <-- extend left   extend right -->
              Stop when cumulative score drops X below maximum
```python

### Phase 3: Gapped Extension (Full Alignment)

The highest-scoring ungapped HSPs are then extended with a **gapped alignment** (Smith-Waterman with dynamic programming). This is the expensive step, but by this point only a tiny fraction of database sequences are being processed.

### Why BLAST Is Approximate

BLAST may **miss** true homologs because:
1. No seed is found (query and target diverged too much in all word-length regions)
2. The X-drop terminates extension early

This is the sensitivity/speed tradeoff. PSI-BLAST and slower tools (e.g., SW search via SSEARCH) are more sensitive but slower.

### Summary of Key Parameters

| Parameter | Default | Effect |
|---|---|---|
| Word size ($W$) | 3 (protein), 11 (DNA) | Larger → faster, less sensitive |
| Word threshold ($T$) | 11 (blastp) | Higher → faster, less sensitive |
| Ungapped X-drop ($X_g$) | 20 | Higher → extends seeds further |
| Gapped X-drop ($X_u$) | 50 | Affects final extension |
| E-value cutoff | 10 | Higher → more hits reported |

```python
# Demonstration: word matching and neighborhood concept

from Bio.Align import substitution_matrices
import numpy as np

blosum62 = substitution_matrices.load("BLOSUM62")

def word_neighborhood(query_word, matrix, threshold=11):
    """
    Find all 3-letter amino acid words whose BLOSUM62 score with
    query_word exceeds the threshold T.
    
    This is what blastp does during the seeding phase.
    """
    amino_acids = "ACDEFGHIKLMNPQRSTVWY"
    neighbors = []
    for a1 in amino_acids:
        for a2 in amino_acids:
            for a3 in amino_acids:
                word = a1 + a2 + a3
                score = (matrix[query_word[0], word[0]] +
                         matrix[query_word[1], word[1]] +
                         matrix[query_word[2], word[2]])
                if score >= threshold:
                    neighbors.append((word, score))
    return sorted(neighbors, key=lambda x: -x[1])


def xdrop_extend(seq1, seq2, match=2, mismatch=-1, x_drop=5):
    """
    Ungapped X-drop extension from the center of two sequences.
    Returns the score of the highest-scoring ungapped alignment.
    
    In real BLAST, this starts from a seed position and extends
    in both directions.
    """
    best_score = 0
    score = 0
    
    for i in range(min(len(seq1), len(seq2))):
        score += match if seq1[i] == seq2[i] else mismatch
        if score > best_score:
            best_score = score
        if best_score - score > x_drop:
            break
    return best_score


# Example 1: show word neighborhoods for "LEW"
print("=== BLAST Word Neighborhood (threshold T=11) ===")
query_word = "LEW"
neighbors = word_neighborhood(query_word, blosum62, threshold=11)
print(f"Query word: {query_word}")
print(f"Neighborhood size: {len(neighbors)} words")
print(f"Top 10 neighbors:")
for word, score in neighbors[:10]:
    print(f"  {word}: score = {score}")

# Example 2: compare neighborhood sizes for different thresholds
print()
print("=== Sensitivity vs Speed: Word Neighborhood Size vs Threshold ===")
print(f"{'Threshold T':>15} {'Neighborhood size':>20} {'Relative sensitivity'}")
print("-" * 55)
for T in [7, 9, 11, 13, 15]:
    n = len(word_neighborhood(query_word, blosum62, threshold=T))
    rel = n / len(word_neighborhood(query_word, blosum62, threshold=7)) * 100
    print(f"{T:>15} {n:>20} {rel:>18.1f}%")

print()
print("Higher T = smaller neighborhood = faster but less sensitive BLAST.")

# Example 3: X-drop extension demonstration
print()
print("=== X-drop Extension Demonstration ===")
seq_ref = "MVLSPADKTNVKAAWGKVGAHAG"
seq_hit = "MVLSGEDKSNIKAAWGKIGGHGAE"
score = xdrop_extend(seq_ref, seq_hit, match=2, mismatch=-1, x_drop=5)
print(f"Query:    {seq_ref}")
print(f"Database: {seq_hit}")
print(f"Ungapped X-drop score (X=5): {score}")
print()
print("The X-drop extension terminates when the running score drops more")
print("than X below the best score seen — avoiding wasted computation")
print("on sequences that diverge after a good seed.")
```python

---

## 3. BLAST Program Variants

Different BLAST programs handle different combinations of query and database types:

| Program | Query | Database | Translation | Use Case |
|---|---|---|---|---|
| **blastn** | Nucleotide | Nucleotide | None | Gene finding, species ID, primer design |
| **blastp** | Protein | Protein | None | Protein homology, function prediction |
| **blastx** | Nucleotide | Protein | Query in 6 frames | EST annotation, finding coding regions |
| **tblastn** | Protein | Nucleotide | DB in 6 frames | Finding unannotated genes in genomes |
| **tblastx** | Nucleotide | Nucleotide | Both in 6 frames | Comparing unannotated genomes |

### Decision Guide

```python
What is your QUERY?                 What is your DATABASE?

   DNA/RNA ----+----> DNA/RNA DB ---------> blastn (or megablast)
               |
               +----> Protein DB ---------> blastx (translate query)

   Protein ----+----> Protein DB ---------> blastp
               |
               +----> DNA/RNA DB ---------> tblastn (translate DB)

   DNA/RNA ----------> DNA/RNA DB ---------> tblastx (translate both)
   (when you suspect     (slow but most sensitive for
    distant homology)     distant nucleotide comparisons)
```python

### Nucleotide BLAST Sub-variants

| Variant | Word size | Best for |
|---|---|---|
| **megablast** | 28 | Highly similar sequences (>95% identity), same species |
| **blastn** | 11 | Somewhat similar sequences, cross-species |
| **discontiguous megablast** | 11-12 | More dissimilar sequences (~80% identity) |

```python
# Helper to recommend BLAST program based on input
def recommend_blast(query_type, db_type, sensitivity="standard"):
    """Recommend the appropriate BLAST program."""
    recommendations = {
        ("nucleotide", "nucleotide"): {
            "high_identity": ("megablast", "Best for >95% identity, same species"),
            "standard": ("blastn", "Standard nucleotide search"),
            "sensitive": ("discontiguous megablast", "Cross-species, ~80% identity"),
            "very_sensitive": ("tblastx", "Both translated -- slowest but most sensitive"),
        },
        ("nucleotide", "protein"): {
            "standard": ("blastx", "Translates query in 6 frames"),
        },
        ("protein", "protein"): {
            "standard": ("blastp", "Standard protein search"),
            "sensitive": ("PSI-BLAST", "Iterative, for remote homologs"),
        },
        ("protein", "nucleotide"): {
            "standard": ("tblastn", "Translates DB in 6 frames"),
        },
    }
    
    key = (query_type, db_type)
    if key not in recommendations:
        return "Invalid combination"
    
    options = recommendations[key]
    if sensitivity in options:
        prog, desc = options[sensitivity]
    else:
        prog, desc = options["standard"]
    
    return f"{prog}: {desc}"


# Test various scenarios
scenarios = [
    ("nucleotide", "nucleotide", "high_identity", "Verify a clone sequence"),
    ("nucleotide", "nucleotide", "standard", "Find gene family members"),
    ("protein", "protein", "standard", "Find protein homologs"),
    ("protein", "protein", "sensitive", "Find distant protein homologs"),
    ("nucleotide", "protein", "standard", "Annotate an EST/mRNA"),
    ("protein", "nucleotide", "standard", "Find unannotated genes"),
]

print(f"{'Task':<35} {'Recommendation'}")
print("=" * 75)
for qt, dt, sens, task in scenarios:
    rec = recommend_blast(qt, dt, sens)
    print(f"{task:<35} {rec}")
```python

## Common Pitfalls

- **Coordinate systems**: BED uses 0-based half-open; VCF/GFF use 1-based inclusive — mixing them causes off-by-one errors
- **Batch effects**: Always check for batch confounding before interpreting biological signal
- **Multiple testing**: Apply FDR correction (Benjamini-Hochberg) when testing thousands of features simultaneously
