# 4. BLAST Searching

**Tier 2: Core Bioinformatics**

How BLAST works and how to use it effectively for database searching. Covers the seed-and-extend heuristic, all five BLAST variants, parameter tuning, and programmatic access via BioPython.

## Topics Covered

- Why BLAST is needed: computational cost of optimal alignment at scale
- Seed-and-extend heuristic: word matching, neighborhood generation, extension
- BLAST variants: `blastn`, `blastp`, `blastx`, `tblastn`, `tblastx`
- Parameter tuning: word size, E-value threshold, scoring matrices
- Running BLAST via NCBI web interface and `NCBIWWW` / `NCBIXML`
- Parsing and interpreting results; taxonomy-filtered searches

## Notebooks

| Notebook | Description |
|----------|-------------|
| [01_blast_searching.ipynb](01_blast_searching.ipynb) | BLAST algorithm, variants, parameter tuning, and programmatic usage |

## Prerequisites

- Tier 2 Module 3: Pairwise Sequence Alignment

---

[<- Previous Module](../03_Pairwise_Sequence_Alignment/) | [Back to Course Overview](../../README.md) | [Next Module ->](../05_Multiple_Sequence_Alignment/)
