---
name: immunogenomics
description: "V(D)J repertoire analysis, HLA typing, and neoantigen prediction pipelines."
tool_type: cli
primary_tool: NumPy
---

## Version Compatibility

Reference examples tested with: muon 0.1+, numpy 1.26+, scipy 1.12+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.


# immunogenomics

V(D)J repertoire analysis, HLA typing, and neoantigen prediction pipelines.

## Quick Reference

| Task | Tool | Output |
|------|------|--------|
| scTCR/BCR analysis | scirpy | Clonotype tables, diversity plots |
| Bulk TCR/BCR from RNA-seq | TRUST4 | CDR3 sequences, V/J assignments |
| HLA typing | OptiType | 4-digit HLA alleles |
| Peptide-MHC binding | NetMHCpan | Rank percentile, binding affinity |
| Neoantigen pipeline | pVACseq | Prioritized epitope list |

## scirpy: 10x V(D)J Analysis

```python
import scirpy as ir
import muon as mu

# Load 10x VDJ + GEX data
mdata = mu.read_10x_vdj('all_contig_annotations.csv')
ir.pp.index_chains(mdata)
ir.tl.chain_qc(mdata)

# Define clonotypes (CDR3 amino acid)
ir.pp.ir_dist(mdata, metric='identity')
ir.tl.define_clonotypes(mdata, within_group='sample')
ir.tl.clonal_expansion(mdata)

# Diversity
ir.tl.alpha_diversity(mdata, groupby='condition', target_col='clone_id')

# V-gene usage
ir.pl.vdj_usage(mdata, full_combination=True)
ir.pl.clonotype_network(mdata, color='condition')
```python

## TRUST4 Bulk Repertoire

```bash
# Extract TCR/BCR from RNA-seq BAM
run-trust4 \
    -b tumor_rna.bam \
    -f hg38_bcrtcr.fa \
    --ref human_IMGT+C.fa \
    --thread 8 \
    -o trust4_output

# Output: trust4_output_report.tsv (CDR3, V, J, abundance)
```python

## OptiType HLA Typing

```bash
# Step 1: Extract reads mapping to HLA loci
bwa mem hla_reference.fa sample_R1.fastq.gz sample_R2.fastq.gz \
    | samtools view -b -F 4 > hla_reads.bam
samtools sort -n hla_reads.bam | samtools fastq -1 hla_R1.fq -2 hla_R2.fq

# Step 2: OptiType
OptiTypePipeline.py \
    -i hla_R1.fq hla_R2.fq \
    --dna --verbose \
    --outdir hla_typing/ \
    --prefix sample
# Output: sample_result.tsv → HLA-A,B,C alleles
```python

## Neoantigen Prediction (Python)

```python
# Simplified pMHC binding prediction workflow
# Full pipeline: use pVACseq

import subprocess

# 1. VEP annotation for somatic variants
# !vep -i somatic.vcf -o annotated.vcf --cache --everything

# 2. Extract mutant peptides (9-11 mers around mutation)
def extract_peptides(mut_aa_seq, position, lengths=[9, 10, 11]):
    peptides = []
    for length in lengths:
        for start in range(max(0, position - length + 1), position + 1):
            end = start + length
            if end <= len(mut_aa_seq):
                peptides.append(mut_aa_seq[start:end])
    return peptides

# 3. NetMHCpan prediction
# !netMHCpan -p peptides.txt -a HLA-A02:01,HLA-B07:02 -l 9,10,11 -BA > binding.txt
# Filter: rank < 0.5 = strong binder, rank < 2.0 = weak binder
```python

## Clonotype Diversity Metrics

```python
import numpy as np
from scipy.stats import entropy

def shannon_entropy(counts):
    """Shannon entropy of clonotype distribution."""
    freq = np.array(counts) / np.sum(counts)
    return entropy(freq, base=2)

def simpson_index(counts):
    """Simpson diversity index (1 - D)."""
    freq = np.array(counts) / np.sum(counts)
    return 1 - np.sum(freq ** 2)

def clonal_expansion_index(counts, threshold=2):
    """Fraction of cells in expanded clones (count > threshold)."""
    return np.sum([c for c in counts if c >= threshold]) / np.sum(counts)
```python

## Key Databases
- **IMGT**: V/J/D/C gene segment reference sequences
- **VDJdb**: antigen-specific TCR/BCR sequences with HLA restrictions
- **McPAS-TCR**: manually curated pathology-associated TCRs
- **IEDB**: immune epitope database for T/B cell epitopes

## Common Pitfalls
- **Chain pairing**: 10x gives paired α/β, but some cells may have 2 TCRα chains
- **HLA resolution**: 2-digit (HLA-A*02) vs 4-digit (HLA-A*02:01) affects binding prediction
- **Neoantigen filtering**: also consider RNA expression and peptide processing (TAP, proteasome)
- **TRUST4 sensitivity**: works best with > 50M reads; tumor purity affects sensitivity

## Module
Tier 3 · Module 34 (Immunogenomics)
