---
name: virology-bioinformatics
description: "Viral genome assembly, intra-host variant calling, phylodynamics, and real-time surveillance."
tool_type: cli
primary_tool: Pandas
---

## Version Compatibility

Reference examples tested with: pandas 2.1+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.


# virology-bioinformatics

Viral genome assembly, intra-host variant calling, phylodynamics, and real-time surveillance.

## Quick Reference

| Task | Tool | Notes |
|------|------|-------|
| Trim ARTIC primers | iVar trim | `-b primer.bed -e` |
| Align to reference | Minimap2 / BWA-MEM | `-ax sr` for Illumina |
| Consensus genome | iVar consensus | `-t 0.5` (majority base) |
| Low-freq variants | LoFreq | SNV at ≥ 1% frequency |
| Clade assignment | Nextclade | SARS-CoV-2, flu, RSV |
| Lineage classification | pangolin | Pango lineage |
| Timed phylogeny | Nextstrain Augur | TreeTime integration |
| Wastewater deconvolution | Freyja | Variant proportions |

## SARS-CoV-2 ARTIC Assembly Pipeline

```bash
# 1. Align reads to Wuhan reference (NC_045512.2)
minimap2 -ax sr NC_045512.2.fa sample_R1.fastq.gz sample_R2.fastq.gz \
    | samtools sort -o sample.bam && samtools index sample.bam

# 2. Trim ARTIC v4.1 primers
ivar trim -i sample.bam -b nCoV-2019.primer.bed -p sample_trimmed -e
samtools sort -o sample_trimmed.sorted.bam sample_trimmed.bam
samtools index sample_trimmed.sorted.bam

# 3. Generate consensus (N masking at < 20x coverage)
samtools mpileup -A -d 0 -Q 0 sample_trimmed.sorted.bam \
    | ivar consensus -p consensus -n N -m 20 -t 0.5

# 4. Quality check
samtools flagstat sample_trimmed.sorted.bam
samtools depth sample_trimmed.sorted.bam | awk '{sum+=$3; count++} END {print sum/count}' # mean coverage
```

## Intra-Host Variant Calling with LoFreq

```bash
# Recalibrate base qualities (optional but recommended)
samtools calmd -b sample_trimmed.sorted.bam NC_045512.2.fa > sample_calmd.bam

# Call variants at >= 1% frequency
lofreq call-parallel --pp-threads 8 \
    -f NC_045512.2.fa \
    -o sample_lofreq.vcf \
    --sig 0.01 --bonf dynamic \
    sample_calmd.bam

# Filter: min AF 1%, min depth 100x
lofreq filter -i sample_lofreq.vcf -o sample_filtered.vcf \
    --af-min 0.01 --cov-min 100
```

## Nextclade Assignment (Python API)

```python
import subprocess
import json

# Command-line Nextclade
result = subprocess.run([
    'nextclade', 'run',
    '--dataset-name', 'sars-cov-2',
    '--output-tsv', 'nextclade.tsv',
    'consensus.fa'
], capture_output=True, text=True)

import pandas as pd
nextclade_results = pd.read_csv('nextclade.tsv', sep='\t')
print(nextclade_results[['seqName', 'clade', 'Nextclade_pango', 'qc.overallStatus']])
```

## Nextstrain Augur Pipeline

```bash
# Filter and subsample
augur filter \
    --sequences sequences.fasta \
    --metadata metadata.tsv \
    --output filtered.fasta \
    --group-by country year month \
    --sequences-per-group 10 \
    --min-date 2020-01-01

# Align to reference
augur align \
    --sequences filtered.fasta \
    --reference-sequence NC_045512.2.gbk \
    --output aligned.fasta --fill-gaps

# Build tree (IQ-TREE)
augur tree \
    --alignment aligned.fasta \
    --output tree_raw.nwk \
    --nthreads 8

# Timed tree with TreeTime
augur refine \
    --tree tree_raw.nwk \
    --alignment aligned.fasta \
    --metadata metadata.tsv \
    --timetree --coalescent skyline \
    --output-tree tree.nwk \
    --output-node-data branch_lengths.json
```

## Freyja Wastewater Deconvolution

```bash
# Call variants from wastewater sample
freyja variants wastewater_sample.bam \
    --variants variants.tsv \
    --depths depths.tsv \
    --ref NC_045512.2.fasta

# Deconvolve variant proportions
freyja demix variants.tsv depths.tsv --output demixed.tsv

# Aggregate multiple samples
freyja aggregate output_dir/ --output aggregated.tsv

# Plot time series
freyja plot aggregated.tsv --output wastewater_variants.pdf --lineages
```

## Common Pitfalls
- **Amplicon dropout**: some ARTIC amplicons fail; expect gaps at pool boundaries
- **Consensus masking**: use N (not random base) at low coverage positions
- **Non-conversion**: for amplicon data, primers must be trimmed BEFORE variant calling
- **Phylogenetic signal**: check root-to-tip regression (TempEst) before clock analysis
- **Recombinants**: SARS-CoV-2 XBB lineages have phylogenetic incongruence

## Key Databases
- **GISAID EpiCoV**: SARS-CoV-2 genomes (> 15M sequences)
- **NCBI Virus**: all viral genomes with metadata
- **Nextstrain**: real-time phylogenies for flu, COVID, mpox
- **Pango designations**: lineage nomenclature repository

## Module
Tier 3 · Module 37 (Virology Bioinformatics)
