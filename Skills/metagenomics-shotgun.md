---
name: metagenomics-shotgun
description: Shotgun metagenomics — host decontamination, Kraken2 taxonomic profiling, Bracken abundance, HUMAnN3 functional annotation, MEGAHIT assembly, MetaBAT2 binning, CheckM, MAGs, QIIME2 16S
---

# Metagenomics (Shotgun & 16S)

## When to Use

Use this skill when:
- Processing whole-metagenome shotgun sequencing (WMS) data
- Performing taxonomic profiling beyond 16S (species/strain resolution)
- Annotating community functional pathways
- Recovering metagenome-assembled genomes (MAGs)
- Running 16S amplicon analysis with QIIME2

## Quick Reference

| Task | Tool | Key Command |
|------|------|-------------|
| Host decontamination | Bowtie2 | `bowtie2 --un-conc-gz decontam` |
| Taxonomic profiling | Kraken2 | `kraken2 --report report.txt` |
| Abundance re-estimation | Bracken | `bracken -d db -i report.txt -l S` |
| Functional pathways | HUMAnN3 | `humann --input reads.fq.gz` |
| Assembly | MEGAHIT | `megahit -1 R1 -2 R2 --min-contig-len 500` |
| Binning | MetaBAT2 | `metabat2 -i contigs.fa -a depths.txt` |
| Bin QC | CheckM | `checkm lineage_wf bins/ out/` |
| MAG annotation | Prokka | `prokka --metagenome bin.fa` |
| 16S analysis | QIIME2 | `qiime dada2 denoise-paired` |

## Key Patterns

**Pattern 1: Taxonomic profiling pipeline**
```bash
# Remove host reads
bowtie2 -x hg38 -1 R1.fq.gz -2 R2.fq.gz \
    --un-conc-gz decontam_%.fq.gz > /dev/null

# Kraken2 classification
kraken2 --db standard/ --paired --gzip-compressed \
    decontam_1.fq.gz decontam_2.fq.gz \
    --report kraken2_report.txt --output kraken2_out.txt

# Bracken species-level re-estimation
bracken -d standard/ -i kraken2_report.txt \
    -o bracken_species.txt -r 150 -l S -t 10
```

**Pattern 2: Functional annotation with HUMAnN3**
```bash
humann --input decontam_merged.fq.gz \
    --output humann3_out/ --threads 8

# Normalize pathways
humann_renorm_table --input humann3_out/sample_pathabundance.tsv \
    --output pathways_cpm.tsv --units cpm

# Join multiple samples
humann_join_tables --input humann3_outputs/ \
    --output all_pathways.tsv --file_name pathabundance
```

**Pattern 3: Assembly and binning**
```bash
# Assembly
megahit -1 decontam_1.fq.gz -2 decontam_2.fq.gz \
    -o megahit/ --min-contig-len 500 -t 16

# Coverage for binning
bowtie2-build megahit/final.contigs.fa contigs_index
bowtie2 -x contigs_index -1 decontam_1.fq.gz -2 decontam_2.fq.gz | \
    samtools sort -o contigs.bam && samtools index contigs.bam

# Binning with MetaBAT2
jgi_summarize_bam_contig_depths --outputDepth depths.txt contigs.bam
metabat2 -i megahit/final.contigs.fa -a depths.txt -o bins/bin --minContig 1500

# CheckM quality
checkm lineage_wf bins/ checkm_out/ -t 8 -x fa
```

**Pattern 4: QIIME2 16S workflow (CLI)**
```bash
# Import paired-end reads
qiime tools import \
    --type 'SampleData[PairedEndSequencesWithQuality]' \
    --input-path manifest.csv \
    --output-path reads.qza \
    --input-format PairedEndFastqManifestPhred33V2

# DADA2 denoising
qiime dada2 denoise-paired \
    --i-demultiplexed-seqs reads.qza \
    --p-trunc-len-f 250 --p-trunc-len-r 200 \
    --o-table table.qza --o-representative-sequences rep_seqs.qza

# Taxonomy classification
qiime feature-classifier classify-sklearn \
    --i-classifier silva138_classifier.qza \
    --i-reads rep_seqs.qza --o-classification taxonomy.qza

# Diversity
qiime diversity core-metrics-phylogenetic \
    --i-table table.qza --i-phylogeny rooted_tree.qza \
    --p-sampling-depth 5000 --m-metadata-file metadata.tsv \
    --output-dir diversity/
```

## MAG Quality Standards (MIMAG)

| Tier | Completeness | Contamination |
|------|-------------|---------------|
| High quality | ≥ 90% | < 5% |
| Medium quality | ≥ 50% | < 10% |
| Low quality | < 50% | — |

## Common Pitfalls

- **Host decontamination is critical** — failure to remove host reads inflates classifications
- **Kraken2 database choice** — standard (archaea + bacteria + viral) vs PlusPF (adds protozoa/fungi) vs custom
- **Bracken threshold** — `-t 10` requires ≥10 reads assigned to a taxon; lower for low-coverage samples
- **HUMAnN3 databases** — must download UniRef90 and ChocoPhlAn databases (~25 GB) separately
- **Binning quality** — MetaBAT2 needs ≥2× coverage and ≥500 bp contigs; CONCOCT or MaxBin2 can be combined (DAS Tool for bin refinement)

## Code Templates

### Parse Kraken2 Report
```python
import pandas as pd

def read_kraken2_report(path):
    cols = ['pct', 'clade_reads', 'direct_reads', 'rank', 'taxid', 'name']
    df = pd.read_csv(path, sep='\t', header=None, names=cols)
    df['name'] = df['name'].str.strip()
    return df

report = read_kraken2_report('kraken2_report.txt')
# Species level
species = report[report['rank'] == 'S'].sort_values('pct', ascending=False)
print(species[['name', 'pct']].head(10).to_string(index=False))
```

### Alpha Diversity from OTU Table
```python
import pandas as pd
import numpy as np

def shannon_diversity(counts):
    """Shannon H' diversity index."""
    counts = np.array(counts, dtype=float)
    counts = counts[counts > 0]
    props = counts / counts.sum()
    return -np.sum(props * np.log(props))

def observed_richness(counts):
    return np.sum(np.array(counts) > 0)

# otu_table: samples × taxa DataFrame
otu_table = pd.read_csv('otu_table.tsv', sep='\t', index_col=0)
alpha = pd.DataFrame({
    'shannon': otu_table.apply(shannon_diversity, axis=1),
    'richness': otu_table.apply(observed_richness, axis=1),
})
print(alpha.describe())
```

### CheckM TSV Parser
```python
import pandas as pd

def read_checkm(path):
    """Parse CheckM lineage_wf output qa file."""
    df = pd.read_csv(path, sep='\t')
    # Normalize column names
    df.columns = [c.strip().lower().replace(' ', '_') for c in df.columns]
    high_q = df[(df['completeness'] >= 90) & (df['contamination'] < 5)]
    med_q  = df[(df['completeness'] >= 50) & (df['contamination'] < 10)]
    print(f"High-quality MAGs: {len(high_q)}")
    print(f"Medium-quality MAGs: {len(med_q)}")
    return df

checkm = read_checkm('checkm_results.tsv')
```

### Bray-Curtis Dissimilarity (Beta Diversity)
```python
from scipy.spatial.distance import braycurtis
from sklearn.manifold import MDS
import pandas as pd
import numpy as np

otu = pd.read_csv('otu_table.tsv', sep='\t', index_col=0)
# Relative abundance
otu_rel = otu.div(otu.sum(axis=1), axis=0)

n = len(otu_rel)
dist = np.zeros((n, n))
for i in range(n):
    for j in range(i+1, n):
        d = braycurtis(otu_rel.iloc[i], otu_rel.iloc[j])
        dist[i, j] = dist[j, i] = d

# PCoA
mds = MDS(n_components=2, dissimilarity='precomputed', random_state=42)
coords = mds.fit_transform(dist)
```

## Related Skills
- `ngs-variant-calling` — FASTQ QC, short-read alignment, samtools
- `rnaseq-metagenomics` — RNA-seq differential expression, DESeq2/edgeR
- `python-core-bio` — FASTA/FASTQ parsing, sequence manipulation
- `biostatistics-r` — statistical testing, diversity analysis in R
