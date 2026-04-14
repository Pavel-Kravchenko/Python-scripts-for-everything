---
name: metagenomics-shotgun
description: Shotgun metagenomics — host decontamination, Kraken2 taxonomic profiling, Bracken abundance, HUMAnN3 functional annotation, MEGAHIT assembly, MetaBAT2 binning, CheckM, MAGs, QIIME2 16S
tool_type: python
primary_tool: NumPy
---

# Metagenomics (Shotgun  16S)

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

## Key Patterns (Diversity Analysis)

### Alpha Diversity Metrics
```python
import numpy as np

def shannon_diversity(counts):
    """Shannon H' — richness + evenness."""
    counts = np.array(counts, dtype=float)
    p = counts[counts > 0] / counts[counts > 0].sum()
    return -np.sum(p * np.log(p))

def simpson_diversity(counts):
    """Simpson 1-D — probability two reads differ; robust to rare taxa."""
    counts = np.array(counts, dtype=float)
    p = counts[counts > 0] / counts[counts > 0].sum()
    return 1 - np.sum(p ** 2)

def observed_richness(counts):
    return int(np.sum(np.array(counts) > 0))
```

| Metric | Formula | Notes |
|--------|---------|-------|
| Observed species | `count(counts > 0)` | Richness only |
| Shannon H' | `-Σ pᵢ ln(pᵢ)` | Richness + evenness; 0 = no diversity |
| Simpson 1-D | `1 - Σ pᵢ²` | Probability two reads differ; robust to rare taxa |
| Faith's PD | `Σ branch lengths (observed taxa)` | Phylogenetic richness |

### Beta Diversity: Bray-Curtis Distance
```python
def bray_curtis(s1, s2):
    s1, s2 = np.array(s1, float), np.array(s2, float)
    return np.sum(np.abs(s1 - s2)) / np.sum(s1 + s2)
```

| Metric | Phylogenetic | Quantitative |
|--------|:-----------:|:------------:|
| Jaccard | No | No |
| Bray-Curtis | No | Yes |
| Unweighted UniFrac | Yes | No |
| Weighted UniFrac | Yes | Yes |

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
species = report[report['rank'] == 'S'].sort_values('pct', ascending=False)
print(species[['name', 'pct']].head(10).to_string(index=False))
```

### Alpha Diversity from OTU Table
```python
import pandas as pd

otu_table = pd.read_csv('otu_table.tsv', sep='\t', index_col=0)
alpha = pd.DataFrame({
    'shannon': otu_table.apply(shannon_diversity, axis=1),
    'simpson': otu_table.apply(simpson_diversity, axis=1),
    'richness': otu_table.apply(observed_richness, axis=1),
})
print(alpha.describe())
```

### Bray-Curtis Distance Matrix & PCoA
```python
from scipy.spatial.distance import braycurtis
import numpy as np

otu = pd.read_csv('otu_table.tsv', sep='\t', index_col=0)
otu_rel = otu.div(otu.sum(axis=1), axis=0)

n = len(otu_rel)
dist = np.zeros((n, n))
for i in range(n):
    for j in range(i + 1, n):
        d = braycurtis(otu_rel.iloc[i], otu_rel.iloc[j])
        dist[i, j] = dist[j, i] = d

def pcoa(dm: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Classical PCoA. Returns (coords n×2, proportion_explained[:2])."""
    n = len(dm)
    H = np.eye(n) - np.ones((n, n)) / n
    B = -0.5 * H @ (dm ** 2) @ H
    eigvals, eigvecs = np.linalg.eigh(B)
    idx = np.argsort(eigvals)[::-1]
    eigvals, eigvecs = eigvals[idx], eigvecs[:, idx]
    pos = eigvals > 0
    coords = eigvecs[:, pos] * np.sqrt(eigvals[pos])
    prop = eigvals[pos] / eigvals[pos].sum()
    return coords[:, :2], prop[:2]

coords, prop = pcoa(dist)
```

### PERMANOVA
```python
def permanova(dm: np.ndarray, grouping: np.ndarray,
              n_perm: int = 999) -> tuple[float, float]:
    """Returns (F-statistic, p-value)."""
    unique = np.unique(grouping)
    n = len(grouping)

    def f_stat(g):
        ss_tot = np.sum(dm ** 2) / n
        ss_w = sum(
            np.sum(dm[np.ix_(g == grp, g == grp)] ** 2) / (2 * (g == grp).sum())
            for grp in unique if (g == grp).sum() > 1
        )
        df_b, df_w = len(unique) - 1, n - len(unique)
        return ((ss_tot - ss_w) / df_b) / (ss_w / df_w) if df_w and ss_w else 0.0

    obs = f_stat(grouping)
    perm_count = sum(f_stat(np.random.permutation(grouping)) >= obs
                     for _ in range(n_perm))
    return obs, (perm_count + 1) / (n_perm + 1)
```

### CheckM TSV Parser
```python
def read_checkm(path):
    """Parse CheckM lineage_wf output qa file."""
    df = pd.read_csv(path, sep='\t')
    df.columns = [c.strip().lower().replace(' ', '_') for c in df.columns]
    high_q = df[(df['completeness'] >= 90) & (df['contamination'] < 5)]
    med_q  = df[(df['completeness'] >= 50) & (df['contamination'] < 10)]
    print(f"High-quality MAGs: {len(high_q)}")
    print(f"Medium-quality MAGs: {len(med_q)}")
    return df
```

## Pitfalls

- **Host decontamination is critical** — failure to remove host reads inflates classifications
- **Kraken2 database choice** — standard (archaea + bacteria + viral) vs PlusPF (adds protozoa/fungi) vs custom
- **Bracken threshold** — `-t 10` requires ≥10 reads assigned to a taxon; lower for low-coverage samples
- **HUMAnN3 databases** — must download UniRef90 and ChocoPhlAn databases (~25 GB) separately
- **Binning quality** — MetaBAT2 needs ≥2× coverage and ≥500 bp contigs; CONCOCT or MaxBin2 can be combined (DAS Tool for bin refinement)
- **Rarefying without checking curves** — rarefy only when rarefaction curves plateau; otherwise you discard real signal
- **Comparing raw alpha diversity without rarefaction** — deeper-sequenced samples will always appear more diverse
- **PERMANOVA sensitivity to dispersion** — PERMANOVA is sensitive to differences in within-group variance (not just centroids). Pair with a dispersion test (`betadisper` in R)

## Related Skills
- `ngs-variant-calling` — FASTQ QC, short-read alignment, samtools
- `rnaseq` — RNA-seq differential expression, DESeq2/edgeR
- `python-core-bio` — FASTA/FASTQ parsing, sequence manipulation
- `biostatistics-r` — statistical testing, diversity analysis in R
