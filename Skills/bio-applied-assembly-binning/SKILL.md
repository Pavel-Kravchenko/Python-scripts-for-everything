---
name: bio-applied-assembly-binning
description: "**Tier 3 — Applied Bioinformatics | Module 26 · Notebook 3**"
tool_type: python
source_notebook: "Tier_3_Applied_Bioinformatics/26_Metagenomics_Shotgun/03_assembly_binning.ipynb"
---

# Metagenomic Assembly, Binning & MAGs

*Source: Course notebook `Tier_3_Applied_Bioinformatics/26_Metagenomics_Shotgun/03_assembly_binning.ipynb`*

# Metagenomic Assembly, Binning & MAGs

**Tier 3 — Applied Bioinformatics | Module 26 · Notebook 3**

*Prerequisites: Notebook 1 (Taxonomic Profiling), Module 17 (Genome Assembly)*

---

**By the end of this notebook you will be able to:**
1. Assemble a metagenome with MEGAHIT
2. Map reads back to contigs to estimate coverage (for binning)
3. Bin contigs into metagenome-assembled genomes (MAGs) with MetaBAT2
4. Assess MAG quality with CheckM (completeness and contamination)
5. Annotate high-quality MAGs with Prokka



**Key resources:**
- [MEGAHIT documentation](https://github.com/voutcn/megahit)
- [MetaBAT2 documentation](https://bitbucket.org/berkeleylab/metabat)
- [CheckM documentation](https://github.com/Ecogenomics/CheckM)

## 1. Metagenomic Assembly with MEGAHIT

**MEGAHIT** assembles metagenomes using a **succinct de Bruijn graph** built from multiple k-mer sizes. Unlike single-genome assemblers that use one k-mer length, MEGAHIT iterates from a small k (e.g., 21) to a large k (e.g., 141), merging assemblies at each step. This multi-k strategy handles the extreme variation in genome coverage present in metagenomes (a highly abundant organism at 1000x coverage vs a rare organism at 1x coverage).

**Key parameters:**
- `--min-contig-len 500`: discard contigs shorter than 500 bp (most contain too little signal for binning)
- `-t 8`: threads
- `--k-list 21,29,39,59,79,99,119,141`: k-mer sizes to iterate (default; can be shortened for speed)

```bash
megahit \
    -1 decontam_1.fastq.gz \
    -2 decontam_2.fastq.gz \
    -o megahit_assembly/ \
    --min-contig-len 500 \
    -t 16 \
    -m 0.5   # use at most 50% of available RAM

# Check assembly statistics
grep '>' megahit_assembly/final.contigs.fa | wc -l   # number of contigs
seqkit stats megahit_assembly/final.contigs.fa         # N50, total length, etc.
```

**Alternative assembler — SPAdes metaSPAdes:** More accurate for low-coverage organisms but requires significantly more RAM (~100–400 GB for complex communities). Use when RAM is available and accuracy is paramount.

**Expected output:** For a typical gut metagenome (5 Gbp data), MEGAHIT typically produces ~100,000–500,000 contigs, with N50 of 1–10 kb. Only contigs ≥1.5 kb are typically used for binning.

```python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

rng = np.random.default_rng(42)

# Simulate MEGAHIT assembly statistics
# Typical gut metagenome (50M paired reads, 5 Gbp data)
n_contigs = 320_000
contig_lengths = rng.lognormal(mean=6.8, sigma=1.1, size=n_contigs).astype(int)
contig_lengths = np.sort(np.clip(contig_lengths, 500, 500_000))[::-1]

# Compute N50
total_length = contig_lengths.sum()
cumsum = np.cumsum(contig_lengths)
n50_idx = np.searchsorted(cumsum, total_length / 2)
n50 = contig_lengths[n50_idx]

print('=== MEGAHIT Assembly Statistics ===')
print(f'Total contigs:          {n_contigs:>10,}')
print(f'Total assembly length:  {total_length/1e9:>10.2f} Gbp')
print(f'N50:                    {n50:>10,} bp')
print(f'Longest contig:         {contig_lengths[0]:>10,} bp')
print(f'Contigs >= 1500 bp:     {(contig_lengths >= 1500).sum():>10,}')
print(f'Contigs >= 10000 bp:    {(contig_lengths >= 10000).sum():>10,}')

fig, axes = plt.subplots(1, 2, figsize=(13, 5))

# Contig length distribution
bins = np.logspace(np.log10(500), np.log10(500_000), 50)
axes[0].hist(contig_lengths, bins=bins, color='steelblue', edgecolor='white', alpha=0.8)
axes[0].axvline(n50, color='red', linestyle='--', lw=1.5, label=f'N50 = {n50:,} bp')
axes[0].axvline(1500, color='orange', linestyle=':', lw=1.5, label='1500 bp threshold')
axes[0].set_xscale('log')
axes[0].set_xlabel('Contig length (bp, log scale)')
axes[0].set_ylabel('Count')
axes[0].set_title('MEGAHIT contig length distribution')
axes[0].legend()

# Cumulative length curve (Nx plot)
nx_pcts = np.arange(1, 101)
nx_vals = [contig_lengths[np.searchsorted(cumsum, total_length * (p / 100))]
           for p in nx_pcts]
axes[1].plot(nx_pcts, np.array(nx_vals) / 1000, color='darkorange', lw=2)
axes[1].axvline(50, color='red', linestyle='--', lw=1, label=f'N50 = {n50/1000:.1f} kb')
axes[1].set_xlabel('x (%)')
axes[1].set_ylabel('Nx (kb)')
axes[1].set_title('Nx curve (contig contiguity)')
axes[1].legend()

plt.tight_layout()
plt.show()
```

## 2. Coverage Estimation for Binning

**Binning** (grouping contigs into bins representing individual genomes) relies on two signals:
1. **Tetranucleotide frequency (TNF)**: each genome has a characteristic k-mer composition independent of coverage
2. **Coverage depth**: contigs from the same organism have similar depth across all samples

Coverage depth is computed by aligning the original reads back to the assembled contigs.

```bash
# Build Bowtie2 index for the assembly
bowtie2-build megahit_assembly/final.contigs.fa contigs_index

# Map reads back (produces coverage signal)
bowtie2 \
    -x contigs_index \
    -1 decontam_1.fastq.gz \
    -2 decontam_2.fastq.gz \
    -p 16 \
    --no-unal \
    | samtools sort -@ 8 -o contigs_mapped.bam
samtools index contigs_mapped.bam

# Generate depth file for MetaBAT2
jgi_summarize_bam_contig_depths \
    --outputDepth depths.txt \
    contigs_mapped.bam
# Output: contigName, contigLen, totalAvgDepth, sampleDepth, sampleDepthVar
```

**Multi-sample binning**: If you have multiple samples from the same environment (e.g., time series, or different individuals with similar community), mapping all samples to the same assembly and providing all BAMs to `jgi_summarize_bam_contig_depths` dramatically improves binning quality by providing differential abundance patterns as an additional signal.

```python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

rng = np.random.default_rng(7)

# Simulate contig depth file (jgi_summarize output)
n_contigs_long = 15000  # contigs >= 1500 bp (subset used for binning)
contig_lengths  = rng.lognormal(mean=8, sigma=0.9, size=n_contigs_long).astype(int)
contig_lengths  = np.clip(contig_lengths, 1500, 300_000)

# Simulate 5 organisms with different mean coverages
org_assignments = rng.integers(0, 5, size=n_contigs_long)
org_coverages   = [45, 120, 8, 25, 3]   # depths per organism
gc_contents     = [0.38, 0.55, 0.42, 0.48, 0.61]  # GC% per organism

depths = np.array([
    org_coverages[org] * rng.lognormal(0, 0.15)
    for org in org_assignments
])
gcs = np.array([
    gc_contents[org] + rng.normal(0, 0.02)
    for org in org_assignments
])
gcs = np.clip(gcs, 0, 1)

depth_df = pd.DataFrame({
    'contigLen':     contig_lengths,
    'totalAvgDepth': depths,
    'GC':            gcs,
    'organism':      [f'Org_{o+1}' for o in org_assignments],
})

print(f'Contigs for binning (>= 1500 bp): {n_contigs_long:,}')
print('\nMean depth per organism:')
print(depth_df.groupby('organism')['totalAvgDepth'].agg(['mean', 'count']).round(1))

fig, axes = plt.subplots(1, 2, figsize=(13, 5))

# Coverage distribution
axes[0].hist(np.log10(depths + 1), bins=50, color='#4C72B0', edgecolor='white', alpha=0.8)
axes[0].set_xlabel('log10(coverage + 1)')
axes[0].set_ylabel('Count')
axes[0].set_title('Contig coverage distribution')

# GC vs coverage colored by organism (key binning signal)
colors = ['#4C72B0', '#DD8452', '#55A868', '#C44E52', '#8172B2']
for i, org in enumerate([f'Org_{j+1}' for j in range(5)]):
    mask = depth_df['organism'] == org
    cov_sample = org_coverages[i]
    axes[1].scatter(
        gcs[mask], depths[mask],
        c=colors[i], alpha=0.4, s=5, label=f'{org} ({cov_sample}x)'
    )
axes[1].set_xlabel('GC content')
axes[1].set_ylabel('Coverage depth')
axes[1].set_yscale('log')
axes[1].set_title('GC content vs Coverage (MetaBAT2 binning signal)')
axes[1].legend(markerscale=3, fontsize=8)

plt.tight_layout()
plt.show()
```

## 3. Contig Binning with MetaBAT2

**MetaBAT2** clusters contigs into bins (putative genome-resolved MAGs) using a **Bayesian mixture model** that jointly models tetranucleotide composition and coverage depth from one or more samples.

```bash
# Run MetaBAT2
metabat2 \
    -i megahit_assembly/final.contigs.fa \
    -a depths.txt \
    -o bins/bin \
    --minContig 1500 \   # minimum contig size for binning
    --minClsSize 100000 \ # minimum bin size (100 kb) to report
    -t 8 \
    --saveCls             # save cluster assignment file
# Output: bins/bin.1.fa, bins/bin.2.fa, ... (one FASTA per bin)
```

**Alternative binners:**
- **CONCOCT**: uses k-mer frequency PCA + coverage; good for low-coverage data
- **MaxBin2**: uses marker gene abundance as signal
- **DAS_Tool**: dereplicates and refines bins from multiple binners — produces the highest quality MAG set when all three binners are run and combined

**Practical workflow:**
```bash
# Run all three binners then combine with DAS_Tool
das_tool \
    -i metabat2_bins,concoct_bins,maxbin2_bins \
    -l MetaBAT2,CONCOCT,MaxBin2 \
    -c megahit_assembly/final.contigs.fa \
    -o dastool_output/ \
    -t 8 \
    --write_bins
```

```python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

rng = np.random.default_rng(10)

# Simulate MetaBAT2 binning output
# Represent contigs in (tetranucleotide + coverage) space
n_contigs_long = 15000
org_assignments = rng.integers(0, 5, size=n_contigs_long)
org_coverages   = [45, 120, 8, 25, 3]
gc_contents     = [0.38, 0.55, 0.42, 0.48, 0.61]
colors          = ['#4C72B0', '#DD8452', '#55A868', '#C44E52', '#8172B2']

# Simulate 4-mer embedding as 2D (PCA-reduced from 256-dim space)
org_centers = np.array([
    [1.2, -0.8], [-0.5, 1.5], [2.0, 0.5], [-1.5, -0.5], [0.3, 2.0]
])
tnf_2d = np.array([
    org_centers[o] + rng.normal(0, 0.4, size=2)
    for o in org_assignments
])

depths    = np.array([org_coverages[o] * rng.lognormal(0, 0.15) for o in org_assignments])
gcs       = np.clip([gc_contents[o] + rng.normal(0, 0.02) for o in org_assignments], 0, 1)
contig_lengths = rng.lognormal(mean=8, sigma=0.9, size=n_contigs_long).astype(int)
contig_lengths = np.clip(contig_lengths, 1500, 300_000)

bin_ids   = [f'Bin_{o+1}' for o in org_assignments]

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# TNF 2D projection colored by bin
for i, org in enumerate(range(5)):
    mask = org_assignments == org
    n_contigs_in_bin = mask.sum()
    total_len = contig_lengths[mask].sum()
    axes[0].scatter(tnf_2d[mask, 0], tnf_2d[mask, 1],
                    c=colors[i], alpha=0.3, s=4,
                    label=f'Bin_{org+1} ({n_contigs_in_bin} contigs, {total_len/1e6:.1f} Mb)')
axes[0].set_xlabel('TNF PC1')
axes[0].set_ylabel('TNF PC2')
axes[0].set_title('MetaBAT2: contig clusters (TNF space)')
axes[0].legend(fontsize=7, markerscale=4)

# Summary: bin sizes
bin_df = pd.DataFrame({'bin': bin_ids, 'length': contig_lengths})
bin_summary = bin_df.groupby('bin')['length'].agg(['sum', 'count']).reset_index()
bin_summary.columns = ['Bin', 'Total_length_bp', 'N_contigs']
bin_summary = bin_summary.sort_values('Total_length_bp', ascending=False)
print('MetaBAT2 bin summary:')
print(bin_summary.to_string(index=False))

axes[1].barh(range(len(bin_summary)), bin_summary['Total_length_bp'] / 1e6,
             color=colors[:len(bin_summary)], alpha=0.85)
axes[1].set_yticks(range(len(bin_summary)))
axes[1].set_yticklabels(bin_summary['Bin'])
axes[1].set_xlabel('Total bin length (Mbp)')
axes[1].set_title('MAG sizes (MetaBAT2 bins)')
axes[1].axvline(0.5, color='red', linestyle=':', lw=1.5, label='500 kb threshold')
axes[1].legend()

plt.tight_layout()
plt.show()
```

## 4. MAG Quality Assessment with CheckM

**CheckM** evaluates MAG completeness and contamination using **single-copy marker genes** — a set of ~100 conserved genes expected to appear exactly once in any complete bacterial genome. It selects the appropriate marker gene set based on phylogenetic placement of each bin.

**Quality thresholds (MIMAG standard — Bowers et al. 2017 *Nature Biotechnology*):**

| Quality tier | Completeness | Contamination |
|---|---|---|
| High quality | ≥ 90% | < 5% |
| Medium quality | ≥ 50% | < 10% |
| Low quality | < 50% | — |

```bash
# Run CheckM lineage workflow
checkm lineage_wf \
    bins/ \            # directory containing .fa bins
    checkm_out/ \
    -x fa \            # extension of bin FASTA files
    -t 16 \
    --pplacer_threads 4

# Generate summary table
checkm qa \
    checkm_out/lineage.ms \
    checkm_out/ \
    -o 2 \             # output format 2 = extended table with strain heterogeneity
    > checkm_summary.txt
```

**Contamination vs Strain Heterogeneity:** CheckM's contamination score counts marker genes present > 1 time. However, multi-copy markers can arise from strain mixture (multiple strains in one bin) vs true contamination from a different organism. The **strain heterogeneity** metric measures the proportion of duplicate markers that are ≥90% identical at the amino acid level — high strain heterogeneity suggests strain mixing rather than contamination.
