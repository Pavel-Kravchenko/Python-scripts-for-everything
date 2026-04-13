---
name: bio-applied-microbial-diversity
description: "Microorganisms represent the vast majority of biological diversity on Earth. Studying microbial communities -- who is there, how many, and how they differ across environments -- is central to fields f"
tool_type: python
source_notebook: "Tier_3_Applied_Bioinformatics/04_Microbial_Diversity/01_microbial_diversity.ipynb"
primary_tool: Matplotlib
---

## Version Compatibility

Reference examples tested with: matplotlib 3.8+, numpy 1.26+, pandas 2.1+, scikit-learn 1.4+, scipy 1.12+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.


# Microbial Diversity Analysis

*Source: Course notebook `Tier_3_Applied_Bioinformatics/04_Microbial_Diversity/01_microbial_diversity.ipynb`*

# Microbial Diversity Analysis

Microorganisms represent the vast majority of biological diversity on Earth. Studying microbial communities -- who is there, how many, and how they differ across environments -- is central to fields from ecology to clinical medicine. This notebook covers the computational tools and concepts used to analyze microbial diversity from amplicon sequencing data.

---

## Learning Objectives

By the end of this notebook, you will be able to:

1. Explain the 16S rRNA gene approach to microbial community profiling
2. Distinguish between OTU clustering and ASV denoising
3. Compute and interpret alpha diversity metrics (observed species, Shannon, Simpson, Faith's PD)
4. Construct rarefaction curves and explain their purpose
5. Compute and interpret beta diversity metrics (Bray-Curtis, UniFrac)
6. Apply ordination (PCoA, NMDS) to visualize community differences
7. Perform statistical testing (PERMANOVA, ANOSIM) on distance matrices
8. Create standard microbiome visualizations (stacked bar plots, diversity boxplots, ordination plots)

## 1. The 16S rRNA Gene and Amplicon Sequencing

### Why 16S rRNA?

The 16S ribosomal RNA gene is the standard marker for bacterial identification because it:

- Is present in **all bacteria and archaea** (universal marker)
- Contains **conserved regions** (good primer binding sites) alternating with **variable regions** (V1-V9, which differ between species)
- Has a **large reference database** (SILVA, Greengenes, RDP)
- Is ~1500 bp long -- suitable for sequencing

```
16S rRNA gene (~1500 bp):

  V1   V2   V3   V4   V5   V6   V7   V8   V9
|----|----|----|----|----|----|----|----|----|----|
  C    C    C    C    C    C    C    C    C    C

C = conserved region    V = variable region
```

### Common Primer Targets

| Region | Primers | Length | Notes |
|--------|---------|--------|-------|
| V4 | 515F/806R | ~250 bp | Earth Microbiome Project standard |
| V3-V4 | 341F/785R | ~460 bp | Good taxonomic resolution |
| V1-V3 | 27F/534R | ~500 bp | Alternative target |

### Other Marker Genes

| Marker | Target | Use |
|--------|--------|-----|
| 18S rRNA | Eukaryotes | Fungi, protists |
| ITS (Internal Transcribed Spacer) | Fungi | Better fungal resolution than 18S |
| Shotgun metagenomics | All DNA | Full functional + taxonomic profiling |

### Amplicon Sequencing Workflow

```
Environmental sample (soil, gut, water, etc.)
      |
      v
DNA extraction
      |
      v
PCR amplification of 16S rRNA gene (target V region)
      |
      v
Library preparation (add barcodes + adapters)
      |
      v
Sequencing (Illumina MiSeq: paired-end 250-300 bp)
      |
      v
Bioinformatics pipeline:
  - Quality filtering
  - Denoising (ASVs) or clustering (OTUs)
  - Taxonomy assignment
  - Diversity analysis
```

## 2. OTU Clustering vs. ASV Denoising

After quality filtering, we need to group sequences into meaningful biological units.

### OTU Clustering (Traditional)

**Operational Taxonomic Units (OTUs)** group sequences by similarity (typically 97%, roughly corresponding to species).

Clustering approaches:
- **De novo**: Cluster sequences against each other (e.g., UCLUST, VSEARCH)
- **Closed reference**: Map sequences to a reference database
- **Open reference**: Closed reference first, then de novo on unmatched reads

The 97% threshold is a convention. The actual relationship between sequence similarity and species boundaries varies across taxonomic groups.

### ASV Denoising (Modern)

**Amplicon Sequence Variants (ASVs)** use statistical models to distinguish true biological sequences from sequencing errors, resolving sequences down to single-nucleotide differences.

| Feature | OTUs (97%) | ASVs |
|---------|-----------|------|
| Resolution | Species-level (approximate) | Single-nucleotide |
| Reproducibility | Depends on clustering algorithm | Consistent across studies |
| Biological meaning | Arbitrary threshold | Exact sequences |
| Database dependency | Some methods need reference | Reference-free |
| Tools | UCLUST, VSEARCH, CD-HIT | DADA2, Deblur, UNOISE3 |

**ASVs are now the recommended approach** (Callahan et al., 2017).

### DADA2 Algorithm (Simplified)

DADA2 builds an error model from the data:
1. Learn the error rates for each possible nucleotide transition at each quality score
2. For each unique sequence, test whether it could be explained as an error from a more abundant sequence
3. Output: a set of corrected (denoised) sequences with exact counts

## 3. Taxonomy Assignment

After obtaining OTUs or ASVs, we assign taxonomy to each sequence.

### Methods

| Method | Approach | Example |
|--------|----------|--------|
| Naive Bayes classifier | Machine learning on k-mer frequencies | QIIME2 classify-sklearn |
| BLAST/VSEARCH | Sequence similarity search | Top hit or consensus |
| Phylogenetic placement | Place on reference tree | pplacer, EPA-ng |

### Taxonomy Ranks

```
Kingdom: Bacteria
  Phylum: Firmicutes
    Class: Bacilli
      Order: Lactobacillales
        Family: Lactobacillaceae
          Genus: Lactobacillus
            Species: Lactobacillus acidophilus
```

### Reference Databases

| Database | Size | Notes |
|----------|------|-------|
| SILVA | ~2M sequences | Most comprehensive, regularly updated |
| Greengenes2 | Standardized taxonomy | Uses genome phylogeny |
| RDP | ~3M sequences | Long-standing, Naive Bayes classifier |
| GTDB | Genome-based | Gold standard for genome taxonomy |

```python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from scipy import stats
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, dendrogram
from sklearn.manifold import MDS
from itertools import combinations

np.random.seed(42)

plt.rcParams['figure.dpi'] = 100
plt.rcParams['font.size'] = 11
```

## 5. The Feature Table

The central data structure in microbial diversity analysis is the **feature table** (also called an OTU table or BIOM table):

- Rows = features (OTUs or ASVs)
- Columns = samples
- Values = counts (number of reads assigned to each feature in each sample)

Let's create a simulated dataset representing a human microbiome study with samples from three body sites.

```python
def simulate_microbiome_data(n_samples_per_site=8, n_taxa=50, total_reads_range=(5000, 20000)):
    """
    Simulate a microbiome feature table with samples from three body sites.
    Each body site has a characteristic community composition.
    """
    sites = ['Gut', 'Skin', 'Oral']
    n_samples = n_samples_per_site * len(sites)
    
    # Define taxa with phylum-level classification
    taxa_info = {
        'Bacteroides_sp1': 'Bacteroidetes', 'Bacteroides_sp2': 'Bacteroidetes',
        'Prevotella_sp1': 'Bacteroidetes', 'Prevotella_sp2': 'Bacteroidetes',
        'Alistipes_sp1': 'Bacteroidetes',
        'Faecalibacterium_sp1': 'Firmicutes', 'Roseburia_sp1': 'Firmicutes',
        'Ruminococcus_sp1': 'Firmicutes', 'Clostridium_sp1': 'Firmicutes',
        'Clostridium_sp2': 'Firmicutes', 'Lactobacillus_sp1': 'Firmicutes',
        'Lactobacillus_sp2': 'Firmicutes', 'Streptococcus_sp1': 'Firmicutes',
        'Streptococcus_sp2': 'Firmicutes', 'Streptococcus_sp3': 'Firmicutes',
        'Staphylococcus_sp1': 'Firmicutes', 'Staphylococcus_sp2': 'Firmicutes',
        'Enterococcus_sp1': 'Firmicutes', 'Veillonella_sp1': 'Firmicutes',
        'Dialister_sp1': 'Firmicutes',
        'Bifidobacterium_sp1': 'Actinobacteria', 'Bifidobacterium_sp2': 'Actinobacteria',
        'Corynebacterium_sp1': 'Actinobacteria', 'Corynebacterium_sp2': 'Actinobacteria',
        'Propionibacterium_sp1': 'Actinobacteria', 'Cutibacterium_sp1': 'Actinobacteria',
        'Escherichia_sp1': 'Proteobacteria', 'Haemophilus_sp1': 'Proteobacteria',
        'Neisseria_sp1': 'Proteobacteria', 'Pseudomonas_sp1': 'Proteobacteria',
        'Acinetobacter_sp1': 'Proteobacteria',
        'Fusobacterium_sp1': 'Fusobacteria', 'Fusobacterium_sp2': 'Fusobacteria',
        'Porphyromonas_sp1': 'Bacteroidetes', 'Tannerella_sp1': 'Bacteroidetes',
        'Treponema_sp1': 'Spirochaetes',
    }
    # Add more taxa to reach n_taxa
    existing = len(taxa_info)
    phyla = ['Firmicutes', 'Bacteroidetes', 'Actinobacteria', 'Proteobacteria',
             'Verrucomicrobia', 'Tenericutes']
    for i in range(n_taxa - existing):
        phylum = phyla[i % len(phyla)]
        taxa_info[f'Taxon_{i+1:03d}'] = phylum
    
    taxa_names = list(taxa_info.keys())[:n_taxa]
    
    # Define characteristic community profiles for each body site
    # (relative abundance templates -- Dirichlet concentration parameters)
    profiles = {
        'Gut': np.array([15, 10, 8, 5, 4, 12, 8, 6, 5, 3, 2, 1, 1, 1, 0.5,
                         0.2, 0.1, 1, 1, 0.5, 6, 4, 0.5, 0.1, 0.1, 0.1, 2,
                         0.5, 0.1, 0.1, 0.1, 0.5, 0.2, 0.5, 0.2, 0.1] +
                        [0.3] * (n_taxa - 36)),
        'Skin': np.array([0.5, 0.3, 0.2, 0.1, 0.1, 0.5, 0.3, 0.2, 0.5, 0.3,
                          0.5, 0.3, 1, 0.5, 0.2, 15, 12, 0.3, 0.1, 0.1,
                          0.5, 0.3, 10, 8, 8, 6, 0.5, 0.3, 0.2, 1, 2,
                          0.2, 0.1, 0.5, 0.2, 0.1] +
                         [0.3] * (n_taxa - 36)),
        'Oral': np.array([1, 0.5, 3, 2, 0.3, 1, 0.5, 0.3, 0.5, 0.2,
                          2, 1, 12, 10, 8, 0.5, 0.3, 0.5, 5, 3,
                          0.5, 0.3, 0.5, 0.3, 0.2, 0.1, 1, 6, 5, 1,
                          0.5, 8, 5, 3, 2, 1] +
                         [0.3] * (n_taxa - 36)),
    }
    
    # Generate samples
    sample_names = []
    sample_sites = []
    sample_subjects = []
    counts_data = []
    
    for site in sites:
        for i in range(n_samples_per_site):
            sample_name = f"{site}_{i+1}"
            sample_names.append(sample_name)
            sample_sites.append(site)
            sample_subjects.append(f"Subject_{(i % 4) + 1}")
            
            total_reads = np.random.randint(*total_reads_range)
            # Dirichlet-multinomial: sample relative abundances, then counts
            rel_abund = np.random.dirichlet(profiles[site][:n_taxa] + 0.1)
            sample_counts = np.random.multinomial(total_reads, rel_abund)
            counts_data.append(sample_counts)
    
    feature_table = pd.DataFrame(
        np.array(counts_data).T,
        index=taxa_names,
        columns=sample_names
    )
    
    sample_metadata = pd.DataFrame({
        'body_site': sample_sites,
        'subject': sample_subjects,
    }, index=sample_names)
    
    taxonomy = pd.DataFrame({
        'Phylum': [taxa_info.get(t, 'Unknown') for t in taxa_names]
    }, index=taxa_names)
    
    return feature_table, sample_metadata, taxonomy


# Generate the dataset
feature_table, sample_metadata, taxonomy = simulate_microbiome_data()

print(f"Feature table: {feature_table.shape[0]} taxa x {feature_table.shape[1]} samples")
print(f"\nSample metadata:")
print(sample_metadata.head(6))
print(f"\nSamples per body site:")
print(sample_metadata['body_site'].value_counts())
print(f"\nLibrary sizes range: {feature_table.sum(axis=0).min()} - {feature_table.sum(axis=0).max()}")
```

```python
# Preview the feature table
print("Feature table (first 8 taxa, first 6 samples):")
feature_table.iloc[:8, :6]
```

## 6. Alpha Diversity: Who Is There? How Many?

**Alpha diversity** measures the diversity *within* a single sample. Two key aspects:

- **Richness**: How many different types are present?
- **Evenness**: How evenly are abundances distributed?

A sample with 100 species all at equal abundance is more diverse (in terms of evenness) than one with 100 species where a single species dominates 99% of reads.

### Qualitative vs. Quantitative

- **Qualitative** metrics consider only presence/absence (e.g., observed species)
- **Quantitative** metrics also consider relative abundance (e.g., Shannon, Simpson)

### Phylogenetic vs. Non-phylogenetic

- **Non-phylogenetic**: treat all taxa as equally related
- **Phylogenetic**: incorporate the evolutionary tree (e.g., Faith's PD)
