---
name: bio-applied-immune-repertoire
description: "**Tier 3 — Applied Bioinformatics | Module 34 · Notebook 2**"
tool_type: python
source_notebook: "Tier_3_Applied_Bioinformatics/34_Immunogenomics/02_immune_repertoire.ipynb"
primary_tool: Matplotlib
---

## Version Compatibility

Reference examples tested with: matplotlib 3.8+, numpy 1.26+, pandas 2.1+, scipy 1.12+, seaborn 0.13+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.


# Immune Repertoire Sequencing and Clonal Evolution

*Source: Course notebook `Tier_3_Applied_Bioinformatics/34_Immunogenomics/02_immune_repertoire.ipynb`*

# Immune Repertoire Sequencing and Clonal Evolution

**Tier 3 — Applied Bioinformatics | Module 34 · Notebook 2**

*Prerequisites: Notebook 1 (V(D)J Biology)*

---

**By the end of this notebook you will be able to:**
1. Process bulk immune repertoire data (TRUST4, MiXCR) from RNA-seq or amplicon data
2. Perform clonal tracking across time points or tissue compartments
3. Infer B-cell lineage trees from somatic hypermutation patterns
4. Identify public clonotypes (shared across individuals) vs private
5. Integrate TCR/BCR clonotypes with scRNA-seq gene expression



**Key resources:**
- [TRUST4 documentation](https://github.com/liulab-dfci/TRUST4)
- [MiXCR documentation](https://docs.milaboratories.com/mixcr/)
- [VDJdb database](https://vdjdb.cdr3.net/)

## 1. Bulk Repertoire Sequencing with TRUST4

TRUST4 (T-cell Receptor and immunoglobulin reconstruction from bulk tumor RNA-seq) extracts immune receptor sequences from standard RNA-seq BAM files — **no targeted amplification required**.

### TRUST4 Workflow

```
Tumor RNA-seq BAM
       ↓
   TRUST4 extracts reads mapping to BCR/TCR loci
       ↓
   de novo assembly of CDR3-containing contigs
       ↓
   V(D)J gene assignment via IMGT alignment
       ↓
   Output: CDR3 amino acid sequence + V/J calls + abundance
```

### TRUST4 Output Format (barcode_report.tsv)
```
#barcode  count_TRA  count_TRB  ...  CDR3_TRA           CDR3_TRB
CELL001   2          3          ...  CAVSAVGGKLTF       CASSGLAGNTGELFF
CELL002   1          2          ...  CILREGQKLIF        CASSYSLAGELF
```

### MiXCR (amplicon-based repertoire)
For targeted amplicon libraries, MiXCR is the gold standard:
```bash
mixcr analyze amplicon \
  --species hsa \
  --starting-material rna \
  --5-end v-primers \
  --3-end j-primers \
  sample_R1.fastq.gz sample_R2.fastq.gz result
```

### Key quality metrics
- **Total productive reads**: reads with valid CDR3 + in-frame
- **Clonotype count**: unique V+J+CDR3 combinations
- **D50 index**: number of top clones covering 50% of reads
- **Productive/total ratio**: typically >80% for good libraries

```python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage, dendrogram

np.random.seed(42)

# Simulate TRUST4-like bulk repertoire output across 4 samples:
# Pre-treatment (2) and Post-treatment (2) tumor biopsies
trbv_genes = ["TRBV2", "TRBV5-1", "TRBV6-5", "TRBV9", "TRBV12-3",
               "TRBV14", "TRBV18", "TRBV20-1", "TRBV25-1", "TRBV28",
               "TRBV29-1", "TRBV4-1", "TRBV7-2", "TRBV30"]
trbj_genes = ["TRBJ1-1", "TRBJ1-2", "TRBJ1-3", "TRBJ1-4", "TRBJ1-5",
               "TRBJ2-1", "TRBJ2-2", "TRBJ2-3", "TRBJ2-5", "TRBJ2-7"]
amino_acids = list("ACDEFGHIKLMNPQRSTVWY")

def make_cdr3(length=None):
    if length is None:
        length = np.random.randint(8, 18)
    mid = "".join(np.random.choice(amino_acids, length-2))
    return "C" + mid + "F"

def make_repertoire(n_clones, bias_vgene=None, bias_strength=1.0, seed_offset=0):
    np.random.seed(42 + seed_offset)
    vgene_probs = np.ones(len(trbv_genes))
    if bias_vgene:
        idx = trbv_genes.index(bias_vgene)
        vgene_probs[idx] *= bias_strength
    vgene_probs /= vgene_probs.sum()
    
    counts = np.sort(np.random.zipf(1.8, n_clones))[::-1].astype(int)
    return pd.DataFrame({
        'clone_id': [f"c{i:04d}" for i in range(n_clones)],
        'v_call': np.random.choice(trbv_genes, n_clones, p=vgene_probs),
        'j_call': np.random.choice(trbj_genes, n_clones),
        'junction_aa': [make_cdr3() for _ in range(n_clones)],
        'count': counts
    })

# 4 samples from 2 patients, pre/post anti-PD1 treatment
samples = {
    'Pt1_Pre':  make_repertoire(150, seed_offset=0),
    'Pt1_Post': make_repertoire(200, bias_vgene="TRBV20-1", bias_strength=6, seed_offset=1),
    'Pt2_Pre':  make_repertoire(130, seed_offset=2),
    'Pt2_Post': make_repertoire(180, bias_vgene="TRBV9",    bias_strength=5, seed_offset=3),
}

for name, df in samples.items():
    total = df['count'].sum()
    df['freq'] = df['count'] / total
    print(f"{name}: {len(df)} clones, {total:,} reads, top clone: {df['freq'].max():.1%}")
```

## 2. Clonal Tracking Across Time Points

A key application of longitudinal repertoire sequencing is tracking clonal persistence, expansion, and contraction across conditions (e.g., pre/post treatment, tumor vs blood).

### Clone Tracking Strategy
1. Match clonotypes by CDR3 amino acid sequence (and optionally V+J gene)
2. Compute relative frequency change: $\log_2(freq_{post}/freq_{pre})$
3. Classify into: **expanded**, **persistent**, **contracted**, **new**, **lost**

### Alluvial/Sankey Plots
Used to visualize how top clones persist or change in frequency between time points. Width of flow = clone frequency.

### D50 Index
The number of clones needed to account for 50% of the total repertoire reads. Lower D50 = more clonal/less diverse.

## 3. Repertoire Overlap Metrics

To quantify shared clonotypes between samples (patients, tissue compartments, time points), several metrics are used:

### Morisita-Horn Index
Most commonly used for immune repertoire overlap. Takes both presence and abundance into account:
$$C_{MH} = \frac{2\sum_i p_{1i} \cdot p_{2i}}{\lambda_1 + \lambda_2}$$

Where $\lambda_k = \sum_i p_{ki}^2$ (Simpson concentration for sample $k$).

- Range: 0 (no overlap) to 1 (identical)
- Robust to differences in sample size
- Sensitive to frequency, not just presence/absence

### Jaccard Index (binary)
$$J = \frac{|S_1 \cap S_2|}{|S_1 \cup S_2|}$$
Ignores clone frequency; only presence/absence.

### Cosine Similarity
$$\cos\theta = \frac{\mathbf{p_1} \cdot \mathbf{p_2}}{|\mathbf{p_1}||\mathbf{p_2}|}$$
Equivalent to Morisita-Horn when both vectors are length-normalized.

### Public vs Private Clonotypes
- **Public clones**: shared across multiple individuals (e.g., CMV-specific TRBV20-1 clones)
- **Private clones**: unique to one individual
- ~1-5% of the repertoire is public for T cells; much less for B cells

## 4. TCR Specificity Prediction and VDJdb

**Public clonotypes** with known antigen specificity are catalogued in databases:
- **VDJdb**: curated TCR sequences with known epitope/HLA context
- **McPAS-TCR**: pathology-associated TCRs
- **IEDB**: all immune epitopes including TCR data

### VDJdb Entry Format
```
CDR3aa       V-gene      J-gene   Species  Epitope     MHC    Source
CASSLAPGATNEKLFF  TRBV6-5  TRBJ1-4  HomoSapiens  GILGFVFTL  HLA-A*02:01  PMID:12345
```

### GLIPH2 / TCRdist for Clustering
Groups TCRs likely to recognize the same antigen:
1. **GLIPH2**: identifies shared CDR3 motifs enriched above background
2. **TCRdist**: CDR3 + germline CDR1/CDR2 distance metric

### CDR3 Hamming Distance
Simple pairwise distance between CDR3s of equal length:
$$d_{Hamming}(s_1, s_2) = \sum_i \mathbf{1}[s_1[i] \neq s_2[i]]$$
Used to group near-identical clones into clonal families.

```python
# Simulate a VDJdb-style lookup + CDR3 pairwise distance clustering
np.random.seed(5)

# Known public TCR sequences for common viral epitopes
vdjdb_public = pd.DataFrame({
    'junction_aa': [
        'CASSLAPGATNEKLFF',   # Influenza M1-specific, HLA-A*02:01
        'CASSLGQGNTEAFF',     # EBV BMLF1-specific, HLA-A*02:01
        'CASSIRSSYEQYF',      # CMV pp65-specific, HLA-A*02:01
        'CASSPGTSGANEQFF',    # SARS-CoV-2 specific
        'CASSVGQGYEQYF',      # CMV variant
    ],
    'epitope': ['GILGFVFTL', 'GLCTLVAML', 'NLVPMVATV', 'YLQPRTFLL', 'NLVPMVATV'],
    'v_call': ['TRBV6-5', 'TRBV20-1', 'TRBV7-2', 'TRBV5-1', 'TRBV7-2'],
    'mhc': ['HLA-A*02:01'] * 5
})

# Add a few of these to our repertoire
til_with_public = samples['Pt1_Post'].copy()
for _, row in vdjdb_public.iterrows():
    if np.random.random() < 0.6:  # 60% chance of having each public clone
        new_row = {
            'clone_id': f"public_{row['epitope'][:4]}",
            'v_call': row['v_call'],
            'j_call': 'TRBJ1-4',
            'junction_aa': row['junction_aa'],
            'count': np.random.randint(50, 500),
            'freq': 0
        }
        til_with_public = pd.concat([til_with_public, pd.DataFrame([new_row])], ignore_index=True)

total = til_with_public['count'].sum()
til_with_public['freq'] = til_with_public['count'] / total

# VDJdb lookup
matched = til_with_public.merge(vdjdb_public[['junction_aa', 'epitope', 'mhc']],
                                  on='junction_aa', how='inner')
print("=== VDJdb Public Clonotype Lookup ===")
print(matched[['junction_aa', 'epitope', 'mhc', 'freq']].to_string(index=False))

# CDR3 Hamming distance on equal-length CDR3s
def hamming(s1, s2):
    if len(s1) != len(s2):
        return np.nan
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))

# Take top 20 clones of length 14 (typical TCRβ CDR3)
cdr3_14 = til_with_public[til_with_public['junction_aa'].str.len() == 14].head(15)
cdr3_seqs = cdr3_14['junction_aa'].values

dist_matrix = np.array([[hamming(a, b) for b in cdr3_seqs] for a in cdr3_seqs], dtype=float)
print(f"\nCDR3 pairwise Hamming distances (top 15 length-14 CDR3s):")
print(f"Median distance: {np.nanmedian(dist_matrix[dist_matrix > 0]):.1f} aa substitutions")

fig, ax = plt.subplots(figsize=(7, 6))
im = ax.imshow(dist_matrix, cmap='YlOrRd')
ax.set_xticks(range(len(cdr3_seqs)))
ax.set_yticks(range(len(cdr3_seqs)))
ax.set_xticklabels(cdr3_seqs, rotation=90, fontsize=7)
ax.set_yticklabels(cdr3_seqs, fontsize=7)
ax.set_title('CDR3 Hamming Distance Matrix\n(Clonal Family Clustering)')
plt.colorbar(im, ax=ax, label='Hamming distance')
plt.tight_layout()
plt.savefig('cdr3_hamming.png', dpi=120, bbox_inches='tight')
plt.show()
```

## 3. B-Cell Lineage Tree Inference

> Germline IGHV sequence as root. SHM distance matrix → phylogenetic tree per clonal family. dowser / IgPhyML for B-cell phylogenetics. Visualize antibody affinity maturation.

## 4. TCR Specificity Prediction

> VDJdb lookup for antigen-specific TCR sequences. GLIPH2/TCRdist for grouping TCRs by specificity. Predict pMHC binding from CDR3 sequence.

## 5. Diversity Comparison and Visualization

A comprehensive repertoire analysis combines multiple metrics into a summary visualization. Below, we compare all four samples across diversity indices.

```python
def shannon_entropy(freqs):
    freqs = np.array(freqs); freqs = freqs[freqs > 0]
    return -np.sum(freqs * np.log(freqs))

def normalized_shannon(freqs):
    n = len([f for f in freqs if f > 0])
    h = shannon_entropy(freqs)
    return h / np.log(n) if n > 1 else 0

def simpson_d(freqs):
    return np.sum(np.array(freqs)**2)

def clonality(freqs):
    return 1 - normalized_shannon(freqs)

# Summary table
summary_rows = []
for name, df in samples.items():
    freqs = df['freq'].values
    summary_rows.append({
        'Sample': name,
        'N_clones': len(df),
        'Shannon_H': f"{shannon_entropy(freqs):.3f}",
        'Norm_H_prime': f"{normalized_shannon(freqs):.3f}",
        'Simpson_D': f"{simpson_d(freqs):.4f}",
        'Clonality': f"{clonality(freqs):.3f}",
        'D50': d50_index(freqs),
        'Top1_pct': f"{freqs.max()*100:.1f}%"
    })

summary_df = pd.DataFrame(summary_rows)
print("=== Repertoire Diversity Summary Table ===")
print(summary_df.to_string(index=False))

# Visualization: multi-metric comparison
fig, axes = plt.subplots(2, 2, figsize=(11, 8))
colors = ['steelblue', 'cornflowerblue', 'firebrick', 'tomato']
metric_titles = [("Normalized Shannon H'", 'Norm_H_prime'),
                  ("Clonality", 'Clonality'),
                  ("Simpson Concentration D", 'Simpson_D'),
                  ("D50 Index", 'D50')]

for ax, (title, col) in zip(axes.flat, metric_titles):
    vals = summary_df[col].astype(float).values
    bars = ax.bar(summary_df['Sample'], vals, color=colors, edgecolor='black', linewidth=0.5)
    ax.set_title(title)
    ax.set_xticklabels(summary_df['Sample'], rotation=20, ha='right')
    for bar, val in zip(bars, vals):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.005,
                f"{val:.3f}", ha='center', va='bottom', fontsize=9)

plt.suptitle('Repertoire Diversity Metrics Across Samples', fontsize=13, fontweight='bold')
plt.tight_layout()
plt.savefig('diversity_summary.png', dpi=120, bbox_inches='tight')
plt.show()
```
