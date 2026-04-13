---
name: bio-applied-mirna-seq-pipeline
description: "**Tier 3 — Applied Bioinformatics | Module 35 · Notebook 1**"
tool_type: python
source_notebook: "Tier_3_Applied_Bioinformatics/35_Small_RNA_and_ncRNA/01_mirna_seq_pipeline.ipynb"
primary_tool: Python
---

## Version Compatibility

Reference examples tested with: scipy 1.12+, statsmodels 0.14+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.


# miRNA-seq Processing and Analysis

*Source: Course notebook `Tier_3_Applied_Bioinformatics/35_Small_RNA_and_ncRNA/01_mirna_seq_pipeline.ipynb`*

# miRNA-seq Processing and Analysis

**Tier 3 — Applied Bioinformatics | Module 35 · Notebook 1**

*Prerequisites: Module 01 (NGS Fundamentals), Module 03 (RNA-seq Analysis)*

---

**By the end of this notebook you will be able to:**
1. Describe miRNA biogenesis: Drosha/DGCR8 → pre-miRNA → Dicer → mature miRNA
2. Process small RNA-seq: adapter trimming, size selection, genome alignment
3. Quantify miRNA expression using miRBase annotation and featureCounts
4. Perform differential miRNA expression analysis with DESeq2
5. Predict miRNA target genes and perform pathway enrichment analysis



**Key resources:**
- [miRBase database](https://www.mirbase.org/)
- [miRDeep2 documentation](https://github.com/rajewsky-lab/mirdeep2)
- [TargetScan](https://www.targetscan.org/)
- [mirTarBase](https://mirtarbase.cuhk.edu.cn/)

## 1. miRNA Biology and Biogenesis

**MicroRNAs (miRNAs)** are ~22 nucleotide non-coding RNAs that post-transcriptionally silence gene expression. They are critical regulators of development, differentiation, and disease.

### Biogenesis Pathway

```
Genomic DNA → RNA Pol II transcription
       ↓
   pri-miRNA (primary transcript, 1-3 kb, capped + polyadenylated)
       ↓  DROSHA/DGCR8 microprocessor complex (nucleus)
   pre-miRNA (precursor, ~60-70 nt stem-loop)
       ↓  Exportin-5 + RanGTP nuclear export
   pre-miRNA (cytoplasm)
       ↓  DICER (RNase III) cleavage
   ~22 nt miRNA duplex (guide + passenger strand)
       ↓  Strand selection into RISC (RNA-Induced Silencing Complex)
   Mature miRNA + AGO2 protein → RISC
       ↓
   mRNA target binding via seed sequence (positions 2-7/8)
       ↓
   Translational repression + mRNA deadenylation/degradation
```

### Key Features
- **Guide strand** (miR-X-5p or miR-X-3p): loaded into RISC, guides target silencing
- **Passenger strand** (miR-X*): usually degraded, occasionally functional
- **Seed region**: positions 2-7 of mature miRNA are critical for target recognition
- **Canonical seed types**: 8mer, 7mer-m8, 7mer-A1, 6mer (decreasing efficacy)

### miRNA Naming Convention
- `miR-21`: mature guide miRNA
- `miR-21-5p`: 5' arm (was miR-21)
- `miR-21-3p`: 3' arm (was miR-21*)
- `hsa-miR-21-5p`: human (Homo sapiens) miR-21, 5' arm
- `mir-21`: the gene/precursor (lowercase)

### miRBase Database
- Current version: miRBase v22 (2018): 1,917 human precursor miRNAs → 2,654 mature miRNAs
- Annotation: chromosome position, sequence, family membership, star/guide designation

## 2. Small RNA-seq Library Design and Processing

### Why Standard RNA-seq Fails for miRNA
- Standard poly-A selection: miRNAs lack poly-A tails → excluded
- Standard fragmentation: miRNAs are already short → further fragmentation → poor alignment
- Solution: **3' adapter ligation-based libraries** that capture all small RNAs

### Library Protocol
1. 3' adapter ligation (Illumina: `TGGAATTCTCGGGTGCCAAGG`) → to 3'-OH of RNA
2. Optional: 5' adapter ligation → for accurate 5' end mapping
3. Reverse transcription with primer complementary to 3' adapter
4. PCR amplification (typically 15-20 cycles)
5. Size selection (PAGE gel): 15-30 nt window → removes unligated adapter
6. Sequencing: 50 bp single-end is sufficient (miRNAs are 22 nt)

### Cutadapt Trimming Parameters
```bash
cutadapt \
  -a TGGAATTCTCGGGTGCCAAGG \  # 3' adapter sequence
  -m 16 \                      # minimum length after trimming
  -M 28 \                      # maximum length (removes longer contaminants)
  --discard-untrimmed \         # discard reads without adapter (unligated)
  -j 8 \                        # threads
  -o trimmed.fastq.gz \
  sample.fastq.gz
```

### Bowtie Alignment to miRBase
```bash
# Build index from miRBase mature.fa
bowtie-build mature.fa mirbase_index

# Align: allow 0-1 mismatch, no reverse complement (miRNAs are strand-specific)
bowtie -x mirbase_index \
  --norc \
  -v 1 \                        # allow 1 mismatch
  -m 5 \                        # discard reads mapping >5 locations
  -p 8 \
  trimmed.fastq.gz \
  -S aligned.sam
```

### Quantification with featureCounts
```bash
featureCounts \
  -a hg38_mirbase_v22.gff3 \
  -o raw_counts.txt \
  -t miRNA \
  -g Name \
  -s 1 \                        # strand-specific
  aligned.bam
```

## 3. Processing Pipeline: Normalization and Quality Control

### Normalization Strategies
1. **TMM (Trimmed Mean of M-values)**: recommended for miRNA-seq (same as mRNA)
2. **DESeq2 median-of-ratios**: robust to outliers
3. **RPM/CPM** (Reads/Counts Per Million): simple but affected by highly expressed miRNAs

### Quality Control Metrics
| Metric | Acceptable threshold |
|--------|---------------------|
| Total raw reads | >5M per sample |
| Adapter trimming rate | >60% (most reads contain adapter) |
| Alignment rate to miRBase | >50% |
| Top expressed miRNA | <50% total (if >50%, consider degradation) |
| CV of technical replicates | <10% |

### Differential Expression: DESeq2 for miRNA
miRNA counts follow a negative binomial distribution (over-dispersed Poisson) — the same statistical model as mRNA-seq. DESeq2 applies a Wald test after size-factor normalization:
$$K_{ij} \sim \text{NB}(\mu_{ij}, \alpha_i)$$

Where $\mu_{ij} = s_j \cdot q_{ij}$ (size factor × expected expression).

```python
# Perform simplified differential expression analysis
# (Mimicking DESeq2 approach with log2FC + t-test as approximation)

from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests

# Normalize counts (CPM normalization for visualization)
lib_sizes = counts_df.sum(axis=0)
cpm_df = (counts_df / lib_sizes * 1e6)

# Log2 normalize for DE analysis
log2_cpm = np.log2(cpm_df + 1)

# Compute fold changes and p-values
normal_cols = [c for c in log2_cpm.columns if c.startswith('Normal')]
tumor_cols = [c for c in log2_cpm.columns if c.startswith('Tumor')]

de_results = []
for mirna in log2_cpm.index:
    normal_vals = log2_cpm.loc[mirna, normal_cols].values
    tumor_vals = log2_cpm.loc[mirna, tumor_cols].values
    
    log2fc = tumor_vals.mean() - normal_vals.mean()
    t_stat, pval = ttest_ind(tumor_vals, normal_vals)
    
    de_results.append({
        'miRNA': mirna,
        'log2FC': log2fc,
        'mean_Normal': normal_vals.mean(),
        'mean_Tumor': tumor_vals.mean(),
        'pvalue': pval,
        'BaseMean': (normal_vals.mean() + tumor_vals.mean()) / 2
    })

de_df = pd.DataFrame(de_results)

# Multiple testing correction
reject, padj, _, _ = multipletests(de_df['pvalue'].fillna(1), method='fdr_bh')
de_df['padj'] = padj
de_df['significant'] = (de_df['padj'] < 0.05) & (de_df['log2FC'].abs() > 1.0)

print("=== Differential miRNA Expression: Tumor vs Normal ===")
print(de_df[['miRNA', 'log2FC', 'pvalue', 'padj', 'significant']].sort_values('log2FC', ascending=False).to_string(index=False, float_format=lambda x: f"{x:.3f}"))
```

## 4. Differential miRNA Expression: Visualization

Volcano plots and heatmaps are the standard visualizations for DE miRNA results.

## 5. Target Prediction and Seed Matching

### Seed Sequence Rules
The miRNA seed region (positions 2-7 of the 5' end) base-pairs with the **3' UTR** of target mRNAs:

| Site type | Definition | Repression |
|-----------|-----------|-----------|
| 8mer | Seed + A at position 1 of 3'UTR | Strongest |
| 7mer-m8 | Seed + match at position 8 | Strong |
| 7mer-A1 | Seed + A at position 1 | Moderate |
| 6mer | Seed only | Weak |

### Target Prediction Databases
- **TargetScan 8.0**: context++ score considers seed type, site accessibility, conservation
- **miRDB**: machine learning-based, provides MirTarget score
- **miRTarBase**: experimentally validated interactions (CLASH, luciferase reporter)
- **DIANA-microT**: combines both sequence and accessibility

### miR-21-5p (oncomiR) Target Examples
Seed: `AGCUUAU` (positions 2-8)

Known validated targets:
- **PTEN** (tumor suppressor): contains 8mer site in 3'UTR → loss of PI3K/Akt regulation
- **PDCD4** (programmed cell death 4): anti-apoptotic → promotes survival
- **RECK**: matrix metalloprotease inhibitor → promotes invasion

### ceRNA Hypothesis (competing endogenous RNA)
lncRNAs and circRNAs can act as "sponges" to sequester miRNAs, indirectly regulating mRNA targets. This creates regulatory networks where non-coding RNAs compete for the same miRNA binding sites.

```python
def find_seed_sites(mirna_seq, utr_seq, seed_start=1, seed_end=8):
    """
    Find seed matches in a 3'UTR sequence.
    mirna_seq: RNA sequence of mature miRNA (5'→3')
    utr_seq: DNA/RNA sequence of 3'UTR
    Returns positions of seed matches.
    """
    # Seed is positions 2-7 (0-indexed: 1-7)
    # The mRNA must be complementary (and reverse) to the miRNA seed
    seed = mirna_seq[seed_start:seed_end]  # positions 2-7
    
    # Complement mapping
    complement = str.maketrans('ACGUTacgut', 'UGCAAgugca'.lower().upper() + 'ugca'.upper())
    seed_rc = seed.translate(complement)[::-1]
    seed_rc_dna = seed_rc.replace('U', 'T')  # convert to DNA for matching
    utr_dna = utr_seq.replace('U', 'T').upper()
    
    positions = []
    for i in range(len(utr_dna) - len(seed_rc_dna) + 1):
        if utr_dna[i:i+len(seed_rc_dna)] == seed_rc_dna:
            site_type = 'unknown'
            # Check for 8mer (check position i-1 = A)
            if i > 0 and utr_dna[i-1] == 'A':
                site_type = '7mer-A1'
            if i + len(seed_rc_dna) < len(utr_dna) and utr_dna[i+len(seed_rc_dna)] == 'A':
                site_type = '7mer-m8'
            if (i > 0 and utr_dna[i-1] == 'A' and
                    i + len(seed_rc_dna) < len(utr_dna)):
                site_type = '8mer'
            positions.append({'position': i, 'site_type': site_type, 'match': utr_dna[i:i+len(seed_rc_dna)]})
    return positions

# Simulate PTEN 3'UTR and look for miR-21-5p seed sites
# miR-21-5p: UAGCUUAUCAGACUGAUGUUGA
mir21_5p = "UAGCUUAUCAGACUGAUGUUGA"
seed_2_7 = mir21_5p[1:7]  # "AGCUUA"
print(f"miR-21-5p sequence: {mir21_5p}")
print(f"Seed region (positions 2-7): {seed_2_7}")

# Simulated PTEN 3'UTR fragment with known miR-21 binding sites
pten_utr_fragment = "AAACCCCCAGAGCCCCCAAGCUUAUAGCUAAGCUUAUCAGAAATAGGTTCCAG"
# embed a real 7mer site
sites = find_seed_sites(mir21_5p, pten_utr_fragment)
print(f"\nSeed match sites in PTEN 3'UTR fragment:")
for site in sites:
    print(f"  Position {site['position']}: {site['match']} ({site['site_type']})")

# Show predicted targets for top upregulated miRNAs
mirna_targets = {
    'hsa-miR-21-5p':  ['PTEN', 'PDCD4', 'RECK', 'SPRY1', 'BTG2', 'NFIB'],
    'hsa-miR-155-5p': ['SHIP1', 'SOCS1', 'AID', 'JARID2', 'C/EBPβ'],
    'hsa-miR-210-3p': ['ISCU', 'COX10', 'GPAM', 'RAD52', 'HIF3A'],
    'hsa-let-7a-5p':  ['KRAS', 'NRAS', 'LIN28A', 'HMGA2', 'MYC', 'IGF1R'],
    'hsa-let-7b-5p':  ['KRAS', 'DICER1', 'HMGA2', 'CDK6', 'CDC25A'],
}

print("\n=== miRNA-Target Network (key examples) ===")
for mirna, targets in mirna_targets.items():
    print(f"{mirna}: targets → {', '.join(targets)}")

print("\nBiological significance:")
print("• miR-21 targets PTEN → activates PI3K/AKT → promotes proliferation/survival")
print("• let-7 targets KRAS/NRAS → inhibits RAS signaling → tumor suppression")
print("• miR-210 targets ISCU → impairs iron-sulfur cluster assembly → metabolic shift")
```

## 3. Processing Pipeline

> Cutadapt for 3' adapter trimming. Length filtering (16-28 nt). Alignment with Bowtie (allow 0-1 mismatch). miRBase GFF annotation. Count with featureCounts.

```python
# Example: Adapter trimming and alignment
# !cutadapt -a TGGAATTCTCGGGTGCCAAGG -m 16 -M 28 -o trimmed.fastq.gz sample.fastq.gz
# !bowtie -x hg38 -p 8 --norc trimmed.fastq.gz -S aligned.sam
# !featureCounts -a miRBase_hg38.gff -o counts.txt aligned.bam
```

## 4. Differential miRNA Expression

> DESeq2 for count-based DE testing (same as RNA-seq). MA plot. Top differentially expressed miRNAs. Heatmap of DE miRNAs. miRNA naming convention (miR-21 vs miR-21-5p).

## 5. Target Prediction and Network Analysis

> TargetScan context++ scores. miRTarBase curated experimental targets. miRNA-target interaction network in NetworkX/Cytoscape. KEGG pathway enrichment of targets. miRNA sponge hypothesis and ceRNA networks.
