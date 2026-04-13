---
name: bio-applied-lncrna-classification
description: "**Tier 3 — Applied Bioinformatics | Module 35 · Notebook 2**"
tool_type: python
source_notebook: "Tier_3_Applied_Bioinformatics/35_Small_RNA_and_ncRNA/02_lncrna_classification.ipynb"
primary_tool: Python
---

## Version Compatibility

Reference examples tested with: Python 3.10+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.


# Long Non-Coding RNA: Discovery and Classification

*Source: Course notebook `Tier_3_Applied_Bioinformatics/35_Small_RNA_and_ncRNA/02_lncrna_classification.ipynb`*


**Tier 3 — Applied Bioinformatics | Module 35 · Notebook 2**

*Prerequisites: Notebook 1 (miRNA-seq), Module 03 (RNA-seq Analysis)*

---

**By the end of this notebook you will be able to:**
1. Classify ncRNA types: lncRNA, circRNA, piRNA, snoRNA, snRNA, and their functions
2. Identify novel lncRNAs from RNA-seq using StringTie / TACO assembly
3. Distinguish lncRNA from protein-coding transcripts using coding potential calculators
4. Predict lncRNA functions via co-expression, guilt-by-association
5. Analyze circular RNA (circRNA) from back-splice junctions with CIRI2



**Key resources:**
- [LNCipedia database](https://lncipedia.org/)
- [GENCODE lncRNA annotation](https://www.gencodegenes.org/)
- [CIRI2 for circRNA](https://github.com/bioinfo-biols/CIRI2)
- [CPC2 coding potential](https://cpc2.gao-lab.org/)

## 1. The ncRNA Universe

Non-coding RNAs make up >98% of the human transcriptome by sequence count. They span a huge range in size, biogenesis, and function.

### ncRNA Classification by Size

| Class | Size | Biogenesis | Function |
|-------|------|-----------|---------|
| **miRNA** | ~22 nt | Drosha/Dicer | mRNA silencing via RISC |
| **siRNA** | 21-23 nt | Dicer from dsRNA | mRNA cleavage (exogenous) |
| **piRNA** | 26-31 nt | PIWI-dependent | Transposon silencing in germline |
| **snoRNA** | 60-300 nt | Intron-derived | rRNA/tRNA modification |
| **snRNA** | ~150 nt | Pol II | Splicing (U1, U2, U4, U5, U6) |
| **lncRNA** | >200 nt | Pol II | Diverse: chromatin, splicing, RISC decoy |
| **circRNA** | Variable | Back-splicing | miRNA sponge, protein interaction |

### lncRNA Subtypes

```python
Gene A  ████████████████████
                 ↕ (opposite strand)
lncRNA  ←←←←←←←←←←←←←←   antisense lncRNA

Exon1   Intron  Exon2
████████        ████████    mRNA
        ████████            intronic lncRNA

Gene A  ████    Gene B  ████
        ████████████████    lincRNA (long intergenic ncRNA)
```python

- **lincRNA**: between protein-coding genes (intergenic); most common type
- **Antisense lncRNA**: transcribed from opposite strand, overlapping a coding gene
- **Intronic lncRNA**: within introns of coding genes
- **Enhancer RNA (eRNA)**: transcribed from active enhancers; usually unstable
- **Bidirectional lncRNA**: divergent transcript from bidirectional promoter

### Known lncRNA Examples
| lncRNA | Function | Mechanism |
|--------|---------|-----------|
| **XIST** | X-chromosome inactivation | Coats inactive X, recruits PRC2 |
| **HOTAIR** | Cancer metastasis | Bridges PRC2 and LSD1 complexes |
| **NEAT1** | Paraspeckle structure | Nuclear body scaffold |
| **MALAT1** | Alternative splicing | Regulates SR proteins |
| **H19** | Imprinting | Regulates IGF2 expression |

## 2. lncRNA Discovery from RNA-seq

### StringTie De Novo Assembly
StringTie assembles transcripts from aligned RNA-seq reads and produces GTF files:
```bash
# Step 1: Align with STAR or HISAT2
STAR --genomeDir hg38_index/ --readFilesIn R1.fastq.gz R2.fastq.gz \
     --outSAMtype BAM SortedByCoordinate --outFileNamePrefix sample_

# Step 2: StringTie assembly (each sample)
stringtie sample_Aligned.sortedByCoord.out.bam \
  -G gencode.v44.annotation.gtf \
  -o sample.gtf \
  -p 8 \
  --conservative  # reduces spurious assembly

# Step 3: Merge all samples
stringtie --merge \
  -G gencode.v44.annotation.gtf \
  -o merged.gtf \
  sample1.gtf sample2.gtf sample3.gtf ...

# Step 4: Compare to reference with gffcompare
gffcompare -r gencode.v44.annotation.gtf \
           -G merged.gtf \
           -o comparison
```python

### gffcompare Class Codes for Novel Transcripts

| Code | Meaning | Action |
|------|---------|--------|
| `=` | Exact match to reference | Known transcript |
| `c` | Contained within reference | Partial match |
| `j` | Novel isoform (shares some junctions) | Check carefully |
| `u` | **Intergenic novel transcript** | Candidate lincRNA |
| `x` | Antisense to reference gene | Candidate antisense lncRNA |
| `i` | Completely intronic | Candidate intronic lncRNA |
| `o` | Overlaps reference (other strand) | Antisense candidate |

### Filtering Criteria for lncRNA Candidates
1. Length ≥ 200 nt
2. ≥ 2 exons (multi-exonic preferred for assembly confidence)
3. gffcompare class code: `u`, `x`, or `i`
4. Not overlapping known protein-coding gene on same strand
5. Expressed in ≥ 2 samples (TPM ≥ 0.1)
6. Pass coding potential assessment (CPC2 score < 0 or CPAT coding probability < 0.36)

```python
np.random.seed(7)

# Simulate the output from a full lncRNA discovery pipeline
n_assembled = 8000

# gffcompare class code distribution (realistic)
class_codes = np.random.choice(
    ['=', 'c', 'j', 'u', 'x', 'i', 'o', 'e', 'p'],
    n_assembled,
    p=[0.30, 0.15, 0.12, 0.20, 0.08, 0.06, 0.04, 0.03, 0.02]
)

assembled_transcripts = pd.DataFrame({
    'transcript_id': [f'TCONS_{i+1:06d}' for i in range(n_assembled)],
    'class_code': class_codes,
    'length': np.random.lognormal(np.log(1200), 0.9, n_assembled).astype(int),
    'exon_count': np.random.choice(range(1, 12), n_assembled,
                                    p=[0.3, 0.25, 0.15, 0.1, 0.07, 0.05, 0.03, 0.02, 0.01, 0.01, 0.01]),
    'tpm_mean': np.random.lognormal(np.log(1.5), 1.2, n_assembled),
})

# Apply lncRNA discovery filters
candidate_mask = (
    (assembled_transcripts['length'] >= 200) &
    (assembled_transcripts['exon_count'] >= 2) &
    (assembled_transcripts['class_code'].isin(['u', 'x', 'i'])) &
    (assembled_transcripts['tpm_mean'] >= 0.1)
)

candidates = assembled_transcripts[candidate_mask].copy()
print(f"Total assembled transcripts: {n_assembled}")
print(f"After lncRNA filters: {len(candidates)}")
print(f"\ngffcompare class code distribution in candidates:")
print(candidates['class_code'].value_counts().to_string())
print(f"\n  'u' = intergenic lincRNA")
print(f"  'x' = antisense lncRNA")
print(f"  'i' = intronic lncRNA")

# Visualize
fig, axes = plt.subplots(1, 3, figsize=(14, 4))

# Panel 1: All transcripts by class code
code_counts = pd.Series(class_codes).value_counts()
axes[0].bar(code_counts.index, code_counts.values, color='steelblue', edgecolor='black')
axes[0].set_xlabel('gffcompare Class Code')
axes[0].set_ylabel('Count')
axes[0].set_title(f'Assembled Transcripts\n(n={n_assembled})')
candidate_codes = ['u', 'x', 'i']
for bar, code in zip(axes[0].patches, code_counts.index):
    if code in candidate_codes:
        bar.set_facecolor('firebrick')
axes[0].text(0.97, 0.97, 'Red = lncRNA candidates', transform=axes[0].transAxes,
             ha='right', va='top', fontsize=8, color='firebrick')

# Panel 2: Length distribution of candidates
axes[1].hist(np.log10(candidates['length']), bins=30, color='darkorange', edgecolor='black', alpha=0.8)
axes[1].set_xlabel('log10(Transcript Length, nt)')
axes[1].set_ylabel('Count')
axes[1].set_title('Candidate lncRNA Length Distribution')
median_len = candidates['length'].median()
axes[1].axvline(np.log10(median_len), color='red', linestyle='--',
                 label=f'Median = {int(median_len)} nt')
axes[1].legend()

# Panel 3: Expression level distribution
axes[2].hist(np.log2(candidates['tpm_mean'] + 0.01), bins=30,
              color='mediumseagreen', edgecolor='black', alpha=0.8)
axes[2].set_xlabel('log2(Mean TPM + 0.01)')
axes[2].set_ylabel('Count')
axes[2].set_title('lncRNA Expression Distribution\n(Note: lower than mRNA)')

plt.suptitle('lncRNA Discovery: Assembly Statistics', fontsize=12, fontweight='bold')
plt.tight_layout()
plt.savefig('lncrna_discovery.png', dpi=120, bbox_inches='tight')
plt.show()
```python

## 3. Coding Potential Assessment

The key challenge in lncRNA identification is confidently distinguishing non-coding transcripts from unannotated or short protein-coding genes.

### CPC2 (Coding Potential Calculator 2)
Uses 6 features to compute a coding probability:
1. **Fickett TESTCODE**: nucleotide composition bias (coding vs random)
2. **Hexamer bias**: codon usage in all 3 frames
3. **ORF length**: longest open reading frame
4. **ORF integrity**: presence of start + stop codon
5. **Isoelectric point** of predicted peptide
6. **PhyloCSF score**: evolutionary conservation in CDSs

Output: `label` (coding/noncoding) + `fickett_score` + `hexamer_score` + `coding_probability`

### CPAT (Coding Potential Assessment Tool)
Logistic regression using:
- **Fickett score**: nucleotide composition
- **Hexamer score**: log-odds of codon triplets vs background
- Trained species-specifically on GENCODE annotations

Threshold: human coding probability **> 0.364** = protein-coding

### PhyloCSF (Phylogenetic Codon Substitution Frequencies)
Tests whether a region evolves under protein-coding constraints by comparing codon substitution rates across 58 mammalian genomes. Negative score = non-coding evolution.

### Ribosome Profiling (gold standard)
Ribo-seq shows ribosome occupancy in 3-nt periodicity for translated ORFs. Most lncRNAs show no 3-nt periodicity → confirm non-coding status.

## 4. Conservation Analysis and circRNA

### lncRNA Conservation
lncRNAs show much lower sequence conservation than protein-coding genes:
- ~60% of lncRNAs are primate-specific (not found in mouse)
- Conserved lncRNAs (e.g., XIST, MALAT1, NEAT1) tend to have important functions
- Conservation of **secondary structure** may be more important than primary sequence
- **PhyloP scores** (0 = neutral, positive = conserved, negative = accelerated): lncRNAs average ~0.8 vs mRNAs ~3.2

### circRNA Detection
**Circular RNAs** are formed by **back-splicing**: the 5' end of a downstream exon joins to the 3' end of an upstream exon, creating a covalently closed loop.

Detection signature: reads spanning the **back-splice junction (BSJ)**:
```python
Linear:  ...Exon2–Exon3–Exon4...
Circular: ...Exon4–Exon2... (back-splice)
```python

CIRI2 algorithm:
```bash
# Align to genome first (BSMAP or BWA with chimeric reads)
bwa mem -T 19 hg38.fa R1.fq R2.fq > genome.sam

# CIRI2 detects BSJ reads
perl CIRI2.pl \
  -I genome.sam \
  -O sample_ciri.txt \
  -F hg38.fa \
  -A gencode.v44.gtf \
  -T 8
```python

circRNA quantification: reads per circRNA / total mapped reads × 10⁶ (= CPM)

### Functional Mechanisms of lncRNAs (4 archetypes)
1. **Signal**: expression itself is the function (marker of cell state)
2. **Decoy**: sequesters miRNAs (ceRNA) or transcription factors
3. **Guide**: directs chromatin remodeling complexes (HOTAIR → PRC2)
4. **Scaffold**: assembles protein complexes (NEAT1 → paraspeckle)

```python
# Example: StringTie assembly
# !stringtie sample.bam -o sample.gtf -p 8
# # Merge samples
# !stringtie --merge -G gencode.v44.gtf sample1.gtf sample2.gtf -o merged.gtf
# # Filter for lncRNAs
# !gffcompare -r gencode.v44.gtf merged.gtf
```python

## 3. Coding Potential Assessment

> CPC2 (Coding Potential Calculator). CPAT (Fickett score + hexamer bias). PhyloCSF for evolutionary constraint. Riboseq footprint analysis to confirm non-coding status.

## 4. circRNA Analysis with CIRI2

> Back-splice junction reads as evidence for circularization. CIRI2 alignment to detect BSJ. Quantification with CIRI-quant. circRNA database (circBase) for known entries. Expression relative to linear counterpart.

## 5. lncRNA Co-expression Networks

> WGCNA co-expression modules including lncRNAs. Guilt-by-association: assign lncRNA function from co-expressed protein-coding genes. Identify lncRNA-disease associations (LncDisease, Lnc2Cancer).

## Common Pitfalls

- **Coordinate systems**: BED uses 0-based half-open; VCF/GFF use 1-based inclusive — mixing them causes off-by-one errors
- **Batch effects**: Always check for batch confounding before interpreting biological signal
- **Multiple testing**: Apply FDR correction (Benjamini-Hochberg) when testing thousands of features simultaneously
