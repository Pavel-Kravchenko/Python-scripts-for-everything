---
name: bio-applied-hla-typing
description: HLA Typing and Antigen Presentation
tool_type: python
primary_tool: Python
---

# HLA Typing and Antigen Presentation

- [IMGT/HLA database](https://www.ebi.ac.uk/ipd/imgt/hla/)
- [OptiType documentation](https://github.com/FRED-2/OptiType)
- [NetMHCpan 4.1](https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/)
- [pVACtools for neoantigen prioritization](https://pvactools.readthedocs.io/)

## MHC/HLA Biology

The **Major Histocompatibility Complex (MHC)** encodes cell-surface proteins that present peptide fragments to T cells. In humans, these are called **Human Leukocyte Antigens (HLA)**.

### MHC Class I (HLA-A, -B, -C)
- Expressed on **all nucleated cells**
- Presents **8–11 aa intracellular peptides** (from protein degradation in proteasome)
- Recognized by **CD8+ cytotoxic T cells** (via TCR + CD8 co-receptor)
- Key for viral immunity and tumor neoantigen presentation
- Structure: α chain (HLA gene) + β₂-microglobulin (β2M) + peptide

### MHC Class II (HLA-DR, -DQ, -DP)
- Expressed on **professional APCs** (dendritic cells, macrophages, B cells)
- Presents **13–25 aa extracellular/endosomal peptides**
- Recognized by **CD4+ helper T cells**
- Structure: α chain + β chain (both HLA genes) + peptide

### HLA Polymorphism
HLA genes are the **most polymorphic loci** in the human genome:
- HLA-A: >7,500 alleles; HLA-B: >9,500 alleles; HLA-C: >7,300 alleles
- HLA-DRB1: >3,700 alleles
- Polymorphism is concentrated in the **peptide-binding groove** (exons 2-3)
- This variation determines which peptides each individual can present → population-level immune diversity

### Nomenclature (WHO/IMGT-HLA)
Format: `Gene*Field1:Field2:Field3:Field4[Suffix]`

| Resolution | Example | Meaning |
|------------|---------|---------|
| 2-digit (serological) | A2 | Broad antigen group |
| 4-digit (allele group) | A*02:01 | Same protein sequence in binding groove |
| 6-digit | A*02:01:01 | Synonymous coding variants |
| 8-digit | A*02:01:01:01 | Non-coding variants |

**A*02:01** is the most common HLA-A allele in Europeans (~28% frequency).

## HLA Typing from NGS Data

### Computational HLA Typing Tools

| Tool | Method | Input | Resolution |
|------|--------|-------|-----------|
| **OptiType** | Integer linear programming | DNA or RNA reads | 4-digit |
| **HLA*LA** | Graph-based alignment | WGS BAM | 4-6 digit |
| **arcasHLA** | RNA-seq based | RNA-seq BAM | 4-digit |
| **HLA-HD** | Bidirectional alignment | WGS | 6-digit |
| **Polysolver** | Bayesian | WES BAM | 4-digit |

### OptiType Workflow
```python
FASTQ reads
    ↓  (map to HLA reference with RazerS3/yara)
HLA-mapped reads only
    ↓  (integer linear programming)
Best-fitting allele pair per locus
    ↓
Output: sample_result.tsv  →  A*02:01,A*03:01,B*07:02,B*44:02,C*05:01,C*07:02
```

### arcasHLA (from RNA-seq)
More practical for studies with existing RNA-seq data:
```bash
arcasHLA extract --unmapped -o hla_reads/ tumor_rna.bam
arcasHLA genotype hla_reads/sample.extracted.1.fq.gz \
                  hla_reads/sample.extracted.2.fq.gz \
                  -g A,B,C,DPB1,DQB1,DQA1,DRB1 \
                  -o hla_typing/
```

### Homozygosity and Loss of Heterozygosity (LOH)
- Tumor cells can lose one HLA allele (LOH) to escape T-cell recognition
- HLA LOH is detected by comparing tumor vs matched normal HLA typing
- ~15% of NSCLC tumors show HLA LOH

np.random.seed(42)

# Simulate HLA typing results for a cohort of 20 patients (TCGA-like)
hla_a_pool = ['A*02:01', 'A*01:01', 'A*03:01', 'A*24:02', 'A*11:01',
               'A*29:02', 'A*23:01', 'A*26:01', 'A*31:01', 'A*32:01']
hla_b_pool = ['B*07:02', 'B*08:01', 'B*44:02', 'B*44:03', 'B*35:01',
               'B*51:01', 'B*40:01', 'B*15:01', 'B*18:01', 'B*57:01']
hla_c_pool = ['C*07:01', 'C*07:02', 'C*03:04', 'C*05:01', 'C*04:01',
               'C*06:02', 'C*01:02', 'C*02:02', 'C*08:02', 'C*16:01']

a_probs = np.array([0.283, 0.161, 0.143, 0.098, 0.072,
                     0.054, 0.041, 0.038, 0.031, 0.026])
a_probs /= a_probs.sum()
b_probs = np.array([0.126, 0.099, 0.093, 0.068, 0.061,
                     0.058, 0.048, 0.043, 0.039, 0.036])
b_probs /= b_probs.sum()
c_probs = np.ones(10) / 10

n_patients = 20
hla_typing_results = []
for i in range(n_patients):
    hla_typing_results.append({
        'Patient': f'P{i+1:02d}',
        'HLA-A1': np.random.choice(hla_a_pool, p=a_probs),
        'HLA-A2': np.random.choice(hla_a_pool, p=a_probs),
        'HLA-B1': np.random.choice(hla_b_pool, p=b_probs),
        'HLA-B2': np.random.choice(hla_b_pool, p=b_probs),
        'HLA-C1': np.random.choice(hla_c_pool),
        'HLA-C2': np.random.choice(hla_c_pool),
    })

hla_df = pd.DataFrame(hla_typing_results)

# Mark homozygous patients
hla_df['A_homozygous'] = hla_df['HLA-A1'] == hla_df['HLA-A2']
hla_df['B_homozygous'] = hla_df['HLA-B1'] == hla_df['HLA-B2']

print("=== HLA Typing Results (first 10 patients) ===")
print(hla_df[['Patient', 'HLA-A1', 'HLA-A2', 'HLA-B1', 'HLA-B2', 'HLA-C1', 'HLA-C2']].head(10).to_string(index=False))
print(f"\nA homozygous: {hla_df['A_homozygous'].sum()}/{n_patients}")
print(f"B homozygous: {hla_df['B_homozygous'].sum()}/{n_patients}")

# Count HLA-A allele frequencies in cohort
all_a = pd.concat([hla_df['HLA-A1'], hla_df['HLA-A2']]).value_counts()
print(f"\nHLA-A allele counts in cohort:\n{all_a.to_string()}")

## Peptide-MHC Binding Prediction

### Peptide Presentation Pathway (Class I)
1. Protein → **proteasomal degradation** → 8-25 aa peptides
2. Peptides enter ER via **TAP transporter**
3. **MHC-I loading**: peptide trimmed by ERAP1/2 to 8-11 aa, loaded onto MHC-I
4. MHC-I + peptide → cell surface → **TCR surveillance**

### Binding Affinity (IC50 nM)
- Strong binder: IC50 < 50 nM (%Rank < 0.5)
- Weak binder: IC50 50–500 nM (%Rank 0.5–2.0)
- Non-binder: IC50 > 500 nM

### NetMHCpan 4.1
State-of-the-art predictor using pan-allele neural network:
```bash
netMHCpan -p peptides.txt -a HLA-A02:01 -l 9 -BA > predictions.txt
```

Output columns: `Pos | Peptide | Allele | 1-log50k(aff) | Affinity(nM) | %Rank_EL | BindLevel`

### Anchor Positions
Each HLA allele has characteristic **anchor residues** at positions 2 and 9 (for 9-mers):
- **HLA-A*02:01**: P2=L/M (leucine/methionine), P9=V/L (aliphatic)
- **HLA-B*07:02**: P2=P (proline), P9=L/R
- **HLA-A*01:01**: P3=D/E (acidic), P9=Y (tyrosine)

### Neoantigen Pipeline
```python
Tumor somatic SNVs (VCF)
    ↓  VEP/ANNOVAR → peptide window extraction
8-11 mer peptides containing each mutation
    ↓  NetMHCpan with patient-specific HLA type
Rank by binding affinity (%Rank_EL < 0.5)
    ↓  Filter for expression, clonal fraction, foreignness
Ranked neoantigen candidates for vaccine/TCR therapy
```

## HLA Disease Associations

HLA genotype is one of the strongest genetic determinants of disease susceptibility and drug response.

### Disease Associations (Selected)

| HLA Allele | Disease/Phenotype | Odds Ratio | Mechanism |
|------------|------------------|------------|-----------|
| HLA-B*57:01 | Abacavir hypersensitivity (HIV) | >1000 | Abacavir alters peptide repertoire |
| HLA-B*15:02 | Carbamazepine-SJS/TEN | ~80 | Drug-peptide complex |
| HLA-DQ2/DQ8 | Celiac disease | 5–10 | Gluten peptide presentation |
| HLA-B*27 | Ankylosing spondylitis | 90 | Molecular mimicry / arthritogenic peptides |
| HLA-DR4 | Rheumatoid arthritis | 4–6 | Shared epitope hypothesis |
| HLA-A*02:01 | Melanoma immunotherapy response | 1.5–2 | Better neoantigen presentation |
| HLA-DR2 | Multiple sclerosis | 3–4 | Myelin peptide presentation |

### Schymanski-like confidence levels for HLA-disease associations
- **Level 1**: Confirmed by genome-wide significant GWAS + functional mechanism
- **Level 2**: GWAS significant, mechanism unclear
- **Level 3**: Population-level association, not genome-wide

### HLA Supertype Classification
Groups alleles with similar binding specificity:
- **A2 supertype**: A*02:01, A*02:02, A*02:06, A*68:02 — all prefer aliphatic/L anchor at P2
- Practical: if patient lacks A*02:01 but has A*02:06, similar neoantigens may still bind

```python
np.random.seed(42)

# Visualize neoantigen predictions and HLA typing summary
fig, axes = plt.subplots(1, 3, figsize=(16, 5))

# Panel 1: Neoantigen waterfall plot (sorted by rank)
colors = {'Strong Binder': 'firebrick', 'Weak Binder': 'orange', 'Non-binder': 'lightgray'}
sorted_df = candidate_neoantigens.sort_values('Rank_EL_pct')
bar_colors = [colors[b] for b in sorted_df['Binding_Level']]
axes[0].barh(sorted_df['mutation'], -np.log10(sorted_df['IC50_nM'] + 1),
             color=bar_colors, edgecolor='black', linewidth=0.5)
axes[0].axvline(-np.log10(51), color='red', linestyle='--', label='IC50=50nM')
axes[0].axvline(-np.log10(501), color='orange', linestyle='--', label='IC50=500nM')
axes[0].set_xlabel('-log10(IC50 nM)')
axes[0].set_title('Neoantigen MHC-I Binding\n(HLA-A*02:01)')
axes[0].legend(fontsize=8)
patches = [mpatches.Patch(color=v, label=k) for k, v in colors.items()]
axes[0].legend(handles=patches, fontsize=8, loc='lower right')

# Panel 2: HLA-A allele distribution in patient cohort
all_a_alleles = pd.concat([hla_df['HLA-A1'], hla_df['HLA-A2']]).value_counts()
axes[1].bar(all_a_alleles.index, all_a_alleles.values, color='steelblue',
            edgecolor='black', linewidth=0.5)
axes[1].set_xticklabels(all_a_alleles.index, rotation=45, ha='right', fontsize=8)
axes[1].set_ylabel('Count (n=20 patients × 2 alleles)')
axes[1].set_title(f'HLA-A Allele Distribution\n(n={n_patients} patients)')

# Panel 3: HLA-B*57:01 prevalence visualization across populations
populations = ['European', 'East Asian', 'South Asian', 'African', 'Latin American']
b5701_freq = [0.036, 0.002, 0.018, 0.001, 0.012]  # approximate from AFND
axes[2].bar(populations, [f*100 for f in b5701_freq],
            color=['steelblue', 'coral', 'mediumseagreen', 'gold', 'mediumpurple'],
            edgecolor='black', linewidth=0.5)
axes[2].set_ylabel('HLA-B*57:01 Frequency (%)')
axes[2].set_title('HLA-B*57:01 Population Frequencies\n(Abacavir Screen Relevance)')
axes[2].set_xticklabels(populations, rotation=30, ha='right', fontsize=9)
axes[2].axhline(1.0, color='red', linestyle='--', label='1% threshold')
axes[2].legend(fontsize=9)

plt.suptitle('HLA Analysis Summary', fontsize=13, fontweight='bold')
plt.tight_layout()
plt.savefig('hla_analysis_summary.png', dpi=120, bbox_inches='tight')
plt.show()

print("\nKey clinical message:")
print("• Pharmacogenomics: Screen for HLA-B*57:01 BEFORE prescribing abacavir")
print("  European patients: ~3.6% carry this allele → pre-treatment HLA testing is standard of care")
```

## Neoantigen Prediction Pipeline

> Somatic SNVs → peptide sequences (VEP annotation). 8–11-mer peptide windows around mutation. NetMHCpan prediction for patient-specific HLA alleles. Rank percentile threshold (< 0.5% strong binder). pVACseq pipeline.

## HLA Disease Associations

> HLA-B*57:01 and abacavir hypersensitivity. HLA-DQ and celiac disease. HLA-A*02:01 enrichment in melanoma response to immunotherapy. GWAS HLA fine-mapping.

## Summary

| Topic | Key Points |
|-------|-----------|
| HLA-A, -B, -C | Class I; present 8-11 aa peptides to CD8+ T cells |
| HLA-DR, -DQ, -DP | Class II; present 13-25 aa peptides to CD4+ T cells |
| Nomenclature | `A*02:01` = gene A, field1=02 (allele group), field2=01 (protein) |
| OptiType / arcasHLA | Computational HLA typing from NGS reads |
| NetMHCpan | %Rank_EL < 0.5% = strong binder; used for neoantigen filtering |
| HLA-B*57:01 | Mandatory screen before abacavir → prevents severe hypersensitivity |
| HLA LOH | Tumor immune escape; ~15% of solid tumors |

### Module 34 Complete — Key Connections
- V(D)J recombination (Nb1) → diverse TCR CDR3 → recognizes pMHC
- HLA genotype (Nb3) → determines which peptides are presented
- Repertoire clonality (Nb2) → reflects antigen-specific T cell expansion
- Together: the immunogenomics triangle connecting sequence, specificity, and clinical response

## Pitfalls

- **Coordinate systems**: BED uses 0-based half-open; VCF/GFF use 1-based inclusive — mixing them causes off-by-one errors
- **Batch effects**: Always check for batch confounding before interpreting biological signal
- **Multiple testing**: Apply FDR correction (Benjamini-Hochberg) when testing thousands of features simultaneously
