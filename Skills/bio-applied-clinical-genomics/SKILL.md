---
name: bio-applied-clinical-genomics
description: "This notebook covers the principles and practice of clinical genomics: how genetic variants are classified, interpreted, and used to guide patient care. Inspired by the Medical Genomics course (V.E. R"
tool_type: python
source_notebook: "Tier_3_Applied_Bioinformatics/11_Clinical_Genomics/01_clinical_genomics.ipynb"
primary_tool: Python
---

## Version Compatibility

Reference examples tested with: Python 3.10+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.


# Clinical Genomics: From Variants to Patient Care

*Source: Course notebook `Tier_3_Applied_Bioinformatics/11_Clinical_Genomics/01_clinical_genomics.ipynb`*

# Clinical Genomics: From Variants to Patient Care

## Tier 3 - Applied Bioinformatics

This notebook covers the principles and practice of clinical genomics: how genetic variants are classified, interpreted, and used to guide patient care. Inspired by the Medical Genomics course (V.E. Ramensky, A.A. Zharikova, MSU).

### Learning Objectives
- Understand the scope of clinical genomics and types of genetic testing
- Apply ACMG/AMP variant classification guidelines
- Query clinical databases (ClinVar, gnomAD) programmatically
- Evaluate in silico variant effect predictors (SIFT, PolyPhen-2, CADD, REVEL, AlphaMissense)
- Understand pharmacogenomics, cancer genomics, and polygenic risk scores
- Recognize ethical considerations in clinical genetic testing

### Prerequisites
- Module 01: NGS Fundamentals
- Module 02: Variant Calling and SNP Analysis
- Basic understanding of Mendelian genetics

---

## 1. Introduction to Clinical Genomics

### 1.1 Precision Medicine and Genomics

**Precision medicine** (sometimes called personalized medicine) uses an individual's genomic information to guide clinical decisions -- from diagnosis to treatment selection to risk assessment.

The journey from the Human Genome Project (completed 2003) to routine clinical sequencing has been driven by:
- **Dramatic cost reduction**: Whole-genome sequencing dropped from ~$100M (2001) to ~$200 (2024)
- **Improved interpretation**: Growing databases of variant-disease associations
- **Clinical validation**: Proven utility in rare disease diagnosis, cancer treatment, pharmacogenomics

```
Precision Medicine Pipeline:

  Patient         Sequencing        Bioinformatics       Clinical
  Sample    --->   (WGS/WES/    --->  Variant       --->  Interpretation
  (blood,          Panel)            Calling &            & Report
   tumor)                            Annotation
                                         |
                                         v
                                    Databases:
                                    ClinVar, OMIM,
                                    gnomAD, PharmGKB
                                         |
                                         v
                                    Treatment
                                    Decision
```

### 1.2 Types of Genetic Testing

| Type | Purpose | Examples | Typical Approach |
|------|---------|---------|------------------|
| **Diagnostic** | Identify cause of existing disease | Rare disease workup, syndromic child | WES/WGS or gene panels |
| **Predictive/Presymptomatic** | Assess future disease risk | BRCA1/2 for breast cancer, Huntington's | Targeted testing |
| **Carrier testing** | Identify heterozygous carriers of recessive conditions | Cystic fibrosis, sickle cell | Carrier panels |
| **Pharmacogenomic** | Guide drug selection/dosing | Warfarin dosing, 5-FU toxicity risk | PGx panels |
| **Prenatal/Newborn** | Screen or diagnose in fetus/newborn | NIPT, newborn screening panels | cfDNA, targeted panels |
| **Somatic/Tumor** | Guide cancer treatment | Actionable mutations in tumors | Tumor panels, WES |

### 1.3 Clinical vs Research Sequencing

Clinical and research sequencing differ in critical ways:

| Aspect | Clinical (CLIA/CAP) | Research |
|--------|-------------------|----------|
| Regulation | CLIA-certified lab required | No specific lab certification |
| Validation | Analytically validated assay | Method may be exploratory |
| Reporting | Standardized clinical report | Varies widely |
| Turnaround | Defined TAT (often 2-4 weeks) | Variable |
| Variants reported | Clinically actionable only | All variants of interest |
| Confirmatory testing | Orthogonal confirmation (Sanger) | Not required |
| Return of results | Duty to report to patient | Not always returned |

```python
# Imports for the entire notebook
import json
import math
from collections import defaultdict
from typing import NamedTuple

# Optional: for ClinVar API queries (Section 3)
try:
    import requests
    HAS_REQUESTS = True
except ImportError:
    HAS_REQUESTS = False
    print("requests not installed -- ClinVar API examples will use mock data")

print("Imports ready.")
```

---

## 2. Variant Classification (ACMG/AMP Guidelines)

The **American College of Medical Genetics and Genomics (ACMG)** and the **Association for Molecular Pathology (AMP)** published landmark guidelines in 2015 (Richards et al., *Genetics in Medicine*) for classifying sequence variants. These guidelines are the global standard for clinical variant interpretation.

### 2.1 The 5-Tier Classification System

Every variant evaluated in a clinical context is assigned to one of five categories:

```
ACMG 5-Tier Classification:

  Pathogenic  >  Likely       >  Variant of    >  Likely   >  Benign
     (P)         Pathogenic      Uncertain         Benign      (B)
                   (LP)          Significance       (LB)
                                   (VUS)

  |-- Reportable as   --|  |-- Typically not  --|  |-- Reportable as  --|
     disease-causing        acted upon clinically     not disease-causing
```

- **Pathogenic (P)**: Strong evidence that the variant causes disease
- **Likely Pathogenic (LP)**: >90% certainty of being disease-causing
- **VUS**: Insufficient evidence to classify -- the most challenging category
- **Likely Benign (LB)**: >90% certainty of being benign
- **Benign (B)**: Strong evidence that the variant does not cause disease

### 2.2 ACMG/AMP Evidence Categories

The classification is based on accumulating evidence from multiple categories:

**1. Population data**
- Allele frequency in large reference populations (gnomAD, ExAC)
- If a variant is common (e.g., MAF > 5%), it is almost certainly benign (BA1 criterion)
- Absence from population databases supports pathogenicity (PM2)

**2. Computational/predictive data**
- In silico predictions (SIFT, PolyPhen-2, CADD, REVEL)
- Conservation across species
- Protein domain/structural impact

**3. Functional data**
- Well-established functional assays showing a damaging effect
- Functional studies in model organisms

**4. Segregation data**
- Co-segregation with disease in affected families
- De novo occurrence in affected individual with confirmed parentage

**5. De novo / allelic data**
- De novo variants in a patient with disease and no family history
- Variant detected in trans with a known pathogenic variant (compound heterozygosity)

**6. Other data**
- Variant reported as pathogenic by reputable source
- Patient phenotype highly specific for gene

### 2.3 ACMG/AMP Criteria Codes

Each piece of evidence is assigned a **code** with a **strength level**:

**Pathogenic criteria:**

| Strength | Codes | Examples |
|----------|-------|----------|
| Very Strong | PVS1 | Null variant in a gene where loss-of-function is a known disease mechanism |
| Strong | PS1-PS4 | Same amino acid change as established pathogenic variant (PS1); de novo with confirmed parentage (PS2); well-established functional study (PS3); prevalence in affecteds vs controls (PS4) |
| Moderate | PM1-PM6 | Located in mutational hotspot (PM1); absent from population databases (PM2); protein length change (PM4); novel missense in gene with low missense rate (PM5); assumed de novo (PM6) |
| Supporting | PP1-PP5 | Co-segregation (PP1); computational evidence (PP3); patient phenotype specific for gene (PP4); reputable source reports pathogenic (PP5) |

**Benign criteria:**

| Strength | Codes | Examples |
|----------|-------|----------|
| Stand-alone | BA1 | Allele frequency >5% in any population |
| Strong | BS1-BS4 | Frequency greater than expected for disorder (BS1); observed in healthy adult (BS2); functional study shows no effect (BS3); non-segregation with disease (BS4) |
| Supporting | BP1-BP7 | Missense in gene where only truncating cause disease (BP1); in trans with pathogenic variant in dominant gene (BP2); in silico predictions suggest benign (BP4); synonymous with no splice impact (BP7) |

```python
# ACMG classification rules (simplified)
# Based on Richards et al. 2015, Table 5

CRITERIA_STRENGTH = {
    # Pathogenic criteria
    'PVS1': 'very_strong',
    'PS1': 'strong', 'PS2': 'strong', 'PS3': 'strong', 'PS4': 'strong',
    'PM1': 'moderate', 'PM2': 'moderate', 'PM3': 'moderate',
    'PM4': 'moderate', 'PM5': 'moderate', 'PM6': 'moderate',
    'PP1': 'supporting', 'PP2': 'supporting', 'PP3': 'supporting',
    'PP4': 'supporting', 'PP5': 'supporting',
    # Benign criteria
    'BA1': 'stand_alone',
    'BS1': 'strong', 'BS2': 'strong', 'BS3': 'strong', 'BS4': 'strong',
    'BP1': 'supporting', 'BP2': 'supporting', 'BP3': 'supporting',
    'BP4': 'supporting', 'BP5': 'supporting', 'BP6': 'supporting',
    'BP7': 'supporting',
}


def classify_variant(criteria: list[str]) -> str:
    """Classify a variant based on ACMG/AMP criteria codes.
    
    Implements the combining rules from Richards et al. 2015, Table 5.
    Returns one of: Pathogenic, Likely Pathogenic, VUS, Likely Benign, Benign.
    """
    path_criteria = [c for c in criteria if c.startswith(('PVS', 'PS', 'PM', 'PP'))]
    benign_criteria = [c for c in criteria if c.startswith(('BA', 'BS', 'BP'))]
    
    # Count by strength
    pvs = sum(1 for c in path_criteria if CRITERIA_STRENGTH.get(c) == 'very_strong')
    ps = sum(1 for c in path_criteria if CRITERIA_STRENGTH.get(c) == 'strong')
    pm = sum(1 for c in path_criteria if CRITERIA_STRENGTH.get(c) == 'moderate')
    pp = sum(1 for c in path_criteria if CRITERIA_STRENGTH.get(c) == 'supporting')
    
    ba = sum(1 for c in benign_criteria if CRITERIA_STRENGTH.get(c) == 'stand_alone')
    bs = sum(1 for c in benign_criteria if CRITERIA_STRENGTH.get(c) == 'strong')
    bp = sum(1 for c in benign_criteria if CRITERIA_STRENGTH.get(c) == 'supporting')
    
    # Benign rules (checked first -- BA1 alone is stand-alone benign)
    if ba >= 1:
        return 'Benign'
    if bs >= 2:
        return 'Benign'
    if bs >= 1 and bp >= 1:
        return 'Likely Benign'
    if bp >= 2:
        return 'Likely Benign'
    
    # Pathogenic rules
    # Pathogenic: (i) 1 Very Strong AND (>=1 Strong OR >=2 Moderate OR 1 Moderate+1 Supporting OR >=2 Supporting)
    if pvs >= 1:
        if (ps >= 1 or pm >= 2 or (pm >= 1 and pp >= 1) or pp >= 2):
            return 'Pathogenic'
    # Pathogenic: (ii) >=2 Strong
    if ps >= 2:
        return 'Pathogenic'
    # Pathogenic: (iii) 1 Strong AND (>=3 Moderate OR 2 Moderate+>=2 Supporting OR 1 Moderate+>=4 Supporting)
    if ps >= 1:
        if (pm >= 3 or (pm >= 2 and pp >= 2) or (pm >= 1 and pp >= 4)):
            return 'Pathogenic'
    
    # Likely Pathogenic
    # (i) 1 Very Strong AND 1 Moderate
    if pvs >= 1 and pm >= 1:
        return 'Likely Pathogenic'
    # (ii) 1 Strong AND 1-2 Moderate
    if ps >= 1 and 1 <= pm <= 2:
        return 'Likely Pathogenic'
    # (iii) 1 Strong AND >=2 Supporting
    if ps >= 1 and pp >= 2:
        return 'Likely Pathogenic'
    # (iv) >=3 Moderate
    if pm >= 3:
        return 'Likely Pathogenic'
    # (v) 2 Moderate AND >=2 Supporting
    if pm >= 2 and pp >= 2:
        return 'Likely Pathogenic'
    # (vi) 1 Moderate AND >=4 Supporting
    if pm >= 1 and pp >= 4:
        return 'Likely Pathogenic'
    # (vii) 1 Very Strong alone -> LP
    if pvs >= 1:
        return 'Likely Pathogenic'
    
    return 'VUS'


# Test with examples
test_cases = [
    (['PVS1', 'PM2'],                'Null variant + absent from populations'),
    (['PS1', 'PM2', 'PP3'],          'Same AA change + absent + computational'),
    (['PS3', 'PS4'],                 'Two strong pathogenic'),
    (['PM2', 'PP3'],                 'Absent + computational only'),
    (['BA1'],                        'Common variant (>5% freq)'),
    (['BS1', 'BP4'],                 'Higher freq than expected + in silico benign'),
    (['PVS1', 'PS2', 'PM2'],         'Null + de novo + absent'),
]

print(f'{"Criteria":<30} {"Classification":<20} Description')
print('-' * 85)
for criteria, desc in test_cases:
    result = classify_variant(criteria)
    print(f'{", ".join(criteria):<30} {result:<20} {desc}')
```
