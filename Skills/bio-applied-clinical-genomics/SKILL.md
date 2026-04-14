---
name: bio-applied-clinical-genomics
description: "Clinical genomics — ACMG/AMP variant classification, ClinVar queries, and clinical reporting workflows"
tool_type: python
primary_tool: Python
---

# Clinical Genomics: Variant Classification

## Types of Genetic Testing

| Type | Purpose | Typical Approach |
|------|---------|------------------|
| Diagnostic | Identify cause of existing disease | WES/WGS or gene panels |
| Predictive | Assess future disease risk | Targeted testing |
| Carrier | Identify heterozygous carriers | Carrier panels |
| Pharmacogenomic | Guide drug selection/dosing | PGx panels |
| Prenatal/Newborn | Screen or diagnose fetus/newborn | cfDNA, targeted panels |
| Somatic/Tumor | Guide cancer treatment | Tumor panels, WES |

## ACMG/AMP 5-Tier Classification

```
Pathogenic > Likely Pathogenic > VUS > Likely Benign > Benign
  (P)            (LP)                      (LB)          (B)
```

- **LP/P**: >90% certainty disease-causing, reportable
- **VUS**: insufficient evidence, not acted upon clinically
- **LB/B**: >90% certainty benign, reportable as not disease-causing

## ACMG/AMP Evidence Criteria

**Pathogenic:**

| Strength | Codes | Examples |
|----------|-------|----------|
| Very Strong | PVS1 | Null variant in LOF-mechanism gene |
| Strong | PS1-PS4 | Same AA change as known pathogenic; de novo confirmed; functional study; prevalence in affected |
| Moderate | PM1-PM6 | Mutational hotspot; absent from populations; protein length change; novel missense in low-missense gene |
| Supporting | PP1-PP5 | Co-segregation; computational evidence; phenotype-specific; reputable source |

**Benign:**

| Strength | Codes | Examples |
|----------|-------|----------|
| Stand-alone | BA1 | Allele frequency >5% in any population |
| Strong | BS1-BS4 | Frequency > expected; healthy adult carrier; functional no-effect; non-segregation |
| Supporting | BP1-BP7 | Missense in truncating-only gene; in silico benign; synonymous no splice |

## Classification Implementation

```python
CRITERIA_STRENGTH = {
    'PVS1': 'very_strong',
    'PS1': 'strong', 'PS2': 'strong', 'PS3': 'strong', 'PS4': 'strong',
    'PM1': 'moderate', 'PM2': 'moderate', 'PM3': 'moderate',
    'PM4': 'moderate', 'PM5': 'moderate', 'PM6': 'moderate',
    'PP1': 'supporting', 'PP2': 'supporting', 'PP3': 'supporting',
    'PP4': 'supporting', 'PP5': 'supporting',
    'BA1': 'stand_alone',
    'BS1': 'strong', 'BS2': 'strong', 'BS3': 'strong', 'BS4': 'strong',
    'BP1': 'supporting', 'BP2': 'supporting', 'BP3': 'supporting',
    'BP4': 'supporting', 'BP5': 'supporting', 'BP6': 'supporting', 'BP7': 'supporting',
}

def classify_variant(criteria: list[str]) -> str:
    """ACMG/AMP combining rules (Richards et al. 2015, Table 5)."""
    path_criteria = [c for c in criteria if c.startswith(('PVS', 'PS', 'PM', 'PP'))]
    benign_criteria = [c for c in criteria if c.startswith(('BA', 'BS', 'BP'))]

    pvs = sum(1 for c in path_criteria if CRITERIA_STRENGTH.get(c) == 'very_strong')
    ps = sum(1 for c in path_criteria if CRITERIA_STRENGTH.get(c) == 'strong')
    pm = sum(1 for c in path_criteria if CRITERIA_STRENGTH.get(c) == 'moderate')
    pp = sum(1 for c in path_criteria if CRITERIA_STRENGTH.get(c) == 'supporting')
    ba = sum(1 for c in benign_criteria if CRITERIA_STRENGTH.get(c) == 'stand_alone')
    bs = sum(1 for c in benign_criteria if CRITERIA_STRENGTH.get(c) == 'strong')
    bp = sum(1 for c in benign_criteria if CRITERIA_STRENGTH.get(c) == 'supporting')

    # Benign rules (BA1 alone = Benign)
    if ba >= 1: return 'Benign'
    if bs >= 2: return 'Benign'
    if bs >= 1 and bp >= 1: return 'Likely Benign'
    if bp >= 2: return 'Likely Benign'

    # Pathogenic rules
    if pvs >= 1 and (ps >= 1 or pm >= 2 or (pm >= 1 and pp >= 1) or pp >= 2):
        return 'Pathogenic'
    if ps >= 2: return 'Pathogenic'
    if ps >= 1 and (pm >= 3 or (pm >= 2 and pp >= 2) or (pm >= 1 and pp >= 4)):
        return 'Pathogenic'

    # Likely Pathogenic rules
    if pvs >= 1 and pm >= 1: return 'Likely Pathogenic'
    if ps >= 1 and 1 <= pm <= 2: return 'Likely Pathogenic'
    if ps >= 1 and pp >= 2: return 'Likely Pathogenic'
    if pm >= 3: return 'Likely Pathogenic'
    if pm >= 2 and pp >= 2: return 'Likely Pathogenic'
    if pm >= 1 and pp >= 4: return 'Likely Pathogenic'
    if pvs >= 1: return 'Likely Pathogenic'

    return 'VUS'
```

## Pitfalls

- **Coordinate systems**: BED uses 0-based half-open; VCF/GFF use 1-based inclusive — mixing them causes off-by-one errors
- **Batch effects**: Always check for batch confounding before interpreting biological signal
- **Multiple testing**: Apply FDR correction (Benjamini-Hochberg) when testing thousands of features simultaneously
