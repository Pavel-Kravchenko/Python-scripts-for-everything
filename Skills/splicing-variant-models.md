---
name: splicing-variant-models
description: SpliceAI-style donor/acceptor delta scoring and splice variant prioritization workflows.
---

# splicing-variant-models

## When to Use

Use this skill when:
- Prioritizing variants for potential splice disruption
- Computing donor/acceptor gain/loss deltas
- Combining splice-specific and multi-task model outputs

## Quick Reference

| Metric | Meaning |
|---|---|
| DS_DG | Donor gain |
| DS_DL | Donor loss |
| DS_AG | Acceptor gain |
| DS_AL | Acceptor loss |

Practical heuristic: rank by `max(DS_DG, DS_DL, DS_AG, DS_AL)`.

## Key Patterns

**Pattern 1: splice delta computation**
```python
def splice_deltas(seq, pos, alt, donor_score, acceptor_score, w=9):
    half = w // 2
    s, e = max(0, pos-half), min(len(seq), pos+half+1)
    ref_w = seq[s:e]
    alt_w = (seq[:pos] + alt + seq[pos+1:])[s:e]
    d_ref, d_alt = donor_score(ref_w), donor_score(alt_w)
    a_ref, a_alt = acceptor_score(ref_w), acceptor_score(alt_w)
    return {
        'DS_DG': max(0.0, d_alt - d_ref),
        'DS_DL': max(0.0, d_ref - d_alt),
        'DS_AG': max(0.0, a_alt - a_ref),
        'DS_AL': max(0.0, a_ref - a_alt),
    }
```

**Pattern 2: variant ranking by splice impact**
```python
def rank_splice_variants(records):
    return sorted(records, key=lambda r: max(r['ds'].values()), reverse=True)
```

## Code Templates

### Toy donor/acceptor scorers
```python
def donor_score(window):
    c = window[len(window)//2 - 1: len(window)//2 + 1]
    return 0.8 if c == 'GT' else 0.1

def acceptor_score(window):
    c = window[len(window)//2 - 1: len(window)//2 + 1]
    return 0.8 if c == 'AG' else 0.1
```

### Integrate splice and broader signals
```python
def integrated_priority(max_ds, expr_delta, chrom_delta):
    return 0.6 * max_ds + 0.25 * abs(expr_delta) + 0.15 * abs(chrom_delta)
```

## Common Pitfalls

- Assuming canonical GT/AG is sufficient for full splice modeling
- Ranking on one DS channel only
- Ignoring transcript context and exon/intron boundaries

## Related Skills

- `enformer-regulatory-prediction`
- `epigenomic-sequence-models`
- `clinical-genomics`

