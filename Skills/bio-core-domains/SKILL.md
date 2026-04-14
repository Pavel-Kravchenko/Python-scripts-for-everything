---
name: bio-core-domains
description: Sequence motifs and protein domains — PWM construction, PWM scanning, information content, sequence logos, PROSITE patterns, and Pfam/HMMER concepts
tool_type: python
primary_tool: NumPy
---

# Sequence Motifs and Protein Domains

## Key Concepts

- **Motif**: short functional pattern (<20 aa), does not fold independently (e.g. NLS, phosphorylation site)
- **Domain**: structurally independent unit (50–300 aa), folds on its own, recurs across proteins (e.g. SH2, kinase domain)
- **Profile HMM (Pfam)**: represents a domain family with match, insert, and delete states; scores in bits vs null model

## PWM Construction and Scanning

```python
import numpy as np
import re

BASES = ['A', 'C', 'G', 'T']

def build_ppm(sequences, pseudocount=0.1):
    """Position Probability Matrix from aligned sequences."""
    n_pos = len(sequences[0])
    counts = np.full((4, n_pos), pseudocount)
    for seq in sequences:
        for pos, base in enumerate(seq.upper()):
            if base in BASES:
                counts[BASES.index(base), pos] += 1
    return counts / counts.sum(axis=0)

def ppm_to_pwm(ppm, background=None):
    """Convert PPM to log-odds PWM."""
    if background is None:
        background = np.array([0.25, 0.25, 0.25, 0.25])
    bg = background[:, np.newaxis]
    return np.log2((ppm + 1e-10) / bg)

def information_content(ppm):
    """IC per position in bits. Max = 2 bits for DNA."""
    ic = np.zeros(ppm.shape[1])
    for pos in range(ppm.shape[1]):
        p = ppm[:, pos]
        entropy = -np.sum(p * np.log2(p + 1e-12))
        ic[pos] = 2.0 - entropy
    return ic

def score_sequence(pwm, sequence):
    """Sum of log-odds scores for a sequence against a PWM."""
    return sum(pwm[BASES.index(b), i]
               for i, b in enumerate(sequence.upper()) if b in BASES)

def scan_sequence(pwm, sequence, threshold_pct=0.6):
    """Return (position, subsequence, score) hits above threshold."""
    motif_len = pwm.shape[1]
    max_score = np.sum(np.max(pwm, axis=0))
    threshold = threshold_pct * max_score
    hits = []
    for i in range(len(sequence) - motif_len + 1):
        subseq = sequence[i:i+motif_len]
        s = score_sequence(pwm, subseq)
        if s >= threshold:
            hits.append((i, subseq, s))
    return sorted(hits, key=lambda x: -x[2])

def reverse_complement(seq):
    comp = str.maketrans('ATGCatgc', 'TACGtacg')
    return seq.translate(comp)[::-1]

def scan_both_strands(pwm, sequence, threshold_pct=0.6):
    """Scan forward and reverse strands; report hits in forward coordinates."""
    motif_len = pwm.shape[1]
    fwd = [(p, s, sc, '+') for p, s, sc in scan_sequence(pwm, sequence, threshold_pct)]
    rev = scan_sequence(pwm, reverse_complement(sequence), threshold_pct)
    rev_fwd = [(len(sequence)-p-motif_len, s, sc, '-') for p, s, sc in rev]
    return sorted(fwd + rev_fwd, key=lambda x: -x[2])
```

## PROSITE Pattern Syntax and Conversion

| Syntax | Meaning |
|--------|---------|
| `C` | Specific amino acid |
| `x` | Any amino acid |
| `[LIVMF]` | One of listed |
| `{P}` | Not P |
| `x(3)` | Exactly 3 of any |
| `x(2,4)` | 2 to 4 of any |

```python
def prosite_to_regex(pattern):
    """Convert PROSITE pattern string to Python regex."""
    pattern = pattern.strip('.')
    regex_parts = []
    for elem in pattern.split('-'):
        m = re.match(r'^(.+?)\((\d+)(?:,(\d+))?\)$', elem)
        if m:
            core, low, high = m.group(1), m.group(2), m.group(3)
        else:
            core, low, high = elem, None, None

        if core == 'x':        r = '.'
        elif core.startswith('['):  r = core
        elif core.startswith('{'):  r = f'[^{core[1:-1]}]'
        else:                  r = core

        if low is not None:
            r += f'{{{low},{high}}}' if high else f'{{{low}}}'
        regex_parts.append(r)
    return ''.join(regex_parts)

# Examples:
# N-glycosylation:          'N-{P}-[ST]-{P}'  → 'N[^P][ST][^P]'
# PKC phosphorylation:      '[ST]-x-[RK]'     → '[ST].[RK]'
# Zinc finger C2H2:         'C-x(2,4)-C-x(3)-[LIVMFYWC]-x(8)-H-x(3,5)-H'
```

## Domain Database Overview

| Database | Content | Best use |
|----------|---------|----------|
| Pfam | Profile HMMs, ~20k families | General domain annotation |
| InterPro | Integrates Pfam, SMART, CDD, PROSITE | First-pass comprehensive scan |
| SMART | Signaling domains | Kinases, receptors |
| CDD | NCBI-curated, integrates Pfam/SMART | Free with BLAST |
| PROSITE | Patterns and profiles | Active sites, short motifs |

```bash
# HMMER search
hmmscan --domtblout results.domtbl Pfam-A.hmm protein.fasta

# InterProScan
interproscan.sh -i protein.fasta -f tsv -o results.tsv
```

## Pitfalls

- **Profile HMM E-values**: hmmscan reports two E-values — sequence-level and domain-level; use **domain E-value** when annotating individual domain hits
- **Gathering threshold (GA)**: Pfam uses family-specific GA thresholds, not a universal bit-score cutoff; respect per-family thresholds
- **Domain boundaries are approximate**: Pfam HMM boundaries reflect consensus, not exact structural extents — check the alignment for precise boundaries
- **PROSITE false positives**: short degenerate patterns (e.g. `N-x-[ST]`) match frequently by chance; always cross-validate with structural evidence
- **Scanning both strands**: TF binding sites occur on either strand — always scan reverse complement too
- **PWM threshold choice**: 60% of max score is a starting point; tune based on known sites in your organism's genome
- **Pseudocount matters**: without pseudocounts, a single missing base in the training set gives −∞ log-odds; use ≥0.1 pseudocount per position
