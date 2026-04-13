---
name: bio-core-sequence-motifs
description: "Split from `01_sequence_motifs_and_domains.ipynb` to keep this topic self-contained."
tool_type: python
source_notebook: "Tier_2_Core_Bioinformatics/10_Sequence_Motifs_and_Domains/01_sequence_motifs.ipynb"
---

# Sequence Motifs

*Source: Course notebook `Tier_2_Core_Bioinformatics/10_Sequence_Motifs_and_Domains/01_sequence_motifs.ipynb`*

# Sequence Motifs

Split from `01_sequence_motifs_and_domains.ipynb` to keep this topic self-contained.

**Navigation:** [Topic overview](./01_sequence_motifs_and_domains.ipynb) · [Next: Domains](./02_domains.ipynb)

## How to use this notebook

1. Run cells top-to-bottom in order — later cells depend on earlier ones
2. Run the environment check cell first to confirm all imports work
3. Read the explanatory text before each code cell — the context matters
4. The exercises at the end are designed to be done after reading each section
5. If a code cell requires internet access (Entrez, PDB download), it is marked — these can be skipped if offline

## Complicated moments explained

- **PFM vs PPM vs PWM**: A Position Frequency Matrix (PFM) contains raw counts. Divide by the number of sequences to get the Position Probability Matrix (PPM). Convert PPM to log-odds relative to background to get the Position Weight Matrix (PWM). Only the PWM is suitable for scoring — raw probabilities are too sensitive to small differences near zero.
- **Pseudocounts are essential**: Without pseudocounts, a position where no sequences have base 'G' would give log(0) = -∞ in the PWM. Adding a small pseudocount (typically 0.1-1 times the background frequency) prevents this while barely affecting well-sampled positions. The appropriate pseudocount is debated; 0.5 is a common default.
- **Information content (IC) and entropy**: IC = log2(4) - H, where H is the Shannon entropy. A fully conserved position has H=0, IC=2 bits. A random position has H=2, IC=0. The total IC of a motif (sum over all positions) indicates how specific it is. A CTCF motif has ~15-16 bits of total IC; a TATA box has ~10-12 bits.
- **Threshold selection is empirical**: There is no single correct PWM score threshold. Common approaches: use the p-value of a score (fraction of random sequences that would score that high), use a fixed fraction of the maximum possible score (e.g., 80% of max), or calibrate empirically from ChIP-seq data. The threshold controls the precision-recall tradeoff.
- **PROSITE patterns are qualitative, PWMs are quantitative**: PROSITE patterns specify exact matches or character classes at each position (e.g., `N-{P}-[ST]-{P}` for N-glycosylation). PWMs allow graded scoring. For highly degenerate binding sites, PWMs are far more sensitive than PROSITE patterns.

## Environment check (run this first)

```python
# Environment check
import numpy as np
import matplotlib.pyplot as plt
from collections import Counter
import re

print("Imports ready.")

# Quick motif demo: build a toy PFM
sequences = ["ACGTACGT", "ACGAACGT", "TCGTACGT", "ACGTACGC"]
L = len(sequences[0])
bases = "ACGT"
pfm = np.zeros((4, L))
for seq in sequences:
    for pos, base in enumerate(seq):
        pfm[bases.index(base), pos] += 1

# Convert to PPM with pseudocount
ppm = (pfm + 0.1) / (len(sequences) + 0.4)

# Compute information content
ic = np.sum(ppm * np.log2(ppm / 0.25), axis=0)
print(f"\nToy motif from {len(sequences)} sequences, length {L}:")
print(f"IC per position: {ic.round(2)}")
print(f"Total IC: {ic.sum():.2f} bits (max = {2*L} bits for fully conserved)")
print(f"Consensus: {''.join(bases[np.argmax(ppm[:,i])] for i in range(L))}")
```

```python
BASES = ["A", "C", "G", "T"]


def build_pfm(sequences):
    """Build a Position Frequency Matrix (count matrix) from aligned sequences."""
    length = len(sequences[0])
    pfm = np.zeros((4, length), dtype=int)

    for seq in sequences:
        for pos, base in enumerate(seq.upper()):
            if base in BASES:
                pfm[BASES.index(base), pos] += 1

    return pfm


def pfm_to_ppm(pfm, pseudocount=0.1):
    """Convert PFM to Position Probability Matrix with pseudocounts."""
    n_seq = pfm.sum(axis=0)[0]  # total sequences
    return (pfm + pseudocount) / (n_seq + 4 * pseudocount)


def ppm_to_pwm(ppm, background=0.25):
    """Convert PPM to log-odds PWM."""
    return np.log2(ppm / background)


# Example: build PWM from known binding sites
binding_sites = [
    "ATGACTCA",
    "ATGACTCA",
    "ATGACTTA",
    "GTGACTCA",
    "ATGACTCG",
]

pfm = build_pfm(binding_sites)
ppm = pfm_to_ppm(pfm)
pwm = ppm_to_pwm(ppm)

print("Position Frequency Matrix (counts):")
header = "     " + "  ".join(f"P{i+1}" for i in range(pfm.shape[1]))
print(header)
for i, base in enumerate(BASES):
    print(f"  {base}: {'  '.join(f'{int(v):2d}' for v in pfm[i])}")

print(f"\nPWM score range: {np.sum(np.min(pwm, axis=0)):.2f} to {np.sum(np.max(pwm, axis=0)):.2f}")
```

```python
def display_ppm_heatmap(ppm, title="Position Probability Matrix"):
    """Display the PPM as a heatmap."""
    fig, ax = plt.subplots(figsize=(max(6, ppm.shape[1] * 0.8), 2.5))
    im = ax.imshow(ppm, aspect="auto", cmap="YlOrRd", vmin=0, vmax=1)

    ax.set_yticks(range(4))
    ax.set_yticklabels(BASES)
    ax.set_xticks(range(ppm.shape[1]))
    ax.set_xticklabels([str(i + 1) for i in range(ppm.shape[1])])
    ax.set_xlabel("Position")
    ax.set_title(title)

    # Annotate cells
    for i in range(4):
        for j in range(ppm.shape[1]):
            color = "white" if ppm[i, j] > 0.6 else "black"
            ax.text(j, i, f"{ppm[i, j]:.2f}", ha="center", va="center",
                    fontsize=8, color=color)

    plt.colorbar(im, ax=ax, label="Probability")
    plt.tight_layout()
    plt.show()


display_ppm_heatmap(ppm, "PPM for TF Binding Sites")
```

### 2.4 Scoring Sequences with a PWM

```python
def score_sequence(pwm, sequence):
    """Score a sequence against a PWM. Returns sum of log-odds scores."""
    score = 0.0
    for pos, base in enumerate(sequence.upper()):
        if base in BASES:
            score += pwm[BASES.index(base), pos]
    return score


def scan_sequence(pwm, sequence, threshold=None):
    """
    Scan a long sequence for PWM matches.
    Returns list of (position, subsequence, score) sorted by score descending.
    """
    motif_len = pwm.shape[1]
    max_score = np.sum(np.max(pwm, axis=0))
    if threshold is None:
        threshold = 0.6 * max_score  # default: 60% of max

    hits = []
    for i in range(len(sequence) - motif_len + 1):
        subseq = sequence[i:i + motif_len]
        s = score_sequence(pwm, subseq)
        if s >= threshold:
            hits.append((i, subseq, s))

    return sorted(hits, key=lambda x: -x[2])


# Test scoring
test_seqs = ["ATGACTCA", "GTGACTCA", "CCCCCCCC", "ATGACTTA"]
max_s = np.sum(np.max(pwm, axis=0))

print(f"{'Sequence':<12} {'Score':>7} {'% of max':>9}")
print("-" * 30)
for seq in test_seqs:
    s = score_sequence(pwm, seq)
    print(f"{seq:<12} {s:7.2f} {100 * s / max_s:8.1f}%")
```

```python
# Scan a longer sequence
target = "GCCTAGATGACTCAGGTTTCCCGTGACTCAATGCAATGACTTACCC"

hits = scan_sequence(pwm, target, threshold=0.5 * np.sum(np.max(pwm, axis=0)))

print(f"Scanning: {target}")
print(f"\n{'Pos':>4} {'Subseq':<10} {'Score':>7}")
print("-" * 25)
for pos, subseq, score in hits:
    print(f"{pos:4d} {subseq:<10} {score:7.2f}")
```

### 2.5 Both-Strand Scanning

Transcription factor binding sites can occur on either DNA strand. To search both
strands, also scan the reverse complement.

```python
def reverse_complement(seq):
    """Return the reverse complement of a DNA sequence."""
    comp = str.maketrans("ATGCatgc", "TACGtacg")
    return seq.translate(comp)[::-1]


def scan_both_strands(pwm, sequence, threshold=None):
    """Scan both strands for PWM matches."""
    fwd_hits = scan_sequence(pwm, sequence, threshold)
    rc_seq = reverse_complement(sequence)
    rev_hits = scan_sequence(pwm, rc_seq, threshold)

    # Convert reverse hit positions to forward-strand coordinates
    motif_len = pwm.shape[1]
    seq_len = len(sequence)
    all_hits = [(pos, subseq, score, "+") for pos, subseq, score in fwd_hits]
    for pos, subseq, score in rev_hits:
        fwd_pos = seq_len - pos - motif_len
        all_hits.append((fwd_pos, subseq, score, "-"))

    return sorted(all_hits, key=lambda x: -x[2])


dual_hits = scan_both_strands(pwm, target,
                               threshold=0.5 * np.sum(np.max(pwm, axis=0)))

print(f"{'Pos':>4} {'Strand':>6} {'Match':<10} {'Score':>7}")
print("-" * 32)
for pos, subseq, score, strand in dual_hits:
    print(f"{pos:4d} {strand:>6} {subseq:<10} {score:7.2f}")
```

---

## 3. Information Content and Sequence Logos

### 3.1 Information Content

The **information content** (IC) at each position quantifies how conserved it is:

$$IC_i = \log_2(K) - H_i = \log_2(K) + \sum_{b} p_{b,i} \log_2(p_{b,i})$$

For DNA ($K=4$): $IC$ ranges from 0 bits (uniform, no conservation) to 2 bits
(perfectly conserved).

The **total information content** of a motif is $\sum_i IC_i$ and reflects how specific
the motif is overall.

```python
def information_content(ppm):
    """Calculate information content per position (DNA: max 2 bits)."""
    ic = np.zeros(ppm.shape[1])
    for pos in range(ppm.shape[1]):
        probs = ppm[:, pos]
        entropy = -np.sum(probs * np.log2(probs + 1e-12))
        ic[pos] = 2.0 - entropy  # 2 = log2(4)
    return ic


ic = information_content(ppm)

print(f"{'Position':>8} {'IC (bits)':>10} {'Bar'}")
print("-" * 35)
for i, val in enumerate(ic):
    bar = '#' * int(val * 15)
    print(f"{i + 1:8d} {val:10.3f} {bar}")
print(f"\nTotal IC: {ic.sum():.2f} bits (max possible: {2 * len(ic):.0f} bits)")
```

### 3.2 Sequence Logos

A **sequence logo** is the standard visualization of a motif. At each position:
- The total height of the stack equals the information content (in bits)
- Each letter's height is proportional to its frequency: $h_{b,i} = p_{b,i} \times IC_i$
- Letters are stacked from most to least frequent (top to bottom)

We can create a basic logo with matplotlib (for publication quality, use
[WebLogo](https://weblogo.berkeley.edu/) or the `logomaker` Python package).

```python
BASE_COLORS = {"A": "green", "C": "blue", "G": "orange", "T": "red"}


def plot_sequence_logo(ppm, title="Sequence Logo"):
    """
    Plot a simple sequence logo using stacked colored bars.

    Each bar's height = p(base, pos) * IC(pos).
    """
    ic = information_content(ppm)
    n_pos = ppm.shape[1]

    fig, ax = plt.subplots(figsize=(max(5, n_pos * 0.7), 3))

    for pos in range(n_pos):
        # Sort bases by probability at this position
        order = np.argsort(ppm[:, pos])  # ascending
        y_bottom = 0
        for idx in order:
            height = ppm[idx, pos] * ic[pos]
            if height > 0.01:
                base = BASES[idx]
                ax.bar(pos + 1, height, bottom=y_bottom, width=0.8,
                       color=BASE_COLORS[base], edgecolor="none")
                if height > 0.15:
                    ax.text(pos + 1, y_bottom + height / 2, base,
                            ha="center", va="center", fontweight="bold",
                            fontsize=10, color="white")
                y_bottom += height

    ax.set_xlim(0.4, n_pos + 0.6)
    ax.set_ylim(0, 2.1)
    ax.set_xlabel("Position")
    ax.set_ylabel("Information content (bits)")
    ax.set_title(title)
    ax.set_xticks(range(1, n_pos + 1))
    plt.tight_layout()
    plt.show()


plot_sequence_logo(ppm, "Sequence Logo: Example TF Binding Motif")
```

---

## 4. PROSITE Patterns

[PROSITE](https://prosite.expasy.org/) is a database of protein families and domains
described as **regular-expression-like patterns** or **profiles**.

### PROSITE Pattern Syntax

| Syntax | Meaning | Example |
|--------|---------|---------|
| Uppercase letter | Specific amino acid | `C` means cysteine |
| `x` | Any amino acid | `x` matches anything |
| `[ABC]` | One of the listed | `[LIVMF]` = hydrophobic |
| `{P}` | Not proline | `{PC}` = not Pro or Cys |
| `(3)` | Repeat 3 times | `x(3)` = any 3 residues |
| `(2,4)` | Repeat 2-4 times | `x(2,4)` = 2 to 4 residues |
| `-` | Separator | between elements |

**Example:** N-glycosylation site: `N-{P}-[ST]-{P}`  
Meaning: Asn, then any residue except Pro, then Ser or Thr, then any except Pro.

### Converting PROSITE to Python Regex

```python
def prosite_to_regex(pattern):
    """
    Convert a PROSITE pattern to a Python regular expression.

    Handles: specific residues, x (any), [set], {exclusion}, (n), (n,m), separators.
    """
    # Remove terminal dots if present
    pattern = pattern.strip(".")
    elements = pattern.split("-")
    regex_parts = []

    for elem in elements:
        # Check for repeat quantifier
        match_rep = re.match(r'^(.+?)\((\d+)(?:,(\d+))?\)$', elem)

        if match_rep:
            core = match_rep.group(1)
            low = match_rep.group(2)
            high = match_rep.group(3)
        else:
            core = elem
            low = high = None

        # Convert core
        if core == "x":
            r = "."
        elif core.startswith("[") and core.endswith("]"):
            r = core  # already valid regex character class
        elif core.startswith("{") and core.endswith("}"):
            # Exclusion set -> negative character class
            excluded = core[1:-1]
            r = f"[^{excluded}]"
        elif len(core) == 1 and core.isalpha():
            r = core
        else:
            r = core  # fallback

        # Add quantifier
        if low is not None:
            if high is not None:
                r += f"{{{low},{high}}}"
            else:
                r += f"{{{low}}}"

        regex_parts.append(r)

    return "".join(regex_parts)


# Examples
patterns = {
    "N-glycosylation": "N-{P}-[ST]-{P}",
    "Protein kinase C phosphorylation": "[ST]-x-[RK]",
    "Zinc finger C2H2": "C-x(2,4)-C-x(3)-[LIVMFYWC]-x(8)-H-x(3,5)-H",
}

for name, pat in patterns.items():
    regex = prosite_to_regex(pat)
    print(f"{name}:")
    print(f"  PROSITE: {pat}")
    print(f"  Regex:   {regex}\n")
```
