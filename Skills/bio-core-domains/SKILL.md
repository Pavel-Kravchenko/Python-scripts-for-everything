---
name: bio-core-domains
description: "- Explain what sequence motifs are and why they matter biologically - Construct and use Position Weight Matrices (PWMs) to score sequences - Understand information content and sequence logos - Work wi"
tool_type: python
source_notebook: "Tier_2_Core_Bioinformatics/10_Sequence_Motifs_and_Domains/02_domains.ipynb"
primary_tool: NumPy
---

## Version Compatibility

Reference examples tested with: matplotlib 3.8+, numpy 1.26+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.


# Sequence Motifs and Protein Domains

*Source: Course notebook `Tier_2_Core_Bioinformatics/10_Sequence_Motifs_and_Domains/02_domains.ipynb`*


---

## Learning Objectives

By the end of this notebook you will be able to:

- Explain what sequence motifs are and why they matter biologically
- Construct and use Position Weight Matrices (PWMs) to score sequences
- Understand information content and sequence logos
- Work with PROSITE patterns
- Describe protein domains and the major domain databases (Pfam, InterPro, CDD)
- Search for domains using hmmscan and InterProScan
- Interpret domain architectures of multi-domain proteins
- Apply these concepts to find transcription factor binding sites and annotate protein domains

## How to use this notebook

1. Run cells top-to-bottom in order — later cells depend on earlier ones
2. Run the environment check cell first to confirm all imports work
3. Read the explanatory text before each code cell — the context matters
4. The exercises at the end are designed to be done after reading each section

## Complicated moments explained

- **Domain vs. motif**: A domain is a structurally and functionally independent unit, typically 50-300 residues, that can fold independently and often recurs in diverse protein contexts. A motif is a shorter pattern (typically <20 residues) that is functionally important but does not necessarily fold independently.
- **Profile HMM scoring**: Pfam represents each domain family as a profile HMM with match states, insert states, and delete states. Scores are in bits relative to a null model; the gathering threshold (GA) is family-specific — do not apply a universal bit-score cutoff.
- **E-value vs. score for domain search**: hmmscan reports two E-values: sequence E-value (is there any domain match in the whole sequence?) and domain E-value (how significant is this specific domain hit?). For annotating individual domains, use the domain E-value.
- **Domain boundaries are approximate**: Pfam domain boundaries represent the HMM consensus. Always check the actual alignment to determine the true structural domain extents for your protein.
- **InterPro integrates multiple databases**: An InterPro entry (IPR000xxx) integrates matching entries from Pfam, SMART, ProSite, PRINTS, CDD, and SUPERFAMILY. InterProScan is the best first-pass tool for comprehensive domain annotation.

## Environment check (run this first)

```python
# Environment check
import numpy as np
import matplotlib.pyplot as plt
import shutil

print("Imports ready.")

hmmer_found = shutil.which('hmmscan')
if hmmer_found:
    print(f"HMMER hmmscan found at: {hmmer_found}")
else:
    print("HMMER not found (install: conda install -c bioconda hmmer)")
    print("The notebook uses simulated HMMER output for demonstrations.")

print("\nDomain databases:")
databases = [
    ("Pfam",     "Profile HMMs, ~20,000 families",      "Best overall coverage"),
    ("InterPro", "Integrates Pfam, SMART, CDD, etc.",   "Most comprehensive"),
    ("SMART",    "Focus on signaling domains",           "Excellent for kinases, receptors"),
    ("CDD",      "NCBI-curated, integrates Pfam/SMART", "Free with local BLAST"),
    ("PROSITE",  "Patterns and profiles",                "Motifs and active sites"),
]
for db, content, notes in databases:
    print(f"  {db:<12} {content:<42} {notes}")
```python

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
```python

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
```python

```python
# Scan a longer sequence
target = "GCCTAGATGACTCAGGTTTCCCGTGACTCAATGCAATGACTTACCC"

hits = scan_sequence(pwm, target, threshold=0.5 * np.sum(np.max(pwm, axis=0)))

print(f"Scanning: {target}")
print(f"\n{'Pos':>4} {'Subseq':<10} {'Score':>7}")
print("-" * 25)
for pos, subseq, score in hits:
    print(f"{pos:4d} {subseq:<10} {score:7.2f}")
```python

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
```python

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
```python

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
```python

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
```python

## Common Pitfalls

- **Coordinate systems**: BED uses 0-based half-open; VCF/GFF use 1-based inclusive — mixing them causes off-by-one errors
- **Batch effects**: Always check for batch confounding before interpreting biological signal
- **Multiple testing**: Apply FDR correction (Benjamini-Hochberg) when testing thousands of features simultaneously
