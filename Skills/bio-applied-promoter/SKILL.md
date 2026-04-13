---
name: bio-applied-promoter
description: "**Tier 3 — Module 05** | [Next: Regulatory Analysis →](./02_regulatory_analysis.ipynb)"
tool_type: python
source_notebook: "Tier_3_Applied_Bioinformatics/05_Promoter_and_Regulatory_Analysis/01_promoter.ipynb"
---

# Promoter and Regulatory Sequence Analysis — Part 1: Core Promoter Elements

*Source: Course notebook `Tier_3_Applied_Bioinformatics/05_Promoter_and_Regulatory_Analysis/01_promoter.ipynb`*

# Promoter and Regulatory Sequence Analysis — Part 1: Core Promoter Elements

**Tier 3 — Module 05** | [Next: Regulatory Analysis →](./02_regulatory_analysis.ipynb)

# Promoter and Regulatory Sequence Analysis

**Tier 3 -- Applied Bioinformatics**

Gene expression is orchestrated by regulatory elements embedded in DNA. This notebook covers how to identify and analyze these elements computationally: promoter regions, CpG islands, transcription factor binding sites (TFBS), and methylation patterns.

**Prerequisites:** Tier 2 (BioPython basics, sequence handling, basic statistics)  
**Libraries:** `numpy`, `pandas`, `matplotlib`, `re`

```python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import re
from collections import Counter

%matplotlib inline
plt.rcParams['figure.figsize'] = (12, 5)
plt.rcParams['font.size'] = 12
```

---
## 1. Gene Regulation: The Big Picture

Every cell in an organism contains the same DNA, yet a neuron and a liver cell behave very differently. The key lies in **gene regulation** -- controlling which genes are transcribed, when, and how much.

### Key regulatory elements

| Element | Location | Function |
|---------|----------|----------|
| **Promoter** | Immediately upstream of gene (typically -1 to -1000 bp from TSS) | Recruits RNA polymerase, initiates transcription |
| **Enhancer** | Can be thousands of bp away (upstream, downstream, or within introns) | Boosts transcription when bound by activators |
| **Silencer** | Variable distance from gene | Represses transcription |
| **Insulator** | Between regulatory elements | Blocks enhancer-promoter interactions |

### Transcription factors (TFs)

Transcription factors are proteins that bind specific short DNA motifs (6-20 bp) to regulate gene expression. They can be:
- **Activators** -- enhance transcription by recruiting RNA polymerase
- **Repressors** -- block transcription

The **transcription start site (TSS)** is position +1, where RNA polymerase begins transcribing. Positions upstream are negative (e.g., -30 for the TATA box), and downstream positions are positive.

---
## 2. Promoter Elements

### 2.1 The TATA box

The TATA box is a core promoter element found approximately 25-35 bp upstream of the TSS. Its consensus sequence is **TATAAA**, though variations exist (TATAWAW where W = A or T).

- Present in roughly 10-20% of human genes (more common in highly regulated, tissue-specific genes)
- Bound by the TATA-binding protein (TBP), part of the TFIID complex
- TATA-less promoters often rely on other elements like the Initiator (Inr) or DPE

```python
def find_tata_boxes(sequence, strict=True):
    """Find TATA box motifs in a DNA sequence.
    
    strict=True: exact TATAAA
    strict=False: TATA[AT]A[AT] (IUPAC TATAWAW)
    """
    sequence = sequence.upper()
    pattern = 'TATAAA' if strict else 'TATA[AT]A[AT]'
    matches = [(m.start(), m.group()) for m in re.finditer(pattern, sequence)]
    return matches

# Example: a synthetic promoter region (1000 bp upstream of TSS at position 1000)
np.random.seed(42)
bases = 'ACGT'
# Generate random sequence, then insert a TATA box near expected position
random_seq = ''.join(np.random.choice(list(bases), size=2000))
# Insert TATA box at position 970 (i.e., -30 relative to TSS at 1000)
promoter_seq = random_seq[:970] + 'TATAAAG' + random_seq[977:]

hits = find_tata_boxes(promoter_seq, strict=True)
print(f"TATA box hits in 2000 bp region (TSS at position 1000):")
for pos, motif in hits:
    relative = pos - 1000
    print(f"  Position {pos} (TSS{relative:+d}): {motif}")
```

### 2.2 CpG islands

**CpG dinucleotides** (cytosine followed by guanine, written CpG to distinguish from C-G base pairing) are statistically underrepresented in vertebrate genomes because methylated cytosines in CpG context tend to mutate to thymine over evolutionary time.

However, clusters of CpGs called **CpG islands** are preserved near promoters of ~70% of human genes. They are defined by:

1. **Length** >= 200 bp
2. **GC content** >= 50%
3. **Observed/Expected CpG ratio** >= 0.6

The observed/expected ratio is:

$$\text{CpG O/E} = \frac{N_{CpG} \times L}{N_C \times N_G}$$

where $N_{CpG}$ is the count of CpG dinucleotides, $L$ is the window length, and $N_C$, $N_G$ are counts of C and G nucleotides.

```python
def cpg_island_scanner(sequence, window=200, step=1, gc_thresh=0.5, oe_thresh=0.6):
    """Scan a sequence for CpG islands using a sliding window.
    
    Returns a list of (start, end, gc_content, cpg_oe_ratio) for qualifying windows.
    """
    sequence = sequence.upper()
    islands = []
    
    for i in range(0, len(sequence) - window + 1, step):
        win = sequence[i:i + window]
        n_c = win.count('C')
        n_g = win.count('G')
        n_cpg = win.count('CG')  # CpG dinucleotide count
        gc_content = (n_c + n_g) / window
        
        # Avoid division by zero
        if n_c > 0 and n_g > 0:
            cpg_oe = (n_cpg * window) / (n_c * n_g)
        else:
            cpg_oe = 0.0
        
        if gc_content >= gc_thresh and cpg_oe >= oe_thresh:
            islands.append((i, i + window, gc_content, cpg_oe))
    
    return islands

print("CpG island scanner defined. We will use it in the practical section below.")
```

### 2.3 Transcription factor binding sites (TFBS)

TFs recognize short, degenerate DNA motifs. These are represented as **position weight matrices (PWMs)** or **position frequency matrices (PFMs)**.

A PWM stores the log-likelihood of each base at each position of the motif relative to background:

$$w_{b,i} = \log_2\frac{p(b, i)}{p_{bg}(b)}$$

where $p(b, i)$ is the probability of base $b$ at position $i$ and $p_{bg}(b)$ is the background probability.

### Key databases

| Database | Description | Access |
|----------|-------------|--------|
| **JASPAR** | Open-access, curated TF binding profiles | https://jaspar.elixir.no |
| **TRANSFAC** | Comprehensive commercial/academic TF database | https://genexplain.com/transfac |
| **HOCOMOCO** | Human and mouse TF motifs from ChIP-seq | https://hocomoco11.autosome.org |

In the source data for this course, TFBS were identified using TRANSFAC matrices aligned to rice promoter regions. Each position in a gene's promoter was annotated with the number of TFBS overlapping it.

```python
# Build a simple PWM for a TATA box motif from example binding sites
tata_sites = [
    'TATAAAG',
    'TATAAAT',
    'TATAAAA',
    'TATATAG',
    'TATAAAC',
    'TATAAAT',
    'TATAAAG',
    'TATAAAT',
    'TATAAAA',
    'TATAAAG',
]

def build_pfm(sites):
    """Build a position frequency matrix from aligned binding sites."""
    length = len(sites[0])
    pfm = {b: [0] * length for b in 'ACGT'}
    for site in sites:
        for i, base in enumerate(site.upper()):
            pfm[base][i] += 1
    return pfm

def pfm_to_pwm(pfm, pseudocount=0.5, bg=None):
    """Convert PFM to PWM (log-odds scores)."""
    if bg is None:
        bg = {'A': 0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25}
    n = sum(pfm[b][0] for b in 'ACGT')  # total sequences
    length = len(pfm['A'])
    pwm = {b: [0.0] * length for b in 'ACGT'}
    for b in 'ACGT':
        for i in range(length):
            freq = (pfm[b][i] + pseudocount) / (n + 4 * pseudocount)
            pwm[b][i] = np.log2(freq / bg[b])
    return pwm

pfm = build_pfm(tata_sites)
pwm = pfm_to_pwm(pfm)

print("Position Frequency Matrix (PFM):")
pfm_df = pd.DataFrame(pfm, index=[f"pos {i+1}" for i in range(len(pfm['A']))])
print(pfm_df)
print("\nPosition Weight Matrix (PWM, log2 odds):")
pwm_df = pd.DataFrame(pwm, index=[f"pos {i+1}" for i in range(len(pwm['A']))])
print(pwm_df.round(3))
```

```python
def score_sequence(sequence, pwm):
    """Score a sequence against a PWM. Returns the sum of log-odds."""
    score = 0.0
    for i, base in enumerate(sequence.upper()):
        if base in pwm:
            score += pwm[base][i]
    return score

def scan_with_pwm(sequence, pwm, threshold=0.0):
    """Scan a sequence with a PWM and return hits above threshold."""
    motif_len = len(pwm['A'])
    hits = []
    for i in range(len(sequence) - motif_len + 1):
        subseq = sequence[i:i + motif_len]
        sc = score_sequence(subseq, pwm)
        if sc >= threshold:
            hits.append((i, subseq, sc))
    return hits

# Scan our synthetic promoter with the TATA PWM
hits = scan_with_pwm(promoter_seq, pwm, threshold=5.0)
print(f"PWM hits (score >= 5.0) in synthetic promoter:")
for pos, seq, sc in sorted(hits, key=lambda x: -x[2])[:10]:
    relative = pos - 1000
    print(f"  Position {pos} (TSS{relative:+d}): {seq}  score={sc:.2f}")
```

---
## 3. TSS Prediction Concepts

Identifying the **transcription start site (TSS)** is fundamental for promoter analysis. Computational approaches include:

1. **Signal-based methods**: Look for known motifs (TATA box at -25 to -30, Inr at +1, DPE at +30)
2. **Content-based methods**: CpG island presence, GC content profiles
3. **Experimental data integration**: CAGE (Cap Analysis of Gene Expression), which captures 5' ends of mRNAs
4. **Machine learning**: Combine multiple features (motifs, chromatin marks, conservation) to predict TSS

In the rice promoter dataset used in this course, TSS positions were known, and each gene had 1000 bp upstream and 1000 bp downstream of the TSS. This allows us to study the distribution of regulatory signals relative to TSS.

```python
# Simulate signal distribution around TSS (based on patterns from the rice promoter study)
positions = np.arange(-1000, 1000)

# TATA box: peak around -30
tata_signal = np.exp(-0.5 * ((positions + 30) / 5) ** 2) * 50
tata_signal += np.random.poisson(1, len(positions))

# CpG density: enriched around TSS
cpg_signal = np.exp(-0.5 * (positions / 200) ** 2) * 30
cpg_signal += np.random.poisson(2, len(positions))

# TFBS density: enriched upstream
tfbs_signal = np.where(positions < 0, 15 + 10 * np.exp(-0.5 * ((positions + 200) / 150) ** 2), 5)
tfbs_signal += np.random.poisson(2, len(positions))

fig, axes = plt.subplots(3, 1, figsize=(14, 10), sharex=True)

axes[0].fill_between(positions, tata_signal, alpha=0.5, color='orange')
axes[0].axvline(0, color='red', linestyle='--', label='TSS')
axes[0].set_ylabel('TATA box signal')
axes[0].set_title('Regulatory Signal Distribution Around the TSS')
axes[0].legend()

axes[1].fill_between(positions, cpg_signal, alpha=0.5, color='blue')
axes[1].axvline(0, color='red', linestyle='--', label='TSS')
axes[1].set_ylabel('CpG density')
axes[1].legend()

axes[2].fill_between(positions, tfbs_signal, alpha=0.5, color='green')
axes[2].axvline(0, color='red', linestyle='--', label='TSS')
axes[2].set_ylabel('TFBS density')
axes[2].set_xlabel('Position relative to TSS (bp)')
axes[2].legend()

plt.tight_layout()
plt.show()
```

---
## 4. CpG Island Detection

Let us build a complete CpG island detection pipeline. We will:
1. Generate a realistic promoter sequence with an embedded CpG island
2. Scan it with our sliding window algorithm
3. Merge overlapping windows into island calls
4. Visualize the results

```python
def generate_promoter_with_cpg_island(length=3000, island_start=1200, island_length=400):
    """Generate a synthetic promoter with an embedded CpG island."""
    np.random.seed(123)
    
    # Background: AT-rich, CpG-depleted (typical vertebrate genome)
    bg_weights = [0.29, 0.21, 0.21, 0.29]  # A, C, G, T
    bg_seq = list(np.random.choice(list('ACGT'), size=length, p=bg_weights))
    
    # CpG island: GC-rich, CpG-enriched
    island_weights = [0.18, 0.32, 0.32, 0.18]
    island_seq = list(np.random.choice(list('ACGT'), size=island_length, p=island_weights))
    
    # Further enrich CpGs in the island
    for i in range(0, island_length - 1, 3):
        if np.random.random() < 0.3:
            island_seq[i] = 'C'
            island_seq[i + 1] = 'G'
    
    # Insert island into background
    bg_seq[island_start:island_start + island_length] = island_seq
    return ''.join(bg_seq)

synth_promoter = generate_promoter_with_cpg_island()
print(f"Synthetic promoter length: {len(synth_promoter)} bp")
print(f"Expected CpG island: positions 1200-1600")
```
