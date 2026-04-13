---
name: bio-applied-ngs-fundamentals
description: "This notebook covers the complete NGS analysis pipeline: sequencing technologies, raw data formats, quality control, read trimming, alignment, and working with SAM/BAM files."
tool_type: python
source_notebook: "Tier_3_Applied_Bioinformatics/01_NGS_Fundamentals/01_ngs_fundamentals.ipynb"
primary_tool: NumPy
---

## Version Compatibility

Reference examples tested with: matplotlib 3.8+, numpy 1.26+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.


# NGS Fundamentals: From Sequencing to Alignment

*Source: Course notebook `Tier_3_Applied_Bioinformatics/01_NGS_Fundamentals/01_ngs_fundamentals.ipynb`*

# NGS Fundamentals: From Sequencing to Alignment

## Tier 3 - Applied Bioinformatics

This notebook covers the complete NGS analysis pipeline: sequencing technologies, raw data formats, quality control, read trimming, alignment, and working with SAM/BAM files.

### Learning Objectives
- Understand how Illumina, PacBio, and Oxford Nanopore sequencing work
- Parse and analyze FASTQ files
- Interpret quality scores and FastQC reports
- Understand read trimming concepts (adapters, quality trimming)
- Work with SAM/BAM format: fields, flags, CIGAR strings
- Use samtools for basic BAM operations
- Calculate genome coverage and depth

---

## 1. Sequencing Technologies

### 1.1 Illumina (Sequencing by Synthesis)

Illumina is the dominant short-read sequencing platform, used in the vast majority of genomics studies.

**How it works:**
1. **Library preparation**: DNA is fragmented, adapters are ligated to both ends
2. **Cluster generation**: Fragments bind to a flow cell, are amplified by bridge PCR into clonal clusters (~1000 copies each)
3. **Sequencing by synthesis**: Fluorescently labeled, reversibly terminated nucleotides are added one at a time. After each incorporation, the flow cell is imaged to identify the base. The terminator is then cleaved, allowing the next cycle.
4. **Base calling**: Images are processed to call bases and assign quality scores

**Key characteristics:**
- Read length: 50-300 bp (paired-end: 2 x 150 bp most common)
- Error rate: ~0.1% (substitution errors dominate)
- Throughput: Up to 6 Tb per run (NovaSeq 6000)
- Cost: Lowest per-base cost

```
Illumina Sequencing by Synthesis:

Flow cell surface
    |  Bridge PCR         Clonal cluster
    |  ~~~>~~~>~~~   ->   |||||||||||||
    |                     |||||||||||||
    
Cycle 1:   A-fluor  ->  Image  ->  Cleave terminator
Cycle 2:   T-fluor  ->  Image  ->  Cleave terminator
Cycle 3:   G-fluor  ->  Image  ->  Cleave terminator
  ...       ...          ...          ...
```

### 1.2 PacBio (Single Molecule Real Time - SMRT)

PacBio produces long reads by observing a single polymerase incorporating nucleotides in real time.

**How it works:**
1. **SMRTbell library**: Circular DNA template with hairpin adapters on both ends
2. **Zero-mode waveguide (ZMW)**: A single polymerase sits at the bottom of a tiny well. As it incorporates fluorescently labeled nucleotides, the fluorescence is detected in real time.
3. **Circular consensus**: The polymerase reads around the circular template multiple times, generating subreads that are combined into a high-accuracy HiFi read.

**Key characteristics:**
- CLR (Continuous Long Read): 10-25 kb, ~10-15% error
- HiFi (Circular Consensus): 10-25 kb, >99.9% accuracy (Q30+)
- Throughput: ~30 Gb per SMRT cell (Revio)
- Errors: Insertions dominate (random, not systematic)

```
PacBio SMRT Sequencing:

  SMRTbell template (circular):
  
  5'---ACGTACGT---3'
  3'---TGCATGCA---5'
      ^         ^
   hairpin   hairpin
   adapter   adapter

  Polymerase reads around and around:
  Pass 1: ------>  (subread 1)
  Pass 2: <------  (subread 2)
  Pass 3: ------>  (subread 3)
  
  Consensus of passes = HiFi read (Q30+)
```

### 1.3 Oxford Nanopore (Nanopore Sequencing)

Oxford Nanopore passes DNA through a protein nanopore and measures changes in ionic current to identify bases.

**How it works:**
1. **Library preparation**: Adapters with motor proteins are ligated to DNA
2. **Nanopore translocation**: A motor protein feeds the DNA strand through a protein nanopore embedded in a membrane. As each base (or k-mer) passes through the pore, it partially blocks the ionic current in a characteristic way.
3. **Base calling**: The current signal is decoded by neural networks (e.g., Guppy, Dorado) into a DNA sequence.

**Key characteristics:**
- Read length: Theoretically unlimited (>1 Mb demonstrated)
- Typical: 10-100 kb, depending on DNA integrity
- Error rate: ~5% raw (R9.4.1), ~1% with R10.4.1 + latest basecallers
- Real-time: Data available as sequencing runs
- Portable: MinION is a USB-powered device

```
Oxford Nanopore Sequencing:

  Membrane
  ==========[PORE]==========
             | |
    Motor -> A T G C A T ...
    protein  | |
             v v
  
  Current signal:
  ___/\___/\_____/\__/\___
    A   T   G  C  A  T
  
  Neural network -> base calls
```

### 1.4 Platform Comparison

| Feature | Illumina | PacBio HiFi | Oxford Nanopore |
|---------|----------|-------------|------------------|
| Read length | 50-300 bp | 10-25 kb | 10 kb - 1 Mb+ |
| Accuracy | ~99.9% | ~99.9% | ~99% (R10.4.1 simplex), ~99.9% (duplex) |
| Error type | Substitutions | Random substitutions (CCS averages out CLR indels) | Homopolymer indels (R9), substitutions (R10.4.1) |
| Throughput | Up to 6 Tb/run | ~30 Gb/cell | 50-200 Gb/cell |
| Time to data | 1-3 days | Hours | Real-time |
| Cost/Gb | Lowest | Medium | Medium |
| Best for | WGS, RNA-seq, ChIP-seq | Genome assembly, SVs | Structural variants, field work |

**Choosing a platform:**
- **High coverage WGS, RNA-seq**: Illumina (cost-effective, high accuracy)
- **De novo assembly, structural variants**: PacBio HiFi (long + accurate)
- **Rapid diagnostics, ultra-long reads**: Oxford Nanopore (portable, real-time)
- **Hybrid approaches**: Combine short and long reads for best results

---

## 2. FASTQ Format

FASTQ is the standard format for storing sequencing reads with per-base quality scores.

### 2.1 Structure

Each read occupies exactly 4 lines:

```
@SEQ_ID                           <- Line 1: Header (starts with @), read identifier
GATCGATCGATCGATCGATC              <- Line 2: DNA sequence
+                                 <- Line 3: Separator (+ optionally followed by ID)
IIIIIIIIIHHHHHGGGFF               <- Line 4: Quality scores (ASCII encoded, same length as seq)
```

### 2.2 Quality Scores (Phred+33)

Quality scores encode the probability that a base call is wrong:

**Phred score** = -10 * log10(P_error)

| Phred Score | Error Probability | Accuracy | ASCII Character |
|-------------|-------------------|----------|------------------|
| 0 | 1.0 (100%) | 0% | ! (33) |
| 10 | 0.1 (10%) | 90% | + (43) |
| 20 | 0.01 (1%) | 99% | 5 (53) |
| 30 | 0.001 (0.1%) | 99.9% | ? (63) |
| 40 | 0.0001 (0.01%) | 99.99% | I (73) |

**Encoding**: ASCII character value - 33 = Phred score (Phred+33 / Sanger encoding, used by all modern platforms)

```python
import math
import numpy as np
from collections import Counter, defaultdict

# Phred score utilities

def phred_to_error_prob(phred):
    """Convert Phred quality score to error probability."""
    return 10 ** (-phred / 10)

def error_prob_to_phred(prob):
    """Convert error probability to Phred quality score."""
    if prob <= 0:
        return 40  # cap at Q40
    return -10 * math.log10(prob)

def ascii_to_phred(char, offset=33):
    """Convert ASCII quality character to Phred score."""
    return ord(char) - offset

def phred_to_ascii(phred, offset=33):
    """Convert Phred score to ASCII quality character."""
    return chr(phred + offset)

# Demonstrate the encoding
print("Phred Score Encoding (Phred+33):")
print(f"{'Phred':>6} {'ASCII':>6} {'Error Prob':>12} {'Accuracy':>10}")
print("-" * 40)
for q in [0, 5, 10, 15, 20, 25, 30, 35, 40]:
    p_err = phred_to_error_prob(q)
    print(f"{q:>6} {phred_to_ascii(q):>6} {p_err:>12.6f} {(1 - p_err)*100:>9.4f}%")

# Decode an example quality string
qual_string = "IIIII?55!!##FFFF"
print(f"\nDecoding quality string: {qual_string}")
scores = [ascii_to_phred(c) for c in qual_string]
print(f"Phred scores: {scores}")
print(f"Mean quality: Q{np.mean(scores):.1f}")
```

### 2.3 Paired-End Reads

Most Illumina experiments generate **paired-end** reads: two reads from opposite ends of the same DNA fragment.

```
DNA fragment:
5'----[====Read 1====>........<====Read 2====]----3'
      |<------------- insert size ------------->|

File naming convention:
  sample_R1.fastq.gz  (forward reads)
  sample_R2.fastq.gz  (reverse reads)
  
Reads are in the same order in both files:
  R1 read 1  <-->  R2 read 1  (same fragment)
  R1 read 2  <-->  R2 read 2  (same fragment)
  ...
```

**Insert size**: Distance between the start of Read 1 and the end of Read 2. Typical: 300-500 bp for standard libraries.

```python
def parse_fastq(filepath, max_reads=None):
    """
    Parse a FASTQ file and yield (header, sequence, quality_string) tuples.
    Handles both plain text and conceptually gzipped files.
    """
    count = 0
    with open(filepath, 'r') as f:
        while True:
            header = f.readline().strip()
            if not header:
                break
            sequence = f.readline().strip()
            f.readline()  # + line
            quality = f.readline().strip()
            yield (header[1:], sequence, quality)  # strip leading @
            count += 1
            if max_reads and count >= max_reads:
                break

# We'll work with simulated data for the examples
# In practice, you would use real FASTQ files

import random
random.seed(42)

def simulate_fastq_reads(n_reads=1000, read_length=150, mean_quality=30):
    """Generate simulated FASTQ reads for demonstration."""
    bases = 'ACGT'
    reads = []
    for i in range(n_reads):
        # Random sequence
        seq = ''.join(random.choice(bases) for _ in range(read_length))
        
        # Quality degrades toward end of read (realistic behavior)
        quals = []
        for pos in range(read_length):
            # Quality decreases with position
            base_q = mean_quality - (pos / read_length) * 10
            q = max(2, min(40, int(random.gauss(base_q, 4))))
            quals.append(q)
        
        qual_str = ''.join(phred_to_ascii(q) for q in quals)
        header = f"SIM:{i+1:06d} length={read_length}"
        reads.append((header, seq, qual_str))
    return reads

sim_reads = simulate_fastq_reads(n_reads=5000)

# Show first 3 reads
print("Simulated FASTQ reads (first 3):\n")
for header, seq, qual in sim_reads[:3]:
    scores = [ascii_to_phred(c) for c in qual]
    print(f"@{header}")
    print(f"{seq[:60]}...")
    print("+")
    print(f"{qual[:60]}...")
    print(f"  Mean Q: {np.mean(scores):.1f}\n")
```

---

## 3. Quality Control

### 3.1 FastQC: The Standard QC Tool

FastQC provides a comprehensive quality assessment of sequencing data. It produces a report with multiple modules:

```bash
# Install
# macOS: brew install fastqc
# Ubuntu: sudo apt install fastqc
# conda: conda install -c bioconda fastqc

# Run FastQC
fastqc sample_R1.fastq.gz sample_R2.fastq.gz -o qc_output/ -t 4

# Output: sample_R1_fastqc.html (visual report)
#         sample_R1_fastqc.zip  (raw data)
```

### 3.2 Key FastQC Metrics

| Module | What to Look For | Common Issues |
|--------|-----------------|---------------|
| Per-base quality | Q28+ across all positions | Quality drops at 3' end |
| Per-sequence quality | Peak at Q30+ | Bimodal = some reads failed |
| Per-base sequence content | Uniform A/T/G/C | Bias in first 10-15 bp (normal for RNA-seq) |
| GC content | Normal distribution matching genome | Shifted = contamination |
| Sequence length | Uniform (before trimming) | Variable after trimming |
| Duplication levels | Low (<20%) | High in PCR-heavy or targeted experiments |
| Overrepresented sequences | None or few | Adapter dimers, rRNA |
| Adapter content | <5% at read ends | >10% needs trimming |

```python
import matplotlib.pyplot as plt

def compute_per_position_quality(reads):
    """Compute per-position quality statistics (like FastQC)."""
    position_quals = defaultdict(list)
    
    for header, seq, qual in reads:
        for pos, qchar in enumerate(qual):
            position_quals[pos].append(ascii_to_phred(qchar))
    
    positions = sorted(position_quals.keys())
    stats = {
        'positions': positions,
        'median': [np.median(position_quals[p]) for p in positions],
        'q25': [np.percentile(position_quals[p], 25) for p in positions],
        'q75': [np.percentile(position_quals[p], 75) for p in positions],
        'p10': [np.percentile(position_quals[p], 10) for p in positions],
        'p90': [np.percentile(position_quals[p], 90) for p in positions],
    }
    return stats

# Compute and plot per-position quality
qual_stats = compute_per_position_quality(sim_reads)

fig, ax = plt.subplots(figsize=(14, 5))

pos = qual_stats['positions']
ax.fill_between(pos, qual_stats['p10'], qual_stats['p90'],
                alpha=0.15, color='blue', label='10th-90th percentile')
ax.fill_between(pos, qual_stats['q25'], qual_stats['q75'],
                alpha=0.3, color='blue', label='IQR (25th-75th)')
ax.plot(pos, qual_stats['median'], 'b-', linewidth=2, label='Median')

# Quality zones
ax.axhspan(30, 42, alpha=0.08, color='green')
ax.axhspan(20, 30, alpha=0.08, color='yellow')
ax.axhspan(0, 20, alpha=0.08, color='red')
ax.axhline(y=30, color='green', linestyle='--', alpha=0.5, label='Q30')
ax.axhline(y=20, color='orange', linestyle='--', alpha=0.5, label='Q20')

ax.set_xlabel('Position in Read (bp)')
ax.set_ylabel('Phred Quality Score')
ax.set_title('Per-Position Quality Scores (FastQC-style)')
ax.set_ylim(0, 42)
ax.legend(loc='lower left')
plt.tight_layout()
plt.show()
```
