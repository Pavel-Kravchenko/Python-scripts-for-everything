---
name: bio-applied-ngs-fundamentals
description: "NGS platform comparison, FASTQ format, Phred quality scores, and QC metrics. Reference for sequencing technology selection and read quality assessment."
tool_type: python
primary_tool: NumPy
---

# NGS Fundamentals

## Platform Comparison

| Feature | Illumina | PacBio HiFi | Oxford Nanopore |
|---------|----------|-------------|-----------------|
| Read length | 50–300 bp | 10–25 kb | 10 kb – 1 Mb+ |
| Accuracy | ~99.9% | ~99.9% | ~99% (R10.4.1 simplex) |
| Error type | Substitutions | Random (CCS averages indels) | Homopolymer indels (R9), substitutions (R10) |
| Throughput | Up to 6 Tb/run | ~30 Gb/cell | 50–200 Gb/cell |
| Best for | WGS, RNA-seq, ChIP-seq | De novo assembly, SVs | Structural variants, field diagnostics |

**Platform selection:**
- High-coverage WGS / RNA-seq: Illumina (cost, accuracy)
- De novo assembly / complex SVs: PacBio HiFi (long + accurate)
- Rapid diagnostics / ultra-long reads: Oxford Nanopore (portable, real-time)
- Best assemblies: hybrid Illumina + long-read

## FASTQ Format

Each read = 4 lines: `@header`, sequence, `+`, quality string (same length as sequence).

### Phred+33 encoding (all modern platforms)

| Phred | Error prob | Accuracy | ASCII char |
|-------|-----------|----------|------------|
| 10 | 10% | 90% | + |
| 20 | 1% | 99% | 5 |
| 30 | 0.1% | 99.9% | ? |
| 40 | 0.01% | 99.99% | I |

```python
def phred_to_error_prob(q): return 10 ** (-q / 10)
def ascii_to_phred(char, offset=33): return ord(char) - offset
def phred_to_ascii(q, offset=33): return chr(q + offset)

# Decode quality string
qual_string = "IIIII?55!!"
scores = [ascii_to_phred(c) for c in qual_string]
mean_q = sum(scores) / len(scores)
```

### Parsing FASTQ

```python
def parse_fastq(filepath):
    with open(filepath, 'r') as f:
        while True:
            header = f.readline().strip()
            if not header:
                break
            seq = f.readline().strip()
            f.readline()  # + line
            qual = f.readline().strip()
            yield header[1:], seq, qual  # strip leading @
```

For gzipped files: `import gzip; open = gzip.open(filepath, 'rt')`.

## FastQC Key Metrics

| Module | Pass criteria | Common failure cause |
|--------|--------------|---------------------|
| Per-base quality | Q28+ across all positions | Quality drop at 3' end (normal) |
| Per-sequence quality | Peak at Q30+ | Bimodal = subset of failed reads |
| GC content | Normal distribution | Shifted = contamination |
| Duplication level | <20% | High in targeted / PCR-heavy libs |
| Adapter content | <5% at ends | >10% → trim with Trimmomatic/fastp |

```bash
fastqc sample_R1.fastq.gz sample_R2.fastq.gz -o qc_output/ -t 4
```

```python
from collections import defaultdict
import numpy as np

def per_position_quality(reads):
    """Compute per-position quality stats (FastQC-style)."""
    pos_quals = defaultdict(list)
    for _, _, qual in reads:
        for pos, qchar in enumerate(qual):
            pos_quals[pos].append(ascii_to_phred(qchar))
    positions = sorted(pos_quals)
    return {
        'positions': positions,
        'median': [np.median(pos_quals[p]) for p in positions],
        'q25':    [np.percentile(pos_quals[p], 25) for p in positions],
        'q75':    [np.percentile(pos_quals[p], 75) for p in positions],
    }
```

## Pitfalls

- **Coordinate systems**: BED = 0-based half-open; VCF/GFF = 1-based inclusive — mixing causes off-by-one errors
- **Phred+33 vs Phred+64**: Old Illumina pipeline (CASAVA <1.8) used +64 offset. Check `fastqc` encoding warning. `ord('@') - 64 = 0`, `ord('!') - 33 = 0`.
- **Quality drop at 3' end**: Normal behavior for SBS. Trim low-quality tails before alignment (fastp `--cut_tail` or Trimmomatic `TRAILING:20`).
- **Paired-end read order**: R1 and R2 must be in the same order. If a read is filtered from R1, remove the same read from R2.
- **Batch effects**: Always check for batch confounding before interpreting biological signal
- **Multiple testing**: Apply FDR correction (Benjamini-Hochberg) when testing thousands of features
