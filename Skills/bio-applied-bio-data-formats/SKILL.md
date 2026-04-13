---
name: bio-applied-bio-data-formats
description: "Modern bioinformatics pipelines exchange data through a rich ecosystem of plain-text and binary file formats. This notebook surveys **every major format** you will encounter in NGS, structural biology"
tool_type: python
source_notebook: "Tier_3_Applied_Bioinformatics/01_NGS_Fundamentals/02_bio_data_formats.ipynb"
primary_tool: Python
---

## Version Compatibility

Reference examples tested with: biopython 1.83+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.


# Bioinformatics Data Formats: A Comprehensive Guide

*Source: Course notebook `Tier_3_Applied_Bioinformatics/01_NGS_Fundamentals/02_bio_data_formats.ipynb`*


## Tier 3 – Applied Bioinformatics

Modern bioinformatics pipelines exchange data through a rich ecosystem of
plain-text and binary file formats. This notebook surveys **every major format**
you will encounter in NGS, structural biology, and phylogenetics work, with
hands-on Python parsing examples for each one.

---

### Contents

1. [FASTA – Sequences](#1.-FASTA)
2. [FASTQ – Sequencing Reads](#2.-FASTQ)
3. [SAM / BAM / CRAM – Alignments](#3.-SAM-/-BAM-/-CRAM)
4. [VCF / BCF – Variants](#4.-VCF-/-BCF)
5. [BED – Genomic Intervals](#5.-BED)
6. [GFF / GTF – Gene Annotations](#6.-GFF-/-GTF)
7. [WIG / BedGraph / BigWig – Coverage Tracks](#7.-WIG-/-BedGraph-/-BigWig)
8. [PDB – Protein Structure](#8.-PDB)
9. [mmCIF / PDBx – Modern Structure Format](#9.-mmCIF-/-PDBx)
10. [FAST5 / POD5 – Oxford Nanopore Raw Signal](#10.-FAST5-/-POD5)
11. [Newick / Nexus / NHX – Phylogenetic Trees](#11.-Newick-/-Nexus-/-NHX)
12. [Format Quick-Reference Table](#12.-Quick-Reference)

```python
# Standard-library and lightweight helpers used throughout this notebook
import io
import re
import gzip
import struct
import textwrap
from pathlib import Path
from collections import defaultdict

# BioPython covers most format I/O we need
try:
    from Bio import SeqIO, AlignIO
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq
    BIOPYTHON = True
except ImportError:
    BIOPYTHON = False
    print("BioPython not installed – pure-Python fallbacks will be used.")

print("Ready.")
```python

---

## 1. FASTA

### 1.1 Format specification

FASTA is the simplest sequence format and the *lingua franca* of bioinformatics.

```python
>sequence_id [optional description]
AGCTAGCTAGCTAGCTAGCT
AGCTAGCTAGCTAGCTAGCT   ← wrap at 60–80 chars (conventional)
```python

Key rules
- The `>` character **must** be the first character on the header line.
- Sequence lines may span multiple lines; there is no canonical line length.
- Blank lines between records are **allowed but discouraged**.
- Sequence characters are case-insensitive; lowercase is used for soft-masked
  (repeat-masked) regions.

### 1.2 Common variants

| Variant | Extension | Notes |
|---------|-----------|-------|
| Nucleotide sequences | `.fa`, `.fna`, `.fasta` | DNA/RNA |
| Protein sequences | `.faa`, `.fasta` | Amino-acid one-letter codes |
| Multi-FASTA | `.fa`, `.fasta` | Multiple records in one file |
| FASTA quality (QUAL) | `.qual` | Phred scores in FASTA-style (legacy) |
| Indexed FASTA | `.fa` + `.fai` | `samtools faidx` index for random access |

```python
# ── Pure-Python FASTA parser ──────────────────────────────────────────────
def parse_fasta(text_or_path):
    """Yield (header, sequence) tuples from a FASTA string or file path."""
    if isinstance(text_or_path, (str, Path)) and Path(text_or_path).exists():
        opener = gzip.open if str(text_or_path).endswith('.gz') else open
        with opener(text_or_path, 'rt') as fh:
            yield from _fasta_iter(fh)
    else:
        yield from _fasta_iter(io.StringIO(str(text_or_path)))

def _fasta_iter(fh):
    header, seqs = None, []
    for line in fh:
        line = line.rstrip('\n')
        if line.startswith('>'):
            if header is not None:
                yield header, ''.join(seqs)
            header, seqs = line[1:], []
        elif header is not None:
            seqs.append(line)
    if header is not None:
        yield header, ''.join(seqs)

# ── Example FASTA data ────────────────────────────────────────────────────
sample_fasta = """\
>sp|P69905|HBA_HUMAN Hemoglobin subunit alpha OS=Homo sapiens
MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHG
KKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTP
AVHASLDKFLASVSTVLTSKYR
>sp|P68871|HBB_HUMAN Hemoglobin subunit beta OS=Homo sapiens
MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPK
VKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFG
KEFTPPVQAAYQKVVAGVANALAHKYH
"""

for header, seq in parse_fasta(sample_fasta):
    seq_id = header.split()[0]
    print(f"ID   : {seq_id}")
    print(f"Desc : {' '.join(header.split()[1:])}")
    print(f"Len  : {len(seq)} aa")
    print()
```python

```python
# ── Writing FASTA ─────────────────────────────────────────────────────────
def write_fasta(records, wrap=60):
    """Return a FASTA-formatted string from an iterable of (header, seq) pairs."""
    lines = []
    for header, seq in records:
        lines.append(f'>{header}')
        for i in range(0, len(seq), wrap):
            lines.append(seq[i:i+wrap])
    return '\n'.join(lines) + '\n'

output = write_fasta([
    ("my_gene|123 A synthetic demo sequence", "ATGCGTACGATCGATCGATCGATCGATCGATCGATCGATCGATCG"),
    ("another_gene|456 Second record",         "ATCGATCGATCGATCGATCG"),
])
print(output)
```python

```python
# ── FASTA index (.fai) ────────────────────────────────────────────────────
# samtools faidx creates a 5-column tab-separated index:
# NAME  LENGTH  OFFSET  LINEBASES  LINEWIDTH

fai_example = """\
chr1\t248956422\t52\t60\t61
chr2\t242193529\t253404903\t60\t61
chrX\t156040895\t2699520056\t60\t61
"""

print("FAI column meanings:")
cols = ["NAME", "LENGTH", "OFFSET (bytes)", "BASES PER LINE", "BYTES PER LINE"]
for col, val in zip(cols, fai_example.strip().split('\n')[0].split('\t')):
    print(f"  {col:<20} = {val}")

# Random-access using the index
def faidx_fetch(fasta_path, fai_path, chrom, start, end):
    """
    Fetch a substring from an indexed FASTA without loading the whole file.
    Coordinates are 0-based half-open [start, end).
    """
    # Load index
    index = {}
    with open(fai_path) as fh:
        for line in fh:
            name, length, offset, bases_per_line, bytes_per_line = line.split()
            index[name] = (int(length), int(offset),
                           int(bases_per_line), int(bytes_per_line))
    length, offset, bases_per_line, bytes_per_line = index[chrom]
    # Byte position of first base of the query region
    n_full_lines, remainder = divmod(start, bases_per_line)
    byte_start = offset + n_full_lines * bytes_per_line + remainder
    # Upper bound (over-estimate; we'll trim)
    bases_needed = end - start
    bytes_needed = (bases_needed // bases_per_line + 2) * bytes_per_line
    with open(fasta_path, 'rb') as fh:
        fh.seek(byte_start)
        raw = fh.read(bytes_needed).decode()
    seq = re.sub(r'\s', '', raw)[:bases_needed]
    return seq

print("\nfaidx_fetch is defined – call it on a real .fa/.fai pair to fetch regions.")
```python

---

## 2. FASTQ

### 2.1 Format specification

FASTQ extends FASTA by embedding per-base quality scores.

```python
@SEQ_ID [optional description]    ← header line (@ not >)
AGCTAGCTAGCTAGCT                  ← raw sequence (may NOT wrap)
+                                 ← separator (optionally repeats the header)
IIIIIIIIIIIIIIII                  ← quality string, same length as sequence
```python

Each ASCII character in the quality string encodes a **Phred score**:

```python
Phred(Q) = -10 × log₁₀(P_error)
ASCII offset 33 (Sanger/Illumina ≥1.8): quality_char = chr(Q + 33)
```python

### 2.2 Quality encoding history

| Encoding | Offset | Range | Used by |
|----------|--------|-------|---------|
| Sanger / Illumina ≥1.8 | 33 | Q0–Q40 | Current standard |
| Illumina 1.3–1.7 | 64 | Q0–40 | Legacy |
| Solexa | 64 | Q-5 to Q40 | Very old data |

> 💡 **Tip**: `file.fastq.gz` gzip-compressed FASTQ is the de-facto standard
> for storing raw reads; always gzip before storing.

```python
# ── FASTQ parser and quality stats ───────────────────────────────────────
import math

def parse_fastq(source, max_reads=None):
    """Yield (header, sequence, quality_string) from a FASTQ source."""
    lines = []
    if isinstance(source, str) and '\n' in source:
        fh = io.StringIO(source)
    else:
        opener = gzip.open if str(source).endswith('.gz') else open
        fh = opener(source, 'rt')
    try:
        count = 0
        it = iter(fh)
        for line in it:
            header = line.rstrip('\n')[1:]   # drop '@'
            seq    = next(it).rstrip('\n')
            next(it)                          # skip '+'
            qual   = next(it).rstrip('\n')
            yield header, seq, qual
            count += 1
            if max_reads and count >= max_reads:
                break
    finally:
        if not isinstance(source, str):
            fh.close()

def phred(char, offset=33):
    return ord(char) - offset

def mean_quality(qual_str):
    scores = [phred(c) for c in qual_str]
    return sum(scores) / len(scores)

# ── Sample FASTQ ──────────────────────────────────────────────────────────
sample_fastq = """\
@read1 length=50
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIA
@read2 length=50
TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
+
!"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQR
"""

for hdr, seq, qual in parse_fastq(sample_fastq):
    q_vals = [phred(c) for c in qual]
    print(f"Read : {hdr}")
    print(f"  Seq len : {len(seq)}")
    print(f"  Mean Q  : {mean_quality(qual):.1f}")
    print(f"  Min Q   : {min(q_vals)}  Max Q: {max(q_vals)}")
    print(f"  Q≥30 %  : {sum(q>=30 for q in q_vals)/len(q_vals)*100:.1f}%")
    print()
```python

```python
# ── FASTQ quality histogram (text-based) ─────────────────────────────────
from collections import Counter

qual_string = 'IIIIIIIIIII888888888!!!!!!!!!!!!!!!IIIIIIIIIIIIIIII'
counts = Counter(phred(c) for c in qual_string)

print(f"Quality score distribution (n={len(qual_string)} bases):")
print(f"{'Q':>4}  {'count':>5}  bar")
for q in sorted(counts):
    bar = '█' * counts[q]
    print(f"{q:>4}  {counts[q]:>5}  {bar}")
```python

---

## 3. SAM / BAM / CRAM

### 3.1 SAM format

**SAM (Sequence Alignment/Map)** is the universal format for read alignments.

```python
@HD  VN:1.6  SO:coordinate          ← header section (lines starting with @)
@SQ  SN:chr1  LN:248956422
@PG  ID:bwa  PN:bwa  VN:0.7.17
read1  0  chr1  100  60  50M  *  0  0  ACGT...  IIII...  NM:i:0
```python

#### SAM FLAG field (column 2)

The FLAG is a bitwise combination:

| Bit  | Hex  | Meaning |
|------|------|---------|
| 1    | 0x1  | Read is paired |
| 2    | 0x2  | Pair properly mapped |
| 4    | 0x4  | Read unmapped |
| 8    | 0x8  | Mate unmapped |
| 16   | 0x10 | Read on reverse strand |
| 32   | 0x20 | Mate on reverse strand |
| 64   | 0x40 | Read 1 of pair |
| 128  | 0x80 | Read 2 of pair |
| 256  | 0x100 | Secondary alignment |
| 512  | 0x200 | Not passing filters |
| 1024 | 0x400 | PCR / optical duplicate |
| 2048 | 0x800 | Supplementary alignment |

#### CIGAR string

Describes how the read aligns to the reference:

| Op | Meaning |
|----|---------|
| M  | Match / mismatch (consumes both) |
| I  | Insertion to reference |
| D  | Deletion from reference |
| N  | Skipped region (intron in RNA-seq) |
| S  | Soft clip (bases present in read) |
| H  | Hard clip (bases removed) |
| =  | Sequence match |
| X  | Sequence mismatch |

### 3.2 BAM – Binary SAM

BAM is the **BGZF-compressed binary encoding** of SAM. It supports O(1) random
access via a `.bai` index created by `samtools index`.

### 3.3 CRAM – Compact Random Access Map

CRAM achieves ~50–60 % compression over BAM by storing differences from a
reference genome rather than raw bases. Requires access to the reference
at decoding time.

| Format | Typical size | Random access | Tools |
|--------|-------------|--------------|-------|
| SAM    | ~5–10 GB    | Sequential   | text editors, awk |
| BAM    | ~1–3 GB     | ✓ (with .bai) | samtools, pysam |
| CRAM   | ~0.5–1.5 GB | ✓ (with .crai) | samtools, htslib |

```python
# ── SAM FLAG decoder ─────────────────────────────────────────────────────
FLAG_BITS = {
    1:    'paired',
    2:    'proper_pair',
    4:    'unmapped',
    8:    'mate_unmapped',
    16:   'reverse_strand',
    32:   'mate_reverse_strand',
    64:   'read1',
    128:  'read2',
    256:  'secondary',
    512:  'qc_fail',
    1024: 'duplicate',
    2048: 'supplementary',
}

def decode_flag(flag):
    """Return a dict of set flag bits."""
    return {name: bool(flag & bit) for bit, name in FLAG_BITS.items()}

def flag_summary(flag):
    bits_set = [name for bit, name in FLAG_BITS.items() if flag & bit]
    return ', '.join(bits_set) if bits_set else 'none'

# Typical Illumina paired-end properly-mapped read1
for flag in [99, 147, 4, 1024, 2064]:
    print(f"FLAG {flag:>5} (0x{flag:04X}): {flag_summary(flag)}")
```python

## Common Pitfalls

- **Coordinate systems**: BED uses 0-based half-open; VCF/GFF use 1-based inclusive — mixing them causes off-by-one errors
- **Batch effects**: Always check for batch confounding before interpreting biological signal
- **Multiple testing**: Apply FDR correction (Benjamini-Hochberg) when testing thousands of features simultaneously
