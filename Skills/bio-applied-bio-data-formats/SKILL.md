---
name: bio-applied-bio-data-formats
description: Quick reference for bioinformatics file formats — FASTA, FASTQ, SAM/BAM/CRAM, VCF, BED, GFF/GTF, BigWig, PDB, Newick — specs, coordinate systems, and parsing patterns.
tool_type: python
primary_tool: Python
---

## Format Quick Reference

| Format | Content | Coord system | Random access | Tools |
|--------|---------|-------------|--------------|-------|
| FASTA | Sequences | N/A | With .fai index | samtools faidx, BioPython |
| FASTQ | Reads + quality | N/A | No (stream) | fastp, cutadapt |
| SAM | Alignments (text) | 1-based | No | samtools, awk |
| BAM | Alignments (binary) | 1-based | With .bai | samtools, pysam |
| CRAM | Alignments (ref-compressed) | 1-based | With .crai | samtools, htslib |
| VCF | Variants | 1-based | With .tbi (tabix) | bcftools, PyVCF |
| BED | Intervals | **0-based half-open** | With .bai | bedtools, pybedtools |
| GFF3/GTF | Gene annotations | **1-based inclusive** | With tabix | gffutils |
| BigWig | Coverage/signal | 0-based | Yes (UCSC) | pyBigWig, deeptools |
| PDB/mmCIF | 3D structure | 1-based residues | No | Bio.PDB, MDAnalysis |
| Newick | Phylogenetic tree | N/A | N/A | Bio.Phylo, ete3 |
| FAST5/POD5 | ONT raw signal | N/A | No | pod5, ont-fast5-api |

## Format Specs

### FASTA
```
>sequence_id [optional description]
AGCTAGCTAGCTAGCTAGCT   ← wrap at 60–80 chars (conventional)
```
- Lowercase = soft-masked (repeat) regions
- `.fai` index: `samtools faidx file.fa` → enables O(1) region fetch

### FASTQ
```
@SEQ_ID [optional description]
AGCTAGCTAGCTAGCT       ← sequence (must NOT wrap)
+                      ← separator
IIIIIIIIIIIIIIII       ← quality, same length as sequence
```

Phred encoding: `Q = -10 × log₁₀(P_error)`, ASCII offset 33 (Illumina ≥1.8 / Sanger standard)

| Q | Error prob | Accuracy |
|---|-----------|----------|
| 10 | 1/10 | 90% |
| 20 | 1/100 | 99% |
| 30 | 1/1,000 | 99.9% |
| 40 | 1/10,000 | 99.99% |

### SAM FLAG field

| Bit | Hex | Meaning |
|-----|-----|---------|
| 1 | 0x1 | Paired |
| 4 | 0x4 | Read unmapped |
| 16 | 0x10 | Reverse strand |
| 64 | 0x40 | Read 1 |
| 128 | 0x80 | Read 2 |
| 256 | 0x100 | Secondary alignment |
| 1024 | 0x400 | PCR duplicate |
| 2048 | 0x800 | Supplementary |

```python
FLAG_BITS = {1:'paired',4:'unmapped',16:'reverse',64:'read1',
             128:'read2',256:'secondary',1024:'duplicate',2048:'supplementary'}
def decode_flag(flag):
    return [name for bit, name in FLAG_BITS.items() if flag & bit]
```

### CIGAR Operations

| Op | Consumes query | Consumes ref | Meaning |
|----|---------------|-------------|---------|
| M | Yes | Yes | Match/mismatch |
| I | Yes | No | Insertion |
| D | No | Yes | Deletion |
| N | No | Yes | Intron skip (RNA-seq) |
| S | Yes | No | Soft clip |
| H | No | No | Hard clip |

### VCF Essential Columns
```
CHROM  POS  ID  REF  ALT  QUAL  FILTER  INFO  FORMAT  sample1...
chr1   925952  .  G   A    .     PASS    AF=0.5  GT:DP  0/1:30
```
- POS is **1-based**; REF/ALT for indels include anchor base

### BED vs GFF coordinate trap
```
BED:  chr1  0   100   → covers bases 1–100 (0-based, half-open)
GFF:  chr1  1   100   → covers bases 1–100 (1-based, inclusive)
VCF:  chr1  925952    → position 925952 (1-based)
```

## Parsing Patterns

### Streaming FASTQ parser
```python
def parse_fastq(source):
    import gzip, io
    opener = gzip.open if str(source).endswith('.gz') else open
    with opener(source, 'rt') as fh:
        for line in fh:
            header = line.rstrip('\n')[1:]
            seq    = next(fh).rstrip('\n')
            next(fh)                       # skip '+'
            qual   = next(fh).rstrip('\n')
            yield header, seq, qual

def phred(char, offset=33):
    return ord(char) - offset
```

### FASTA random access with .fai
```python
def faidx_fetch(fasta_path, fai_path, chrom, start, end):
    """0-based half-open [start, end)."""
    import re
    index = {}
    with open(fai_path) as fh:
        for line in fh:
            name, length, offset, bases_per_line, bytes_per_line = line.split()
            index[name] = (int(offset), int(bases_per_line), int(bytes_per_line))
    offset, bpl, bypl = index[chrom]
    n_lines, rem = divmod(start, bpl)
    byte_start = offset + n_lines * bypl + rem
    bases_needed = end - start
    with open(fasta_path, 'rb') as fh:
        fh.seek(byte_start)
        raw = fh.read((bases_needed // bpl + 2) * bypl).decode()
    return re.sub(r'\s', '', raw)[:bases_needed]
```

### pysam for BAM
```python
import pysam
bam = pysam.AlignmentFile("sample.bam", "rb")
for read in bam.fetch("chr1", 1000, 2000):
    if read.is_unmapped or read.is_duplicate:
        continue
    print(read.query_name, read.reference_start, read.cigarstring)
bam.close()
```

## Pitfalls

- **Coordinate systems**: BED is 0-based half-open; VCF/GFF/SAM are 1-based inclusive. Mixing them causes silent off-by-one errors — the most common bioinformatics bug.
- **FASTQ must not wrap sequence lines**: unlike FASTA, the sequence and quality strings must each be exactly one line.
- **BAM requires sorted + indexed file for `fetch()`**: `samtools sort` then `samtools index`. Unsorted BAM raises an exception or returns nothing.
- **CRAM requires the reference genome at decode time**: store the reference path in `@HD UR:` header or set `REF_PATH` env var.
- **VCF indels include anchor base**: `REF=GACT ALT=G` is a 3-bp deletion, not 4-bp. POS points to the anchor base.
- **Always gzip + index VCFs with bgzf, not gzip**: tabix requires bgzf compression; standard gzip blocks random access.
- **GTF attribute parsing**: GTF uses key `"value"` semicolon-delimited format; do not split on `;` naively — values can contain semicolons inside quotes.
