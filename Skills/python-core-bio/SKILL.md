---
name: python-core-bio
description: Core Python for bioinformatics — types, bio-string manipulation (codons, reverse complement), functions, and file I/O for FASTA/GenBank
primary_tool: Python
---

## Version Compatibility

Reference examples tested with: Python 3.10+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.


# Python Core for Bioinformatics

## When to Use
- Writing DNA/RNA/protein sequence manipulation code from scratch
- Parsing FASTA, FASTQ, CSV, TSV, JSON, GenBank files
- Implementing GC content, codon tables, ORF finding, reverse complement
- Designing reusable functions for sequence pipelines
- Handling large files without loading them into memory

## Quick Reference

### DNA/RNA Alphabet and Validation
```python
VALID_DNA = set("ATGC")
VALID_RNA = set("AUGC")
is_valid = all(c in VALID_DNA for c in seq.upper())  # use all() + generator

def detect_seq_type(seq):
    chars = set(seq.upper())
    if chars <= set("ATGC"): return "DNA"
    if chars <= set("AUGC"): return "RNA"
    if chars <= set("ACDEFGHIKLMNPQRSTVWY"): return "Protein"
    return "Unknown"
```

### Reverse Complement
```python
RC_TABLE = str.maketrans("ATGC", "TACG")
rc = seq.upper().translate(RC_TABLE)[::-1]
complement = seq.upper().translate(RC_TABLE)          # 3'->5', no reversal
```

### GC Content
```python
def gc_content(seq: str) -> float:
    s = seq.upper()
    return (s.count("G") + s.count("C")) / len(s) * 100
```

### Transcription / Reverse Transcription
```python
mrna = dna.upper().replace("T", "U")   # DNA coding strand -> mRNA
dna  = rna.upper().replace("U", "T")   # mRNA -> DNA coding strand
```

### Codon Extraction
```python
codons = [seq[i:i+3] for i in range(0, len(seq), 3) if len(seq[i:i+3]) == 3]
# reading frame offset: position % 3 gives frame (0, 1, 2)
```

### Codon Table (Standard Genetic Code)
```python
CODON_TABLE = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L',
    'CTT':'L','CTC':'L','CTA':'L','CTG':'L',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M',
    'GTT':'V','GTC':'V','GTA':'V','GTG':'V',
    'TCT':'S','TCC':'S','TCA':'S','TCG':'S',
    'CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'ACT':'T','ACC':'T','ACA':'T','ACG':'T',
    'GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*',
    'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K',
    'GAT':'D','GAC':'D','GAA':'E','GAG':'E',
    'TGT':'C','TGC':'C','TGA':'*','TGG':'W',
    'CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'AGT':'S','AGC':'S','AGA':'R','AGG':'R',
    'GGT':'G','GGC':'G','GGA':'G','GGG':'G',
}
STOP_CODONS = {'TAA', 'TAG', 'TGA'}
```

### Translation
```python
def translate(dna: str) -> str:
    dna = dna.upper()
    protein = []
    for i in range(0, len(dna) - 2, 3):
        aa = CODON_TABLE.get(dna[i:i+3], '?')
        if aa == '*': break
        protein.append(aa)
    return ''.join(protein)
```

## Key Patterns

### Melting Temperature
```python
def tm(primer: str) -> float:
    s = primer.upper()
    a, t, g, c = s.count('A'), s.count('T'), s.count('G'), s.count('C')
    if len(s) < 14:
        return 2 * (a + t) + 4 * (g + c)          # Wallace rule
    return 64.9 + 41 * (g + c - 16.4) / len(s)    # salt-adjusted
```

### DNA Molecular Weight (ssDNA)
```python
NUC_MW = {'A': 331.2, 'T': 322.2, 'G': 347.2, 'C': 307.2}
def dna_mw(seq: str) -> float:
    s = seq.upper()
    return sum(NUC_MW[n] for n in s) - (len(s) - 1) * 18.02
```

### ORF Finding (all 3 forward frames)
```python
def find_orfs(seq: str, min_len: int = 30) -> list[dict]:
    seq = seq.upper()
    orfs = []
    for frame in range(3):
        i = frame
        while i <= len(seq) - 3:
            if seq[i:i+3] == 'ATG':
                j = i + 3
                while j <= len(seq) - 3:
                    if seq[j:j+3] in STOP_CODONS:
                        length = j + 3 - i
                        if length >= min_len:
                            orfs.append({'start': i, 'end': j+3,
                                         'length': length, 'frame': frame+1,
                                         'seq': seq[i:j+3]})
                        break
                    j += 3
            i += 3
    return orfs

# Both strands: run find_orfs on seq and on rc separately, tag strand +/-
```

### Sliding Window GC
```python
def sliding_gc(seq: str, window: int = 100, step: int = 1) -> list[tuple]:
    s = seq.upper()
    return [(i, (s[i:i+window].count('G') + s[i:i+window].count('C')) / window * 100)
            for i in range(0, len(s) - window + 1, step)]
```

### Motif / Restriction Site Finder
```python
def find_motif(seq: str, motif: str) -> list[int]:
    positions, pos = [], seq.find(motif)
    while pos != -1:
        positions.append(pos)
        pos = seq.find(motif, pos + 1)
    return positions
```

### Palindrome Check (restriction enzyme sites)
```python
def is_palindrome(seq: str) -> bool:
    s = seq.upper()
    return s == s.translate(RC_TABLE)[::-1]
```

### Sequence Identity (Hamming)
```python
def identity(s1: str, s2: str) -> float:
    matches = sum(a == b for a, b in zip(s1, s2))
    return matches / len(s1) * 100
```

### CDS Validation
```python
def validate_cds(seq: str) -> list[tuple]:
    s = seq.upper().strip()
    checks = []
    invalid = set(s) - set('ATGC')
    checks.append(("Valid chars", not invalid, str(invalid) if invalid else "OK"))
    checks.append(("Length % 3", len(s) % 3 == 0, f"{len(s)} bp"))
    checks.append(("Start codon", s[:3] == 'ATG', s[:3]))
    checks.append(("Stop codon", s[-3:] in STOP_CODONS, s[-3:]))
    return checks
```

## Code Templates

### FASTA Parser (memory-efficient generator)
```python
def parse_fasta(filename: str):
    """Yield (header, sequence) tuples."""
    header, parts = None, []
    with open(filename) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if header is not None:
                    yield header, ''.join(parts)
                header, parts = line[1:], []
            else:
                parts.append(line)
    if header is not None:
        yield header, ''.join(parts)

# Usage:
for header, seq in parse_fasta('seqs.fasta'):
    seq_id = header.split()[0]
    print(seq_id, len(seq))
```

### FASTA Writer (with line wrapping)
```python
def write_fasta(seqs: dict, filename: str, width: int = 60):
    with open(filename, 'w') as f:
        for header, seq in seqs.items():
            f.write(f">{header}\n")
            for i in range(0, len(seq), width):
                f.write(seq[i:i+width] + '\n')
```

### FASTQ Parser
```python
def parse_fastq(filename: str):
    """Yield (id, sequence, quality_scores) tuples."""
    with open(filename) as f:
        while True:
            header = f.readline().strip()
            if not header: break
            seq   = f.readline().strip()
            f.readline()                            # '+' line
            qual  = f.readline().strip()
            scores = [ord(c) - 33 for c in qual]   # Phred+33 encoding
            yield header[1:], seq, scores
```

### CSV / TSV Reading and Writing
```python
import csv

# Read (DictReader gives column-name access)
with open('genes.csv') as f:
    for row in csv.DictReader(f):
        gc = float(row['gc_content'])

# Read TSV
with open('genes.bed') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        ...

# Write with header
with open('results.csv', 'w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=['gene', 'fold_change', 'p_value'])
    writer.writeheader()
    writer.writerows(results)   # results: list of dicts
```

### GenBank Minimal Sequence Extractor
```python
def extract_genbank_seq(filename: str) -> str:
    """Extract raw sequence from GenBank flat file."""
    seq_parts, in_origin = [], False
    with open(filename) as f:
        for line in f:
            if line.startswith('ORIGIN'):
                in_origin = True; continue
            if line.startswith('//'):
                break
            if in_origin:
                seq_parts.append(''.join(c for c in line if c.isalpha()))
    return ''.join(seq_parts).upper()
# Full parsing: use Bio.SeqIO.read(filename, "genbank")
```

### Parse TSV Data Line (common in bioinformatics)
```python
# e.g., "BRCA1\t17\t43044295\t43170245\t0.423"
gene, chrom, start, end, gc = line.split('\t')
chrom, start, end = int(chrom), int(start), int(end)
gc = float(gc)
```

### UniProt FASTA Header Parsing
```python
# ">sp|P04637|P53_HUMAN Cellular tumor antigen p53 OS=Homo sapiens GN=TP53 PE=1 SV=4"
header = line[1:]                          # strip '>'
parts = header.split('|')
accession  = parts[1]                      # P04637
entry_name = parts[2].split()[0]           # P53_HUMAN
gn = next((p[3:] for p in header.split() if p.startswith('GN=')), None)  # TP53
```

### Functional Patterns for Sequence Collections
```python
sequences = [("seq1", "ATGCGCGC"), ("seq2", "ATATAT"), ("seq3", "GCGCGCGC")]

# Sort by GC content descending
sorted_seqs = sorted(sequences, key=lambda t: gc_content(t[1]), reverse=True)

# Filter GC-rich
gc_rich = list(filter(lambda t: gc_content(t[1]) > 50, sequences))

# Map: uppercase all
upper = list(map(lambda t: (t[0], t[1].upper()), sequences))

# List comprehension (preferred over map/filter for readability)
lengths = {sid: len(seq) for sid, seq in sequences}
```

### File Path Handling
```python
from pathlib import Path

p = Path('data/sequences/gene.fasta')
p.stem    # 'gene'
p.suffix  # '.fasta'
p.parent  # Path('data/sequences')

output = Path('results') / 'analysis' / 'output.csv'
output.parent.mkdir(parents=True, exist_ok=True)

for fasta_file in Path('data').glob('*.fasta'):
    ...
```

## Common Pitfalls

| Pitfall | Fix |
|---|---|
| `seq[0] = 'G'` raises TypeError | Strings are immutable: `seq = 'G' + seq[1:]` |
| Mixed case breaks counting | Always `.upper()` before processing |
| `0.1 + 0.2 == 0.3` is False | Use `abs(a - b) < 1e-9` for float comparison |
| `int(3.7)` truncates to 3, not 4 | Use `round()` when rounding is intended |
| Mutable default arg `def f(items=[])` | Use `items=None`, then `if items is None: items = []` |
| `int()` on file values — forgetting to convert | All values from file reads are `str`; cast explicitly |
| `seq.count('ATG')` counts overlaps | Returns overlapping matches; fine for simple motifs |
| `is` vs `==` for None | Always use `if x is None`, never `if x == None` |
| `f.close()` not called on error | Always use `with open(...) as f:` |
| Loading entire genome into RAM | Use generator-based parsers (`yield`) for large files |
| GC formula without parentheses: `g + c / total * 100` | Always `(g + c) / total * 100` |
| `3 % 3 == 0` reading frame confusion | Frame 0 starts at index 0; frame = `pos % 3` |

## Related Skills
- `biopython-databases` — SeqIO, Entrez, BLAST wrappers, SeqRecord objects
- `numpy-pandas-wrangling` — vectorised sequence stats, expression data, DataFrames
- `data-visualization-bio` — GC profile plots, sequence logos, coverage plots
- `python-advanced-sql` — when pipeline results go into SQLite/PostgreSQL
