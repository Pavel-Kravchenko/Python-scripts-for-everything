---
name: python-collections-regex
description: Python collections (namedtuple, defaultdict, Counter, set) and regex for bioinformatics — k-mer counting, FASTA/FASTQ streaming, restriction maps, PROSITE patterns.
primary_tool: Python
---

# Python Collections & Regex for Bioinformatics

## When to Use

| Type | Use case |
|---|---|
| `list` | ordered sequences, sliding windows, k-mer extraction |
| `tuple` | fixed records: gene coords `(name, chr, start, end)`, SNP `(chr, pos, ref, alt)` |
| `namedtuple` | readable fixed records with field names; dict keys |
| `dict` | codon table, gene lengths, annotations; O(1) access |
| `defaultdict(list/set/int)` | grouping by chromosome, reverse codon table |
| `Counter` | k-mer frequencies, nucleotide counts, amino acid composition |
| `set` | gene set ops (union/intersection/difference), deduplication, O(1) membership |
| `generator` | streaming FASTA/FASTQ; genome-scale pipelines that won't fit in RAM |
| `regex` | restriction sites, FASTA/GenBank header parsing, PROSITE motifs, homopolymers |

## Pitfalls

| Pitfall | Fix |
|---|---|
| `alias = my_list` shares reference | `.copy()` or `[:]` for shallow copy |
| `{}` creates empty dict, not set | Use `set()` for an empty set |
| `re.findall` with capture group returns group text, not full match | Use `(?:...)` non-capturing or `re.finditer` |
| Greedy `ATG.*TAG` spans multiple ORFs | Use lazy `ATG.*?TAG` |
| `re.match` only checks string start | Use `re.search` to match anywhere |
| Dict iteration while modifying | Iterate over `list(d.items())` |
| Generator exhausted silently after one pass | Re-create the generator; generators are single-pass |
| `defaultdict` creates key on `d[key]` even for reads | Use `.get()` to avoid side-effects |

## Quick Reference

### Collections
```python
from collections import namedtuple, defaultdict, Counter

Gene = namedtuple('Gene', ['name', 'chrom', 'start', 'end', 'strand'])
g = Gene('BRCA1', 'chr17', 43044295, 43125483, '-')
length = g.end - g.start

by_chrom = defaultdict(list)       # grouping
kmer_counts = Counter(kmers)       # counting; .most_common(n)
kmers_a & kmers_b                  # Counter intersection (min counts)
kmers_a + kmers_b                  # Counter union (sum counts)
```

### Set operations on gene lists
```python
cancer = {"BRCA1", "TP53", "EGFR", "MYC", "KRAS"}
repair = {"BRCA1", "BRCA2", "ATM", "MLH1", "TP53"}
cancer | repair          # union
cancer & repair          # intersection -> {"BRCA1", "TP53"}
cancer - repair          # cancer-only
cancer ^ repair          # symmetric difference
```

### re module API
```python
import re
re.search(pat, s)        # first match anywhere; None if no match
re.findall(pat, s)       # list of all match strings
re.finditer(pat, s)      # iterator of match objects (use for positions)
re.sub(pat, repl, s)     # replace; repl can be callable
re.compile(pat, flags)   # pre-compile for repeated use

m.group()  m.group(1)  m.group('name')   # text
m.start()  m.end()     m.span()           # positions
```

## Key Patterns

### Codon table + reverse
```python
GENETIC_CODE = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V',
    'TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E',
    'TGT':'C','TGC':'C','TGA':'*','TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'AGT':'S','AGC':'S','AGA':'R','AGG':'R','GGT':'G','GGC':'G','GGA':'G','GGG':'G',
}
reverse_code = defaultdict(set)
for codon, aa in GENETIC_CODE.items():
    reverse_code[aa].add(codon)
```

### K-mer counting
```python
from collections import Counter

def kmer_counts(seq, k):
    return Counter(seq[i:i+k] for i in range(len(seq) - k + 1))

# Jaccard similarity between two sequences
profile_a = kmer_counts(seq_a, 3)
profile_b = kmer_counts(seq_b, 3)
jaccard = len(profile_a & profile_b) / len(profile_a | profile_b)
```

### Comprehension idioms
```python
# Codon splitting
codons = [cds[i:i+3] for i in range(0, len(cds) - len(cds) % 3, 3)]

# Translation (up to first stop)
protein = ''.join(aa for aa in (GENETIC_CODE.get(c, 'X') for c in codons) if aa != '*')

# Reverse complement
rc = {'A':'T','T':'A','G':'C','C':'G'}
rev_comp = ''.join(rc[n] for n in seq[::-1])

# GC content
gc = sum(1 for n in seq if n in 'GC') / len(seq)

# All reading frames
frames = [[dna[i:i+3] for i in range(f, len(dna)-2, 3)] for f in range(3)]

# All k-mers of length k
from itertools import product
all_kmers = [''.join(c) for c in product('ATGC', repeat=k)]
```

### Streaming FASTA generator
```python
def read_fasta(filename):
    """Yield (header, sequence) tuples. Memory: O(single record)."""
    header, parts = None, []
    with open(filename) as f:
        for line in f:
            line = line.strip()
            if not line: continue
            if line.startswith('>'):
                if header is not None:
                    yield header, ''.join(parts)
                header, parts = line[1:], []
            else:
                parts.append(line)
    if header is not None:
        yield header, ''.join(parts)
```

### Streaming FASTQ pipeline
```python
def read_fastq(filename):
    with open(filename) as f:
        while True:
            header = f.readline().strip()
            if not header: break
            seq  = f.readline().strip()
            f.readline()    # '+' line
            qual = f.readline().strip()
            yield {'id': header[1:].split()[0], 'sequence': seq, 'quality': qual}

def quality_filter(records, min_q=30):
    for r in records:
        if sum(ord(c) - 33 for c in r['quality']) / len(r['quality']) >= min_q:
            yield r

# Chain generators — nothing computed until consumed
pipeline = quality_filter(read_fastq('reads.fastq'))
```

### Regex: biological sequence patterns
```python
import re

re.findall(r'TAA|TAG|TGA', dna)                        # stop codons
re.finditer(r'(.)\1{4,}', dna)                         # homopolymer runs ≥5
re.finditer(r'([ATGC]{3})\1+', dna)                    # tandem repeats
re.findall(r'[^ATGCatgcNn]', seq)                       # non-standard bases
re.findall(r'(?=(ATG))', dna)                           # overlapping matches via lookahead
re.findall(r'ATG.*?TAG', dna)                           # lazy ORF (vs greedy ATG.*TAG)
```

### Regex: FASTA/GenBank header parsing
```python
# UniProt SwissProt: >sp|P04637|P53_HUMAN ...
pat = re.compile(
    r'>(?:sp|tr)\|(?P<accession>[^|]+)\|(?P<entry>\S+)\s+'
    r'(?P<desc>.+?)\s+OS=(?P<organism>[^=]+?)(?:\s+\w+=|$)'
)
# GenBank feature table
genes    = re.findall(r'/gene="([^"]+)"', gb_text)
prot_ids = re.findall(r'/protein_id="([^"]+)"', gb_text)
coords   = re.findall(r'(\d+)\.\.(\d+)', gb_text)
```

### PROSITE pattern → Python regex
```python
def prosite_to_regex(pattern):
    """N-{P}-[ST]-{P}  ->  N[^P][ST][^P]"""
    result, parts = [], pattern.replace('-', '')
    i = 0
    while i < len(parts):
        if parts[i] == '[':
            j = parts.index(']', i); result.append(parts[i:j+1]); i = j+1
        elif parts[i] == '{':
            j = parts.index('}', i); result.append(f'[^{parts[i+1:j]}]'); i = j+1
        elif parts[i] == '(':
            j = parts.index(')', i); result.append(f'{{{parts[i+1:j]}}}'); i = j+1
        elif parts[i] == 'x':
            result.append('.'); i += 1
        else:
            result.append(parts[i]); i += 1
    return ''.join(result)

PROSITE = {
    'N-glycosylation':     'N-{P}-[ST]-{P}',
    'PKC_phosphorylation': '[ST]-x-[RK]',
    'CK2_phosphorylation': '[ST]-x(2)-[DE]',
}
```

### Restriction enzyme map
```python
RESTRICTION_ENZYMES = {
    'EcoRI':'GAATTC', 'BamHI':'GGATCC', 'HindIII':'AAGCTT',
    'NotI':'GCGGCCGC', 'XhoI':'CTCGAG', 'NdeI':'CATATG',
}
def restriction_map(seq, enzymes):
    return {e: [m.start() for m in re.finditer(site, seq)]
            for e, site in enzymes.items() if re.search(site, seq)}

def predict_fragments(seq, site):
    return [len(f) for f in re.split(site, seq) if f]
```
