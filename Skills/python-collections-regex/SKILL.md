---
name: python-collections-regex
description: Python collections (lists, dicts, sets), comprehensions, generators for large bio files, and regex for sequence motifs
primary_tool: Python
---

## Version Compatibility

Reference examples tested with: Python 3.10+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.


# Python Collections & Regex for Bioinformatics

## When to Use
- **list**: ordered, mutable collection of sequences/records; sorting, sliding windows, k-mer extraction
- **tuple**: fixed immutable records — gene coordinates `(name, chr, start, end)`, SNP `(chr, pos, ref, alt)`; dict keys
- **namedtuple**: readable fixed records with field names; drop-in for plain tuples
- **dict**: lookup tables (codon table, gene lengths, annotations); O(1) access
- **defaultdict(list/set/int)**: grouping genes by chromosome, building reverse codon table
- **Counter**: k-mer frequencies, nucleotide counts, amino acid composition
- **set**: gene set operations (union/intersection/difference), deduplication, O(1) membership
- **generator**: streaming FASTA/FASTQ; pipelines over genome-scale data that won't fit in RAM
- **regex**: restriction sites, FASTA/GenBank header parsing, PROSITE motifs, homopolymer detection

## Quick Reference

### Collections
```python
from collections import namedtuple, defaultdict, Counter

Gene = namedtuple('Gene', ['name', 'chrom', 'start', 'end', 'strand'])
g = Gene('BRCA1', 'chr17', 43044295, 43125483, '-')
length = g.end - g.start

by_chrom = defaultdict(list)          # grouping
kmer_counts = Counter(kmers)          # counting; .most_common(n)
kmers_a & kmers_b                     # Counter intersection (min counts)
kmers_a + kmers_b                     # Counter union (sum counts)
```python

### Set operations on gene lists
```python
cancer = {"BRCA1", "TP53", "EGFR", "MYC", "KRAS"}
repair = {"BRCA1", "BRCA2", "ATM", "MLH1", "TP53"}
cancer | repair          # union
cancer & repair          # intersection -> {"BRCA1", "TP53"}
cancer - repair          # difference (cancer-only)
cancer ^ repair          # symmetric difference
cancer.issubset(repair)  # subset test
```python

### re module
```python
import re
re.search(pat, s)        # first match anywhere; returns match object or None
re.findall(pat, s)       # list of all match strings
re.finditer(pat, s)      # iterator of match objects (use for positions)
re.sub(pat, repl, s)     # replace; repl can be string or callable
re.split(pat, s)         # split at matches
re.compile(pat, flags)   # pre-compile for repeated use

m.group()   m.group(1)  m.group('name')   # text
m.start()   m.end()     m.span()           # positions
```python

## Key Patterns

### Codon table (dict)
```python
GENETIC_CODE = {
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

# Reverse: amino acid -> set of codons
from collections import defaultdict
reverse_code = defaultdict(set)
for codon, aa in GENETIC_CODE.items():
    reverse_code[aa].add(codon)
```python

### Pathway / gene set analysis
```python
pathways = {
    'DNA_repair':  {'BRCA1','BRCA2','ATM','TP53','MLH1'},
    'Cell_cycle':  {'TP53','RB1','CDK4','CCND1','MYC'},
    'Apoptosis':   {'TP53','BCL2','BAX','CASP3','FAS'},
}
tp53_pathways = [name for name, genes in pathways.items() if 'TP53' in genes]
all_genes = set().union(*pathways.values())
```python

### Restriction enzyme map (dict + regex)
```python
RESTRICTION_ENZYMES = {
    'EcoRI':'GAATTC', 'BamHI':'GGATCC', 'HindIII':'AAGCTT',
    'NotI':'GCGGCCGC', 'XhoI':'CTCGAG', 'NdeI':'CATATG',
}
def restriction_map(seq, enzymes):
    return {e: [m.start() for m in re.finditer(site, seq)]
            for e, site in enzymes.items()
            if re.search(site, seq)}

def predict_fragments(seq, site):
    return [len(f) for f in re.split(site, seq) if f]
```python

## Code Templates

### K-mer counting with Counter
```python
from collections import Counter

def kmer_counts(seq, k):
    return Counter(seq[i:i+k] for i in range(len(seq) - k + 1))

# Compare two sequences
profile_a = Counter(seq_a[i:i+3] for i in range(len(seq_a) - 2))
profile_b = Counter(seq_b[i:i+3] for i in range(len(seq_b) - 2))
shared = profile_a & profile_b          # min counts
jaccard = len(profile_a & profile_b) / len(profile_a | profile_b)
```python

### Comprehension idioms
```python
# Codon splitting
codons = [cds[i:i+3] for i in range(0, len(cds) - len(cds) % 3, 3)]

# Translation (up to first stop)
protein = ''.join(
    aa for aa in (GENETIC_CODE.get(c, 'X') for c in codons)
    if aa != '*'
)

# Reverse complement (one-liner)
rc = {'A':'T','T':'A','G':'C','C':'G'}
rev_comp = ''.join(rc[n] for n in seq[::-1])

# GC content
gc = sum(1 for n in seq if n in 'GC') / len(seq)

# Filter sequences by GC and length
filtered = {name: seq for name, seq in fasta.items()
            if len(seq) > 100 and 0.4 <= gc_content(seq) <= 0.6}

# Sliding window GC profile
gc_profile = [(i, (w.count('G') + w.count('C')) / window)
              for i in range(len(seq) - window + 1)
              for w in [seq[i:i+window]]]

# All reading frames (list of lists)
frames = [[dna[i:i+3] for i in range(f, len(dna) - 2, 3)] for f in range(3)]

# Generate all k-mers of length k (4^k)
from itertools import product
all_kmers = [''.join(c) for c in product('ATGC', repeat=k)]

# Unique k-mer spectrum
spectrum = sorted({seq[i:i+k] for i in range(len(seq) - k + 1)})
```python

### Streaming FASTA generator
```python
def read_fasta(filename):
    """Yield (header, sequence) tuples. Memory: O(single record)."""
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

# Usage: process without loading all sequences
gc_stats = {hdr.split()[0]: gc_content(seq) for hdr, seq in read_fasta('genome.fa')}
```python

### Streaming FASTQ generator + quality filter pipeline
```python
def read_fastq(filename):
    """Yield {'id', 'sequence', 'quality'} dicts. One record at a time."""
    with open(filename) as f:
        while True:
            header = f.readline().strip()
            if not header:
                break
            seq  = f.readline().strip()
            f.readline()                      # '+' line
            qual = f.readline().strip()
            yield {'id': header[1:].split()[0], 'sequence': seq, 'quality': qual}

def quality_filter(records, min_q=30):
    for r in records:
        if sum(ord(c) - 33 for c in r['quality']) / len(r['quality']) >= min_q:
            yield r

def trim_ns(records):
    for r in records:
        s = r['sequence'].rstrip('N')
        yield {**r, 'sequence': s, 'quality': r['quality'][:len(s)]}

# Chain generators — nothing computed until consumed
pipeline = trim_ns(quality_filter(read_fastq('reads.fastq')))
passing = sum(1 for _ in pipeline)
```python

### Regex: biological sequence patterns
```python
import re

# Stop codons
re.findall(r'TAA|TAG|TGA', dna)

# Find all ORFs (ATG ... stop, in-frame)
for start_m in re.finditer(r'ATG', dna):
    s = start_m.start()
    stop_m = re.search(r'TAA|TAG|TGA', dna[s+3:])   # rough; check frame separately
    if stop_m and stop_m.start() % 3 == 0:
        orf = dna[s : s+3+stop_m.end()]

# Homopolymer runs >= 5 bases
for m in re.finditer(r'(.)\1{4,}', dna):
    print(m.group(), m.start())

# Tandem repeats (any triplet repeated 2+ times)
re.finditer(r'([ATGC]{3})\1+', dna)

# Non-standard bases
invalid = re.findall(r'[^ATGCatgcNn]', seq)

# Overlapping k-mer positions using lookahead
positions = re.findall(r'(?=(ATG))', dna)

# Greedy vs lazy ORF matching
re.findall(r'ATG.*TAG', dna)     # greedy: longest span
re.findall(r'ATG.*?TAG', dna)    # lazy: shortest span
```python

### Regex: FASTA/GenBank header parsing
```python
# UniProt SwissProt header
# >sp|P04637|P53_HUMAN Cellular tumor antigen p53 OS=Homo sapiens
pat = re.compile(
    r'>(?:sp|tr)\|(?P<accession>[^|]+)\|(?P<entry>\S+)\s+'
    r'(?P<desc>.+?)\s+OS=(?P<organism>[^=]+?)(?:\s+\w+=|$)'
)

# NCBI RefSeq
# >NM_001234.5 description [Homo sapiens]
pat = re.compile(r'>(?P<acc>[A-Z]{2}_\d+\.?\d*)\s+(?P<desc>.+?)\s*\[(?P<org>[^\]]+)\]')

# GenBank feature table
genes    = re.findall(r'/gene="([^"]+)"', gb_text)
prot_ids = re.findall(r'/protein_id="([^"]+)"', gb_text)
coords   = re.findall(r'(\d+)\.\.(\d+)', gb_text)       # list of (start, end) strings
```python

### PROSITE pattern -> Python regex converter
```python
def prosite_to_regex(pattern):
    """Convert PROSITE pattern to Python regex.
    N-{P}-[ST]-{P}  ->  N[^P][ST][^P]
    x -> .   (n) -> {n}   (n,m) -> {n,m}
    """
    result, parts = [], pattern.replace('-', '')
    i = 0
    while i < len(parts):
        if parts[i] == '[':
            j = parts.index(']', i)
            result.append(parts[i:j+1]); i = j+1
        elif parts[i] == '{':
            j = parts.index('}', i)
            result.append(f'[^{parts[i+1:j]}]'); i = j+1
        elif parts[i] == '(':
            j = parts.index(')', i)
            result.append(f'{{{parts[i+1:j]}}}'); i = j+1
        elif parts[i] == 'x':
            result.append('.'); i += 1
        else:
            result.append(parts[i]); i += 1
    return ''.join(result)

# Common PROSITE patterns
PROSITE = {
    'N-glycosylation':     'N-{P}-[ST]-{P}',
    'PKC_phosphorylation': '[ST]-x-[RK]',
    'CK2_phosphorylation': '[ST]-x(2)-[DE]',
}
for name, pat in PROSITE.items():
    for m in re.finditer(prosite_to_regex(pat), protein):
        print(f"{name} at {m.start()+1}: {m.group()}")
```python

## Common Pitfalls

| Pitfall | Fix |
|---|---|
| `alias = my_list` shares reference | Use `.copy()` or `[:]` for a shallow copy |
| `{}` creates empty dict, not set | Use `set()` for an empty set |
| `re.findall` with capture group returns group text, not full match | Use `(?:...)` non-capturing group or `re.finditer` |
| Greedy `ATG.*TAG` spans multiple ORFs | Use lazy `ATG.*?TAG` |
| `re.match` only checks string start | Use `re.search` to match anywhere |
| Dict iteration while modifying | Iterate over `list(d.items())` |
| `namedtuple` field access by index still works — silent wrong field | Use field names exclusively |
| Generator exhausted silently after one pass | Re-create the generator; generators are single-pass |
| `defaultdict` creates key on `d[key]` even for reads | Use `.get()` if you don't want side-effects |

## Related Skills
- `numpy-pandas-wrangling` — vectorized sequence/count operations, DataFrame from Counter
- `python-advanced-sql` — storing gene annotations in SQLite, querying with dicts
