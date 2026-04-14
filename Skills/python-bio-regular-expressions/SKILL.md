---
name: python-bio-regular-expressions
description: Regular expressions for biological sequence pattern matching — ORF finding, motifs, restriction sites, FASTA headers, overlapping matches
tool_type: python
primary_tool: Python
---

# Regular Expressions for Bioinformatics

## Quick Reference

| Function | Returns | Use |
|----------|---------|-----|
| `re.search(pat, s)` | First match object or None | Check if pattern exists |
| `re.findall(pat, s)` | List of strings (or tuples with groups) | All non-overlapping matches |
| `re.finditer(pat, s)` | Iterator of match objects | All matches with positions |
| `re.sub(pat, repl, s)` | Modified string | Replace matches |
| `re.split(pat, s)` | List of strings | Split on pattern |
| `re.compile(pat)` | Compiled regex | Reuse expensive patterns |

## Bio-Specific Patterns

```python
import re

# Stop codons
stop = re.compile(r'TAA|TAG|TGA')

# Overlapping ATG starts — requires lookahead
for m in re.finditer(r'(?=(ATG))', dna):
    print(m.start(), m.group(1))

# ORF: greedy (longest) vs lazy (shortest)
greedy = re.findall(r'ATG.*TAG', dna)    # first ATG to LAST TAG
lazy   = re.findall(r'ATG.*?TAG', dna)  # first ATG to nearest TAG

# Homopolymer runs (4+ bases)
for m in re.finditer(r'(.)\1{3,}', dna):
    print(m.group(), m.start(), len(m.group()))

# Tandem repeats
for m in re.finditer(r'([ATGC]{3})\1+', dna):
    unit = m.group(1)
    copies = len(m.group()) // len(unit)

# IUPAC ambiguity: translate before matching
IUPAC = {'N':'[ATGC]','R':'[AG]','Y':'[CT]','W':'[AT]','S':'[GC]',
         'M':'[AC]','K':'[GT]','B':'[CGT]','D':'[AGT]','H':'[ACT]','V':'[ACG]'}
def iupac_to_regex(seq):
    return ''.join(IUPAC.get(b, b) for b in seq.upper())

# Restriction digest (overlapping sites)
def find_cut_positions(dna, site, cut_offset):
    return [m.start() + cut_offset for m in re.finditer(f'(?={site})', dna.upper())]

# FASTA header parsing (named groups)
header = ">sp|P04637|P53_HUMAN Cellular tumor antigen OS=Homo sapiens"
pat = r'>sp\|(?P<acc>[^|]+)\|(?P<entry>\S+)\s+(?P<desc>.+?)\s+OS=(?P<org>.+)'
m = re.search(pat, header)
# m.group('acc'), m.group('org') etc.

# Transcription: T→U
rna = re.sub('T', 'U', dna)

# Split at restriction sites
fragments = re.split(r'GAATTC|GGATCC|AAGCTT', dna)

# Multi-line FASTA headers
headers = re.findall(r'^>.*', fasta, re.MULTILINE)

# Non-standard nucleotides
non_std = re.findall(r'[^ATGCN]', seq)

# Mask homopolymers with N
def mask_repeat(m):
    return 'N' * len(m.group()) if len(m.group()) >= 5 else m.group()
masked = re.sub(r'(.)\1{2,}', mask_repeat, dna)
```

## Lookaround Quick Reference

| Syntax | Meaning | Bio use |
|--------|---------|---------|
| `(?=pat)` | Positive lookahead | Overlapping matches; ATG not consuming |
| `(?!pat)` | Negative lookahead | ATG not followed by stop |
| `(?<=pat)` | Positive lookbehind | Bases after restriction site |
| `(?<!pat)` | Negative lookbehind | Context-excluded matches |

## Pitfalls

- **Greedy vs non-greedy**: `ATG.*TAG` matches ATG to the *last* TAG in the string; `ATG.*?TAG` stops at the *nearest* — always choose deliberately for ORF searches
- **Overlapping matches**: `re.findall('ATG', dna)` skips overlaps; use `re.finditer(r'(?=(ATG))', dna)` to find all
- **IUPAC codes are not regex**: `N` in a sequence means any base, but `N` in a regex matches the literal `N` — always translate with `iupac_to_regex()` first
- **`re.MULTILINE` for FASTA**: `^>` only anchors to the start of the whole string by default; add `re.MULTILINE` to match each line
- **`findall` with groups**: with one capturing group returns group contents (not full match); with multiple groups returns list of tuples
- **Off-by-one with `match.start()`**: positions are 0-based; bioinformatics coordinates are often 1-based — add 1 when reporting
