---
name: python-bio-regular-expressions
description: "**Why this matters:** Biological data is text-heavy -- FASTA headers, GenBank annotations, protein motifs, and sequence patterns all require sophisticated text matching. Regular expressions are the st"
tool_type: python
source_notebook: "Tier_1_Python_for_Bioinformatics/12_Regular_Expressions/01_regular_expressions.ipynb"
primary_tool: Python
---

## Version Compatibility

Reference examples tested with: Python 3.10+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.


# Regular Expressions for Bioinformatics

*Source: Course notebook `Tier_1_Python_for_Bioinformatics/12_Regular_Expressions/01_regular_expressions.ipynb`*


## Learning Objectives
- Master the `re` module: `match`, `search`, `findall`, `finditer`, `sub`, `compile`
- Understand regex patterns: character classes, quantifiers, groups, lookahead/lookbehind
- Use flags: `IGNORECASE`, `MULTILINE`, `DOTALL`, `VERBOSE`
- Apply regex to real bioinformatics tasks: restriction sites, gene IDs, protein motifs, PROSITE patterns

**Why this matters:** Biological data is text-heavy -- FASTA headers, GenBank annotations, protein motifs, and sequence patterns all require sophisticated text matching. Regular expressions are the standard tool for this work.

## How to use this notebook

Run cells from top to bottom on the first pass — later cells depend on functions and variables defined earlier. Once you have run through the notebook, feel free to modify parameters and re-run individual sections.

Each section has runnable examples first, followed by exercises. Try the exercise before looking at the solution cell below it.

## Complicated moments explained

**1. Greedy vs non-greedy**
`ATG.*TAG` matches from the first ATG to the **last** TAG (greedy). `ATG.*?TAG` stops at the **nearest** TAG (non-greedy). In ORF finding this distinction is critical.

**2. `re.findall` with capturing groups**
With no groups, `findall` returns matched strings. With one group, it returns the group contents. With multiple groups, it returns a list of tuples.

**3. Overlapping matches**
`re.findall('ATG', dna)` skips overlapping matches. Use a lookahead: `re.findall('(?=(ATG))', dna)` to find all overlapping occurrences.

**4. IUPAC ambiguity codes are not regex**
`N` in a DNA sequence means any base, but `N` in a regex pattern matches the literal character N. Always translate IUPAC codes to character classes (`[ATGC]`) before matching.

**5. Anchors and `re.MULTILINE`**
`^>` anchors to the start of the entire string by default. To match at the start of each line in a multi-line FASTA string, compile with `re.MULTILINE`.

```python
# re.search() -- find first match anywhere in the string
dna = "GCTATGCGATCGATCGTAA"

match = re.search('ATG', dna)
if match:
    print(f"Found '{match.group()}' at position {match.start()}-{match.end()}")
```python

```python
# re.match() -- matches only at the BEGINNING of the string
dna = "ATGCGATCGATCG"

print(f"match('ATG', ...): {re.match('ATG', dna)}")   # Match object (starts at 0)
print(f"match('CGA', ...): {re.match('CGA', dna)}")   # None (CGA not at start)
```python

```python
# re.findall() -- return all non-overlapping matches as a list
dna = "ATGCGATCGATGCGATGCTAA"

starts = re.findall('ATG', dna)
print(f"Start codons found: {starts}")
print(f"Count: {len(starts)}")
```python

```python
# re.finditer() -- returns an iterator of match objects (preferred for positions)
dna = "ATGCGATCGATGCGATGCTAA"

print("All ATG positions:")
for match in re.finditer('ATG', dna):
    print(f"  Position {match.start()}: {match.group()}")
```python

```python
# Match objects carry useful information
dna = "GCTATGAAACGATCGTAA"
match = re.search('ATG', dna)

print(f"Matched text:  {match.group()}")
print(f"Start position: {match.start()}")
print(f"End position:   {match.end()}")
print(f"Span:           {match.span()}")
```python

---
## 2. Metacharacters and Character Classes

Special characters that have meaning in regex patterns:

| Symbol | Meaning |
|---|---|
| `.` | Any single character (except newline by default) |
| `^` | Start of string |
| `$` | End of string |
| `[ABC]` | Any one character from the set |
| `[^ABC]` | Any character NOT in the set |
| `[A-Z]` | Any character in the range |
| `\d` | Any digit (0-9) |
| `\w` | Any word character (letter, digit, underscore) |
| `\s` | Any whitespace |
| `\b` | Word boundary |
| `\|` | Alternation (OR) |

```python
# . (dot) -- matches any single character
dna = "ATGCGATCG"

# Find A followed by any character, then G
print(f"Pattern A.G matches: {re.findall('A.G', dna)}")
```python

```python
# ^ and $ -- anchors for start and end of string
dna = "ATGCGATCGTAA"

print(f"Starts with ATG: {bool(re.match('^ATG', dna))}")
print(f"Ends with TAA:   {bool(re.search('TAA$', dna))}")
print(f"Ends with ATG:   {bool(re.search('ATG$', dna))}")
```python

```python
# [] -- character classes
dna = "ATGCGATCGATCG"

# Purines (A or G)
purines = re.findall('[AG]', dna)
print(f"Purines: {purines} (count: {len(purines)})")

# Pyrimidines (C or T)
pyrimidines = re.findall('[CT]', dna)
print(f"Pyrimidines: {pyrimidines} (count: {len(pyrimidines)})")
```python

```python
# [^] -- negated character class
dna = "ATGCGATNNCGATXCG"

# Find non-standard nucleotides (not A, T, G, C)
non_standard = re.findall('[^ATGC]', dna)
print(f"Non-standard characters: {non_standard}")

# Find their positions
for m in re.finditer('[^ATGC]', dna):
    print(f"  Position {m.start()}: '{m.group()}'")
```python

```python
# | (pipe) -- alternation (OR)
dna = "ATGCGATCGATCGTAGTGATAA"

# Find all stop codons
stop_codons = re.findall('TAA|TAG|TGA', dna)
print(f"Stop codons found: {stop_codons}")
```python

---
## 3. Quantifiers

| Quantifier | Meaning |
|---|---|
| `*` | Zero or more |
| `+` | One or more |
| `?` | Zero or one |
| `{n}` | Exactly n times |
| `{n,m}` | Between n and m times |
| `{n,}` | n or more times |

By default, quantifiers are **greedy** (match as much as possible). Append `?` to make them **lazy** (match as little as possible).

```python
dna = "AAATTTAAAAAATTTTTGCGCGC"

print(f"A+ (one or more A): {re.findall('A+', dna)}")
print(f"T+ (one or more T): {re.findall('T+', dna)}")
print(f"[GC]+ (GC runs):    {re.findall('[GC]+', dna)}")
```python

```python
# Exact repetition counts
dna = "ATGATGATGATGATG"

print(f"Triplets (.{{3}}): {re.findall('.{3}', dna)}")
print(f"4-6 chars ([ATGC]{{4,6}}): {re.findall('[ATGC]{4,6}', dna)}")
```python

```python
# Greedy vs lazy matching -- critical for biological pattern search
dna = "ATGCGATCGTAAATGCCCTAG"

# Greedy: matches the LONGEST possible string
greedy = re.findall('ATG.*TAG', dna)
print(f"Greedy  (ATG.*TAG):  {greedy}")

# Lazy: matches the SHORTEST possible string
lazy = re.findall('ATG.*?TAG', dna)
print(f"Lazy    (ATG.*?TAG): {lazy}")
```python

```python
# Practical example: find homopolymer runs of 4+ bases
dna = "ATGCAAAAGCCCCCGATTTTTCG"

print(f"Sequence: {dna}")
print("\nHomopolymer runs (4+ bases):")
for m in re.finditer(r'(.)\1{3,}', dna):
    print(f"  '{m.group()}' at position {m.start()} ({len(m.group())} bases)")
```python

---
## 4. Groups

Parentheses `()` create **capturing groups**. They allow you to extract specific parts of a match.

| Syntax | Meaning |
|---|---|
| `(pattern)` | Capturing group |
| `(?:pattern)` | Non-capturing group |
| `(?P<name>pattern)` | Named group |
| `\1`, `\2` | Backreference to group 1, 2, etc. |

```python
# Capturing groups -- extract parts of a match
header = ">gene1|BRCA1|Homo sapiens breast cancer gene"

match = re.search(r'>([^|]+)\|([^|]+)\|(.+)', header)
if match:
    print(f"Full match:  {match.group(0)}")
    print(f"Gene ID:     {match.group(1)}")
    print(f"Gene name:   {match.group(2)}")
    print(f"Description: {match.group(3)}")
```python

```python
# Named groups -- (?P<name>pattern)
header = ">sp|P04637|P53_HUMAN Cellular tumor antigen p53 OS=Homo sapiens"

pattern = r'>sp\|(?P<accession>[^|]+)\|(?P<entry_name>\S+)\s+(?P<description>.+?)\s+OS=(?P<organism>.+)'
match = re.search(pattern, header)

if match:
    print(f"Accession:   {match.group('accession')}")
    print(f"Entry name:  {match.group('entry_name')}")
    print(f"Description: {match.group('description')}")
    print(f"Organism:    {match.group('organism')}")
```python

```python
# Non-capturing groups -- (?:pattern)
# Useful for alternation without capturing
dna = "ATGCGATCGATCGTAGTGATAA"

# With capturing group: findall returns the group, not the full match
with_capture = re.findall('(TAA|TAG|TGA)', dna)
print(f"Capturing group:     {with_capture}")

# With non-capturing group: findall returns the full match
without_capture = re.findall('(?:TAA|TAG|TGA)', dna)
print(f"Non-capturing group: {without_capture}")
```python

```python
# Backreferences -- find repeated codons (tandem repeats)
dna = "ATGATGATGCATGCATGCATGC"

# Find triplets that repeat at least twice
repeated = re.findall(r'([ATGC]{3})\1+', dna)
print(f"Repeated codons: {repeated}")

# With full context
for m in re.finditer(r'([ATGC]{3})\1+', dna):
    unit = m.group(1)
    copies = len(m.group(0)) // len(unit)
    print(f"  '{unit}' x {copies} at position {m.start()}")
```python

---
## 5. Lookahead and Lookbehind

Lookaround assertions match a position without consuming characters.

| Syntax | Name | Meaning |
|---|---|---|
| `(?=pattern)` | Positive lookahead | Followed by pattern |
| `(?!pattern)` | Negative lookahead | NOT followed by pattern |
| `(?<=pattern)` | Positive lookbehind | Preceded by pattern |
| `(?<!pattern)` | Negative lookbehind | NOT preceded by pattern |

```python
# Positive lookahead: find ATG only when followed by at least 30 bases before a stop
dna = "ATGCGATCGATCGATCGATCGATCGATCGATCGTAA"

# Find 'G' only when followed by 'C' (G in GC dinucleotide)
gc_dinucs = re.findall('G(?=C)', dna)
print(f"G followed by C: {len(gc_dinucs)} occurrences")

# Find all overlapping 3-mers starting with AT using lookahead
overlapping = re.findall('(?=(AT.))', dna)
print(f"Overlapping 3-mers starting with AT: {overlapping}")
```python

```python
# Negative lookahead: find a codon that is NOT a stop codon
dna = "ATGAAATTTCCCGGG"
codons = re.findall('.{3}', dna)
print(f"All codons: {codons}")

# Find positions where ATG is NOT followed by a stop codon
non_stop_starts = re.findall('ATG(?!TAA|TAG|TGA)', dna)
print(f"ATG not followed by stop: {non_stop_starts}")
```python

```python
# Lookbehind: find bases that come after a specific context
dna = "GAATTCGGATCCAAGCTT"

# Find 'C' only when preceded by 'AT' (ATC codons)
after_at = re.findall('(?<=AT)C', dna)
print(f"C preceded by AT: {len(after_at)} occurrences")

# Find bases after EcoRI site (GAATTC)
after_ecori = re.findall('(?<=GAATTC).{3}', dna)
print(f"3 bases after EcoRI site: {after_ecori}")
```python

---
## 6. Substitution and Splitting

```python
# re.sub() -- find and replace
dna = "ATGCGATCGATCG"

# Transcription: T -> U
rna = re.sub('T', 'U', dna)
print(f"DNA: {dna}")
print(f"RNA: {rna}")
```python

```python
# Using groups in substitution
# Reformat FASTA headers: add a prefix to gene IDs
fasta = ">gene1 description one\nATGCGATC\n>gene2 description two\nGCTAGCTA"

modified = re.sub(r'>(\w+)', r'>Hsapiens_\1', fasta)
print(modified)
```python

```python
# Substitution with a function
# Mask long homopolymer runs in a sequence
def mask_repeat(match):
    """Replace long homopolymer runs with N's."""
    seq = match.group(0)
    if len(seq) >= 5:
        return 'N' * len(seq)
    return seq


dna = "ATGAAAAAACGATTTTTTTCGATCG"
masked = re.sub(r'(.)\1{2,}', mask_repeat, dna)
print(f"Original: {dna}")
print(f"Masked:   {masked}")
```python

```python
# re.split() -- split by pattern
# Split a sequence at restriction sites
dna = "ATGCGAATTCGATCGGATCCATCGAAGCTTATCG"

# Split at any common restriction site
fragments = re.split('GAATTC|GGATCC|AAGCTT', dna)
print(f"Fragments after triple digestion:")
for i, frag in enumerate(fragments, 1):
    print(f"  Fragment {i}: {frag} ({len(frag)} bp)")
```python

---
## 7. Compiled Patterns and Flags

Pre-compile patterns with `re.compile()` for efficiency when reusing the same pattern.

```python
# Compile a pattern for repeated use
start_codon = re.compile('ATG')
stop_pattern = re.compile('TAA|TAG|TGA')

dna = "ATGCGATCGATGCGTAAATGCCCTAG"

print(f"Start codons: {start_codon.findall(dna)}")
print(f"Stop codons:  {stop_pattern.findall(dna)}")

# Compiled patterns have the same methods: search, match, findall, finditer, sub, split
for m in start_codon.finditer(dna):
    print(f"  ATG at position {m.start()}")
```python

```python
# re.IGNORECASE (re.I) -- case-insensitive matching
dna_mixed = "atgCGAtcgATGccc"

# Without flag: misses lowercase ATG
print(f"Case sensitive:   {re.findall('ATG', dna_mixed)}")

# With flag: finds both
print(f"Case insensitive: {re.findall('ATG', dna_mixed, re.IGNORECASE)}")
```python

```python
# re.MULTILINE (re.M) -- ^ and $ match at line boundaries
fasta = """>gene1
ATGCGATCG
>gene2
GCTAGCTAG"""

# Without MULTILINE: ^ only matches start of entire string
print(f"Without MULTILINE: {re.findall('^>.*', fasta)}")

# With MULTILINE: ^ matches start of each line
print(f"With MULTILINE:    {re.findall('^>.*', fasta, re.MULTILINE)}")
```python

## Common Pitfalls

- **Mutable default arguments**: Never use `def f(x=[])` — use `def f(x=None)` and set inside the function
- **Off-by-one errors**: Python ranges are half-open `[start, stop)` — bioinformatics coordinates are often 1-based
- **Deep vs shallow copy**: Nested data structures require `copy.deepcopy()` — `list.copy()` only copies the top level
