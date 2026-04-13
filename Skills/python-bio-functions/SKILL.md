---
name: python-bio-functions
description: "By the end of this module, you will be able to: - Define and call functions with various parameter types - Use default arguments, `*args`, and `**kwargs` - Write lambda functions and use `map`/`filter"
tool_type: python
source_notebook: "Tier_1_Python_for_Bioinformatics/06_Functions/01_functions.ipynb"
---

# Module 6: Functions in Python

*Source: Course notebook `Tier_1_Python_for_Bioinformatics/06_Functions/01_functions.ipynb`*

# Module 6: Functions in Python

## From Basic Definitions to Functional Programming

---

### Learning Objectives

By the end of this module, you will be able to:
- Define and call functions with various parameter types
- Use default arguments, `*args`, and `**kwargs`
- Write lambda functions and use `map`/`filter`/`reduce`
- Understand variable scope (LEGB rule)
- Implement recursive algorithms
- Write proper docstrings and type hints
- Build reusable bioinformatics utility functions

---

**Why functions matter in bioinformatics:** Every analysis pipeline is built from functions -- GC content calculators, sequence translators, file parsers, statistical tests. Writing clean, reusable functions is what separates a quick script from a maintainable tool.

## How to use this notebook

1. Run cells in order. Sections 1–4 cover essentials (definition, parameters, defaults, docstrings). Sections 5–9 cover advanced topics. Sections 10–11 provide complete examples and exercises.
2. Each exercise has a blank cell (write your solution) followed by a solution cell. **Do the blank one first** — the exercises test real bioinformatics tasks.
3. The scope section (Section 7) is conceptually important: understand it before writing functions that modify global state.

## Common stumbling points

- **`return` vs `print()`:** A function that uses `print()` sends output to the screen but returns `None`. If you want to use the result elsewhere, use `return`. Beginners often call `result = my_function()` and get `None` because they forgot `return`.
- **Mutable default arguments are a trap:** `def f(x, data=[])` is a classic Python bug. The default list is created *once* when the function is defined, not each time it is called. On the second call, the list from the first call is still there. Always use `None` as default for mutable arguments: `def f(x, data=None): if data is None: data = []`.
- **Variable scope (LEGB):** A variable inside a function is local — it does not exist outside the function. A function can *read* a global variable but cannot reassign it without `global`. This is intentional: it prevents functions from silently breaking each other.
- **`*args` and `**kwargs` order matters:** The order must be: required args, default args, `*args`, keyword-only args, `**kwargs`. Getting this wrong causes a `SyntaxError`.
- **Lambda functions cannot contain statements:** A `lambda` can only contain a single *expression* — no `if` blocks, no `for` loops, no assignment. Use `def` for anything more complex.

---

[← Previous: Module 5: Control Flow](../05_Control_Flow/01_control_flow.ipynb) | [Next: Module 7: File Operations →](../07_File_Operations/01_file_operations.ipynb)

---

## 1. Function Basics

A function is a reusable block of code that you can call from anywhere in your program.

```python
def function_name(parameter1, parameter2):
    """Docstring: describes what the function does."""
    # function body
    return result
```

**Why use functions?**
- Write code once, call it many times (DRY: Don't Repeat Yourself)
- Give a name to a block of logic — easier to read
- Test in isolation
- Share with others (as a module)

```python
# Simplest function: no parameters, no return value
def greet():
    print("Hello, Bioinformatician!")

greet()
greet()  # call it as many times as you like
```

```python
# Function with a parameter and a return value
def calculate_gc_content(sequence):
    """Calculate GC content of a DNA sequence.
    
    Args:
        sequence: DNA sequence string (ATGC characters)
    
    Returns:
        GC content as a percentage (float)
    """
    sequence = sequence.upper()
    gc_count = sequence.count('G') + sequence.count('C')
    return (gc_count / len(sequence)) * 100

# Call the function
seq = "ATGCGATCGATCGTAGC"
gc = calculate_gc_content(seq)
print(f"GC content of {seq}: {gc:.2f}%")
```

```python
# Functions can return multiple values as a tuple
def analyze_nucleotides(sequence):
    """Count all nucleotides in a DNA sequence."""
    seq = sequence.upper()
    return seq.count('A'), seq.count('T'), seq.count('G'), seq.count('C')

# Unpack the returned tuple
a, t, g, c = analyze_nucleotides("ATGCGATCGATCGATCG")
print(f"A={a}, T={t}, G={g}, C={c}")
print(f"Total: {a + t + g + c} nucleotides")
```

```python
# Returning a dictionary for richer results
def analyze_sequence(sequence):
    """Comprehensive sequence analysis.
    
    Returns:
        dict with length, gc_content, and nucleotide counts
    """
    seq = sequence.upper()
    length = len(seq)
    gc = (seq.count('G') + seq.count('C')) / length * 100 if length > 0 else 0
    counts = {nuc: seq.count(nuc) for nuc in 'ATGC'}
    
    return {'length': length, 'gc_content': round(gc, 2), 'counts': counts}

result = analyze_sequence("ATGCGATCGATCGTAGCGATCGATCG")
print(f"Length: {result['length']} bp")
print(f"GC content: {result['gc_content']}%")
print(f"Nucleotide counts: {result['counts']}")
```

---
## 2. Default Arguments

Parameters can have default values. Required arguments must come before defaults.

```
Positional:  def func(a, b)       -> func(1, 2)
Keyword:     def func(a, b)       -> func(a=1, b=2)
Default:     def func(a, b=10)    -> func(1)  # b defaults to 10
```

```python
import random

def make_random_sequence(length=200, alphabet="ATGC"):
    """Generate a random nucleotide sequence.
    
    Args:
        length: Number of nucleotides (default: 200)
        alphabet: Characters to choose from (default: 'ATGC' for DNA)
    """
    return ''.join(random.choice(alphabet) for _ in range(length))

# Use defaults
dna = make_random_sequence()
print(f"Random DNA (200 bp): {dna[:50]}...")

# Override length
short_dna = make_random_sequence(length=20)
print(f"Short DNA (20 bp): {short_dna}")

# Generate RNA instead
rna = make_random_sequence(length=20, alphabet="AUGC")
print(f"Random RNA (20 bp): {rna}")
```

```python
# Codon translation with default table
def translate_codon(codon, table='standard'):
    """Translate a DNA codon to its amino acid."""
    standard_table = {
        'ATG': 'M', 'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAA': '*', 'TAG': '*', 'TGA': '*',
        'GGG': 'G', 'GGA': 'G', 'GGC': 'G', 'GGT': 'G',
        'AAA': 'K', 'AAG': 'K', 'GCG': 'A', 'TGG': 'W',
    }
    return standard_table.get(codon.upper(), 'X')  # X for unknown

print(f"ATG -> {translate_codon('ATG')}  (Methionine / Start)")
print(f"TAA -> {translate_codon('TAA')}  (Stop)")
print(f"GGC -> {translate_codon('GGC')}  (Glycine)")
```

```python
# WARNING: Mutable default arguments are a classic Python trap!

# BAD -- the default list is created ONCE and shared across calls
def bad_collect(item, items=[]):
    items.append(item)
    return items

print(bad_collect('A'))  # ['A']
print(bad_collect('T'))  # ['A', 'T'] -- Surprise! The list persists!

print()

# GOOD -- use None and create a fresh list each call
def good_collect(item, items=None):
    if items is None:
        items = []
    items.append(item)
    return items

print(good_collect('A'))  # ['A']
print(good_collect('T'))  # ['T'] -- Correct!
```

---
## 3. `*args` and `**kwargs`

Handle variable numbers of arguments:
- `*args` collects extra **positional** arguments into a **tuple**
- `**kwargs` collects extra **keyword** arguments into a **dictionary**

```python
# *args: concatenate any number of sequences
def concatenate_sequences(*sequences):
    """Join any number of DNA/RNA sequences together."""
    return ''.join(sequences)

print(concatenate_sequences('ATG', 'CCC', 'GGG'))  # 'ATGCCCGGG'
print(concatenate_sequences('ATG'))                 # 'ATG'
```

```python
# **kwargs: build a FASTA header with flexible metadata
def create_fasta_header(sequence_id, **metadata):
    """Create a FASTA header line from an ID and arbitrary metadata."""
    header = f">{sequence_id}"
    for key, value in metadata.items():
        header += f" [{key}={value}]"
    return header

header = create_fasta_header(
    "BRCA1_human",
    organism="Homo sapiens",
    gene="BRCA1",
    length=5500,
    chromosome=17
)
print(header)
```

```python
# Combining all argument types
# Order: required, default, *args, keyword-only, **kwargs

def sequence_stats(*sequences, normalize=False, **options):
    """Calculate statistics for multiple sequences.
    
    Args:
        *sequences: Variable number of DNA sequences
        normalize: If True, report nucleotide fractions instead of counts
        **options: Additional options (e.g., include_gc=True)
    """
    include_gc = options.get('include_gc', True)
    results = []
    
    for seq in sequences:
        seq = seq.upper()
        stats = {'length': len(seq)}
        for nuc in 'ATGC':
            count = seq.count(nuc)
            stats[nuc] = count / len(seq) if normalize else count
        if include_gc:
            stats['gc_content'] = (seq.count('G') + seq.count('C')) / len(seq) * 100
        results.append(stats)
    
    return results

stats = sequence_stats('ATGCATGC', 'GGGGCCCC', normalize=True, include_gc=True)
for i, s in enumerate(stats):
    print(f"Sequence {i+1}: {s}")
```

---
## 4. Docstrings and Type Hints

Good documentation makes your functions usable by others (and by future you).

```python
def validate_dna_sequence(sequence: str, allow_n: bool = False) -> bool:
    """Check whether a string is a valid DNA sequence.
    
    Args:
        sequence: The string to validate.
        allow_n: If True, allow 'N' as an ambiguous nucleotide.
    
    Returns:
        True if the sequence contains only valid DNA characters.
    
    Examples:
        >>> validate_dna_sequence('ATGC')
        True
        >>> validate_dna_sequence('ATXG')
        False
        >>> validate_dna_sequence('ATNG', allow_n=True)
        True
    """
    valid_chars = set('ATGCN') if allow_n else set('ATGC')
    return set(sequence.upper()) <= valid_chars

# Test
print(validate_dna_sequence('ATGCGATCG'))          # True
print(validate_dna_sequence('ATGXYZ'))              # False
print(validate_dna_sequence('ATGNNNC', allow_n=True))  # True
```

```python
def find_orfs(sequence: str, min_length: int = 30) -> list[dict]:
    """Find all Open Reading Frames (ORFs) in a DNA sequence.
    
    An ORF starts with ATG and ends with a stop codon (TAA, TAG, TGA).
    
    Args:
        sequence: DNA sequence to search.
        min_length: Minimum ORF length in nucleotides (default: 30).
    
    Returns:
        List of dicts with 'start', 'end', 'length', and 'sequence' keys.
    """
    sequence = sequence.upper()
    stop_codons = {'TAA', 'TAG', 'TGA'}
    orfs = []
    
    for frame in range(3):  # Three reading frames
        i = frame
        while i < len(sequence) - 2:
            codon = sequence[i:i+3]
            if codon == 'ATG':
                # Found a start codon -- scan for stop
                for j in range(i + 3, len(sequence) - 2, 3):
                    next_codon = sequence[j:j+3]
                    if next_codon in stop_codons:
                        orf_seq = sequence[i:j+3]
                        if len(orf_seq) >= min_length:
                            orfs.append({
                                'start': i,
                                'end': j + 3,
                                'length': len(orf_seq),
                                'sequence': orf_seq
                            })
                        i = j + 3  # Continue after this stop codon
                        break
                else:
                    i += 3
                    continue
                continue
            i += 3
    
    return sorted(orfs, key=lambda x: x['start'])

# Test: a sequence containing two ORFs
test_seq = "CCCATGGCCGATCGATAGCCCATGAAAGGGCCCTTTAAATTT"
orfs = find_orfs(test_seq, min_length=9)
for orf in orfs:
    print(f"ORF at {orf['start']}-{orf['end']}: {orf['sequence']} ({orf['length']} bp)")
```

---
## 5. Lambda Functions

Anonymous, single-expression functions. Useful as short callbacks for sorting/filtering.

```python
lambda arguments: expression
```

**Important:** Lambdas provide no speed advantage over regular functions. They are purely a convenience for short, throwaway operations.
