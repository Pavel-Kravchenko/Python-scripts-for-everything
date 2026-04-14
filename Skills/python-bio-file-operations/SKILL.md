---
name: python-bio-file-operations
description: FASTA/FASTQ/CSV parsing recipes, streaming file I/O patterns, and format-specific pitfalls for bioinformatics data files
tool_type: python
primary_tool: Python
---

# File Operations for Bioinformatics

## Pitfalls

- **`strip()` when reading lines:** every line from a file has a trailing `\n`. Forgetting this causes sequences to have invisible newlines that break comparisons.
- **`"w"` mode destroys existing content.** Use `"a"` to append. Use `"x"` to fail-safe (error if file exists).
- **Binary mode for BAM/compressed files:** `open(path, "rb")`. Text mode decodes bytes and corrupts binary data.
- **`f.read()` loads entire file.** For large FASTA/FASTQ, always iterate line-by-line with `for line in f:`.

## Reading Methods

| Method | Memory |
|--------|--------|
| `f.read()` | Loads all |
| `f.readline()` | One line |
| `f.readlines()` | Loads all |
| `for line in f:` | One line (preferred) |

## FASTA Parser (streaming)

```python
def read_fasta(filename):
    """Parse FASTA file -> dict of {header: sequence}. Streams line-by-line."""
    sequences = {}
    current_header = None
    current_seq = []
    with open(filename) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if current_header is not None:
                    sequences[current_header] = ''.join(current_seq)
                current_header = line[1:]
                current_seq = []
            else:
                current_seq.append(line)
        if current_header is not None:           # don't forget the last sequence!
            sequences[current_header] = ''.join(current_seq)
    return sequences
```

## FASTA Writer (with line wrapping)

```python
def write_fasta(sequences, filename, line_width=60):
    """Write sequences (dict or list of tuples) to FASTA with wrapped lines."""
    items = sequences.items() if isinstance(sequences, dict) else sequences
    with open(filename, 'w') as f:
        for header, seq in items:
            f.write(f">{header}\n")
            for i in range(0, len(seq), line_width):
                f.write(seq[i:i+line_width] + '\n')
```

## CSV/TSV with `csv` module

```python
import csv

# Read CSV with column-name access
with open('gene_expression.csv') as f:
    for row in csv.DictReader(f):
        print(row['gene_name'], float(row['expression']))

# Read TSV (BED files, etc.)
with open('genes.bed') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        ...

# Write with header
with open('results.csv', 'w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=['gene', 'fold_change', 'p_value'])
    writer.writeheader()
    writer.writerows(results)
```
