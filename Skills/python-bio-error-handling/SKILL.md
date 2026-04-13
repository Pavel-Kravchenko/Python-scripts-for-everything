---
name: python-bio-error-handling
description: "Bioinformatics data is notoriously messy: FASTA files with empty sequences, VCF lines with the wrong number of fields, gene IDs that don't match between databases, quality scores outside valid range. "
tool_type: python
source_notebook: "Tier_1_Python_for_Bioinformatics/15_Error_Handling/01_error_handling.ipynb"
primary_tool: Python
---

## Version Compatibility

Reference examples tested with: Python 3.10+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.


# Module 15: Error Handling and Debugging

*Source: Course notebook `Tier_1_Python_for_Bioinformatics/15_Error_Handling/01_error_handling.ipynb`*

# Module 15: Error Handling and Debugging

## Writing Robust Bioinformatics Code

---

### Learning Objectives
- Master `try`/`except`/`else`/`finally` for structured error handling
- Raise exceptions and create custom exception classes
- Understand the Python exception hierarchy
- Debug effectively with print debugging, `assert`, and `pdb`/`breakpoint()`
- Use the `logging` module for production code
- Handle real-world bioinformatics errors: malformed FASTA, invalid sequences, missing GFF data

---

## How to use this notebook

Run cells from top to bottom on the first pass — later cells depend on functions and variables defined earlier. Once you have run through the notebook, feel free to modify parameters and re-run individual sections.

Each section has runnable examples first, followed by exercises. Try the exercise before looking at the solution cell below it.

## Complicated moments explained

**1. Always catch specific exceptions**
`except:` catches *everything* including `KeyboardInterrupt` and `SystemExit`, making it impossible to stop the script with Ctrl+C. Always name the exception: `except ValueError`, `except FileNotFoundError`.

**2. `else` on a `try` block**
The `else` clause runs only if the `try` block raised *no* exception. Use it to separate "code that might fail" from "code that should only run on success":
```python
try:
    seq = parse_fasta(file)
except ParseError:
    log.warning(f"bad file: {file}")
else:
    analyze(seq)   # only runs if parse succeeded
```

**3. `raise ... from ...` preserves the chain**
When catching one exception and raising another, use `raise NewError(...) from original` to keep the full traceback chain. Without `from`, the original cause is hidden.

**4. `assert` is for development invariants only**
Python disables `assert` statements when run with `python -O`. Never use `assert` to validate external data. Use `if not condition: raise ValueError(...)` instead.

```python
def gc_content(seq):
    """Calculate GC content with basic error handling."""
    try:
        seq = seq.upper()
        gc = seq.count('G') + seq.count('C')
        return gc / len(seq) * 100
    except ZeroDivisionError:
        print("Warning: empty sequence")
        return 0.0
    except AttributeError:
        print(f"Warning: expected string, got {type(seq).__name__}")
        return None


print(f"Normal: {gc_content('ATGCGC'):.1f}%")
print(f"Empty:  {gc_content('')}%")
print(f"Wrong type: {gc_content(12345)}")
```

**Rule:** Always catch specific exceptions. A bare `except:` hides bugs
by catching everything, including `KeyboardInterrupt` and `SystemExit`.

## 3. The Complete `try`/`except`/`else`/`finally` Block

- `try` -- code that might raise an exception
- `except` -- handle specific exceptions
- `else` -- runs only if **no** exception occurred (the "happy path")
- `finally` -- **always** runs (cleanup code)

```python
def read_fasta(filename):
    """Read a FASTA file with full error handling."""
    sequences = {}
    current_id = None
    file = None

    try:
        file = open(filename, 'r')

    except FileNotFoundError:
        print(f"Error: File '{filename}' not found")
        return {}

    except PermissionError:
        print(f"Error: No permission to read '{filename}'")
        return {}

    else:
        # Only runs if the file opened successfully
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                current_id = line[1:].split()[0]
                sequences[current_id] = []
            elif current_id:
                sequences[current_id].append(line)
        sequences = {k: ''.join(v) for k, v in sequences.items()}
        print(f"Read {len(sequences)} sequences from {filename}")

    finally:
        # Always close the file if it was opened
        if file:
            file.close()
            print("File closed")

    return sequences


# Create a test file
with open('test.fasta', 'w') as f:
    f.write(">seq1 test sequence\nATGCATGC\n>seq2 another\nGGCCAATT\n")

print("--- Existing file ---")
seqs = read_fasta('test.fasta')
print(f"Result: {seqs}\n")

print("--- Non-existent file ---")
seqs = read_fasta('does_not_exist.fasta')
```

## 4. Raising Exceptions

Use `raise` when your code detects an error condition. Provide a clear message
that helps the caller understand what went wrong.

```python
def validate_dna(sequence, min_length=3):
    """Validate a DNA sequence, raising clear errors on problems."""
    if not isinstance(sequence, str):
        raise TypeError(f"Expected string, got {type(sequence).__name__}")

    if len(sequence) < min_length:
        raise ValueError(
            f"Sequence too short: {len(sequence)} bp (minimum {min_length})"
        )

    invalid = set(sequence.upper()) - set('ATGCN')
    if invalid:
        raise ValueError(f"Invalid DNA characters: {invalid}")

    return sequence.upper()


test_cases = [
    ("ATGCATGC", "Valid DNA"),
    ("AT", "Too short"),
    ("ATGXYZ", "Invalid characters"),
    (42, "Wrong type"),
]

for seq, description in test_cases:
    try:
        result = validate_dna(seq)
        print(f"  {description}: OK -> {result}")
    except (TypeError, ValueError) as e:
        print(f"  {description}: FAILED -> {e}")
```

## 5. Custom Exceptions

Custom exceptions make your code's errors self-documenting and allow
callers to catch specific categories of errors.

```python
class BioinformaticsError(Exception):
    """Base exception for all bioinformatics errors."""
    pass


class InvalidSequenceError(BioinformaticsError):
    """Raised when a sequence contains invalid characters."""
    def __init__(self, invalid_chars, seq_type="DNA"):
        self.invalid_chars = invalid_chars
        self.seq_type = seq_type
        super().__init__(f"Invalid {seq_type} characters: {invalid_chars}")


class SequenceLengthError(BioinformaticsError):
    """Raised when a sequence does not meet length requirements."""
    def __init__(self, actual, minimum=None, maximum=None):
        self.actual = actual
        if minimum and actual < minimum:
            msg = f"Sequence too short: {actual} < {minimum}"
        elif maximum and actual > maximum:
            msg = f"Sequence too long: {actual} > {maximum}"
        else:
            msg = f"Invalid sequence length: {actual}"
        super().__init__(msg)


class FastaParseError(BioinformaticsError):
    """Raised on FASTA format violations, with context."""
    def __init__(self, message, filename=None, line_number=None, line_content=None):
        self.filename = filename
        self.line_number = line_number
        self.line_content = line_content
        parts = [message]
        if filename:
            parts.append(f"file={filename}")
        if line_number:
            parts.append(f"line={line_number}")
        if line_content:
            parts.append(f"content='{line_content[:50]}'")
        super().__init__(" | ".join(parts))


class TranslationError(BioinformaticsError):
    """Raised when DNA/RNA translation fails."""
    pass
```

```python
# Using custom exceptions
def validate_and_translate(sequence):
    """Validate DNA and translate to protein."""
    seq = sequence.upper()

    # Check for invalid characters
    invalid = set(seq) - set('ATGC')
    if invalid:
        raise InvalidSequenceError(invalid, seq_type="DNA")

    # Check length
    if len(seq) < 3:
        raise SequenceLengthError(len(seq), minimum=3)

    if len(seq) % 3 != 0:
        raise TranslationError(
            f"Sequence length {len(seq)} is not divisible by 3"
        )

    # Simple translation
    codon_table = {
        'ATG': 'M', 'GCT': 'A', 'GCC': 'A', 'TAA': '*', 'TAG': '*', 'TGA': '*',
        'TTT': 'F', 'TTC': 'F', 'TGG': 'W', 'GAA': 'E', 'GAG': 'E',
    }
    protein = []
    for i in range(0, len(seq), 3):
        codon = seq[i:i+3]
        aa = codon_table.get(codon, 'X')
        if aa == '*':
            break
        protein.append(aa)
    return ''.join(protein)


# Test various error conditions
test_cases = [
    "ATGGCTTAG",     # valid
    "ATGXYZ",        # invalid characters
    "AT",            # too short
    "ATGGC",         # not divisible by 3
]

for seq in test_cases:
    try:
        protein = validate_and_translate(seq)
        print(f"  '{seq}' -> '{protein}'")
    except InvalidSequenceError as e:
        print(f"  '{seq}' -> Invalid: {e}")
    except SequenceLengthError as e:
        print(f"  '{seq}' -> Length: {e}")
    except TranslationError as e:
        print(f"  '{seq}' -> Translation: {e}")
    except BioinformaticsError as e:
        # Catches any BioinformaticsError subclass not caught above
        print(f"  '{seq}' -> General bio error: {e}")
```

## 6. Exception Chaining

Use `raise ... from ...` to preserve the original error context.
This is essential for debugging multi-layered pipelines.

```python
def parse_gff_line(line, line_number=None):
    """Parse a single GFF3 line into a dict."""
    fields = line.strip().split('\t')
    if len(fields) != 9:
        raise FastaParseError(
            f"Expected 9 tab-separated fields, got {len(fields)}",
            line_number=line_number,
            line_content=line.strip()
        )

    try:
        return {
            'seqid': fields[0],
            'source': fields[1],
            'type': fields[2],
            'start': int(fields[3]),
            'end': int(fields[4]),
            'score': fields[5],
            'strand': fields[6],
            'phase': fields[7],
            'attributes': fields[8],
        }
    except ValueError as e:
        # Chain the original ValueError into our custom exception
        raise FastaParseError(
            f"Cannot parse coordinates",
            line_number=line_number,
            line_content=line.strip()
        ) from e


# Valid GFF line
good_line = "chr1\tRefSeq\tgene\t11874\t14409\t.\t+\t.\tID=gene1;Name=DDX11L1"
print(parse_gff_line(good_line))

# Bad coordinate
bad_line = "chr1\tRefSeq\tgene\tNOT_A_NUMBER\t14409\t.\t+\t.\tID=gene1"
try:
    parse_gff_line(bad_line, line_number=42)
except FastaParseError as e:
    print(f"\nParse error: {e}")
    if e.__cause__:
        print(f"  Caused by: {e.__cause__}")
```

## 7. Handling Malformed FASTA Files

A strict FASTA parser that gives useful error messages.

```python
def strict_fasta_parser(filename):
    """Parse FASTA with strict validation and detailed errors."""
    sequences = {}
    current_id = None
    current_seq = []

    with open(filename, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line:
                continue

            if line.startswith('>'):
                # Save previous record
                if current_id:
                    seq = ''.join(current_seq)
                    if not seq:
                        raise FastaParseError(
                            "Empty sequence",
                            filename, line_num, f">{current_id}"
                        )
                    sequences[current_id] = seq

                # Parse header
                if len(line) == 1:
                    raise FastaParseError(
                        "Empty sequence ID (header is just '>')",
                        filename, line_num, line
                    )
                current_id = line[1:].split()[0]
                current_seq = []

            else:
                if current_id is None:
                    raise FastaParseError(
                        "Sequence data before first header",
                        filename, line_num, line
                    )
                # Validate sequence characters
                invalid = set(line.upper()) - set('ATGCNRYSWKMBDHV.-')
                if invalid:
                    raise FastaParseError(
                        f"Invalid characters: {invalid}",
                        filename, line_num, line
                    )
                current_seq.append(line.upper())

    # Save last record
    if current_id:
        seq = ''.join(current_seq)
        if not seq:
            raise FastaParseError("Empty sequence for last record", filename)
        sequences[current_id] = seq

    return sequences


# Test with a malformed file
with open('bad.fasta', 'w') as f:
    f.write(">seq1\nATGC\n>\nGGCC\n")  # second header is empty

try:
    seqs = strict_fasta_parser('bad.fasta')
except FastaParseError as e:
    print(f"Parse error: {e}")

# Test with data before header
with open('bad2.fasta', 'w') as f:
    f.write("ATGCATGC\n>seq1\nGGCC\n")  # sequence before header

try:
    seqs = strict_fasta_parser('bad2.fasta')
except FastaParseError as e:
    print(f"Parse error: {e}")
```
