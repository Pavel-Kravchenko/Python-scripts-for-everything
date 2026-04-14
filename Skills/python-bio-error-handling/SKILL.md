---
name: python-bio-error-handling
description: Robust error handling for bioinformatics pipelines — custom exception hierarchy, try/except/else/finally, exception chaining, strict FASTA parsing.
tool_type: python
primary_tool: Python
---

# Error Handling for Bioinformatics

## Custom Exception Hierarchy

```python
class BioinformaticsError(Exception):
    """Base for all pipeline errors."""

class InvalidSequenceError(BioinformaticsError):
    def __init__(self, invalid_chars, seq_type="DNA"):
        self.invalid_chars = invalid_chars
        super().__init__(f"Invalid {seq_type} characters: {invalid_chars}")

class SequenceLengthError(BioinformaticsError):
    def __init__(self, actual, minimum=None, maximum=None):
        if minimum and actual < minimum:
            msg = f"Sequence too short: {actual} < {minimum}"
        elif maximum and actual > maximum:
            msg = f"Sequence too long: {actual} > {maximum}"
        else:
            msg = f"Invalid sequence length: {actual}"
        super().__init__(msg)

class FastaParseError(BioinformaticsError):
    def __init__(self, message, filename=None, line_number=None, line_content=None):
        parts = [message]
        if filename:    parts.append(f"file={filename}")
        if line_number: parts.append(f"line={line_number}")
        if line_content: parts.append(f"content='{line_content[:50]}'")
        super().__init__(" | ".join(parts))

class TranslationError(BioinformaticsError): pass
```

## try/except/else/finally Pattern

```python
def read_fasta(filename):
    sequences = {}
    try:
        f = open(filename, 'r')
    except FileNotFoundError:
        raise  # let caller handle missing file
    except PermissionError as e:
        raise FastaParseError(f"Cannot read '{filename}'") from e
    else:
        # Only runs if open() succeeded
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                current_id = line[1:].split()[0]
                sequences[current_id] = []
            elif current_id:
                sequences[current_id].append(line)
        sequences = {k: ''.join(v) for k, v in sequences.items()}
    finally:
        f.close()  # always runs
    return sequences
```

## Exception Chaining

```python
def parse_gff_line(line, line_number=None):
    fields = line.strip().split('\t')
    if len(fields) != 9:
        raise FastaParseError(
            f"Expected 9 fields, got {len(fields)}",
            line_number=line_number, line_content=line.strip()
        )
    try:
        return {'start': int(fields[3]), 'end': int(fields[4]), ...}
    except ValueError as e:
        raise FastaParseError(
            "Cannot parse coordinates",
            line_number=line_number, line_content=line.strip()
        ) from e   # preserves original traceback
```

## Strict FASTA Parser

```python
def strict_fasta_parser(filename: str) -> dict[str, str]:
    sequences: dict[str, str] = {}
    current_id: str | None = None
    current_seq: list[str] = []
    VALID = set('ATGCNRYSWKMBDHV.-')

    with open(filename, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if current_id is not None:
                    seq = ''.join(current_seq)
                    if not seq:
                        raise FastaParseError("Empty sequence", filename, line_num, f">{current_id}")
                    sequences[current_id] = seq
                if len(line) == 1:
                    raise FastaParseError("Empty sequence ID", filename, line_num, line)
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                if current_id is None:
                    raise FastaParseError("Sequence data before first header", filename, line_num, line)
                invalid = set(line.upper()) - VALID
                if invalid:
                    raise FastaParseError(f"Invalid characters: {invalid}", filename, line_num, line)
                current_seq.append(line.upper())

    if current_id is not None:
        seq = ''.join(current_seq)
        if not seq:
            raise FastaParseError("Empty sequence for last record", filename)
        sequences[current_id] = seq

    return sequences
```

## Batch Processing with Graceful Degradation

```python
def batch_gc_analysis(sequences: dict[str, str], strict: bool = False):
    results, errors = {}, []
    for seq_id, seq in sequences.items():
        try:
            if not seq:
                raise SequenceLengthError(0, minimum=1)
            invalid = set(seq.upper()) - set('ATGCN')
            if invalid:
                raise InvalidSequenceError(invalid)
            s = seq.upper()
            results[seq_id] = round((s.count('G') + s.count('C')) / len(s) * 100, 2)
        except BioinformaticsError as e:
            if strict:
                raise
            errors.append((seq_id, str(e)))
    return results, errors
```

## Pitfalls

- **Bare `except:`**: catches `KeyboardInterrupt` and `SystemExit` — impossible to stop the script with Ctrl+C. Always name the exception.
- **`assert` for input validation**: disabled with `python -O`. Use `if not condition: raise ValueError(...)` instead.
- **`raise NewError(...) from original`**: omitting `from original` hides the root cause in tracebacks.
- **`else` on `try`**: runs only when *no* exception occurred — separates "might fail" from "should only run on success". Commonly forgotten.
- **`finally` runs even on `return`**: useful for cleanup (file close, temp removal), but mutating state in `finally` can mask real errors.
- **Catching `BioinformaticsError` too broadly**: catch the most specific subclass first; put the base class catch last as a fallback.
