---
name: foundations-character-encodings
description: "Every bioinformatics pipeline starts with reading data from files. FASTA sequences, GenBank records, PDB coordinates, GFF annotations -- they all live as bytes on disk. When your Python script reads g"
tool_type: python
source_notebook: "Tier_0_Computational_Foundations/04_File_Encodings/01_character_encodings.ipynb"
primary_tool: Python
---

## Version Compatibility

Reference examples tested with: Python 3.10+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.


# Character Encodings and Binary Data in Bioinformatics

*Source: Course notebook `Tier_0_Computational_Foundations/04_File_Encodings/01_character_encodings.ipynb`*


**Tier 0 -- Computational Foundations | Module 4**

---

## Why This Matters

Every bioinformatics pipeline starts with reading data from files. FASTA sequences, GenBank records, PDB coordinates, GFF annotations -- they all live as bytes on disk. When your Python script reads garbled text, when special characters in gene names vanish, or when a collaborator's CSV from Windows looks wrong on your Linux server, the root cause is almost always an **encoding mismatch**.

This module gives you a deep, practical understanding of character encodings so you can confidently handle any file you encounter.

### Learning Objectives

By the end of this module you will be able to:
- Explain how text is stored as bytes (bits, bytes, hexadecimal)
- Describe ASCII and why it matters for biological sequence data
- Contrast UTF-8, Latin-1, and legacy Cyrillic encodings
- Diagnose and fix common encoding problems in bioinformatics files
- Write robust Python code that handles encodings correctly
- Deal with special characters in gene names, FASTA headers, and database IDs

## How to use this notebook
1. Run all cells in order — later sections use the `hexdump` function defined earlier.
2. The hex inspection cells are the core skill: when you see unexpected bytes, you need to read them in hex to understand what is happening.
3. Try modifying the examples: change the encoding names, see what different bytes decode to.
4. The exercises cover the most common real-world encoding problems you will encounter with biological data.

## Common stumbling points

- **Latin-1 never fails — but that does not mean it is right**: Latin-1 decodes every byte 0–255 successfully, which is why it is used as a "safety" fallback. But it will silently produce wrong characters when the file is actually UTF-8 with multi-byte sequences.
- **UTF-8 BOM confusion**: Some Windows editors save UTF-8 files with a 3-byte BOM prefix (`\xef\xbb\xbf`). Decoders that do not strip the BOM produce a garbage invisible character at position 0, which silently breaks FASTA header parsing.
- **Windows line endings in sequences**: A `\r` at the end of a DNA sequence line is a legal base in no alphabet, but it is invisible in most text editors and terminal output.
- **Non-breaking space vs. regular space**: They have the same visual width, different byte values. Database lookups, regex matches, and split operations treat them as different characters.

---

## 1. From Bits to Bytes

Computers store everything as sequences of 0s and 1s. Eight binary digits form one **byte**, which can represent 256 distinct values (0--255).

| Notation | Example | Meaning |
|----------|---------|-------|
| Binary | `01000001` | 8 bits = 1 byte |
| Decimal | `65` | Human-friendly |
| Hexadecimal | `0x41` | Compact (2 hex digits = 1 byte) |

Hexadecimal (base-16) uses digits `0-9` and letters `A-F`. It is the standard way to display raw file bytes because each byte maps to exactly two hex digits.

```python
# Explore the relationship between decimal, binary, and hex
print(f"{'Decimal':>7}  {'Binary':>10}  {'Hex':>5}  {'Character'}")
print("-" * 40)

for value in [0, 10, 32, 48, 65, 90, 97, 122, 127, 200, 255]:
    binary = format(value, '08b')
    hexadecimal = format(value, '02X')
    # Only show printable ASCII characters
    char = chr(value) if 32 <= value < 127 else '--'
    print(f"{value:7d}  {binary:>10}  0x{hexadecimal}  '{char}'")
```python

```python
# Python's built-in functions for conversions
print("ord('A') =", ord('A'))        # Character -> integer code point
print("chr(65)  =", chr(65))         # Integer code point -> character
print("hex(65)  =", hex(65))         # Integer -> hex string
print("bin(65)  =", bin(65))         # Integer -> binary string
print("int('41', 16) =", int('41', 16))  # Hex string -> integer
```python

---

## 2. ASCII -- The Universal Foundation

ASCII (American Standard Code for Information Interchange, 1963) defines meanings for byte values 0--127. Every modern encoding is backward-compatible with ASCII.

| Range | Content |
|-------|---------|
| 0--31 | Control characters (newline `\n`=10, tab `\t`=9, carriage return `\r`=13) |
| 32 | Space |
| 48--57 | Digits `0`--`9` |
| 65--90 | Uppercase `A`--`Z` |
| 97--122 | Lowercase `a`--`z` |
| Other | Punctuation, symbols (`>`, `|`, `_`, etc.) |

### Why ASCII matters for bioinformatics

DNA (`ATGCN`), RNA (`AUGC`), and amino acid one-letter codes (`ACDEFGHIKLMNPQRSTVWY`) are all plain ASCII. FASTA headers, FASTQ quality scores, GFF columns -- all ASCII. This means **sequence data itself never has encoding issues**. The problems arise in *annotations*, *author names*, *organism names*, and *free-text fields*.

```python
# DNA nucleotides are pure ASCII
nucleotides = 'ATGC'
amino_acids = 'ACDEFGHIKLMNPQRSTVWY'

print("DNA bases in ASCII:")
for base in nucleotides:
    print(f"  {base} -> code {ord(base)} -> binary {format(ord(base), '08b')}")

print("\nAll amino acid codes are ASCII:")
all_ascii = all(ord(aa) < 128 for aa in amino_acids)
print(f"  {amino_acids} -> all < 128? {all_ascii}")
```python

```python
# FASTQ quality scores: also ASCII
# Phred+33 encoding maps quality 0-93 to ASCII 33-126
quality_string = "IIIIIIIIIIFFFFFFFDDDDDBBBBBB"

print("FASTQ quality decoding (Phred+33):")
for char in quality_string[:10]:  # Show first 10
    phred = ord(char) - 33
    error_prob = 10 ** (-phred / 10)
    print(f"  '{char}' -> ASCII {ord(char)} -> Phred {phred} -> error rate {error_prob:.6f}")
```python

---

## 3. Beyond ASCII: The Encoding Zoo

Byte values 128--255 are where encodings diverge. Different systems historically assigned different characters to these values.

### 3.1 Latin-1 (ISO 8859-1)

Covers Western European languages. Used in many older databases and some BioPython parsers as a fallback. Key property: **every byte value 0--255 is valid**, so decoding *never fails*. This makes it a useful "safety net" encoding, but it may produce wrong characters.

### 3.2 Legacy Cyrillic Encodings

If you work with Russian-language annotations or databases:
- **CP1251** (Windows-1251): Windows standard for Cyrillic
- **CP866**: DOS Cyrillic
- **KOI8-R**: Unix/Linux Cyrillic (older systems)

These are all single-byte, meaning one byte = one character. But the same byte maps to different characters in each!

### 3.3 UTF-8 -- The Modern Standard

UTF-8 is a variable-length encoding for Unicode:
- ASCII characters (0--127): 1 byte (identical to ASCII)
- Latin/Cyrillic/Greek: 2 bytes
- CJK characters: 3 bytes
- Emoji and rare symbols: 4 bytes

**UTF-8 is the default for virtually all modern bioinformatics tools, databases, and file formats.**

```python
# Same text, different encodings = different bytes
text = "alpha-helix"
greek_text = "\u03b1-helix"  # Using actual Greek alpha: alpha-helix

print(f"Text: {greek_text}")
print(f"  UTF-8 bytes:   {greek_text.encode('utf-8')}  ({len(greek_text.encode('utf-8'))} bytes)")
print(f"  Latin-1 fails for Greek alpha!")

try:
    greek_text.encode('latin-1')
except UnicodeEncodeError as e:
    print(f"  Error: {e}")

# But plain ASCII text is identical across all encodings
plain = "ATGCGATCGA"
print(f"\n'{plain}' is identical in UTF-8 and Latin-1:")
print(f"  UTF-8:   {plain.encode('utf-8')}")
print(f"  Latin-1: {plain.encode('latin-1')}")
print(f"  Equal?   {plain.encode('utf-8') == plain.encode('latin-1')}")
```python

```python
# Demonstration: same byte, different interpretation
byte_value = bytes([0xC0])  # Byte 192

print(f"Byte 0xC0 (192) decoded as:")
print(f"  Latin-1: '{byte_value.decode('latin-1')}'  (A with grave accent)")
print(f"  CP1251:  '{byte_value.decode('cp1251')}'   (Cyrillic A)")
print(f"  CP866:   '{byte_value.decode('cp866')}'    (Box-drawing character)")

try:
    byte_value.decode('utf-8')
except UnicodeDecodeError:
    print(f"  UTF-8:   ERROR -- 0xC0 alone is not valid UTF-8")
```python

```python
# UTF-8 variable-length encoding in action
examples = [
    ('A', 'ASCII letter'),
    ('\u03b1', 'Greek alpha (used in protein structure notation)'),
    ('\u00e9', 'e-acute (common in French author names)'),
    ('\u0410', 'Cyrillic A (Russian annotations)'),
    ('\u4e00', 'CJK character (Chinese/Japanese databases)'),
    ('\U0001F9EC', 'DNA emoji'),
]

print(f"{'Char':<5} {'UTF-8 bytes':<20} {'# bytes':<8} {'Context'}")
print("-" * 70)
for char, context in examples:
    encoded = char.encode('utf-8')
    hex_repr = ' '.join(f'{b:02X}' for b in encoded)
    print(f"{char:<5} {hex_repr:<20} {len(encoded):<8} {context}")
```python

---

## 4. Examining Raw File Bytes

When you suspect an encoding issue, the first step is to look at the raw bytes. Open the file in **binary mode** (`'rb'`) to see exactly what is stored on disk, without any decoding.

```python
def hexdump(data, width=16):
    """Display bytes in classic hexdump format."""
    for i in range(0, len(data), width):
        chunk = data[i:i + width]
        hex_part = ' '.join(f'{b:02X}' for b in chunk)
        ascii_part = ''.join(chr(b) if 32 <= b < 127 else '.' for b in chunk)
        print(f"{i:04X}  {hex_part:<{width * 3}}  {ascii_part}")


# Create a realistic FASTA header and examine its bytes
fasta = ">sp|P04637|P53_HUMAN Cellular tumor antigen p53 OS=Homo sapiens\n"
fasta += "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP\n"

print("FASTA content:")
print(fasta)
print("\nHexdump:")
hexdump(fasta.encode('utf-8'))
```python

```python
# What a BOM (Byte Order Mark) looks like
# Some Windows editors prepend UTF-8 BOM: EF BB BF
bom_content = b'\xef\xbb\xbf>sp|P04637|P53_HUMAN\nMEEPQ\n'

print("File with UTF-8 BOM:")
hexdump(bom_content)

# Decoding with utf-8-sig automatically strips the BOM
decoded_with_bom = bom_content.decode('utf-8')
decoded_without_bom = bom_content.decode('utf-8-sig')

print(f"\nWith 'utf-8':     starts with '{decoded_with_bom[0]}' (BOM character, invisible)")
print(f"With 'utf-8-sig': starts with '{decoded_without_bom[0]}' (BOM stripped)")
print(f"\nThe BOM character has code point: U+{ord(decoded_with_bom[0]):04X} (FEFF = ZERO WIDTH NO-BREAK SPACE)")
```python

---

## 5. Common Encoding Problems in Bioinformatics

### 5.1 Windows vs. Unix Line Endings

| System | Line ending | Bytes |
|--------|------------|-------|
| Unix/macOS | `\n` (LF) | `0A` |
| Windows | `\r\n` (CR+LF) | `0D 0A` |
| Old Mac (pre-OS X) | `\r` (CR) | `0D` |

Many bioinformatics tools expect Unix line endings. Windows-style `\r\n` can cause subtle bugs: extra `\r` at the end of sequence lines, off-by-one in coordinates, or FASTA parsers that include `\r` in the sequence.

```python
# Line ending problems
unix_fasta = b">gene1\nATGCGATCGA\nTTTAAAGGGC\n"
windows_fasta = b">gene1\r\nATGCGATCGA\r\nTTTAAAGGGC\r\n"

print("Unix FASTA (correct):")
hexdump(unix_fasta)

print("\nWindows FASTA (problematic):")
hexdump(windows_fasta)

# Naive parsing breaks with Windows line endings
lines = windows_fasta.decode('ascii').split('\n')
for i, line in enumerate(lines):
    print(f"\nLine {i}: {repr(line)}")
    if not line.startswith('>') and line.strip():
        # The \r at the end will be included in the sequence!
        print(f"  WARNING: sequence has trailing carriage return!")
```python

```python
# Fix: always strip lines, or use universal newline mode
def parse_fasta_robust(raw_bytes):
    """Parse FASTA from raw bytes, handling any line ending style."""
    # Decode as UTF-8, replace \r\n with \n, then also replace lone \r
    text = raw_bytes.decode('utf-8').replace('\r\n', '\n').replace('\r', '\n')
    
    sequences = {}
    current_header = None
    current_seq = []
    
    for line in text.strip().split('\n'):
        line = line.strip()
        if line.startswith('>'):
            if current_header is not None:
                sequences[current_header] = ''.join(current_seq)
            current_header = line[1:]
            current_seq = []
        elif line:
            current_seq.append(line)
    
    if current_header is not None:
        sequences[current_header] = ''.join(current_seq)
    
    return sequences


# Works with both line ending styles
result = parse_fasta_robust(windows_fasta)
for header, seq in result.items():
    print(f">{header}")
    print(f"{seq}")
    print(f"Length: {len(seq)} (no extra characters)")
```python

## Common Pitfalls

- **Coordinate systems**: BED uses 0-based half-open; VCF/GFF use 1-based inclusive — mixing them causes off-by-one errors
- **Batch effects**: Always check for batch confounding before interpreting biological signal
- **Multiple testing**: Apply FDR correction (Benjamini-Hochberg) when testing thousands of features simultaneously
