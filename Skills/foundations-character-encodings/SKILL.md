---
name: foundations-character-encodings
description: Every bioinformatics pipeline starts with reading data from files. FASTA sequences, GenBank records, PDB coordinates, GFF annotations -- they all live as bytes on disk. When your Python script reads g
tool_type: python
primary_tool: Python
---

# Character Encodings and Binary Data in Bioinformatics

## Pitfalls

- **Latin-1 never fails — but that does not mean it is right**: Latin-1 decodes every byte 0–255 successfully, silently producing wrong characters when the file is actually UTF-8 with multi-byte sequences.
- **UTF-8 BOM confusion**: Windows editors may prepend `\xef\xbb\xbf`. Decoders that don't strip the BOM insert an invisible garbage character at position 0, silently breaking FASTA header parsing. Use `'utf-8-sig'` to auto-strip.
- **Windows line endings in sequences**: `\r` at the end of a DNA line is invisible in most editors but breaks sequence length, regex, and split operations.
- **Non-breaking space vs. regular space**: Same visual width, different byte values — database lookups and split() treat them differently.

## From Bits to Bytes

| Notation | Example | Meaning |
|----------|---------|---------|
| Binary | `01000001` | 8 bits = 1 byte |
| Decimal | `65` | Human-friendly |
| Hexadecimal | `0x41` | 2 hex digits = 1 byte |

```python
for value in [0, 10, 32, 48, 65, 90, 97, 122, 127, 200, 255]:
    char = chr(value) if 32 <= value < 127 else '--'
    print(f"{value:7d}  {format(value,'08b'):>10}  0x{format(value,'02X')}  '{char}'")
```

## ASCII

Byte values 0–127. All modern encodings are backward-compatible with ASCII.

| Range | Content |
|-------|---------|
| 0–31 | Control chars (`\n`=10, `\t`=9, `\r`=13) |
| 32 | Space |
| 48–57 | `0`–`9` |
| 65–90 | `A`–`Z` |
| 97–122 | `a`–`z` |

**Bioinformatics relevance:** DNA (`ATGCN`), RNA (`AUGC`), and amino acid codes are pure ASCII — sequence data itself never has encoding issues. Problems arise in annotations, author names, organism names, and free-text fields.

```python
# FASTQ quality scores: Phred+33 maps quality 0-93 to ASCII 33-126
quality_string = "IIIIIIIIFFFFFFDDDDDBBBB"
for char in quality_string[:6]:
    phred = ord(char) - 33
    print(f"'{char}' -> Phred {phred} -> error rate {10**(-phred/10):.6f}")
```

## Encoding Zoo

Byte values 128–255 are where encodings diverge.

| Encoding | Scope | Key property |
|----------|-------|--------------|
| **Latin-1** (ISO 8859-1) | Western European | Every byte 0–255 valid; never raises error |
| **CP1251** | Cyrillic (Windows) | Byte 0xC0 = Cyrillic А |
| **CP866** | Cyrillic (DOS) | Byte 0xC0 = box-drawing char |
| **UTF-8** | Universal Unicode | Variable-length: ASCII=1B, Latin/Cyrillic=2B, CJK=3B, emoji=4B |

**UTF-8 is the default for all modern bioinformatics tools and databases.**

```python
# Same byte, three interpretations
byte_value = bytes([0xC0])
print(f"Latin-1: '{byte_value.decode('latin-1')}'")   # À
print(f"CP1251:  '{byte_value.decode('cp1251')}'")    # А (Cyrillic)
print(f"CP866:   '{byte_value.decode('cp866')}'")     # box-drawing
try:
    byte_value.decode('utf-8')
except UnicodeDecodeError:
    print("UTF-8:   ERROR — 0xC0 alone is not valid UTF-8")
```

```python
# UTF-8 variable-length in action
examples = [
    ('A', 'ASCII letter'),
    ('\u03b1', 'Greek alpha (protein structure)'),
    ('\u00e9', 'e-acute (French author names)'),
    ('\u0410', 'Cyrillic A (Russian annotations)'),
    ('\U0001F9EC', 'DNA emoji'),
]
for char, context in examples:
    enc = char.encode('utf-8')
    print(f"{char}  {' '.join(f'{b:02X}' for b in enc):<20} {len(enc)}B  {context}")
```

## Examining Raw File Bytes

```python
def hexdump(data, width=16):
    for i in range(0, len(data), width):
        chunk = data[i:i+width]
        hex_part = ' '.join(f'{b:02X}' for b in chunk)
        ascii_part = ''.join(chr(b) if 32 <= b < 127 else '.' for b in chunk)
        print(f"{i:04X}  {hex_part:<{width*3}}  {ascii_part}")

# BOM detection
bom_content = b'\xef\xbb\xbf>sp|P04637|P53_HUMAN\nMEEPQ\n'
decoded_with_bom    = bom_content.decode('utf-8')
decoded_without_bom = bom_content.decode('utf-8-sig')  # strips BOM
print(f"utf-8:     starts with U+{ord(decoded_with_bom[0]):04X}")   # FEFF (invisible)
print(f"utf-8-sig: starts with '{decoded_without_bom[0]}'")         # >
```

## Windows vs. Unix Line Endings

| System | Ending | Bytes |
|--------|--------|-------|
| Unix/macOS | `\n` (LF) | `0A` |
| Windows | `\r\n` (CR+LF) | `0D 0A` |
| Old Mac | `\r` (CR) | `0D` |

Many bioinformatics tools expect Unix line endings. `\r\n` causes extra `\r` in sequence lines, off-by-one coordinates, or FASTA parsers that include `\r` in the sequence.

```python
def parse_fasta_robust(raw_bytes):
    """Parse FASTA from raw bytes, handling any line ending style."""
    text = raw_bytes.decode('utf-8').replace('\r\n', '\n').replace('\r', '\n')
    sequences = {}; current_header = None; current_seq = []
    for line in text.strip().split('\n'):
        line = line.strip()
        if line.startswith('>'):
            if current_header is not None:
                sequences[current_header] = ''.join(current_seq)
            current_header = line[1:]; current_seq = []
        elif line:
            current_seq.append(line)
    if current_header is not None:
        sequences[current_header] = ''.join(current_seq)
    return sequences
```

## Pitfalls (General)

- **Coordinate systems**: BED 0-based half-open; VCF/GFF 1-based inclusive — mixing causes off-by-one errors
- **Batch effects**: Always check for batch confounding before interpreting biological signal
- **Multiple testing**: Apply FDR correction (Benjamini-Hochberg) for thousands of features
