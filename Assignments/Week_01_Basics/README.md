# Week 1: Python Basics - Assignments

## Assignment 1.1: DNA Sequence Analyzer

**Difficulty**: ⭐ Beginner

### Problem Statement
Write a Python script that analyzes a DNA sequence and provides basic statistics.

### Requirements
1. Accept a DNA sequence as input
2. Validate the sequence (only A, T, G, C allowed)
3. Calculate and display:
   - Total length
   - Count of each nucleotide
   - GC content percentage
   - AT/GC ratio

### Example Input
```
ATGCGATCGATCGTAGC
```

### Expected Output
```
DNA Sequence Analysis
=====================
Sequence: ATGCGATCGATCGTAGC
Length: 17 bp

Nucleotide Counts:
  A: 3 (17.6%)
  T: 4 (23.5%)
  G: 5 (29.4%)
  C: 5 (29.4%)

GC Content: 58.8%
AT/GC Ratio: 0.70
```

### Starter Code
```python
# Assignment 1.1: DNA Sequence Analyzer

def analyze_dna(sequence):
    """
    Analyze a DNA sequence and return statistics.
    
    Parameters:
        sequence (str): DNA sequence (ATGC only)
    
    Returns:
        dict: Dictionary with analysis results
    """
    # Your code here
    pass

def validate_sequence(sequence):
    """
    Validate that sequence contains only ATGC.
    
    Returns:
        bool: True if valid, False otherwise
    """
    # Your code here
    pass

def main():
    sequence = input("Enter DNA sequence: ").upper().strip()
    
    if not validate_sequence(sequence):
        print("Error: Invalid sequence. Use only A, T, G, C.")
        return
    
    results = analyze_dna(sequence)
    # Display results
    # Your code here

if __name__ == "__main__":
    main()
```

---

## Assignment 1.2: Codon Translator

**Difficulty**: ⭐⭐ Intermediate

### Problem Statement
Create a program that translates a DNA coding sequence into a protein sequence.

### Requirements
1. Read a DNA sequence (must be divisible by 3)
2. Split into codons (3-nucleotide chunks)
3. Translate each codon to its amino acid
4. Stop at stop codons (TAA, TAG, TGA)
5. Display the result with codon-by-codon breakdown

### Example Input
```
ATGGCCGATCGATAG
```

### Expected Output
```
DNA Translation
===============
DNA: ATGGCCGATCGATAG

Codon breakdown:
  ATG → M (Methionine) - START
  GCC → A (Alanine)
  GAT → D (Aspartic acid)
  CGA → R (Arginine)
  TAG → * - STOP

Protein: MADR
Length: 4 amino acids
```

---

## Assignment 1.3: Complement Generator

**Difficulty**: ⭐ Beginner

### Problem Statement
Generate the complement and reverse complement of a DNA sequence.

### Requirements
1. Calculate the complementary strand (A↔T, G↔C)
2. Calculate the reverse complement
3. Display with proper 5' and 3' orientation

### Example Input
```
ATGCGATC
```

### Expected Output
```
DNA Strand Viewer
=================

Original:         5'-ATGCGATC-3'
                     ||||||||
Complement:       3'-TACGCTAG-5'

Reverse Complement: 5'-GATCGCAT-3'
```

---

## Assignment 1.4: GC Content Calculator (File-based)

**Difficulty**: ⭐⭐ Intermediate

### Problem Statement
Read sequences from a FASTA file and calculate GC content for each.

### Requirements
1. Parse a FASTA file (handle multiple sequences)
2. Calculate GC content for each sequence
3. Identify the sequence with highest/lowest GC
4. Write results to a CSV file

### Sample FASTA Input (sequences.fasta)
```
>seq1 | Human gene
ATGCGATCGATCGTAGCGATCGATCG
>seq2 | Mouse gene
GCGCGCGCGCGCATATATAT
>seq3 | Rat gene
ATATATATATATGCGCGCGC
```

### Expected Output (results.csv)
```
id,length,gc_content,classification
seq1,26,53.85,Moderate GC
seq2,20,60.00,High GC
seq3,20,40.00,Moderate GC
```

---

## Assignment 1.5: Reading Frame Finder

**Difficulty**: ⭐⭐⭐ Advanced

### Problem Statement
Find all possible Open Reading Frames (ORFs) in a DNA sequence across all 6 reading frames (3 forward, 3 reverse).

### Requirements
1. Identify all ORFs (start with ATG, end with TAA/TAG/TGA)
2. Search all 6 reading frames
3. Report position, length, and sequence
4. Filter ORFs by minimum length (default: 30 bp)

### Example Output
```
ORF Finder Results
==================
Sequence length: 150 bp

Found 3 ORF(s):

Frame +1:
  ORF 1: pos 12-48 (36 bp)
  DNA: ATGGCCGATCGATAGCCAATGGATCGATAGCCA...
  Protein: MADR...

Frame +2:
  No ORFs found

Frame +3:
  ORF 2: pos 45-90 (45 bp)
  ...
```

---

## Submission Guidelines

1. Save each solution as a separate Python file:
   - `assignment_1_1.py`
   - `assignment_1_2.py`
   - etc.

2. Include docstrings for all functions

3. Add comments explaining your logic

4. Test with multiple inputs before submitting

## Grading Criteria

- **Correctness** (40%): Does the code produce correct results?
- **Code Quality** (30%): Is the code readable, well-organized?
- **Error Handling** (20%): Does it handle edge cases?
- **Documentation** (10%): Are functions documented?
