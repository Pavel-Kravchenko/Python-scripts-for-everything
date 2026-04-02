#!/usr/bin/env python3
"""
Alignment Conservation Calculator

Calculate identity and similarity percentages for multiple sequence alignments.
Uses BLOSUM62 matrix for similarity scoring.

Usage:
    python alignment_stats.py -i alignment.fasta [-m matrix_file]
"""

import sys
import argparse
from pathlib import Path


def read_matrix(filepath):
    """
    Read BLOSUM or PAM matrix file.
    
    Expected format:
    # Comment lines start with #
       A  R  N  D  C  Q  E  G  H  I  ...
    A  4 -1 -2 -2  0 -1 -1  0 -2 -1  ...
    R -1  5  0 -2 -3  1  0 -2  0 -3  ...
    ...
    
    Returns: dict of {(aa1, aa2): score}
    """
    matrix = {}
    header = None
    
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            
            parts = line.split()
            
            if header is None:
                # First non-comment line is the header
                header = parts
                continue
            
            # Data rows: first element is row label
            row_aa = parts[0]
            scores = [int(x) for x in parts[1:]]
            
            for i, col_aa in enumerate(header):
                if i < len(scores):
                    matrix[(row_aa, col_aa)] = scores[i]
    
    return matrix


def read_fasta(filepath):
    """
    Read FASTA alignment file.
    
    Returns: list of [accession, description, sequence]
    """
    alignment = []
    
    with open(filepath) as f:
        current = None
        for line in f:
            line = line.strip()
            if not line:
                continue
            
            if line.startswith('>'):
                parts = line[1:].split(None, 1)
                acc = parts[0]
                desc = parts[1] if len(parts) > 1 else ""
                current = [acc, desc, ""]
                alignment.append(current)
            elif current is not None:
                current[2] += line
    
    return alignment


def check_alignment(alignment):
    """Verify all sequences have equal length."""
    if not alignment:
        print("Error: Empty alignment")
        sys.exit(1)
    
    length = len(alignment[0][2])
    for seq in alignment:
        if len(seq[2]) != length:
            print(f"Error: Sequence {seq[0]} has different length")
            sys.exit(1)
    
    return length


def analyze_position(alignment, pos, matrix, gaps={'-', '.', '~'}):
    """
    Analyze a single column of the alignment.
    
    Returns: (is_identical, is_similar)
    - identical: all non-gap characters are the same
    - similar: all pairwise scores are positive (BLOSUM)
    """
    column = [seq[2][pos] for seq in alignment]
    
    # Check for gaps
    has_gaps = any(c in gaps for c in column)
    non_gap_chars = [c for c in column if c not in gaps]
    
    # Identity check
    is_identical = len(set(non_gap_chars)) == 1 and not has_gaps
    
    # Similarity check (all positive pairwise scores)
    is_similar = False
    if not has_gaps and len(non_gap_chars) == len(column):
        is_similar = True
        for i in range(len(column)):
            for j in range(i, len(column)):
                pair = (column[i].upper(), column[j].upper())
                score = matrix.get(pair, matrix.get((pair[1], pair[0]), -999))
                if score <= 0:
                    is_similar = False
                    break
            if not is_similar:
                break
    
    return is_identical, is_similar


def count_positions(alignment, matrix):
    """Count identical and similar positions in alignment."""
    length = len(alignment[0][2])
    identity = 0
    similarity = 0
    
    for i in range(length):
        ident, sim = analyze_position(alignment, i, matrix)
        if ident:
            identity += 1
        if sim:
            similarity += 1
    
    return identity, similarity, length


def main():
    parser = argparse.ArgumentParser(
        description="Calculate alignment conservation statistics"
    )
    parser.add_argument(
        "-i", "--fasta", required=True,
        help="Input alignment in FASTA format"
    )
    parser.add_argument(
        "-m", "--matrix", default=None,
        help="Substitution matrix file (default: built-in BLOSUM62)"
    )
    
    args = parser.parse_args()
    
    # Load or use default matrix
    if args.matrix:
        matrix = read_matrix(args.matrix)
    else:
        # Minimal built-in BLOSUM62 for common amino acids
        # For full matrix, use -m option with EBLOSUM62.txt
        matrix = {
            ('A', 'A'): 4, ('A', 'S'): 1,
            ('V', 'I'): 3, ('V', 'L'): 1, ('V', 'V'): 4,
            ('I', 'I'): 4, ('I', 'L'): 2, ('L', 'L'): 4,
            ('F', 'Y'): 3, ('F', 'F'): 6, ('Y', 'Y'): 7,
            ('F', 'W'): 1, ('W', 'W'): 11, ('W', 'Y'): 2,
            ('S', 'T'): 1, ('S', 'S'): 4, ('T', 'T'): 5,
            ('K', 'R'): 2, ('K', 'K'): 5, ('R', 'R'): 5,
            ('D', 'E'): 2, ('D', 'D'): 6, ('E', 'E'): 5,
            ('N', 'D'): 1, ('N', 'N'): 6, ('Q', 'E'): 2, ('Q', 'Q'): 5,
        }
        print("Note: Using minimal built-in matrix. For full BLOSUM62, use -m option.")
    
    # Read alignment
    alignment = read_fasta(args.fasta)
    length = check_alignment(alignment)
    
    print(f"\nAlignment: {len(alignment)} sequences, {length} positions")
    print("-" * 50)
    
    # Calculate statistics
    identity, similarity, total = count_positions(alignment, matrix)
    
    print(f"Identical positions:  {identity:4d} ({100*identity/total:5.1f}%)")
    print(f"Similar positions:    {similarity:4d} ({100*similarity/total:5.1f}%)")
    print(f"Total positions:      {total:4d}")


if __name__ == "__main__":
    main()
