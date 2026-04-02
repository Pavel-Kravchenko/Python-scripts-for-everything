"""
Assignment 1.1: DNA Sequence Analyzer - SOLUTION
=================================================

This script analyzes a DNA sequence and provides basic statistics.
"""

def validate_sequence(sequence):
    """
    Validate that sequence contains only ATGC.
    
    Parameters:
        sequence (str): Input DNA sequence
    
    Returns:
        bool: True if valid, False otherwise
    """
    valid_nucleotides = set("ATGC")
    sequence = sequence.upper()
    return all(nuc in valid_nucleotides for nuc in sequence)


def analyze_dna(sequence):
    """
    Analyze a DNA sequence and return statistics.
    
    Parameters:
        sequence (str): DNA sequence (ATGC only)
    
    Returns:
        dict: Dictionary with analysis results including:
            - length: sequence length
            - counts: nucleotide counts
            - percentages: nucleotide percentages
            - gc_content: GC percentage
            - at_gc_ratio: AT to GC ratio
    """
    sequence = sequence.upper()
    length = len(sequence)
    
    # Count nucleotides
    counts = {
        'A': sequence.count('A'),
        'T': sequence.count('T'),
        'G': sequence.count('G'),
        'C': sequence.count('C')
    }
    
    # Calculate percentages
    percentages = {nuc: (count / length) * 100 for nuc, count in counts.items()}
    
    # Calculate GC content
    gc_count = counts['G'] + counts['C']
    gc_content = (gc_count / length) * 100
    
    # Calculate AT/GC ratio
    at_count = counts['A'] + counts['T']
    at_gc_ratio = at_count / gc_count if gc_count > 0 else float('inf')
    
    return {
        'sequence': sequence,
        'length': length,
        'counts': counts,
        'percentages': percentages,
        'gc_content': gc_content,
        'at_gc_ratio': at_gc_ratio
    }


def display_results(results):
    """
    Display analysis results in a formatted way.
    
    Parameters:
        results (dict): Analysis results from analyze_dna()
    """
    print("\nDNA Sequence Analysis")
    print("=" * 40)
    print(f"Sequence: {results['sequence']}")
    print(f"Length: {results['length']} bp")
    
    print("\nNucleotide Counts:")
    for nuc in 'ATGC':
        count = results['counts'][nuc]
        percent = results['percentages'][nuc]
        print(f"  {nuc}: {count} ({percent:.1f}%)")
    
    print(f"\nGC Content: {results['gc_content']:.1f}%")
    print(f"AT/GC Ratio: {results['at_gc_ratio']:.2f}")
    
    # Classification
    gc = results['gc_content']
    if gc < 40:
        classification = "AT-rich"
    elif gc <= 60:
        classification = "Balanced"
    else:
        classification = "GC-rich"
    print(f"Classification: {classification}")


def main():
    """Main function to run the DNA analyzer."""
    print("DNA Sequence Analyzer")
    print("-" * 40)
    
    sequence = input("Enter DNA sequence: ").upper().strip()
    
    # Validate input
    if not sequence:
        print("Error: No sequence provided.")
        return
    
    if not validate_sequence(sequence):
        # Find invalid characters
        invalid = set(sequence) - set("ATGC")
        print(f"Error: Invalid characters found: {invalid}")
        print("Use only A, T, G, C.")
        return
    
    # Analyze and display
    results = analyze_dna(sequence)
    display_results(results)


if __name__ == "__main__":
    main()


# =============================================================================
# Test Cases
# =============================================================================

def test_analyze_dna():
    """Test the analyze_dna function."""
    # Test case 1: Balanced sequence
    seq1 = "ATGC"
    result1 = analyze_dna(seq1)
    assert result1['length'] == 4
    assert result1['gc_content'] == 50.0
    assert result1['at_gc_ratio'] == 1.0
    
    # Test case 2: AT-rich sequence
    seq2 = "AATTAATT"
    result2 = analyze_dna(seq2)
    assert result2['gc_content'] == 0.0
    
    # Test case 3: GC-rich sequence
    seq3 = "GCGCGCGC"
    result3 = analyze_dna(seq3)
    assert result3['gc_content'] == 100.0
    
    print("All tests passed! ✓")


def test_validate_sequence():
    """Test sequence validation."""
    assert validate_sequence("ATGC") == True
    assert validate_sequence("atgc") == True
    assert validate_sequence("ATGCX") == False
    assert validate_sequence("ATGCU") == False  # RNA
    assert validate_sequence("") == True  # Empty is technically valid
    
    print("Validation tests passed! ✓")


# Uncomment to run tests:
# test_analyze_dna()
# test_validate_sequence()
