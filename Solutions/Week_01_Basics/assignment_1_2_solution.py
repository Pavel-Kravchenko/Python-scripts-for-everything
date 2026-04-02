"""
Assignment 1.2: Codon Translator - SOLUTION
============================================

This script translates a DNA coding sequence into a protein sequence.
"""

# Complete genetic code dictionary
GENETIC_CODE = {
    'TTT': ('F', 'Phenylalanine'), 'TTC': ('F', 'Phenylalanine'),
    'TTA': ('L', 'Leucine'), 'TTG': ('L', 'Leucine'),
    'CTT': ('L', 'Leucine'), 'CTC': ('L', 'Leucine'),
    'CTA': ('L', 'Leucine'), 'CTG': ('L', 'Leucine'),
    'ATT': ('I', 'Isoleucine'), 'ATC': ('I', 'Isoleucine'),
    'ATA': ('I', 'Isoleucine'), 'ATG': ('M', 'Methionine'),  # START
    'GTT': ('V', 'Valine'), 'GTC': ('V', 'Valine'),
    'GTA': ('V', 'Valine'), 'GTG': ('V', 'Valine'),
    'TCT': ('S', 'Serine'), 'TCC': ('S', 'Serine'),
    'TCA': ('S', 'Serine'), 'TCG': ('S', 'Serine'),
    'CCT': ('P', 'Proline'), 'CCC': ('P', 'Proline'),
    'CCA': ('P', 'Proline'), 'CCG': ('P', 'Proline'),
    'ACT': ('T', 'Threonine'), 'ACC': ('T', 'Threonine'),
    'ACA': ('T', 'Threonine'), 'ACG': ('T', 'Threonine'),
    'GCT': ('A', 'Alanine'), 'GCC': ('A', 'Alanine'),
    'GCA': ('A', 'Alanine'), 'GCG': ('A', 'Alanine'),
    'TAT': ('Y', 'Tyrosine'), 'TAC': ('Y', 'Tyrosine'),
    'TAA': ('*', 'Stop'), 'TAG': ('*', 'Stop'),  # STOP codons
    'CAT': ('H', 'Histidine'), 'CAC': ('H', 'Histidine'),
    'CAA': ('Q', 'Glutamine'), 'CAG': ('Q', 'Glutamine'),
    'AAT': ('N', 'Asparagine'), 'AAC': ('N', 'Asparagine'),
    'AAA': ('K', 'Lysine'), 'AAG': ('K', 'Lysine'),
    'GAT': ('D', 'Aspartic acid'), 'GAC': ('D', 'Aspartic acid'),
    'GAA': ('E', 'Glutamic acid'), 'GAG': ('E', 'Glutamic acid'),
    'TGT': ('C', 'Cysteine'), 'TGC': ('C', 'Cysteine'),
    'TGA': ('*', 'Stop'),  # STOP codon
    'TGG': ('W', 'Tryptophan'),
    'CGT': ('R', 'Arginine'), 'CGC': ('R', 'Arginine'),
    'CGA': ('R', 'Arginine'), 'CGG': ('R', 'Arginine'),
    'AGT': ('S', 'Serine'), 'AGC': ('S', 'Serine'),
    'AGA': ('R', 'Arginine'), 'AGG': ('R', 'Arginine'),
    'GGT': ('G', 'Glycine'), 'GGC': ('G', 'Glycine'),
    'GGA': ('G', 'Glycine'), 'GGG': ('G', 'Glycine'),
}


def extract_codons(sequence):
    """
    Extract codons (triplets) from a DNA sequence.
    
    Parameters:
        sequence (str): DNA sequence
    
    Returns:
        list: List of codon strings
    """
    sequence = sequence.upper().replace(' ', '')
    codons = []
    for i in range(0, len(sequence) - 2, 3):
        codons.append(sequence[i:i+3])
    return codons


def translate_codon(codon):
    """
    Translate a single codon to amino acid.
    
    Parameters:
        codon (str): 3-nucleotide codon
    
    Returns:
        tuple: (single_letter_code, full_name) or ('X', 'Unknown')
    """
    codon = codon.upper()
    return GENETIC_CODE.get(codon, ('X', 'Unknown'))


def translate_sequence(sequence, stop_at_stop=True):
    """
    Translate a DNA sequence to protein.
    
    Parameters:
        sequence (str): DNA coding sequence
        stop_at_stop (bool): Stop translation at stop codon
    
    Returns:
        dict: Translation results including:
            - codons: list of codons
            - amino_acids: list of (codon, aa_code, aa_name) tuples
            - protein: protein sequence string
    """
    codons = extract_codons(sequence)
    amino_acids = []
    protein = []
    
    for codon in codons:
        aa_code, aa_name = translate_codon(codon)
        amino_acids.append((codon, aa_code, aa_name))
        
        if aa_code == '*':  # Stop codon
            if stop_at_stop:
                break
        else:
            protein.append(aa_code)
    
    return {
        'codons': codons,
        'amino_acids': amino_acids,
        'protein': ''.join(protein)
    }


def display_translation(sequence, results):
    """Display translation results in a formatted way."""
    print("\nDNA Translation")
    print("=" * 50)
    print(f"DNA: {sequence}")
    print(f"Length: {len(sequence)} bp ({len(sequence)//3} codons)")
    
    print("\nCodon breakdown:")
    for codon, aa_code, aa_name in results['amino_acids']:
        # Add special labels
        if codon == 'ATG':
            label = " - START"
        elif aa_code == '*':
            label = " - STOP"
        else:
            label = ""
        
        print(f"  {codon} → {aa_code} ({aa_name}){label}")
    
    print(f"\nProtein: {results['protein']}")
    print(f"Length: {len(results['protein'])} amino acids")


def validate_coding_sequence(sequence):
    """
    Validate a coding sequence.
    
    Returns:
        tuple: (is_valid, error_message)
    """
    sequence = sequence.upper().replace(' ', '')
    
    # Check for valid nucleotides
    valid_nucs = set("ATGC")
    if not all(nuc in valid_nucs for nuc in sequence):
        invalid = set(sequence) - valid_nucs
        return False, f"Invalid nucleotides: {invalid}"
    
    # Check length divisible by 3
    if len(sequence) % 3 != 0:
        remainder = len(sequence) % 3
        return False, f"Length ({len(sequence)}) not divisible by 3 (remainder: {remainder})"
    
    # Check for start codon
    if not sequence.startswith('ATG'):
        return False, "Warning: Sequence doesn't start with ATG (start codon)"
    
    return True, "Valid coding sequence"


def main():
    """Main function to run the codon translator."""
    print("Codon Translator")
    print("-" * 50)
    
    sequence = input("Enter DNA coding sequence: ").upper().strip()
    
    if not sequence:
        print("Error: No sequence provided.")
        return
    
    # Validate
    is_valid, message = validate_coding_sequence(sequence)
    if not is_valid:
        print(f"Error: {message}")
        return
    elif "Warning" in message:
        print(message)
    
    # Translate and display
    results = translate_sequence(sequence)
    display_translation(sequence, results)


if __name__ == "__main__":
    main()


# =============================================================================
# Test Cases
# =============================================================================

def test_translate():
    """Test translation function."""
    # Test 1: Simple translation
    seq1 = "ATGGCC"
    result1 = translate_sequence(seq1)
    assert result1['protein'] == "MA", f"Expected 'MA', got '{result1['protein']}'"
    
    # Test 2: Translation with stop codon
    seq2 = "ATGGCCTAG"
    result2 = translate_sequence(seq2)
    assert result2['protein'] == "MA", f"Expected 'MA', got '{result2['protein']}'"
    
    # Test 3: All stop codons
    for stop in ["TAA", "TAG", "TGA"]:
        seq = f"ATG{stop}"
        result = translate_sequence(seq)
        assert result['protein'] == "M"
    
    print("All translation tests passed! ✓")


def test_codon_extraction():
    """Test codon extraction."""
    assert extract_codons("ATGCAT") == ["ATG", "CAT"]
    assert extract_codons("ATG") == ["ATG"]
    assert extract_codons("AT") == []  # Too short
    
    print("Codon extraction tests passed! ✓")


# Uncomment to run tests:
# test_translate()
# test_codon_extraction()
