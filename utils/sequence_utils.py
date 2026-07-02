"""
Sequence utility functions.
Consolidates sequence-related operations like identity calculation.
"""

from Bio import pairwise2


def calculate_identity(seq1, seq2):
    """
    Calculate percentage identity between two sequences using global alignment.
    
    Args:
        seq1: First sequence string
        seq2: Second sequence string
        
    Returns:
        float: Percentage identity (0-100)
    """
    if not seq1 or not seq2:
        return 0.0
    
    alignments = pairwise2.align.globalxx(seq1, seq2, one_alignment_only=True)
    
    if alignments:
        aln1, aln2, score, start, end = alignments[0]
        matches = sum(1 for a, b in zip(aln1, aln2) if a == b and a != '-')
        length = max(len(seq1), len(seq2))
        return (matches / length * 100) if length > 0 else 0.0
    
    return 0.0


def extract_uniprot_id(subject_acc):
    """
    Extract UniProt ID from various formats.
    
    Handles:
    - AlphaFold format: AF-{uniprot_id}-F1
    - UniProt format: {uniprot_id} or {uniprot_id}|...
    
    Args:
        subject_acc: Subject accession string
        
    Returns:
        str: Clean UniProt ID
    """
    if subject_acc.startswith('AF-'):
        # Format: AF-{uniprot_id}-F1
        parts = subject_acc.split('-')
        if len(parts) >= 2:
            return parts[1]
        return subject_acc
    else:
        # Format: {uniprot_id} or {uniprot_id}|...
        return subject_acc.split('|')[0] if '|' in subject_acc else subject_acc


def clean_sequence_id(seq_id):
    """
    Clean sequence ID for use in file names and identifiers.
    
    Args:
        seq_id: Raw sequence ID
        
    Returns:
        str: Cleaned sequence ID
    """
    return seq_id.replace('|', '_').replace('.', '_').replace(' ', '_')


def validate_sequence(sequence):
    """
    Validate that a sequence contains only valid amino acid characters.
    
    Args:
        sequence: Protein sequence string
        
    Returns:
        tuple: (is_valid, error_message)
    """
    if not sequence or len(sequence.strip()) < 20:
        return False, "Sequence must be at least 20 amino acids long"
    
    valid_aa = set('ACDEFGHIKLMNPQRSTVWY')
    clean_seq = ''.join(sequence.split()).upper()
    
    invalid_chars = set(clean_seq) - valid_aa
    if invalid_chars:
        return False, f"Sequence contains invalid characters: {', '.join(invalid_chars)}"
    
    return True, None
