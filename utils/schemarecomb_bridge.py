"""
Bridge utilities to convert between our data formats and schemarecomb formats.

This module provides conversion functions to integrate schemarecomb's
ParentSequences, PDBStructure, and Library objects with our Streamlit app.
"""

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pathlib import Path
import tempfile

try:
    import schemarecomb
    from schemarecomb import ParentSequences, PDBStructure
    from schemarecomb.energy_functions import SCHEMA
    SCHEMARECOMB_AVAILABLE = True
except (ImportError, AttributeError):
    SCHEMARECOMB_AVAILABLE = False
    # Create dummy classes for type hints
    ParentSequences = None
    PDBStructure = None
    SCHEMA = None


def sequences_to_parent_sequences(sequences, msa_file=None, progress_callback=None):
    """
    Convert our (id, sequence) tuple format to schemarecomb ParentSequences.
    
    Args:
        sequences: List of (id, sequence) tuples, where first is ("query", query_seq)
        msa_file: Optional path to MSA file (if already aligned)
        progress_callback: Optional callback for progress updates
        
    Returns:
        schemarecomb.ParentSequences object
    """
    if not SCHEMARECOMB_AVAILABLE:
        raise ImportError("schemarecomb is required. Install with: pip install schemarecomb")
    if progress_callback:
        progress_callback(0.1, "Converting sequences to ParentSequences format...")
    
    # Convert to BioPython SeqRecords
    records = []
    for seq_id, seq in sequences:
        # Clean up sequence ID for BioPython
        clean_id = seq_id.replace(' ', '_').replace('|', '_')
        if clean_id == "query":
            clean_id = "query_sequence"
        
        seq_record = SeqRecord(
            Seq(seq),
            id=clean_id,
            name=clean_id,
            description=f"Sequence {clean_id}"
        )
        records.append(seq_record)
    
    if progress_callback:
        progress_callback(0.5, "Creating ParentSequences object...")
    
    # If MSA file exists, load as prealigned
    if msa_file and Path(msa_file).exists():
        parents = ParentSequences.from_fasta(msa_file, prealigned=True)
        # Verify records match
        if len(parents.records) != len(records):
            # If mismatch, create from records and align
            parents = ParentSequences(records)
            if progress_callback:
                progress_callback(0.7, "Aligning sequences...")
            parents.align()
    else:
        # Create from records and align
        parents = ParentSequences(records)
        if progress_callback:
            progress_callback(0.7, "Aligning sequences...")
        parents.align()
    
    if progress_callback:
        progress_callback(1.0, "ParentSequences created")
    
    return parents


def msa_to_parent_sequences(msa_file, pdb_structure=None, progress_callback=None):
    """
    Load MSA file into schemarecomb ParentSequences.
    
    Handles both FASTA format and SCHEMA format (space-separated).
    
    Args:
        msa_file: Path to MSA file (FASTA or SCHEMA format)
        pdb_structure: Optional schemarecomb.PDBStructure object
        progress_callback: Optional callback for progress updates
        
    Returns:
        schemarecomb.ParentSequences object
    """
    if not SCHEMARECOMB_AVAILABLE:
        raise ImportError("schemarecomb is required. Install with: pip install schemarecomb")
    if msa_file is None:
        raise ValueError("msa_file cannot be None")
    msa_path = Path(msa_file)
    if not msa_path.exists():
        raise FileNotFoundError(f"MSA file not found: {msa_file}")
    if progress_callback:
        progress_callback(0.1, f"Loading MSA from {msa_file}...")
    
    # Try to read the file to determine format
    with open(msa_file, 'r') as f:
        first_line = f.readline().strip()
        f.seek(0)
        content = f.read()
    
    # Check if it's FASTA format (starts with >)
    if first_line.startswith('>'):
        # FASTA format - use from_fasta directly
        try:
            parents = ParentSequences.from_fasta(msa_file, prealigned=True, pdb_structure=pdb_structure)
        except ValueError as e:
            if "records must not be empty" in str(e):
                # Try parsing manually to see what's wrong
                from Bio import SeqIO
                records = list(SeqIO.parse(msa_file, 'fasta'))
                if len(records) == 0:
                    raise ValueError(f"MSA file {msa_file} contains no sequences. File content preview: {content[:200]}")
                raise
            raise
    else:
        # SCHEMA format (space-separated: ID SEQUENCE)
        # Convert to FASTA format first
        if progress_callback:
            progress_callback(0.3, "Converting SCHEMA format to FASTA...")
        
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        import tempfile
        
        records = []
        for line in content.split('\n'):
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            
            # SCHEMA format: ID SEQUENCE (space-separated)
            parts = line.split(None, 1)  # Split on whitespace, max 1 split
            if len(parts) == 2:
                seq_id, sequence = parts
                seq_record = SeqRecord(
                    Seq(sequence),
                    id=seq_id,
                    name=seq_id,
                    description=f"Sequence {seq_id}"
                )
                records.append(seq_record)
        
        if len(records) == 0:
            raise ValueError(f"MSA file {msa_file} contains no sequences. File content preview: {content[:200]}")
        
        # Create ParentSequences from records
        parents = ParentSequences(records, prealigned=True, pdb_structure=pdb_structure)
    
    if progress_callback:
        progress_callback(1.0, f"ParentSequences loaded ({len(parents.records)} sequences)")
    
    return parents


def load_pdb_structure(pdb_file, chain='A', progress_callback=None):
    """
    Load PDB file into schemarecomb PDBStructure.
    
    Args:
        pdb_file: Path to PDB file
        chain: Chain identifier (default: 'A')
        progress_callback: Optional callback for progress updates
        
    Returns:
        schemarecomb.PDBStructure object
    """
    if not SCHEMARECOMB_AVAILABLE:
        raise ImportError("schemarecomb is required. Install with: pip install schemarecomb")
    if progress_callback:
        progress_callback(0.1, f"Loading PDB structure from {pdb_file}...")
    
    pdb_structure = PDBStructure.from_pdb_file(pdb_file, chain=chain)
    
    if progress_callback:
        progress_callback(1.0, "PDBStructure loaded")
    
    return pdb_structure


def align_pdb_to_parents(pdb_structure, parents, progress_callback=None, use_alignment_numbering=True):
    """
    Align PDB structure to parent sequences.
    
    This function:
    1. Trims the alignment to remove terminal gaps from the query (first parent)
    2. Renumbers the PDB structure to match alignment positions
    3. Only includes structure residues that align to query positions
    
    Args:
        pdb_structure: schemarecomb.PDBStructure object
        parents: schemarecomb.ParentSequences object (must be aligned)
        progress_callback: Optional callback for progress updates
        use_alignment_numbering: If True (default), use alignment column positions for
                                 residue numbering. If False, use sequential numbering.
        
    Returns:
        tuple: (schemarecomb.PDBStructure object (renumbered), trimmed ParentSequences)
    """
    if not SCHEMARECOMB_AVAILABLE:
        raise ImportError("schemarecomb is required. Install with: pip install schemarecomb")
    
    if progress_callback:
        progress_callback(0.1, "Trimming alignment to query boundaries...")
    
    # Step 1: Trim alignment to query boundaries (no terminal gaps in query)
    try:
        trimmed_parents = trim_alignment_to_query(parents, query_index=0)
    except Exception as e:
        # Fall back to original if trimming fails
        trimmed_parents = parents
        if progress_callback:
            progress_callback(0.2, f"Note: Could not trim alignment ({e}), using original")
    
    if progress_callback:
        progress_callback(0.3, "Getting trimmed query sequence...")
    
    # Get first parent's aligned sequence from trimmed alignment
    try:
        p0_aligned = trimmed_parents.p0_aligned
    except AttributeError:
        # Fallback: construct from alignment
        if hasattr(trimmed_parents, 'alignment'):
            p0_aligned = ''.join([aminos[0] for aminos in trimmed_parents.alignment])
        else:
            raise ValueError("ParentSequences must be aligned before renumbering PDB structure")
    
    if progress_callback:
        progress_callback(0.5, "Renumbering PDB structure to alignment positions...")
    
    # Step 2: Renumber PDB to alignment positions
    if use_alignment_numbering:
        pdb_structure = renumber_pdb_to_alignment(pdb_structure, p0_aligned, progress_callback)
    else:
        # Use standard schemarecomb renumbering (sequential)
        pdb_structure.renumber(p0_aligned)
    
    if progress_callback:
        progress_callback(1.0, "PDB structure aligned to trimmed domain")
    
    return pdb_structure, trimmed_parents


def contacts_from_pdb_structure(pdb_structure, contact_distance=4.5):
    """
    Extract contacts from renumbered PDBStructure.
    
    Args:
        pdb_structure: schemarecomb.PDBStructure object (must be renumbered)
        contact_distance: Distance threshold in Angstroms (schemarecomb uses 4.5)
        
    Returns:
        list: List of contact tuples (i, j) where i and j are residue indices
              in parent-aligned coordinates
    """
    if not SCHEMARECOMB_AVAILABLE:
        raise ImportError("schemarecomb is required. Install with: pip install schemarecomb")
    # schemarecomb calculates contacts automatically and caches them
    # The contacts are already in parent-aligned coordinates after renumbering
    contacts = pdb_structure.contacts
    
    # Convert to format compatible with our existing code
    # schemarecomb returns [(i, j), (j, i)] for each contact pair
    # We'll return just the unique pairs
    unique_contacts = []
    seen = set()
    for i, j in contacts:
        if (i, j) not in seen and (j, i) not in seen:
            unique_contacts.append((i, j))
            seen.add((i, j))
    
    return unique_contacts


def library_to_results(library):
    """
    Convert schemarecomb Library object to our result format for display.
    
    Args:
        library: schemarecomb.Library object
        
    Returns:
        dict: Dictionary with results in our format
    """
    return {
        'energy': float(library.energy),
        'mutation_rate': float(library.mutation_rate),
        'breakpoints': [bp.position for bp in library.breakpoints],
        'block_indices': library.block_indices,
        'min_block_len': library.min_block_len,
        'max_block_len': library.max_block_len,
        'gg_prob': float(library.gg_prob) if hasattr(library, 'gg_prob') else None,
        'gg_overhangs': library.gg_overhangs if hasattr(library, 'gg_overhangs') else None,
        'dna_blocks': library.dna_blocks if hasattr(library, 'dna_blocks') else None,
    }


def parent_sequences_to_json(parents):
    """
    Serialize ParentSequences to JSON for session state storage.
    
    Args:
        parents: schemarecomb.ParentSequences object
        
    Returns:
        str: JSON string
    """
    if not SCHEMARECOMB_AVAILABLE:
        raise ImportError("schemarecomb is required. Install with: pip install schemarecomb")
    return parents.to_json()


def parent_sequences_from_json(json_str):
    """
    Deserialize ParentSequences from JSON.
    
    Args:
        json_str: JSON string from to_json()
        
    Returns:
        schemarecomb.ParentSequences object
    """
    if not SCHEMARECOMB_AVAILABLE:
        raise ImportError("schemarecomb is required. Install with: pip install schemarecomb")
    return ParentSequences.from_json(json_str)


def save_parents_to_fasta(parents, output_file):
    """
    Save ParentSequences to a FASTA file.
    
    Args:
        parents: schemarecomb.ParentSequences object
        output_file: Path to output FASTA file
    """
    from Bio import SeqIO
    
    # Write records to FASTA
    with open(output_file, 'w') as f:
        SeqIO.write(parents.records, f, 'fasta')
    
    return output_file


def trim_alignment_to_query(parents, query_index=0):
    """
    Trim alignment to remove terminal gaps from the query sequence.
    
    The query sequence defines the domain boundaries. This function:
    1. Finds the first and last non-gap positions in the query
    2. Trims all sequences to those boundaries
    3. Returns a new ParentSequences with trimmed sequences
    
    Args:
        parents: schemarecomb.ParentSequences object (must be aligned)
        query_index: Index of the query sequence (default: 0, the first parent)
        
    Returns:
        schemarecomb.ParentSequences object with trimmed alignment
    """
    if not SCHEMARECOMB_AVAILABLE:
        raise ImportError("schemarecomb is required. Install with: pip install schemarecomb")
    
    # Get the alignment
    if not hasattr(parents, 'alignment') or parents.alignment is None:
        raise ValueError("ParentSequences must be aligned before trimming")
    
    alignment = parents.alignment  # List of tuples, one per column
    
    if not alignment:
        return parents
    
    # Find first and last non-gap positions in query sequence
    first_pos = None
    last_pos = None
    
    for i, column in enumerate(alignment):
        query_char = column[query_index]
        if query_char != '-' and query_char != '.':
            if first_pos is None:
                first_pos = i
            last_pos = i
    
    if first_pos is None:
        raise ValueError("Query sequence has no non-gap characters")
    
    # If no trimming needed, return original
    if first_pos == 0 and last_pos == len(alignment) - 1:
        return parents
    
    # Trim the alignment
    trimmed_alignment = alignment[first_pos:last_pos + 1]
    
    # Build new sequences from trimmed alignment
    num_parents = len(parents.records)
    trimmed_seqs = [''] * num_parents
    
    for column in trimmed_alignment:
        for p_idx in range(num_parents):
            trimmed_seqs[p_idx] += column[p_idx]
    
    # Create new SeqRecords with trimmed sequences
    new_records = []
    for i, record in enumerate(parents.records):
        new_record = SeqRecord(
            Seq(trimmed_seqs[i]),
            id=record.id,
            name=record.name,
            description=record.description
        )
        new_records.append(new_record)
    
    # Create new ParentSequences with trimmed alignment
    trimmed_parents = ParentSequences(new_records, prealigned=True)
    
    return trimmed_parents


def renumber_pdb_to_alignment(pdb_structure, p0_aligned, progress_callback=None):
    """
    Renumber PDB structure to match alignment column positions.
    
    This function:
    - Uses alignment column positions as residue numbers
    - If structure has a gap at column N, there is no residue N
    - If query has a gap at column N (internal gap), structure residue there is excluded
    - Residue numbers correspond to positions in the trimmed alignment
    
    Args:
        pdb_structure: schemarecomb.PDBStructure object
        p0_aligned: First parent's aligned sequence (trimmed, may have internal gaps)
        progress_callback: Optional callback for progress updates
        
    Returns:
        schemarecomb.PDBStructure object (renumbered using alignment column positions)
    """
    if not SCHEMARECOMB_AVAILABLE:
        raise ImportError("schemarecomb is required. Install with: pip install schemarecomb")
    
    from Bio import pairwise2
    
    if progress_callback:
        progress_callback(0.1, "Preparing alignment...")
    
    # Derenumber if already renumbered
    if hasattr(pdb_structure, 'unrenumbered_amino_acids'):
        pdb_structure.derenumber()
    
    # Order amino acids by index
    index_to_aa = {aa.index: aa for aa in pdb_structure.amino_acids}
    ordered_aas = [index_to_aa[i] for i in sorted(index_to_aa)]
    
    # Get PDB sequence
    pdb_seq = ''.join([aa.letter for aa in ordered_aas])
    
    # Build mapping from ungapped query position to alignment column position
    # p0_aligned may have '-' for internal gaps
    ungapped_to_column = {}  # ungapped_pos -> column_pos
    column_to_is_gap = {}    # column_pos -> True if query has gap
    ungapped_pos = 0
    for col, char in enumerate(p0_aligned):
        column_to_is_gap[col] = (char == '-')
        if char != '-':
            ungapped_to_column[ungapped_pos] = col
            ungapped_pos += 1
    
    # Get ungapped query sequence for pairwise alignment
    p0_ungapped = p0_aligned.replace('-', '')
    
    if progress_callback:
        progress_callback(0.3, "Running pairwise alignment...")
    
    # Align PDB sequence to UNGAPPED query sequence
    aln = pairwise2.align.globalxx(p0_ungapped, pdb_seq, gap_char='.',
                                   one_alignment_only=True)[0]
    aligned_query = aln.seqA   # May have '.' for gaps
    aligned_pdb = aln.seqB     # May have '.' for gaps
    
    if progress_callback:
        progress_callback(0.5, "Renumbering residues to alignment column positions...")
    
    # Renumber using alignment column positions
    # 
    # Strategy:
    # 1. Iterate through the pairwise alignment
    # 2. Track position in the ungapped query sequence
    # 3. Use ungapped_to_column to map to alignment column position
    # 4. That column position becomes the residue number
    
    renumbered_aas = []
    unrenumbered_aas = []
    
    pdb_iter = iter(ordered_aas)
    ungapped_query_pos = 0  # Position in ungapped query sequence
    
    for query_char, pdb_char in zip(aligned_query, aligned_pdb):
        if pdb_char != '.':
            # PDB has a residue here
            candidate_aa = next(pdb_iter)
            
            if query_char != '.':
                # Both have residues - get the alignment column position
                column_pos = ungapped_to_column[ungapped_query_pos]
                candidate_aa.renumber(column_pos)
                renumbered_aas.append(candidate_aa)
            else:
                # PDB has residue but query has gap in pairwise alignment
                # This is an insertion in structure - exclude
                unrenumbered_aas.append(candidate_aa)
        
        # Increment ungapped_query_pos for each non-gap query character
        if query_char != '.':
            ungapped_query_pos += 1
    
    if progress_callback:
        progress_callback(0.9, "Finalizing structure...")
    
    # Update the structure
    pdb_structure.amino_acids = renumbered_aas
    pdb_structure.unrenumbered_amino_acids = unrenumbered_aas
    pdb_structure.renumbering_seq = p0_aligned
    
    # Clear cached contacts
    if hasattr(pdb_structure, '_contacts'):
        del pdb_structure._contacts
    
    if progress_callback:
        progress_callback(1.0, "Structure renumbered to alignment positions")
    
    return pdb_structure
