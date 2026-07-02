"""
Wrapper functions for SCHEMA energy calculations.
Provides a Python API for the bundled schema_raspp modules and schemarecomb integration.
"""

from pathlib import Path

try:
    from schema_raspp import schema, schemacontacts, schemaenergy, pdb, pdbseq
except ImportError:
    import schema, schemacontacts, schemaenergy, pdb, pdbseq

try:
    import schemarecomb
    from schemarecomb import ParentSequences, PDBStructure
    from schemarecomb.energy_functions import SCHEMA
    SCHEMARECOMB_AVAILABLE = True
except (ImportError, AttributeError) as e:
    SCHEMARECOMB_AVAILABLE = False

from utils.schemarecomb_bridge import (
    load_pdb_structure, align_pdb_to_parents, contacts_from_pdb_structure,
    msa_to_parent_sequences, sequences_to_parent_sequences, save_parents_to_fasta
)


def calculate_contacts(pdb_file, msa_file, chains=None, contact_distance=5.0, 
                      progress_callback=None):
    """
    Calculate SCHEMA contacts from PDB structure and MSA using schemarecomb.
    
    PDB alignment is handled automatically by schemarecomb's PDBStructure.renumber() method.
    
    Args:
        pdb_file: Path to PDB structure file
        msa_file: Path to multiple sequence alignment file
        chains: List of chain identifiers (default: ['A', ' '])
        contact_distance: Distance threshold for contacts in Angstroms (default: 5.0, schemarecomb uses 4.5)
        progress_callback: Optional callback for progress updates
        
    Returns:
        dict: Dictionary containing contacts and metadata
        
    Raises:
        ImportError: If schemarecomb is not available
        Exception: If contact calculation fails
    """
    if not SCHEMARECOMB_AVAILABLE:
        raise ImportError("schemarecomb is required for contact calculation. Install with: pip install schemarecomb")
    
    pdb_path = Path(pdb_file)
    msa_path = Path(msa_file)
    if not pdb_path.exists():
        raise FileNotFoundError(f"PDB file not found: {pdb_file}")
    if not msa_path.exists():
        raise FileNotFoundError(f"MSA file not found: {msa_file}")
    
    try:
        if progress_callback:
            progress_callback(0.1, "Loading PDB structure with schemarecomb...")
        
        # Determine chain
        chain = 'A'
        if chains:
            if isinstance(chains, list) and len(chains) > 0:
                chain = chains[0] if chains[0] != ' ' else 'A'
            elif isinstance(chains, str) and chains != ' ':
                chain = chains
        
        # Load PDB structure
        pdb_structure = load_pdb_structure(pdb_file, chain=chain, progress_callback=progress_callback)
        
        if progress_callback:
            progress_callback(0.4, "Loading parent sequences from MSA...")
        
        # Load parent sequences from MSA
        parents = msa_to_parent_sequences(msa_file, progress_callback=progress_callback)
        
        if progress_callback:
            progress_callback(0.6, "Aligning PDB structure to parent sequences...")
        
        # Align PDB to parents (trims alignment and renumbers to alignment positions)
        pdb_structure, trimmed_parents = align_pdb_to_parents(
            pdb_structure, parents, progress_callback=progress_callback
        )
        
        if progress_callback:
            progress_callback(0.8, "Extracting contacts...")
        
        # Extract contacts (schemarecomb uses 4.5A by default, but we'll use our threshold)
        # Note: schemarecomb calculates contacts at 4.5A, but we can filter if needed
        contacts_list = contacts_from_pdb_structure(pdb_structure)
        
        # Convert to our format: (i, j, ri, rj) where ri and rj are residue objects
        # For schemarecomb, we'll use indices only since we don't have residue objects
        contacts = []
        for i, j in contacts_list:
            contacts.append((i, j, None, None))  # ri, rj not available from schemarecomb
        
        # Get parent list from TRIMMED alignment (domain boundaries)
        parent_list = [(rec.id, str(rec.seq)) for rec in trimmed_parents.records]
        parent_dict = dict(parent_list)
        
        # Get pdb_key (first parent)
        pdb_key = trimmed_parents.records[0].id if trimmed_parents.records else None
        
        # Save trimmed alignment to a new file
        trimmed_msa_file = None
        try:
            msa_path = Path(msa_file)
            trimmed_msa_file = str(msa_path.parent / f"{msa_path.stem}_trimmed{msa_path.suffix}")
            save_parents_to_fasta(trimmed_parents, trimmed_msa_file)
        except Exception as e:
            # Fall back to original file if save fails
            trimmed_msa_file = msa_file
        
        if progress_callback:
            progress_callback(1.0, "Contacts calculated")
        
        return {
            'contacts': contacts,
            'residues': None,  # Not available from schemarecomb
            'parents': parent_list,
            'parent_dict': parent_dict,
            'pdb_key': pdb_key,
            'chain_identifiers': [chain],
            'pdb_structure': pdb_structure,  # Store for later use
            'parents_object': trimmed_parents,  # Store TRIMMED ParentSequences object
            'original_parents': parents,  # Keep original for reference if needed
            'trimmed_msa_file': trimmed_msa_file  # Path to trimmed alignment
        }
        
    except Exception as e:
        raise Exception(f"Error calculating contacts: {str(e)}")


def calculate_energies(msa_file, contact_file, crossover_file, 
                       chimeras=None, print_E=True, print_m=False,
                       parents_object=None, use_schemarecomb=True):
    """
    Calculate SCHEMA energies for chimeras.
    
    Uses schemarecomb if available, otherwise falls back to SCHEMA-RASPP.
    
    Args:
        msa_file: Path to multiple sequence alignment file
        contact_file: Path to contact file (from calculate_contacts) - optional if use_schemarecomb=True
        crossover_file: Path to crossover points file
        chimeras: Optional list of chimera block patterns (e.g., ['11213312'])
        print_E: Whether to calculate disruption energies (default: True)
        print_m: Whether to calculate mutation distances (default: False)
        parents_object: Optional schemarecomb.ParentSequences object (if already loaded)
        use_schemarecomb: Whether to use schemarecomb (default: True)
        
    Returns:
        dict: Dictionary containing energy results
    """
    # Use schemarecomb if available and requested
    if use_schemarecomb and SCHEMARECOMB_AVAILABLE:
        try:
            # Load parents if not provided
            if parents_object is None:
                parents_object = msa_to_parent_sequences(msa_file)
            
            # Create SCHEMA energy function (requires PDB structure)
            if not hasattr(parents_object, 'pdb_structure'):
                raise ValueError("ParentSequences must have pdb_structure for SCHEMA energy calculation")
            
            schema_energy = SCHEMA(parents_object)
            
            # Read crossovers to determine breakpoints
            with open(crossover_file, 'r') as f:
                crossover_lines = f.readlines()
            
            # Parse crossover positions (assuming format: one position per line)
            breakpoint_positions = []
            for line in crossover_lines:
                line = line.strip()
                if line and not line.startswith('#'):
                    try:
                        pos = int(line)
                        breakpoint_positions.append(pos)
                    except ValueError:
                        continue
            
            # Calculate energy for each chimera pattern
            results = []
            alignment_len = len(parents_object.alignment)
            
            if chimeras is None:
                # Generate all possible chimeras (if small enough)
                num_fragments = len(breakpoint_positions) + 1
                num_parents = len(parents_object.records)
                total_chimeras = num_parents ** num_fragments
                
                if total_chimeras > 10000:
                    raise ValueError(f"Too many possible chimeras ({total_chimeras}). Please specify chimeras explicitly.")
                
                # Generate all chimera patterns
                chimeras = []
                for i in range(total_chimeras):
                    pattern = ''
                    num = i
                    for _ in range(num_fragments):
                        pattern = str((num % num_parents)) + pattern
                        num //= num_parents
                    chimeras.append(pattern)
            
            for chimera_blocks in chimeras:
                # Calculate energy using schemarecomb's block_energy method
                # Convert chimera pattern to breakpoints
                total_energy = 0.0
                prev_pos = 0
                
                for i, parent_idx in enumerate(chimera_blocks):
                    parent_idx = int(parent_idx)
                    if i < len(breakpoint_positions):
                        curr_pos = breakpoint_positions[i]
                    else:
                        curr_pos = alignment_len
                    
                    # Energy of block from prev_pos to curr_pos
                    if prev_pos < curr_pos:
                        block_e = schema_energy.block_energy(prev_pos, curr_pos)
                        total_energy += float(block_e)
                    
                    prev_pos = curr_pos
                
                # Calculate mutation distance if requested
                m = None
                if print_m:
                    # Calculate Hamming distance to nearest parent
                    chimera_seq = ''
                    prev_pos = 0
                    for i, parent_idx in enumerate(chimera_blocks):
                        parent_idx = int(parent_idx)
                        if i < len(breakpoint_positions):
                            curr_pos = breakpoint_positions[i]
                        else:
                            curr_pos = alignment_len
                        
                        # Get block from parent
                        parent_aligned = ''.join([aminos[parent_idx] for aminos in parents_object.alignment[prev_pos:curr_pos]])
                        chimera_seq += parent_aligned
                        prev_pos = curr_pos
                    
                    # Find minimum distance to any parent
                    min_dist = float('inf')
                    for rec in parents_object.records:
                        parent_seq = ''.join([aminos[0] for aminos in parents_object.alignment])
                        dist = sum(1 for a, b in zip(chimera_seq, parent_seq) if a != b)
                        min_dist = min(min_dist, dist)
                    m = min_dist
                
                results.append({
                    'chimera': chimera_blocks,
                    'disruption': total_energy if print_E else None,
                    'mutations': m,
                    'sequence': None  # Can be generated if needed
                })
            
            return {
                'results': results,
                'fragments': None,  # Not used in schemarecomb
                'crossovers': breakpoint_positions,
                'parents': [(rec.id, str(rec.seq)) for rec in parents_object.records],
                'contacts': None,  # Not needed with schemarecomb
                'schema_energy': schema_energy  # Store for later use
            }
            
        except Exception as e:
            # Fall through to old implementation
            if not use_schemarecomb:
                raise
            # Try old implementation as fallback
            pass
    
    # Fallback to old SCHEMA-RASPP implementation
    try:
        # Read MSA
        with open(msa_file, 'r') as f:
            parent_list = schema.readMultipleSequenceAlignmentFile(f)
        parents = [seq for _, seq in parent_list]
        
        # Read contacts
        with open(contact_file, 'r') as f:
            contacts = schema.readContactFile(f)
        
        # Read crossovers
        with open(crossover_file, 'r') as f:
            crossovers = schema.readCrossoverFile(f)
        
        # Generate fragments from crossovers
        fragments = schema.getFragments(crossovers, parents[0])
        
        # Calculate energies
        results = []
        
        if chimeras is None:
            # Generate all possible chimeras (if small enough)
            num_fragments = len(fragments)
            num_parents = len(parents)
            total_chimeras = num_parents ** num_fragments
            
            if total_chimeras > 10000:
                raise ValueError(f"Too many possible chimeras ({total_chimeras}). Please specify chimeras explicitly.")
            
            # Generate all chimera patterns
            chimeras = []
            for i in range(total_chimeras):
                pattern = ''
                num = i
                for _ in range(num_fragments):
                    pattern = str((num % num_parents) + 1) + pattern
                    num //= num_parents
                chimeras.append(pattern)
        
        for chimera_blocks in chimeras:
            if not schema.checkChimera(chimera_blocks, fragments, parents):
                continue
            
            E = None
            m = None
            
            if print_E:
                E = schema.getChimeraDisruption(chimera_blocks, contacts, fragments, parents)
            
            if print_m:
                m = schema.getChimeraShortestDistance(chimera_blocks, fragments, parents)
            
            results.append({
                'chimera': chimera_blocks,
                'disruption': E,
                'mutations': m,
                'sequence': schema.getChimeraSequence(chimera_blocks, fragments, parents) if print_E or print_m else None
            })
        
        return {
            'results': results,
            'fragments': fragments,
            'crossovers': crossovers,
            'parents': parent_list,
            'contacts': contacts
        }
        
    except Exception as e:
        raise Exception(f"Error calculating energies: {str(e)}")


def save_contacts(contacts, output_file):
    """
    Save contacts to a file.
    
    Handles both old format (i, j, resi, resj) with residue objects
    and new format (i, j, None, None) from schemarecomb.
    
    Args:
        contacts: List of contact tuples (i, j) or (i, j, resi, resj)
        output_file: Path to output file
    """
    with open(output_file, 'w') as f:
        f.write("# Fields are number, contacting residues i & j (in library coordinates), and residues i & j (in PDB coordinates)\n")
        f.write("# n\ti\tj\tpdbi\tpdbj\n")
        for k, contact in enumerate(contacts):
            i = contact[0]
            j = contact[1]
            # Handle both old format (with residue objects) and new format (None)
            if len(contact) >= 4 and contact[2] is not None:
                # Old format with residue objects
                resi = contact[2]
                resj = contact[3]
                pdbi = resi.res_seq if hasattr(resi, 'res_seq') else str(i)
                pdbj = resj.res_seq if hasattr(resj, 'res_seq') else str(j)
            else:
                # New format - use indices as PDB coordinates
                pdbi = str(i)
                pdbj = str(j)
            f.write(f"{k}\t{i}\t{j}\t{pdbi}\t{pdbj}\n")
