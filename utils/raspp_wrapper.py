"""
Wrapper functions for RASPP algorithm.
Provides a Python API for the bundled raspp module and schemarecomb integration.
"""

import streamlit as st
import sys

try:
    from schema_raspp import raspp, schema
except ImportError:
    import raspp, schema

try:
    from schemarecomb import generate_libraries, ParentSequences
    from schemarecomb.energy_functions import SCHEMA
    from schemarecomb.optimizers import RASPP
    SCHEMARECOMB_AVAILABLE = True
except (ImportError, AttributeError):
    SCHEMARECOMB_AVAILABLE = False

from utils.schemarecomb_bridge import library_to_results, parent_sequences_from_json
from utils.config import SESSION_KEYS
from pathlib import Path


def _get_diagnostic_info(parents_object, parents, contacts, use_schemarecomb):
    """
    Generate diagnostic information for debugging RASPP errors.
    
    Args:
        parents_object: schemarecomb.ParentSequences object or None
        parents: List of parent sequences or None
        contacts: List of contact tuples or None
        use_schemarecomb: Whether schemarecomb is being used
        
    Returns:
        str: Diagnostic information string
    """
    lines = ["=== RASPP Diagnostic Information ==="]
    
    lines.append(f"use_schemarecomb: {use_schemarecomb}")
    lines.append(f"SCHEMARECOMB_AVAILABLE: {SCHEMARECOMB_AVAILABLE}")
    
    if parents_object is not None:
        lines.append(f"parents_object type: {type(parents_object)}")
        if hasattr(parents_object, 'records'):
            lines.append(f"parents_object.records length: {len(parents_object.records)}")
        if hasattr(parents_object, 'alignment'):
            try:
                lines.append(f"parents_object.alignment length: {len(parents_object.alignment)}")
            except AttributeError:
                lines.append("parents_object.alignment: Not aligned")
        if hasattr(parents_object, 'pdb_structure'):
            lines.append(f"parents_object.pdb_structure: {parents_object.pdb_structure is not None}")
        else:
            lines.append("parents_object.pdb_structure: Attribute not present")
    else:
        lines.append("parents_object: None")
    
    if parents is not None:
        lines.append(f"parents type: {type(parents)}")
        if isinstance(parents, (list, tuple)):
            lines.append(f"parents length: {len(parents)}")
            if len(parents) > 0:
                lines.append(f"parents[0] type: {type(parents[0])}")
        elif isinstance(parents, ParentSequences):
            lines.append(f"parents (ParentSequences) records length: {len(parents.records)}")
    else:
        lines.append("parents: None")
    
    if contacts is not None:
        lines.append(f"contacts type: {type(contacts)}")
        if isinstance(contacts, (list, tuple)):
            lines.append(f"contacts length: {len(contacts)}")
    else:
        lines.append("contacts: None")
    
    return "\n".join(lines)


def _parents_to_list(parents):
    """
    Normalize parents to a list of sequence strings.
    Accepts: dict (id -> seq), list of (name, seq) tuples, or list of strings.
    Returns: list of sequences. Raises ValueError for invalid format.
    """
    if parents is None:
        raise ValueError("parents cannot be None")
    if isinstance(parents, dict):
        return [seq for _, seq in parents.items()]
    if isinstance(parents, list):
        if not parents:
            raise ValueError("parents list cannot be empty")
        if isinstance(parents[0], tuple):
            return [seq for _, seq in parents]
        if isinstance(parents[0], str):
            return parents
    raise ValueError("parents must be dict, list of (name, seq) tuples, or list of sequence strings")


def run_raspp(contacts, parents, num_crossovers, min_fragment_diversity, progress_callback=None,
              parents_object=None, use_schemarecomb=False, allowed_fragment_starts=None):
    """
    Run RASPP algorithm to find optimal crossover points.
    
    Uses SCHEMA-RASPP by default (no Golden Gate calculations).
    Can optionally use schemarecomb if requested.
    
    Args:
        contacts: List of contact tuples (i, j, ri, rj) - required for SCHEMA-RASPP
        parents: List of parent sequences (strings) - required for SCHEMA-RASPP
        num_crossovers: Number of crossover points desired
        min_fragment_diversity: Minimum fragment diversity (length)
        progress_callback: Optional callback function for progress updates
        parents_object: Optional schemarecomb.ParentSequences object (only used if use_schemarecomb=True)
        use_schemarecomb: Whether to use schemarecomb (default: False - uses SCHEMA-RASPP)
        
    Returns:
        list: List of tuples (SCHEMA-RASPP) or Library objects (schemarecomb)
    """
    # Use schemarecomb if available and requested
    if use_schemarecomb and SCHEMARECOMB_AVAILABLE:
        try:
            if progress_callback:
                progress_callback(0.1, "Preparing RASPP calculation with schemarecomb...")
            
            # Get ParentSequences object
            if parents_object is None:
                # Try to load from parents if it's a ParentSequences object
                if isinstance(parents, ParentSequences):
                    parents_object = parents
                elif isinstance(parents, str):
                    # Assume it's a JSON string
                    parents_object = parent_sequences_from_json(parents)
                else:
                    diag_info = _get_diagnostic_info(parents_object, parents, contacts, use_schemarecomb)
                    raise ValueError(
                        f"parents_object or valid parents must be provided for schemarecomb.\n\n"
                        f"{diag_info}"
                    )
            
            # Ensure parents have PDB structure
            if not hasattr(parents_object, 'pdb_structure') or parents_object.pdb_structure is None:
                diag_info = _get_diagnostic_info(parents_object, parents, contacts, use_schemarecomb)
                raise ValueError(
                    f"ParentSequences must have pdb_structure for RASPP optimization.\n"
                    f"The pdb_structure attribute is missing or None.\n\n"
                    f"{diag_info}\n\n"
                    f"SOLUTION: Ensure that:\n"
                    f"1. A PDB file was used when calculating SCHEMA contacts\n"
                    f"2. The PDB file path is available in session state\n"
                    f"3. The ParentSequences object was created with the PDB structure attached\n"
                    f"4. If loading from a checkpoint, the PDB structure may need to be re-attached"
                )
            
            if progress_callback:
                progress_callback(0.3, f"Generating libraries with {num_crossovers} crossovers...")
            
            # Calculate number of blocks (crossovers + 1)
            num_blocks = num_crossovers + 1
            
            # Validate num_blocks (schemarecomb requires num_blocks >= 2)
            if num_blocks < 2:
                raise ValueError(
                    f"Number of blocks must be at least 2 (got {num_blocks} from {num_crossovers} crossovers). "
                    f"Please use at least 1 crossover."
                )
            
            # Calculate block length constraints
            alignment_len = len(parents_object.alignment)
            min_block_len = min_fragment_diversity
            
            # Validate that the alignment is long enough for the requested configuration
            min_required_length = num_blocks * min_block_len
            if alignment_len < min_required_length:
                raise ValueError(
                    f"Cannot generate {num_blocks} blocks with min_block_len={min_block_len}. "
                    f"Alignment length ({alignment_len}) is too short. "
                    f"Minimum required length: {min_required_length}. "
                    f"Try reducing num_crossovers or min_fragment_diversity."
                )
            
            # Don't set max_block_len - let schemarecomb use its default calculation
            # which is: alignment_len // (num_blocks - 1)
            # This allows the optimize() method to explore a wider range of block lengths
            # and find multiple library options (top 10-20), rather than being constrained
            # to a narrow range that only finds 1 solution.
            # 
            # The default calculation is designed to allow schemarecomb to find multiple
            # optimal libraries by testing different min/max block length combinations.
            max_block_len = None
            
            if progress_callback:
                progress_callback(0.5, f"Optimizing libraries (min_block_len={min_fragment_diversity}, max_block_len=default, num_blocks={num_blocks})...")
            
            # Generate libraries using schemarecomb
            # Let schemarecomb use its default max_block_len calculation for better results
            library_kwargs = {
                'parents': parents_object,
                'num_blocks': num_blocks,
                'min_block_len': min_fragment_diversity,
                'algorithm': 'SCHEMA-RASPP'
                # Note: Not setting max_block_len - let schemarecomb use default
            }
            
            try:
                libraries = generate_libraries(**library_kwargs)
            except ValueError as e:
                # Catch specific errors from schemarecomb and provide better context
                error_msg = str(e)
                if "Length of breakpoints must be greater than 1" in error_msg:
                    # This error occurs when schemarecomb's optimizer finds a library, but after removing
                    # the first breakpoint (if it's at position 0 with no overhangs), there are <= 1 breakpoints left.
                    # This can happen when:
                    # 1. Constraints are too restrictive and the optimizer finds a degenerate solution
                    # 2. The optimizer finds only 2 breakpoints (for 3 blocks) but one gets removed
                    # 3. The optimizer finds a solution that doesn't satisfy the num_blocks requirement
                    raise ValueError(
                        f"RASPP optimization failed: schemarecomb generated a library with too few breakpoints.\n"
                        f"This happens when the optimizer finds a solution with fewer breakpoints than required,\n"
                        f"and then the first breakpoint (at position 0) is removed, leaving <= 1 breakpoint.\n\n"
                        f"Requested: {num_blocks} blocks (requires {num_blocks-1} internal breakpoints)\n"
                        f"Parameters: min_block_len={min_fragment_diversity}, max_block_len={max_block_len}, "
                        f"alignment_len={alignment_len}\n\n"
                        f"SOLUTION: The constraints are likely too restrictive. Try:\n"
                        f"  - Reducing num_crossovers from {num_crossovers} to {max(1, num_crossovers-1)} or lower\n"
                        f"  - Reducing min_fragment_diversity from {min_fragment_diversity} to {max(1, min_fragment_diversity-1)} or lower\n"
                        f"  - The alignment may be too short for {num_blocks} blocks with these constraints"
                    ) from e
                else:
                    # Re-raise other ValueError exceptions as-is
                    raise
            
            if progress_callback:
                progress_callback(1.0, f"RASPP calculation complete! Found {len(libraries)} libraries")
            
            return libraries
            
        except Exception as e:
            if progress_callback:
                diag_info = _get_diagnostic_info(parents_object, parents, contacts, use_schemarecomb)
                progress_callback(0.0, f"Error with schemarecomb: {str(e)}\n\nDiagnostics:\n{diag_info}")
            # Re-raise with diagnostic info if it's not already included
            if "=== RASPP Diagnostic Information ===" not in str(e):
                diag_info = _get_diagnostic_info(parents_object, parents, contacts, use_schemarecomb)
                raise Exception(f"Error with schemarecomb: {str(e)}\n\n{diag_info}")
            raise
    
    # Fallback to old SCHEMA-RASPP implementation
    if not use_schemarecomb or not SCHEMARECOMB_AVAILABLE:
        try:
            if progress_callback:
                progress_callback(0.1, "Preparing RASPP calculation...")
            
            # Convert parents to list if needed
            if isinstance(parents, dict):
                parents_list = [seq for _, seq in parents.items()]
            elif isinstance(parents, list) and isinstance(parents[0], tuple):
                # If it's a list of (name, sequence) tuples
                parents_list = [seq for _, seq in parents]
            elif isinstance(parents, list):
                parents_list = parents
            else:
                raise ValueError("Invalid parents format")
            
            # Ensure all parents have the same length
            parent_lengths = [len(p) for p in parents_list]
            if len(set(parent_lengths)) > 1:
                raise ValueError("All parent sequences must have the same length")
            
            if progress_callback:
                progress_callback(0.3, "Calculating average energies...")
            
            # Run RASPP_SCHEMA which handles the full workflow
            results = raspp.RASPP_SCHEMA(
                contacts,
                parents_list,
                num_crossovers,
                min_fragment_diversity,
                allowed_fragment_starts=allowed_fragment_starts,
            )
            
            if progress_callback:
                progress_callback(1.0, "RASPP calculation complete!")
            
            return results
            
        except Exception as e:
            raise Exception(f"Error running RASPP: {str(e)}")


def calculate_average_energies_from_contacts(contacts, parents):
    """
    Calculate average energies from contacts.
    
    Args:
        contacts: List of contact tuples
        parents: List of parent sequences
        
    Returns:
        list: List of average energy tuples (i, j, avg)
    """
    try:
        parents = _parents_to_list(parents)
        # Calculate average energies
        avg_energy_matrix = raspp.calc_average_energies_from_contacts(contacts, parents)
        
        # Convert matrix to list format
        num_residues = len(parents[0])
        avg_energies = []
        for i in range(num_residues - 1):
            for j in range(i + 1, num_residues):
                if avg_energy_matrix[i][j] > 0:
                    avg_energies.append((i, j, avg_energy_matrix[i][j]))
        
        return avg_energies
        
    except Exception as e:
        raise Exception(f"Error calculating average energies: {str(e)}")


def format_raspp_results(results):
    """
    Format RASPP results for display.
    
    Args:
        results: List of RASPP result tuples
        
    Returns:
        list: List of dictionaries with formatted results
    """
    formatted = []
    for result in results:
        avg_energy, crossovers, l_min, l_max = result
        # Ensure avg_energy is a valid number (convert to float if needed)
        if avg_energy is not None:
            try:
                avg_energy = float(avg_energy)
            except (ValueError, TypeError):
                # If conversion fails, set to None
                avg_energy = None
        formatted.append({
            'average_energy': avg_energy,
            'crossovers': crossovers,
            'min_fragment_length': l_min,
            'max_fragment_length': l_max,
            'num_fragments': len(crossovers) + 1
        })
    return formatted


def streamlit_progress_callback(progress_bar, status_text):
    """
    Create a progress callback function for Streamlit.
    
    Args:
        progress_bar: Streamlit progress bar object
        status_text: Streamlit status text object
        
    Returns:
        function: Progress callback function
    """
    def callback(progress, message):
        if progress_bar:
            progress_bar.progress(progress)
        if status_text:
            status_text.text(message)
    
    return callback
