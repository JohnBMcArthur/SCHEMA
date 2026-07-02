"""
Utilities for testing multiple fragment configurations.
"""

from utils.raspp_wrapper import run_raspp, format_raspp_results
import os
import streamlit as st
import random
import copy

try:
    from schemarecomb import ParentSequences
    from utils.schemarecomb_bridge import library_to_results
    SCHEMARECOMB_AVAILABLE = True
except (ImportError, AttributeError):
    SCHEMARECOMB_AVAILABLE = False
    ParentSequences = None
    library_to_results = None

try:
    from schema_raspp.raspp import RASPP_prepare, RASPP_SCHEMA_from_prepared, collapse_parents
except ImportError:
    try:
        from raspp import RASPP_prepare, RASPP_SCHEMA_from_prepared, collapse_parents
    except ImportError:
        RASPP_prepare = None
        RASPP_SCHEMA_from_prepared = None
        collapse_parents = None

try:
    from utils.visualization import allowed_collapsed_starts_from_graph_segment_marks
except ImportError:
    allowed_collapsed_starts_from_graph_segment_marks = None


def _create_parent_subset(parents_object, indices):
    """
    Create a new ParentSequences object from a subset of parent records.
    
    Args:
        parents_object: schemarecomb.ParentSequences object
        indices: List of indices to include in the subset
        
    Returns:
        schemarecomb.ParentSequences object with subset of records
    """
    if not SCHEMARECOMB_AVAILABLE:
        raise ImportError("schemarecomb is required")
    
    # Get subset of records
    subset_records = [copy.deepcopy(parents_object.records[i]) for i in indices]
    
    # Get PDB structure if available (shared reference is fine)
    pdb_structure = None
    if hasattr(parents_object, 'pdb_structure') and parents_object.pdb_structure:
        pdb_structure = parents_object.pdb_structure
    
    # Create new ParentSequences with subset
    subset_parents = ParentSequences(
        subset_records,
        pdb_structure=pdb_structure,
        prealigned=False  # We'll set alignment manually
    )
    
    # Set alignment to subset of original alignment
    if hasattr(parents_object, 'alignment') and parents_object.alignment:
        # alignment is list of tuples, one per position: [(aa1, aa2, ...), (aa1, aa2, ...), ...]
        # Extract subset for each position
        subset_alignment = [tuple(parents_object.alignment[pos][i] for i in indices) 
                           for pos in range(len(parents_object.alignment))]
        # Set the alignment
        subset_parents.new_alignment([''.join(pos) for pos in subset_alignment])
    
    return subset_parents


def _average_results(results_list):
    """
    Average results from multiple RASPP runs.
    
    For each set, takes the best energy, then averages those best energies.
    Also combines all results for display.
    
    Args:
        results_list: List of result dictionaries from format_raspp_results
        
    Returns:
        dict: Averaged results
    """
    if not results_list:
        return None
    
    # Collect best energies from each set
    best_energies_per_set = []
    all_formatted = []
    total_designs = 0
    successful_sets = 0
    
    for result_dict in results_list:
        if 'error' in result_dict:
            continue
        
        successful_sets += 1
        formatted = result_dict.get('results', [])
        if formatted:
            all_formatted.extend(formatted)
            # Get best energy from this set
            energies = []
            for r in formatted:
                # Check both keys explicitly
                energy = r.get('energy') if 'energy' in r else r.get('average_energy')
                if energy is not None:
                    try:
                        energies.append(float(energy))
                    except (ValueError, TypeError):
                        pass
            if energies:
                best_energies_per_set.append(min(energies))
            total_designs += len(formatted)
    
    if not best_energies_per_set:
        return None
    
    # Average the best energies from each set
    avg_best_energy = sum(best_energies_per_set) / len(best_energies_per_set)
    
    # Also calculate overall statistics from all combined results
    all_energies = []
    for r in all_formatted:
        energy = r.get('energy') if 'energy' in r else r.get('average_energy')
        if energy is not None:
            try:
                all_energies.append(float(energy))
            except (ValueError, TypeError):
                pass
    overall_best_energy = min(all_energies) if all_energies else None
    overall_avg_energy = sum(all_energies) / len(all_energies) if all_energies else None
    
    # Return averaged results
    return {
        'results': all_formatted,  # Combine all results for display
        'raw_results': [r.get('raw_results', []) for r in results_list],
        'best_energy': avg_best_energy,  # Average of best energies from each set
        'overall_best_energy': overall_best_energy,  # Best from all combined
        'avg_energy': overall_avg_energy,
        'num_designs': total_designs,
        'num_sets_tested': successful_sets,
        'best_energies_per_set': best_energies_per_set  # For reference
    }


def test_multiple_fragments(contacts, parents, min_fragments=5, max_fragments=20, 
                            min_fragment_diversity=8, progress_callback=None,
                            parents_object=None, use_schemarecomb=False,
                            results_update_callback=None,
                            restrict_crossovers_at_graph_marks=False,
                            segment_conservation_min_pct=80.0,
                            segment_diversity_max=3):
    """
    Test RASPP with multiple fragment counts.
    
    If there are exactly 10 parents, randomly creates 3 sets of 3 parents each
    and averages results across sets for each fragment count.
    
    Args:
        contacts: List of contact tuples (optional if use_schemarecomb=True)
        parents: List of parent sequences or ParentSequences object
        min_fragments: Minimum number of fragments to test (default: 5)
        max_fragments: Maximum number of fragments to test (default: 20)
        min_fragment_diversity: Minimum fragment diversity
        progress_callback: Optional callback for progress updates
        parents_object: Optional schemarecomb.ParentSequences object
        use_schemarecomb: Whether to use schemarecomb (default: True)
        results_update_callback: Optional callback to update results after each fragment
                                Called with (all_results_dict, num_fragments)
        restrict_crossovers_at_graph_marks: If True, only allow fragment boundaries whose
            full SCHEMA crossover matches a segment mark (same rules as the conservation/diversity
            graph: consecutive columns both above conservation A and below diversity B).
        segment_conservation_min_pct: A (%) — both columns must exceed this (strict).
        segment_diversity_max: B — both columns must have diversity strictly less than this.
        
    Returns:
        dict: Dictionary with fragment count as key, results as value
    """
    all_results = {}
    total_tests = max_fragments - min_fragments + 1
    
    # Check if we have 10 parents and should use subset averaging
    # Note: For now, subset averaging is only enabled for schemarecomb-based runs.
    use_subset_averaging = False
    parent_subsets = None
    
    if use_schemarecomb and parents_object and SCHEMARECOMB_AVAILABLE:
        num_parents = len(parents_object.records)
        if num_parents == 10:
            use_subset_averaging = True
            # Randomly create 3 sets of 3 parents each
            all_indices = list(range(10))
            random.shuffle(all_indices)
            parent_subsets = [
                all_indices[0:3],
                all_indices[3:6],
                all_indices[6:9]
            ]
            if progress_callback:
                progress_callback(0.0, f"Detected 10 parents. Using 3 random subsets of 3 parents each: {parent_subsets}")
    # For SCHEMA-RASPP (use_schemarecomb=False), always use all parents for now.
    # This avoids non-deterministic behavior and simplifies interpretation of results.
    
    # Debug: log the range being tested
    if progress_callback and not use_subset_averaging:
        progress_callback(0.0, f"Starting multi-fragment test: {min_fragments} to {max_fragments} fragments (total: {total_tests} tests)")

    # One-time prep for SCHEMA-RASPP: contacts → arc lengths (reused for each fragment count)
    multifrag_prepared = None
    if (
        not use_schemarecomb
        and not use_subset_averaging
        and RASPP_prepare is not None
    ):
        multifrag_prepared = RASPP_prepare(contacts, parents)

    segment_mark_allowed_starts = None
    if (
        restrict_crossovers_at_graph_marks
        and not use_schemarecomb
        and multifrag_prepared is not None
        and allowed_collapsed_starts_from_graph_segment_marks is not None
    ):
        segment_mark_allowed_starts = allowed_collapsed_starts_from_graph_segment_marks(
            parents,
            multifrag_prepared.collapsed_parents,
            multifrag_prepared.identical_sites,
            segment_conservation_min_pct,
            segment_diversity_max,
        )
        if progress_callback and not segment_mark_allowed_starts:
            progress_callback(
                0.0,
                "Graph segment marks: no eligible crossover positions for these "
                "conservation/diversity thresholds (see MSA plot sliders).",
            )
    
    for num_fragments in range(min_fragments, max_fragments + 1):
        num_crossovers = num_fragments - 1
        num_blocks = num_fragments  # fragments = blocks
        
        if progress_callback:
            progress = (num_fragments - min_fragments) / total_tests
            if use_subset_averaging:
                progress_callback(progress, f"Testing {num_fragments} fragments (averaging across 3 parent sets)...")
            else:
                progress_callback(progress, f"Testing {num_fragments} fragments ({num_crossovers} crossovers, {num_blocks} blocks)...")
        
        try:
            if use_subset_averaging:
                # Test each of the 3 parent subsets and average results
                subset_results = []
                for subset_idx, subset_indices in enumerate(parent_subsets):
                    if progress_callback:
                        progress_callback(progress, f"Testing {num_fragments} fragments with parent set {subset_idx + 1}/3 (parents: {subset_indices})...")
                    
                    # Run RASPP on this subset
                    try:
                        if use_schemarecomb and parents_object:
                            # Create subset ParentSequences for schemarecomb
                            subset_parents_obj = _create_parent_subset(parents_object, subset_indices)
                            subset_raspp_results = run_raspp(
                                contacts=None,
                                parents=None,
                                num_crossovers=num_crossovers,
                                min_fragment_diversity=min_fragment_diversity,
                                progress_callback=None,
                                parents_object=subset_parents_obj,
                                use_schemarecomb=True
                            )
                            # Format results (schemarecomb returns Library objects)
                            formatted = []
                            if subset_raspp_results:
                                for lib in subset_raspp_results:
                                    lib_result = library_to_results(lib)
                                    formatted.append({
                                        'energy': lib_result['energy'],
                                        'average_energy': lib_result['energy'],
                                        'crossovers': lib_result['breakpoints'],
                                        'min_fragment_length': lib_result['min_block_len'],
                                        'max_fragment_length': lib_result['max_block_len'],
                                        'num_fragments': len(lib_result['breakpoints']) + 1
                                    })
                        else:
                            # SCHEMA-RASPP: create subset from list of parent sequences
                            subset_parents_list = [parents[i] for i in subset_indices]
                            subset_allowed = None
                            if (
                                restrict_crossovers_at_graph_marks
                                and allowed_collapsed_starts_from_graph_segment_marks is not None
                                and collapse_parents is not None
                            ):
                                collapsed_sub, id_sub = collapse_parents(subset_parents_list)
                                subset_allowed = allowed_collapsed_starts_from_graph_segment_marks(
                                    subset_parents_list,
                                    collapsed_sub,
                                    id_sub,
                                    segment_conservation_min_pct,
                                    segment_diversity_max,
                                )
                            subset_raspp_results = run_raspp(
                                contacts=contacts,  # Required for SCHEMA-RASPP
                                parents=subset_parents_list,  # Subset of parent sequences
                                num_crossovers=num_crossovers,
                                min_fragment_diversity=min_fragment_diversity,
                                progress_callback=None,
                                parents_object=None,
                                use_schemarecomb=False,
                                allowed_fragment_starts=subset_allowed,
                            )
                            # Format results (SCHEMA-RASPP returns tuples)
                            formatted = format_raspp_results(subset_raspp_results)
                        
                        # Get best energy
                        best_energy = None
                        if formatted:
                            # Extract energies - check both 'energy' and 'average_energy' keys
                            energies = []
                            for r in formatted:
                                energy = r.get('energy') if 'energy' in r else r.get('average_energy')
                                # Only filter out None, not 0 (0 is a valid energy value)
                                if energy is not None:
                                    try:
                                        energies.append(float(energy))
                                    except (ValueError, TypeError):
                                        # Skip invalid energy values
                                        pass
                            if energies:
                                best_energy = min(energies)
                        
                        subset_results.append({
                            'num_crossovers': num_crossovers,
                            'num_blocks': num_blocks,
                            'results': formatted,
                            'raw_results': subset_raspp_results,
                            'best_energy': best_energy,
                            'num_designs': len(formatted),
                            'parent_indices': subset_indices
                        })
                    except Exception as subset_error:
                        error_msg = str(subset_error)
                        subset_results.append({
                            'error': error_msg,
                            'num_crossovers': num_crossovers,
                            'num_designs': 0,
                            'num_blocks': num_blocks,
                            'parent_indices': subset_indices
                        })
                
                # Average results across subsets
                averaged = _average_results(subset_results)
                if averaged:
                    all_results[num_fragments] = {
                        'num_crossovers': num_crossovers,
                        'num_blocks': num_blocks,
                        'results': averaged['results'],
                        'raw_results': averaged['raw_results'],
                        'best_energy': averaged['best_energy'],
                        'avg_energy': averaged.get('avg_energy'),
                        'num_designs': averaged['num_designs'],
                        'num_sets_tested': averaged['num_sets_tested'],
                        'subset_results': subset_results  # Keep individual results for reference
                    }
                else:
                    # All subsets failed
                    all_results[num_fragments] = {
                        'error': 'All 3 parent subsets failed',
                        'num_crossovers': num_crossovers,
                        'num_designs': 0,
                        'num_blocks': num_blocks,
                        'subset_results': subset_results
                    }
            else:
                # SCHEMA-RASPP: reuse prepared arc data across fragment counts
                if multifrag_prepared is not None:
                    results = RASPP_SCHEMA_from_prepared(
                        multifrag_prepared,
                        num_crossovers,
                        min_fragment_diversity,
                        allowed_fragment_starts=segment_mark_allowed_starts,
                    )
                else:
                    results = run_raspp(
                        contacts=contacts,
                        parents=parents,
                        num_crossovers=num_crossovers,
                        min_fragment_diversity=min_fragment_diversity,
                        progress_callback=None,
                        parents_object=None,
                        use_schemarecomb=use_schemarecomb,
                        allowed_fragment_starts=segment_mark_allowed_starts,
                    )
                
                # Format results (schemarecomb returns Library objects, fallback returns tuples)
                if use_schemarecomb and SCHEMARECOMB_AVAILABLE and results and hasattr(results[0], 'energy'):
                    # Convert Library objects to our format
                    formatted = []
                    for lib in results:
                        lib_result = library_to_results(lib)
                        formatted.append({
                            'energy': lib_result['energy'],
                            'average_energy': lib_result['energy'],
                            'crossovers': lib_result['breakpoints'],
                            'min_fragment_length': lib_result['min_block_len'],
                            'max_fragment_length': lib_result['max_block_len'],
                            'num_fragments': len(lib_result['breakpoints']) + 1
                        })
                else:
                    # Use standard formatting for tuples
                    formatted = format_raspp_results(results)

                # Get best energy (handle both schemarecomb and old format)
                best_energy = None
                if formatted:
                    # Extract energies - check both 'energy' and 'average_energy' keys
                    energies = []
                    for r in formatted:
                        # Try 'average_energy' first (SCHEMA-RASPP), then 'energy' (schemarecomb)
                        energy = r.get('average_energy') or r.get('energy')
                        # Only filter out None, not 0 (0 is a valid energy value)
                        if energy is not None:
                            try:
                                energy_val = float(energy)
                                energies.append(energy_val)
                            except (ValueError, TypeError):
                                # Skip invalid energy values (None, non-numeric strings, etc.)
                                pass
                    if energies:
                        best_energy = min(energies)
                    else:
                        # Fallback: try to extract energy from raw results if formatted has no valid energies
                        if len(formatted) > 0 and results and len(results) > 0:
                            raw_sample = results[0]
                            if isinstance(raw_sample, tuple) and len(raw_sample) >= 1:
                                raw_energy = raw_sample[0]
                                if raw_energy is not None:
                                    try:
                                        best_energy = float(raw_energy)
                                    except (ValueError, TypeError):
                                        pass
                        if os.environ.get("SCHEMA_DEBUG") and best_energy is None and formatted:
                            import sys
                            sample = formatted[0] if formatted else {}
                            raw_sample = results[0] if results else None
                            print("DEBUG: No valid energies found.", file=sys.stderr)
                            print(f"  Sample keys: {list(sample.keys())}", file=sys.stderr)
                            print(f"  Raw tuple (first 2): {raw_sample[:2] if isinstance(raw_sample, tuple) else raw_sample}", file=sys.stderr)
                
                # Fallback: If best_energy is None but we have formatted results, try to calculate it again
                if best_energy is None and formatted:
                    # Try one more time with a different approach
                    all_energies = []
                    for r in formatted:
                        # Direct access to average_energy (SCHEMA-RASPP format)
                        if 'average_energy' in r:
                            val = r['average_energy']
                            if val is not None:
                                try:
                                    all_energies.append(float(val))
                                except (ValueError, TypeError):
                                    pass
                        # Also check energy key (schemarecomb format)
                        elif 'energy' in r:
                            val = r['energy']
                            if val is not None:
                                try:
                                    all_energies.append(float(val))
                                except (ValueError, TypeError):
                                    pass
                    if all_energies:
                        best_energy = min(all_energies)
                
                all_results[num_fragments] = {
                    'num_crossovers': num_crossovers,
                    'num_blocks': num_blocks,
                    'results': formatted,
                    'raw_results': results,
                    'best_energy': best_energy,
                    'num_designs': len(formatted)
                }
                
                if len(formatted) > 0 and best_energy is None:
                    try:
                        st.warning(f"⚠️ {num_fragments} fragments: Found {len(formatted)} designs but couldn't extract energy.")
                    except Exception:
                        pass

            if progress_callback:
                if use_subset_averaging:
                    num_designs = all_results[num_fragments].get('num_designs', 0)
                    num_sets = all_results[num_fragments].get('num_sets_tested', 0)
                    progress_callback(progress, f"✓ {num_fragments} fragments: {num_designs} designs found (averaged from {num_sets}/3 parent sets)")
                else:
                    num_designs = all_results[num_fragments].get('num_designs', 0)
                    progress_callback(progress, f"✓ {num_fragments} fragments: {num_designs} designs found")
            
            # Update results incrementally after each successful fragment
            if results_update_callback:
                results_update_callback(all_results.copy(), num_fragments)
            
        except Exception as e:
            error_msg = str(e)
            # Include more detail in error message
            if progress_callback:
                progress_callback(progress, f"⚠️ Failed {num_fragments} fragments: {error_msg[:50]}...")
            all_results[num_fragments] = {
                'error': error_msg,
                'num_crossovers': num_crossovers,
                'num_designs': 0,
                'num_blocks': num_blocks
            }
            
            # Update results incrementally after each failed fragment too
            if results_update_callback:
                results_update_callback(all_results.copy(), num_fragments)
    
    if progress_callback:
        progress_callback(1.0, "All fragment tests complete!")

    return all_results
