"""
Page 2: RASPP Library Design
Design optimal recombination libraries using the RASPP algorithm.
"""

import streamlit as st
import tempfile
from pathlib import Path

# Import utilities
import sys
sys.path.insert(0, str(Path(__file__).parent.parent))

from utils.session_manager import init_session_state, has_required_data
from utils.config import DEFAULTS, SESSION_KEYS

from utils.file_handlers import validate_msa_file, save_uploaded_file
from utils.raspp_wrapper import streamlit_progress_callback
from utils.visualization import (
    plot_crossover_comparison,
    plot_fragment_lengths,
    plot_msa_conservation_diversity,
    compute_msa_column_conservation_diversity,
    compute_consensus_segment_midpoints,
)
from schema_raspp import schema

st.set_page_config(
    page_title="RASPP Design",
    page_icon="🔬",
    layout="wide"
)

# Initialize session state
init_session_state()

st.title("🔬 RASPP Library Design")

st.markdown("""
Design optimal protein recombination libraries using the RASPP algorithm. 
Provide SCHEMA energy results and library constraints to generate optimal 
crossover designs.

**Workflow:**
1. Use SCHEMA contacts from Page 1, or upload MSA and contact files
2. Set fragment constraints for multi-fragment testing
3. Run multi-fragment analysis
4. View optimal crossover designs
""")

# Check if we have SCHEMA contacts from previous page
has_data, error_msg = has_required_data('raspp')
has_schema_contacts = has_data  # Alias for clarity

# Initialize parents_for_fallback
parents_for_fallback = []

if has_schema_contacts:
    st.success("✓ Using SCHEMA contacts from previous calculation")
    contacts = st.session_state[SESSION_KEYS['schema_contacts']]['contacts']
    parents_list = st.session_state[SESSION_KEYS['schema_contacts']].get('parents', [])
    msa_path = st.session_state[SESSION_KEYS['msa_path']]
    
    # Extract parents in correct format
    if parents_list and isinstance(parents_list, list) and len(parents_list) > 0:
        if isinstance(parents_list[0], tuple) and len(parents_list[0]) == 2:
            # List of (name, sequence) tuples
            parents_for_fallback = [seq for _, seq in parents_list]
        elif isinstance(parents_list[0], str):
            # Already a list of sequences
            parents_for_fallback = parents_list
        else:
            # Try to extract sequences
            parents_for_fallback = [seq if isinstance(seq, str) else str(seq) for seq in parents_list]
    else:
        parents_for_fallback = []
    
    # Try to get ParentSequences object from session state
    parents_object = None
    if 'parents_object_json' in st.session_state and st.session_state['parents_object_json']:
        from utils.schemarecomb_bridge import parent_sequences_from_json
        try:
            parents_object = parent_sequences_from_json(st.session_state['parents_object_json'])
            
            # Check if PDB structure is missing and attach it if available
            if parents_object and (not hasattr(parents_object, 'pdb_structure') or parents_object.pdb_structure is None):
                if SESSION_KEYS['pdb_path'] in st.session_state and st.session_state[SESSION_KEYS['pdb_path']]:
                    from utils.schemarecomb_bridge import load_pdb_structure, align_pdb_to_parents
                    try:
                        pdb_file = st.session_state[SESSION_KEYS['pdb_path']]
                        pdb_path_obj = Path(pdb_file)
                        if pdb_path_obj.exists():
                            pdb_structure = load_pdb_structure(pdb_file, chain='A')
                            # align_pdb_to_parents returns (pdb_structure, trimmed_parents) tuple
                            aligned_pdb, _ = align_pdb_to_parents(pdb_structure, parents_object)
                            parents_object.pdb_structure = aligned_pdb
                            st.success("✓ Attached PDB structure to loaded ParentSequences")
                        else:
                            st.warning(f"⚠ PDB file path exists in session but file not found: {pdb_file}")
                    except Exception as pdb_error:
                        st.warning(f"⚠ Could not attach PDB structure to loaded ParentSequences: {str(pdb_error)}")
                        import traceback
                        with st.expander("PDB attachment error details"):
                            st.code(traceback.format_exc())
                else:
                    st.warning("⚠ ParentSequences loaded from JSON but no PDB path in session state. PDB structure may be missing.")
        except Exception as e:
            st.warning(f"Could not load ParentSequences object: {str(e)}")
            import traceback
            with st.expander("ParentSequences loading error details"):
                st.code(traceback.format_exc())
else:
    st.info(error_msg or "Upload MSA and contact files, or go to Page 1 to calculate SCHEMA contacts first.")
    parents_object = None
    parents_for_fallback = []

# Sidebar for parameters
with st.sidebar:
    st.header("RASPP Parameters")

    min_fragment_diversity = st.number_input(
        "Minimum Fragment Length",
        min_value=1,
        max_value=50,
        value=DEFAULTS['min_fragment_diversity'],
        step=1,
        help="Minimum block length (in residues) allowed in generated libraries"
    )
    
    min_fragments = st.slider(
        "Minimum fragments",
        min_value=2,
        max_value=100,
        value=DEFAULTS['min_fragments']
    )

    max_fragments = st.slider(
        "Maximum fragments",
        min_value=2,
        max_value=100,
        value=DEFAULTS['max_fragments']
    )

    st.markdown("---")
    st.markdown("### Project Management")
    
    # Show current project if loaded
    if 'current_project' in st.session_state and st.session_state.get('current_project'):
        st.info(f"**📁 {st.session_state['current_project']}**")
    
    # Store parameters in session state for saving
    st.session_state['min_fragment_diversity'] = min_fragment_diversity
    st.session_state['min_fragments'] = min_fragments
    st.session_state['max_fragments'] = max_fragments
    # Load project option
    from utils.session_manager import list_checkpoints, load_checkpoint
    
    checkpoints = list_checkpoints()
    if checkpoints:
        checkpoint_names = [f"{meta.get('project_name', 'Unknown')}" 
                           for path, meta in checkpoints[:5]]
        checkpoint_paths = [path for path, meta in checkpoints[:5]]
        
        selected_idx = st.selectbox(
            "Load Project",
            range(len(checkpoint_names)),
            format_func=lambda x: checkpoint_names[x],
            key='page2_load_project'
        )
        
        if st.button("📂 Load", use_container_width=True, key='page2_load_btn'):
            try:
                metadata = load_checkpoint(checkpoint_paths[selected_idx])
                st.success(f"✓ Loaded: {metadata.get('project_name', 'Unknown')}")
                st.rerun()
            except Exception as e:
                st.error(f"Error: {str(e)}")
    
    st.markdown("---")
    st.markdown("### Input Files")
    
    if not has_schema_contacts:
        msa_file = st.file_uploader(
            "MSA File",
            type=['txt', 'msa', 'fasta', 'fa', 'aln'],
            help="Multiple sequence alignment file"
        )
        
        contact_file = st.file_uploader(
            "Contact File",
            type=['txt'],
            help="SCHEMA contacts file"
        )
    else:
        msa_file = None
        contact_file = None

st.markdown("---")
st.subheader("MSA column conservation and diversity")
st.caption(
    "Conservation % is the frequency of the most abundant residue in each column; "
    "diversity score is the number of distinct residues (including gap characters)."
)

_msa_plot_sequences = None
if has_schema_contacts and parents_for_fallback:
    _msa_plot_sequences = parents_for_fallback
elif not has_schema_contacts and msa_file:
    try:
        import io

        msa_file.seek(0)
        parsed = schema.readMultipleSequenceAlignmentFile(io.BytesIO(msa_file.read()))
        msa_file.seek(0)
        _msa_plot_sequences = [seq for _, seq in parsed]
    except Exception:
        _msa_plot_sequences = None

if _msa_plot_sequences:
    _nseq = len(_msa_plot_sequences)
    _max_div_slider = max(2, min(25, _nseq))
    _col_thr1, _col_thr2 = st.columns(2)
    with _col_thr1:
        _thr_cons = st.slider(
            "Segment markers: min conservation % (both consecutive columns)",
            min_value=0.0,
            max_value=100.0,
            value=80.0,
            step=0.5,
            help=(
                "A consecutive pair (x1, x2) gets a dot only if both columns have conservation "
                "strictly greater than this value."
            ),
            key="msa_seg_conservation_min",
        )
    with _col_thr2:
        _thr_div = st.slider(
            "Segment markers: max diversity (both consecutive columns)",
            min_value=2,
            max_value=max(2, _max_div_slider),
            value=min(3, max(2, _max_div_slider)),
            step=1,
            help=(
                "Both columns must have diversity (distinct residues) strictly less than this "
                "(e.g. 2 allows only columns with diversity 1)."
            ),
            key="msa_seg_diversity_max",
        )

    _pos_tmp, _cons_tmp, _div_tmp = compute_msa_column_conservation_diversity(_msa_plot_sequences)
    _mid_tmp, _ = compute_consensus_segment_midpoints(
        _pos_tmp, _cons_tmp, _div_tmp, _thr_cons, _thr_div
    )
    st.caption(
        f"Segment markers: **{len(_mid_tmp)}** midpoint(s) at y = 100% conservation "
        f"(pairs of consecutive positions both with conservation > {_thr_cons:.1f}% "
        f"and diversity < {_thr_div})."
    )

    _fig_msa = plot_msa_conservation_diversity(
        _msa_plot_sequences,
        segment_conservation_min_pct=_thr_cons,
        segment_diversity_max=_thr_div,
    )
    if _fig_msa:
        st.plotly_chart(_fig_msa, use_container_width=True, key="msa_conservation_diversity_chart")
else:
    if has_schema_contacts:
        st.caption(
            "Parent sequences were not found in session; run SCHEMA contacts on Page 1 with valid parents."
        )
    else:
        st.caption(
            "Upload an MSA in the sidebar or complete Page 1 with parent sequences to plot conservation and diversity."
        )

restrict_graph_segment_marks = st.checkbox(
    "Crossover only at graph segment marks",
    value=DEFAULTS.get("restrict_crossovers_at_graph_marks", False),
    help=(
        "When enabled, RASPP may place fragment boundaries only at SCHEMA crossovers that "
        "match the green dots on the plot above: use the same conservation and diversity "
        "sliders to define which consecutive column pairs get a mark."
    ),
    key="restrict_crossovers_at_graph_marks",
)

# File uploads (if not using previous results)
if not has_schema_contacts:
    col1, col2 = st.columns(2)
    
    with col1:
        st.subheader("Multiple Sequence Alignment")
        if msa_file:
            is_valid, error_msg = validate_msa_file(msa_file)
            if is_valid:
                st.success(f"✓ Valid MSA file: {msa_file.name}")
            else:
                st.error(f"✗ {error_msg}")
    
    with col2:
        st.subheader("SCHEMA Contacts File")
        if contact_file:
            st.info(f"Uploaded: {contact_file.name}")

# Multi-fragment mode only
if not has_schema_contacts:
    st.warning("Multi-fragment analysis requires SCHEMA contacts from Page 1.")

run_button = st.button(
    "Run Multi-Fragment Analysis",
    type="primary",
    disabled=not has_schema_contacts
)

if run_button:
    # Validate inputs
    if not has_data and (not msa_file or not contact_file):
        st.error("Please upload MSA and contact files, or calculate SCHEMA contacts on Page 1")
    else:
        # Prepare inputs
        if has_data:
            # Use existing contacts
            contacts = st.session_state[SESSION_KEYS['schema_contacts']]['contacts']
            parents_list = st.session_state[SESSION_KEYS['schema_contacts']]['parents']
            
            # Extract parents - handle different formats
            if parents_list is None:
                parents = []
            elif isinstance(parents_list, list) and len(parents_list) > 0:
                # Check if it's a list of tuples (name, sequence) or list of sequences
                if isinstance(parents_list[0], tuple) and len(parents_list[0]) == 2:
                    # List of (name, sequence) tuples
                    parents = [seq for _, seq in parents_list]
                elif isinstance(parents_list[0], str):
                    # Already a list of sequences
                    parents = parents_list
                else:
                    # Try to extract sequences anyway
                    parents = [seq if isinstance(seq, str) else str(seq) for seq in parents_list]
            else:
                parents = []
            
            msa_path = st.session_state[SESSION_KEYS['msa_path']]
        else:
            # Load from uploaded files (temp dir cleaned up when block exits)
            from schema_raspp import schema
            with tempfile.TemporaryDirectory() as temp_dir:
                msa_path = save_uploaded_file(msa_file, temp_dir)
                contact_path = save_uploaded_file(contact_file, temp_dir)
                with open(contact_path, 'r') as f:
                    contacts = schema.readContactFile(f)
                with open(msa_path, 'r') as f:
                    parents_list = schema.readMultipleSequenceAlignmentFile(f)
                parents = [seq for _, seq in parents_list]
        
        # Validate fragment range
        if min_fragments is not None and max_fragments is not None and min_fragments > max_fragments:
            st.error("Minimum fragments must be ≤ maximum fragments. Please adjust the sliders.")
            st.stop()

        run_input_summary = {
            "Mode": "Multi-fragment analysis",
            "Contact source": "Page 1 SCHEMA contacts" if has_schema_contacts else "Uploaded files",
            "MSA path": str(msa_path) if msa_path else "Not available",
            "Parent count": len(parents) if 'parents' in locals() and parents else len(parents_for_fallback),
            "Contact count": len(contacts) if contacts else 0,
            "Minimum fragment length": min_fragment_diversity,
            "Minimum fragments": min_fragments,
            "Maximum fragments": max_fragments,
            "Crossover only at graph segment marks": restrict_graph_segment_marks,
            "Segment conservation min %": st.session_state.get(
                "msa_seg_conservation_min",
                DEFAULTS.get("msa_segment_conservation_min", 80.0),
            ),
            "Segment diversity max": st.session_state.get(
                "msa_seg_diversity_max",
                DEFAULTS.get("msa_segment_diversity_max", 3),
            ),
        }
        st.markdown("### Multi-Fragment Run Inputs")
        st.info("Confirming the exact inputs/constraints being used for this run.")
        st.json(run_input_summary)

        # Initialize session state with empty results for incremental updates
        # Always clear stale results before a new multi-fragment run
        st.session_state[SESSION_KEYS['multi_fragment_results']] = {}

        # Create empty containers for live results display
        st.markdown("---")
        st.subheader("📊 Live Results (Updating as tests complete)")
        results_summary_container = st.empty()
        results_table_container = st.empty()

        # Helper function to render current results
        def render_current_results(all_results):
            """Render current results into the containers"""
            import pandas as pd

            if not all_results:
                results_summary_container.info("⏳ Waiting for first results...")
                results_table_container.empty()
                return

            # Summary
            successful = sum(1 for r in all_results.values() if 'error' not in r)
            failed = len(all_results) - successful
            tested = sorted(all_results.keys())

            if failed > 0:
                results_summary_container.warning(
                    f"✓ Tested {len(all_results)} fragment configurations "
                    f"({successful} successful, {failed} failed) | "
                    f"Completed: {min(tested)}-{max(tested)} fragments"
                )
            else:
                results_summary_container.success(
                    f"✓ Tested {len(all_results)} fragment configurations (all successful!) | "
                    f"Completed: {min(tested)}-{max(tested)} fragments"
                )

            # Table
            summary_data = []
            for num_frag, result in sorted(all_results.items()):
                if 'error' not in result:
                    # Check if best_energy exists (not None), 0 is a valid energy value
                    best_energy = result.get('best_energy')
                    if best_energy is not None:
                        energy_str = f"{best_energy:.2f}"
                    else:
                        energy_str = "N/A"
                    summary_data.append({
                        'Fragments': num_frag,
                        'Crossovers': result['num_crossovers'],
                        'Best Energy': energy_str,
                        'Num Designs': result['num_designs']
                    })

            if summary_data:
                df = pd.DataFrame(summary_data)
                results_table_container.dataframe(df, use_container_width=True)
            else:
                results_table_container.info("No successful results yet...")

        # Initial render
        render_current_results(st.session_state.get(SESSION_KEYS['multi_fragment_results'], {}))

        # Create callback to update results incrementally
        def update_results_incrementally(all_results, num_fragments):
            """Update session state and render results without rerunning"""
            st.session_state[SESSION_KEYS['multi_fragment_results']] = all_results
            render_current_results(all_results)

        with st.spinner("Testing multiple fragment configurations..."):
            try:
                from utils.multi_fragment_testing import test_multiple_fragments
                from utils.raspp_wrapper import streamlit_progress_callback

                # Progress tracking
                progress_bar = st.progress(0)
                status_text = st.empty()
                callback = streamlit_progress_callback(progress_bar, status_text)

                # Test multiple fragments with schemarecomb
                # Use parents_for_fallback if available, otherwise use parents
                parents_to_use = parents if 'parents' in locals() and parents else (parents_for_fallback if 'parents_for_fallback' in locals() else [])

                all_results = test_multiple_fragments(
                    contacts=contacts,  # Needed for SCHEMA-RASPP
                    parents=parents_to_use,  # Needed for SCHEMA-RASPP
                    min_fragments=min_fragments,
                    max_fragments=max_fragments,
                    min_fragment_diversity=min_fragment_diversity,
                    progress_callback=callback,
                    parents_object=None,  # Not using schemarecomb
                    use_schemarecomb=False,  # Use SCHEMA-RASPP instead
                    results_update_callback=update_results_incrementally,
                    restrict_crossovers_at_graph_marks=restrict_graph_segment_marks,
                    segment_conservation_min_pct=st.session_state.get(
                        "msa_seg_conservation_min",
                        DEFAULTS.get("msa_segment_conservation_min", 80.0),
                    ),
                    segment_diversity_max=st.session_state.get(
                        "msa_seg_diversity_max",
                        DEFAULTS.get("msa_segment_diversity_max", 3),
                    ),
                )

                # Final update (in case callback wasn't called for last fragment)
                st.session_state[SESSION_KEYS['multi_fragment_results']] = all_results
                # Store parents for later use on Crossover Analysis page
                if parents_to_use:
                    st.session_state[SESSION_KEYS['raspp_parents']] = parents_to_use
                render_current_results(all_results)

                progress_bar.empty()
                status_text.empty()

                # Show summary of what was tested
                successful = sum(1 for r in all_results.values() if 'error' not in r)
                failed = len(all_results) - successful
                if failed > 0:
                    st.warning(f"✓ Tested {len(all_results)} fragment configurations ({successful} successful, {failed} failed)")
                else:
                    st.success(f"✓ Tested {len(all_results)} fragment configurations (all successful)!")

                # Autosave after multi-fragment test completes
                from utils.session_manager import auto_save
                auto_save("multi_fragment_complete")

                from utils.session_manager import render_save_project_ui
                render_save_project_ui('save_multifrag_btn')

                # Display summary table
                import pandas as pd
                summary_data = []
                for num_frag, result in all_results.items():
                    if 'error' not in result:
                        summary_data.append({
                            'Fragments': num_frag,
                            'Crossovers': result['num_crossovers'],
                            'Best Energy': f"{result['best_energy']:.2f}" if result['best_energy'] else "N/A",
                            'Num Designs': result['num_designs']
                        })

                if summary_data:
                    df = pd.DataFrame(summary_data)
                    st.dataframe(df, use_container_width=True)

                    # Plot comparison
                    from utils.visualization import plot_crossover_comparison
                    # Create combined results for visualization
                    combined_results = []
                    for result in all_results.values():
                        if 'results' in result:
                            combined_results.extend(result['results'])

                    if combined_results:
                        fig = plot_crossover_comparison(combined_results)
                        if fig:
                            st.plotly_chart(fig, use_container_width=True)

            except Exception as e:
                st.error(f"Error in multi-fragment testing: {str(e)}")
                import traceback
                st.code(traceback.format_exc())

# Display multi-fragment results
if (SESSION_KEYS['multi_fragment_results'] in st.session_state and
    isinstance(st.session_state[SESSION_KEYS['multi_fragment_results']], dict) and
    len(st.session_state[SESSION_KEYS['multi_fragment_results']]) > 0):
    st.header("Multi-Fragment Analysis Results")
    
    all_results = st.session_state[SESSION_KEYS['multi_fragment_results']]
    
    # Check if all_results is valid
    if all_results is None:
        st.warning("Multi-fragment results are None. Please run the multi-fragment analysis again.")
    elif not isinstance(all_results, dict):
        st.error(f"Multi-fragment results are in unexpected format: {type(all_results)}. Expected dict.")
    elif len(all_results) == 0:
        st.info("No multi-fragment results available. Please run the multi-fragment analysis.")
    else:
        # Summary statistics
        st.subheader("Summary by Fragment Count")
        import pandas as pd

        summary_data = []
        error_data = []
        for num_frag, result in sorted(all_results.items()):
            if 'error' in result:
                error_data.append({
                    'Fragments': num_frag,
                    'Blocks': result.get('num_blocks', num_frag),
                    'Crossovers': result.get('num_crossovers', 'N/A'),
                    'Status': 'Error',
                    'Error': result['error'][:100] + '...' if len(result['error']) > 100 else result['error']
                })
            elif result.get('best_energy') is not None:
                best_energy = result['best_energy']
                summary_data.append({
                    'Fragments': num_frag,
                    'Blocks': result.get('num_blocks', num_frag),
                    'Crossovers': result['num_crossovers'],
                    'Best Energy': f"{best_energy:.2f}" if isinstance(best_energy, (int, float)) else str(best_energy),
                    'Num Designs': result['num_designs']
                })
            else:
                # No error but no energy either - might be empty results
                error_data.append({
                    'Fragments': num_frag,
                    'Blocks': result.get('num_blocks', num_frag),
                    'Crossovers': result.get('num_crossovers', 'N/A'),
                    'Status': 'No Results',
                    'Error': 'No designs found (empty results)'
                })
        
        # Show errors if any
        if error_data:
            st.warning(f"⚠️ {len(error_data)} fragment configurations failed to generate libraries:")
            error_df = pd.DataFrame(error_data)
            with st.expander("View Failed Configurations", expanded=True):
                st.dataframe(error_df, use_container_width=True)
                
                # Show diagnostic info
                if 'schema_contacts' in st.session_state and st.session_state[SESSION_KEYS['schema_contacts']]:
                    contacts_data = st.session_state[SESSION_KEYS['schema_contacts']]
                    if 'parents_object' in contacts_data and contacts_data['parents_object']:
                        parents_obj = contacts_data['parents_object']
                        if hasattr(parents_obj, 'alignment'):
                            alignment_len = len(parents_obj.alignment)
                            st.info(f"""
                            **Diagnostic Information:**
                            - Alignment length: {alignment_len} residues
                            - Min fragment diversity: {min_fragment_diversity}
                            - Minimum total length needed for N fragments: N × {min_fragment_diversity}
                            
                            **Failed fragment counts require:**
                            {chr(10).join([f'  - {row["Fragments"]} fragments: {row["Fragments"]} × {min_fragment_diversity} = {row["Fragments"] * min_fragment_diversity} residues needed' for row in error_data[:5]])}
                            
                            💡 **Tips:**
                            - Try reducing min_fragment_diversity if alignment length is sufficient
                            - Check if alignment length ({alignment_len}) is long enough for {max([int(r['Fragments']) for r in error_data])} fragments
                            - Some fragment counts may be geometrically impossible with current constraints
                            """)
                
                st.info("💡 Tip: Try reducing min_fragment_diversity or check if alignment length is sufficient for the requested number of fragments.")
        
        if summary_data:
            df = pd.DataFrame(summary_data)
            st.dataframe(df, use_container_width=True)
            
            # Find best overall (skip entries with N/A energy)
            valid_results = [r for r in summary_data if r['Best Energy'] != 'N/A']
            if valid_results:
                best_frag = min(valid_results, key=lambda x: float(x['Best Energy']) if isinstance(x['Best Energy'], str) and x['Best Energy'] != 'N/A' else x['Best Energy'])
                energy_display = best_frag['Best Energy'] if isinstance(best_frag['Best Energy'], str) else f"{best_frag['Best Energy']:.2f}"
                st.success(f"Best configuration: {best_frag['Fragments']} fragments with energy {energy_display}")
            else:
                st.warning("No valid energy values found in results.")
            
            # Allow user to select a fragment count to view details
            selected_frag = st.selectbox(
                "Select fragment count to view detailed results",
                [r['Fragments'] for r in summary_data],
                format_func=lambda x: f"{x} fragments - Energy: {all_results[x]['best_energy']:.2f}"
            )
            
            if selected_frag in all_results and 'results' in all_results[selected_frag]:
                # Use selected fragment's results for inline display only; by default do not
                # overwrite the global RASPP results in session state.
                selected_results = all_results[selected_frag]['results']
                selected_raw_results = all_results[selected_frag]['raw_results']

                # Show a brief summary table for the selected fragment count
                if selected_results:
                    import pandas as pd
                    display_rows = []
                    for r in selected_results:
                        display_rows.append({
                            'Average Energy': r.get('average_energy'),
                            'Crossovers': ', '.join(map(str, r.get('crossovers', []))),
                            'Num Fragments': r.get('num_fragments'),
                            'Min Length': r.get('min_fragment_length'),
                            'Max Length': r.get('max_fragment_length'),
                        })
                    detail_df = pd.DataFrame(display_rows)
                    st.subheader(f"Detailed Designs for {selected_frag} Fragments")
                    st.dataframe(detail_df, use_container_width=True)

                    # Optional: allow promoting this fragment count to be the main RASPP result
                    if st.button("Use these designs as main RASPP results", key=f"use_raspp_results_{selected_frag}"):
                        st.session_state[SESSION_KEYS['raspp_results']] = selected_results
                        st.session_state[SESSION_KEYS['raspp_raw_results']] = selected_raw_results

# Display single RASPP results
if SESSION_KEYS['raspp_results'] in st.session_state and \
   st.session_state[SESSION_KEYS['raspp_results']] is not None:
    st.header("RASPP Results")
    
    results = st.session_state[SESSION_KEYS['raspp_results']]
    raw_results = st.session_state[SESSION_KEYS['raspp_raw_results']]
    
    # Summary statistics
    col1, col2, col3, col4 = st.columns(4)
    
    with col1:
        st.metric("Number of Designs", len(results))
    
    with col2:
        avg_energies = [r['average_energy'] for r in results]
        st.metric("Best Energy", f"{min(avg_energies):.2f}")
    
    with col3:
        st.metric("Average Energy", f"{sum(avg_energies)/len(avg_energies):.2f}")
    
    with col4:
        st.metric("Worst Energy", f"{max(avg_energies):.2f}")
    
    # Sort by energy
    sorted_results = sorted(results, key=lambda x: x['average_energy'])
    
    # Display top designs
    st.subheader("Top Designs")
    
    # Visualization
    fig = plot_crossover_comparison(results)
    if fig:
        st.plotly_chart(fig, use_container_width=True)
    
    # Table of results
    import pandas as pd
    display_data = []
    for i, r in enumerate(sorted_results[:20]):  # Top 20
        display_data.append({
            'Rank': i + 1,
            'Average Energy': f"{r['average_energy']:.2f}",
            'Crossovers': ', '.join(map(str, r['crossovers'])),
            'Num Fragments': r['num_fragments'],
            'Min Length': r['min_fragment_length'],
            'Max Length': r['max_fragment_length']
        })
    
    df = pd.DataFrame(display_data)
    st.dataframe(df, use_container_width=True)
    
    # Detailed view of selected design
    st.subheader("Design Details")
    
    selected_idx = st.selectbox(
        "Select design to view details",
        range(min(10, len(sorted_results))),
        format_func=lambda x: f"Design {x+1} - Energy: {sorted_results[x]['average_energy']:.2f}"
    )
    
    selected_design = sorted_results[selected_idx]
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.write("**Crossover Positions:**")
        st.write(selected_design['crossovers'])
        
        st.write("**Fragment Lengths:**")
        parents = st.session_state.get(SESSION_KEYS['raspp_parents'])
        if parents:
            fragments = schema.getFragments(selected_design['crossovers'], parents[0])
            fragment_lengths = [end - start for start, end in fragments]
            st.write(fragment_lengths)
            
            # Plot fragment lengths
            fig_frag = plot_fragment_lengths(
                selected_design['crossovers'],
                len(parents[0])
            )
            if fig_frag:
                st.plotly_chart(fig_frag, use_container_width=True)
    
    with col2:
        st.write("**Design Parameters:**")
        st.json({
            'Average Energy': selected_design['average_energy'],
            'Number of Crossovers': len(selected_design['crossovers']),
            'Number of Fragments': selected_design['num_fragments'],
            'Min Fragment Length': selected_design['min_fragment_length'],
            'Max Fragment Length': selected_design['max_fragment_length']
        })
    
    # Export results
    st.header("Export Results")
    
    # Create CSV
    export_df = pd.DataFrame(sorted_results)
    csv = export_df.to_csv(index=False)
    
    st.download_button(
        label="Download Results as CSV",
        data=csv,
        file_name="raspp_results.csv",
        mime="text/csv"
    )
    
    # Create text file with crossover positions
    crossover_text = ""
    for i, r in enumerate(sorted_results):
        crossover_text += f"Design {i+1}: Energy={r['average_energy']:.2f}, Crossovers={r['crossovers']}\n"
    
    st.download_button(
        label="Download Crossovers as Text",
        data=crossover_text,
        file_name="raspp_crossovers.txt",
        mime="text/plain"
    )
    
    # Save checkpoint
    st.markdown("---")
    from utils.session_manager import render_save_project_ui
    render_save_project_ui('save_raspp_results_btn')