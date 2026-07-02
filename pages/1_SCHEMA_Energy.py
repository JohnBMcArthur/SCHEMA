"""
Page 1: SCHEMA Energy Calculation
Calculate SCHEMA energies from protein structures and sequence alignments.
"""

import streamlit as st
import tempfile
import os
from pathlib import Path

# Import utilities
import sys
sys.path.insert(0, str(Path(__file__).parent.parent))

from utils.file_handlers import (
    validate_pdb_file, validate_msa_file, validate_crossover_file,
    save_uploaded_file
)
from utils.schema_wrapper import (
    calculate_contacts, calculate_energies, save_contacts
)
from utils.visualization import plot_contact_map, plot_energy_distribution
from utils.session_manager import init_session_state, has_required_data
from utils.temp_file_manager import temp_file_manager
from utils.config import DEFAULTS, SESSION_KEYS

st.set_page_config(
    page_title="SCHEMA Energy",
    page_icon="⚡",
    layout="wide"
)

# Initialize session state
init_session_state()

st.title("⚡ SCHEMA Energy Calculation")

st.markdown("""
Calculate SCHEMA disruption energies for chimeric proteins. Upload your protein 
structure (PDB file) and multiple sequence alignment (MSA) to begin.

**Workflow:**
1. Upload PDB structure file and MSA file
2. Optionally upload PDB-parent alignment file
3. Set contact distance threshold
4. Calculate SCHEMA contacts
5. Upload crossover points file to calculate energies
""")

# Default parameters for contacts (prefer existing session values when available)
if 'contact_distance' in st.session_state and st.session_state['contact_distance'] is not None:
    contact_distance = st.session_state['contact_distance']
else:
    contact_distance = DEFAULTS.get('contact_distance', 5.0)
    st.session_state['contact_distance'] = contact_distance

if 'chains' in st.session_state and st.session_state['chains']:
    chains = st.session_state['chains']
else:
    chains = (DEFAULTS.get('chains') or ['A'])[0] if DEFAULTS.get('chains') else "A"
    st.session_state['chains'] = chains

# Sidebar for project management
with st.sidebar:
    st.markdown("### Project Management")
    
    # Show current project if loaded
    if 'current_project' in st.session_state and st.session_state.get('current_project'):
        st.info(f"**📁 {st.session_state['current_project']}**")
    
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
            key='page1_load_project'
        )
        
        if st.button("📂 Load", use_container_width=True, key='page1_load_btn'):
            try:
                metadata = load_checkpoint(checkpoint_paths[selected_idx])
                st.success(f"✓ Loaded: {metadata.get('project_name', 'Unknown')}")
                st.rerun()
            except Exception as e:
                st.error(f"Error: {str(e)}")

# Workflow mode selection
st.header("Input Method")

workflow_mode = st.radio(
    "Choose workflow mode",
    ["Automated (Enter Sequence)", "Manual (Upload Files)"],
    help="Automated mode: Enter a sequence and we'll find similar sequences, align them, and find PDB structure. Manual mode: Upload your own files."
)

if workflow_mode == "Automated (Enter Sequence)":
    # Project name input
    st.subheader("Name Project")
    project_name_input = st.text_input(
        "Project Name",
        value=st.session_state.get('current_project', ''),
        help="Give your project a name. This will be used for autosaves.",
        key="project_name_input"
    )
    # Store in session state for autosave to use
    if project_name_input:
        st.session_state['project_name_for_autosave'] = project_name_input
    
    st.subheader("Enter Protein Sequence")
    
    sequence_input = st.text_area(
        "Protein Sequence",
        value=st.session_state.get('query_sequence', ''),
        height=150,
        placeholder="Enter a single protein sequence in one-letter code (e.g., MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFGLAPFLPDQIHFVHSQELLSRYPDLDAKGRERAIAKDLGAVFLVGIGGKLSDGHRHDVRAPDYDDWSTPSELGHAGLNGDILVWNPVLEDAFELSSMGIRVDADTLKHQLALTGDEDRLELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQAGVWPAAVRESVPSLL",
        help="Enter a single protein sequence in one-letter code (at least 20 amino acids)"
    )
    
    # Parameters for automated workflow
    with st.expander("Automated Workflow Parameters"):
        num_sequences = st.slider(
            "Number of sequences to retrieve",
            min_value=5,
            max_value=20,
            value=10,
            help="Number of similar sequences to find via BLAST"
        )
        
        min_identity = st.slider(
            "Minimum identity for BLAST hits",
            min_value=0.3,
            max_value=0.8,
            value=0.5,
            step=0.05,
            help="Minimum sequence identity for BLAST hits (sequences below this will be excluded)"
        )
        
        max_identity = st.slider(
            "Maximum identity for selection",
            min_value=0.7,
            max_value=0.99,
            value=0.9,
            step=0.01,
            help="Maximum sequence identity for selection (ensures diversity - sequences above this will be excluded)"
        )
    
    if st.button("Run Automated Workflow", type="primary"):
        if not sequence_input or len(sequence_input.strip()) < 20:
            st.error("Please enter a valid protein sequence (at least 20 amino acids)")
        else:
            # Clean sequence (remove whitespace, newlines)
            clean_sequence = ''.join(sequence_input.split()).upper()
            
            # Validate sequence
            valid_aa = set('ACDEFGHIKLMNPQRSTVWY')
            if not all(c in valid_aa for c in clean_sequence):
                st.error("Sequence contains invalid amino acid characters. Please use standard one-letter codes.")
            else:
                with st.spinner("Running automated workflow..."):
                    try:
                        from utils.sequence_automation import (
                            blast_search_sequences, align_sequences_ebi_muscle,
                            find_pdb_structure
                        )
                        from utils.raspp_wrapper import streamlit_progress_callback
                        
                        # Persist query sequence in session so it can be restored with a project
                        st.session_state['query_sequence'] = clean_sequence

                        # Progress tracking
                        progress_bar = st.progress(0)
                        status_text = st.empty()
                        callback = streamlit_progress_callback(progress_bar, status_text)
                        
                        # Step 1: BLAST search
                        callback(0.05, "Step 1/4: Searching for similar sequences...")
                        sequences, tsv_data, tsv_error, best_hit = blast_search_sequences(
                            clean_sequence,
                            num_sequences=num_sequences,
                            min_identity=min_identity,
                            max_identity=max_identity,
                            progress_callback=lambda p, m: callback(0.05 + p * 0.2, f"BLAST: {m}")
                        )
                        
                        # Display BLAST TSV results for debugging
                        st.markdown("---")
                        st.subheader("🔍 BLAST TSV Results (Debug)")
                        if tsv_data:
                            st.info(f"TSV data length: {len(tsv_data)} characters, {len(tsv_data.split(chr(10)))} lines")
                            # Parse and analyze TSV
                            try:
                                import pandas as pd
                                lines = tsv_data.strip().split('\n')
                                
                                # Detect format by checking header
                                header_line = None
                                for line in lines:
                                    if line.strip() and not line.startswith('#'):
                                        header_line = line
                                        break
                                
                                use_alt_format = False
                                if header_line and 'Hit' in header_line and 'Accession' in header_line:
                                    use_alt_format = True
                                
                                # Skip comment lines and header
                                data_lines = [line for line in lines if line and not line.startswith('#') and line != header_line]
                                
                                if data_lines:
                                    # Parse TSV based on detected format
                                    rows = []
                                    filtered_by_identity = []
                                    filtered_by_evalue = []
                                    passed_initial = []
                                    
                                    for line in data_lines:
                                        fields = line.split('\t')
                                        
                                        if use_alt_format:
                                            # Alternative format: Hit, DB, Accession, Description, Organism, Length, Score(Bits), Identities(%), Positives(%), E()
                                            # Note: Description may contain tabs, so parse from the right
                                            if len(fields) < 10:
                                                continue
                                            
                                            try:
                                                # Parse from the right since last columns are fixed
                                                evalue_str = fields[-1].strip()  # Last column: E()
                                                identity_str = fields[-3].rstrip('%').strip()  # Third to last: Identities(%)
                                                score_str = fields[-4].strip()  # Score(Bits)
                                                length_str = fields[-5].strip()  # Length
                                                
                                                # Accession is always column 2
                                                subject_acc = fields[2]  # Accession
                                                
                                                identity_pct = float(identity_str)
                                                evalue = float(evalue_str) if evalue_str else 1e5
                                                
                                                row_data = {
                                                    'Hit': fields[0],
                                                    'Accession': subject_acc,
                                                    'Identity %': f"{identity_pct:.2f}",
                                                    'E-value': f"{evalue:.2e}",
                                                    'Score': score_str,
                                                    'Length': length_str,
                                                    'Description': ' '.join(fields[3:-5])[:50] if len(fields) > 8 else ''  # Description is everything between Accession and Organism
                                                }
                                                rows.append(row_data)
                                                
                                                # Check filters (outside try/except so it runs for successfully parsed rows)
                                                min_identity_pct = min_identity * 100
                                                max_identity_pct = max_identity * 100
                                                
                                                if evalue > 1e-5:
                                                    filtered_by_evalue.append(row_data)
                                                elif identity_pct < min_identity_pct or identity_pct >= max_identity_pct:
                                                    filtered_by_identity.append(row_data)
                                                else:
                                                    passed_initial.append(row_data)
                                            except (ValueError, IndexError):
                                                continue  # Skip if can't parse
                                        else:
                                            # Standard format: query acc, subject acc, identity, alignment length, 
                                            # mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
                                            if len(fields) >= 12:
                                                identity_pct = float(fields[2])
                                                evalue = float(fields[10])
                                                
                                                row_data = {
                                                    'Query': fields[0],
                                                    'Subject': fields[1],
                                                    'Identity %': f"{identity_pct:.2f}",
                                                    'E-value': f"{evalue:.2e}",
                                                    'Bit Score': fields[11],
                                                    'Alignment Length': fields[3]
                                                }
                                                rows.append(row_data)
                                                
                                                # Check filters
                                                min_identity_pct = min_identity * 100
                                                max_identity_pct = max_identity * 100
                                                
                                                if evalue > 1e-5:
                                                    filtered_by_evalue.append(row_data)
                                                elif identity_pct < min_identity_pct or identity_pct >= max_identity_pct:
                                                    filtered_by_identity.append(row_data)
                                                else:
                                                    passed_initial.append(row_data)
                                    
                                    # Show summary statistics
                                    col1, col2, col3, col4 = st.columns(4)
                                    with col1:
                                        st.metric("Total Hits", len(rows))
                                    with col2:
                                        st.metric("Passed E-value", len(rows) - len(filtered_by_evalue))
                                    with col3:
                                        st.metric("Passed Identity Filter", len(passed_initial))
                                    with col4:
                                        st.metric("Selected Sequences", len(sequences) - 1)  # -1 for query
                                    
                                    # Show filtered out by identity
                                    if filtered_by_identity:
                                        st.warning(f"⚠️ {len(filtered_by_identity)} hits filtered by identity range ({min_identity:.0%}-{max_identity:.0%})")
                                        with st.expander(f"View {len(filtered_by_identity)} hits filtered by identity", expanded=False):
                                            df_filtered = pd.DataFrame(filtered_by_identity[:20])  # Show first 20
                                            st.dataframe(df_filtered, use_container_width=True)
                                    
                                    # Show filtered out by e-value
                                    if filtered_by_evalue:
                                        st.info(f"ℹ️ {len(filtered_by_evalue)} hits filtered by e-value (>{1e-5})")
                                    
                                    # Show hits that passed initial filters
                                    if passed_initial:
                                        st.success(f"✓ {len(passed_initial)} hits passed initial filters (identity: {min_identity:.0%}-{max_identity:.0%}, e-value: ≤1e-5)")
                                        with st.expander(f"View {len(passed_initial)} hits that passed initial filters", expanded=True):
                                            df_passed = pd.DataFrame(passed_initial[:50])  # Show first 50
                                            st.dataframe(df_passed, use_container_width=True)
                                    
                                    # Show all hits table
                                    with st.expander("View All BLAST Hits", expanded=False):
                                        df = pd.DataFrame(rows)
                                        st.dataframe(df, use_container_width=True)
                                    
                                    # Show full TSV as code
                                    with st.expander("View Full TSV Data (Raw)", expanded=False):
                                        st.code(tsv_data, language='text')
                            except Exception as e:
                                st.warning(f"Could not parse TSV as table: {str(e)}")
                                # Fallback: show raw TSV
                                with st.expander("View Full TSV Data", expanded=True):
                                    st.code(tsv_data, language='text')
                        else:
                            if tsv_error:
                                st.error(f"❌ TSV retrieval failed: {tsv_error}")
                                st.info("ℹ️ Using XML fallback for BLAST results parsing")
                                
                                # Show what we got from TSV attempt (tsv_response may not exist when TSV failed)
                                with st.expander("TSV Response Details", expanded=True):
                                    st.code(f"TSV request failed. Error: {tsv_error}\n\nUsing XML fallback for BLAST results.", language='text')
                            else:
                                st.warning("No TSV data available (may have used XML fallback)")
                            
                            # Show debug info about sequences found
                            st.info(f"**Debug Info:**")
                            st.write(f"- Sequences found: {len(sequences)}")
                            st.write(f"- Candidate hits processed: Check BLAST job status")
                            st.write(f"- Identity filters: {min_identity:.0%} ≤ identity < {max_identity:.0%}")
                            
                            # Show what sequences were found
                            if len(sequences) > 1:
                                st.success(f"✓ Found {len(sequences)-1} additional sequences beyond query")
                                with st.expander("View Selected Sequences"):
                                    for seq_id, seq in sequences:
                                        st.text(f"{seq_id}: {seq[:100]}...")
                            else:
                                st.error("⚠️ Only query sequence found - no additional sequences passed filters")
                        
                        # Validate that we have enough sequences for alignment
                        if not sequences or len(sequences) < 2:
                            st.error(
                                f"❌ Not enough sequences found for alignment. "
                                f"Found {len(sequences) if sequences else 0} sequence(s), but at least 2 are required. "
                                f"\n\n**Suggestions:**\n"
                                f"- Try adjusting the identity thresholds (min: {min_identity:.0%}, max: {max_identity:.0%})\n"
                                f"- Try a different query sequence\n"
                                f"- Increase the number of sequences to retrieve"
                            )
                            raise ValueError("Not enough sequences found for alignment")
                        
                        st.success(f"✓ Found {len(sequences)} sequence(s) (including query)")
                        
                        # Display found sequences
                        with st.expander("View Found Sequences"):
                            for seq_id, seq in sequences[:10]:  # Show first 10
                                st.text(f"{seq_id}: {seq[:50]}...")
                        
                        # Step 2: Alignment
                        callback(0.3, f"Step 2/4: Aligning {len(sequences)} sequences using EBI MUSCLE API...")
                        
                        # Debug: Show FASTA format before alignment
                        with st.expander("🔍 Debug: FASTA Sequences (for troubleshooting)", expanded=False):
                            st.markdown("**Sequences in FASTA format that will be sent to EBI MUSCLE API:**")
                            fasta_preview = ""
                            for seq_id, seq in sequences:
                                clean_id = seq_id.replace(' ', '_').replace('|', '_')
                                fasta_preview += f">{clean_id}\n{seq}\n"
                            st.code(fasta_preview, language='fasta')
                            st.info(f"Total sequences: {len(sequences)} | Total characters: {len(fasta_preview)}")
                        
                        msa_file = align_sequences_ebi_muscle(
                            sequences,
                            progress_callback=lambda p, m: callback(0.3 + p * 0.2, f"Alignment: {m}")
                        )
                        
                        st.success(f"✓ {len(sequences)} sequences aligned successfully")
                        
                        # Step 3: Download AlphaFold structure (find best match for query)
                        callback(0.6, "Step 3/4: Finding AlphaFold structure...")
                        from utils.sequence_automation import find_best_structure_for_query, download_alphafold_structure
                        
                        # Use best hit if available (highest identity, regardless of filtering)
                        if best_hit and best_hit[0]:
                            uniprot_id, identity_pct = best_hit
                            callback(0.65, f"Using best BLAST hit: {uniprot_id} ({identity_pct:.1f}% identity)")
                            pdb_id, pdb_file = download_alphafold_structure(
                                uniprot_id,
                                progress_callback=lambda p, m: callback(0.65 + p * 0.15, f"AlphaFold: {m}")
                            )
                        else:
                            # Fall back to finding best structure from sequences
                            pdb_id, pdb_file = find_best_structure_for_query(
                                clean_sequence,
                                sequences,
                                progress_callback=lambda p, m: callback(0.6 + p * 0.2, f"AlphaFold: {m}")
                            )
                        
                        if pdb_file:
                            st.success(f"✓ Found AlphaFold structure: {pdb_id}")
                        else:
                            st.warning("⚠ No AlphaFold structure found. You can upload one manually below.")
                            pdb_id = None
                        
                        # Step 4: Calculate contacts (if PDB found)
                        if pdb_file:
                            callback(0.9, "Step 4/4: Calculating SCHEMA contacts...")
                            from utils.schema_wrapper import calculate_contacts
                            from utils.schemarecomb_bridge import parent_sequences_to_json
                            
                            # Calculate contacts using schemarecomb (automatic PDB alignment)
                            contacts_result = calculate_contacts(
                                pdb_file,
                                msa_file,
                                chains=['A', ' '],
                                contact_distance=contact_distance,
                                progress_callback=lambda p, m: callback(0.9 + p * 0.1, f"Contacts: {m}")
                            )
                            
                            # Store in session state
                            st.session_state.schema_contacts = contacts_result
                            # Use trimmed MSA file if available (domain boundaries applied)
                            st.session_state.msa_path = contacts_result.get('trimmed_msa_file', msa_file)
                            st.session_state.pdb_path = pdb_file
                            st.session_state.pdb_id = pdb_id
                            st.session_state.sequences = sequences
                            
                            # Store ParentSequences object as JSON for later use
                            if 'parents_object' in contacts_result:
                                st.session_state.parents_object_json = parent_sequences_to_json(contacts_result['parents_object'])
                            
                            callback(1.0, "✓ Workflow complete!")
                            st.success("✓ SCHEMA contacts calculated successfully!")
                            
                            # Display summary
                            num_contacts = len(contacts_result['contacts']) if contacts_result.get('contacts') else 0
                            num_residues = len(contacts_result.get('residues', [])) if contacts_result.get('residues') else "N/A"
                            st.info(f"""
                            **Workflow Summary:**
                            - Sequences found: {len(sequences)}
                            - PDB structure: {pdb_id if pdb_id else 'Not found'}
                            - Contacts calculated: {num_contacts}
                            - Number of residues: {num_residues}
                            """)
                            
                            # Autosave after contacts calculated
                            from utils.session_manager import auto_save
                            auto_save("contacts_calculated")
                        else:
                            # Store alignment for manual PDB upload
                            st.session_state.msa_path = msa_file
                            st.session_state.sequences = sequences
                            st.info("Alignment complete. Please upload a PDB file below to calculate contacts.")
                            
                            # Autosave after MSA alignment (even without PDB)
                            from utils.session_manager import auto_save
                            auto_save("msa_aligned")
                        
                        progress_bar.empty()
                        status_text.empty()
                        
                    except Exception as e:
                        st.error(f"Error in automated workflow: {str(e)}")
                        import traceback
                        with st.expander("Error Details"):
                            st.code(traceback.format_exc())
    
    # Manual PDB upload option (if automated workflow didn't find one)
    if 'msa_path' in st.session_state and ('pdb_path' not in st.session_state or st.session_state.get('pdb_path') is None):
        st.markdown("---")
        st.subheader("Upload PDB Structure (if not found automatically)")
        manual_pdb_file = st.file_uploader(
            "Upload PDB file",
            type=['pdb'],
            help="If no PDB was found automatically, upload one here to calculate contacts",
            key="manual_pdb_upload"
        )
        
        if manual_pdb_file and st.button("Calculate Contacts with Uploaded PDB", key="calc_manual_pdb"):
            from utils.schema_wrapper import calculate_contacts
            import tempfile
            
            temp_dir = tempfile.mkdtemp()
            pdb_path = save_uploaded_file(manual_pdb_file, temp_dir)
            
            try:
                contacts_result = calculate_contacts(
                    pdb_path,
                    st.session_state[SESSION_KEYS['msa_path']],
                    None,
                    ['A', ' '],
                    contact_distance
                )
                
                st.session_state[SESSION_KEYS['schema_contacts']] = contacts_result
                st.session_state[SESSION_KEYS['pdb_path']] = pdb_path
                # Update MSA path to trimmed version if available
                if contacts_result.get('trimmed_msa_file'):
                    st.session_state[SESSION_KEYS['msa_path']] = contacts_result['trimmed_msa_file']
                st.success("✓ SCHEMA contacts calculated!")
                
                # Display contact statistics
                num_contacts = len(contacts_result['contacts'])
                # Handle case where residues might be None
                num_residues = len(contacts_result['residues']) if contacts_result.get('residues') else "N/A"
                num_parents = len(contacts_result['parents'])
                
                st.info(f"""
                **Contact Statistics:**
                - Number of contacts: {num_contacts}
                - Number of residues: {num_residues}
                - Number of parent sequences: {num_parents}
                """)
            except Exception as e:
                st.error(f"Error calculating contacts: {str(e)}")

else:
    # Manual file upload mode
    # Project name input
    st.subheader("Name Project")
    project_name_input_manual = st.text_input(
        "Project Name",
        value=st.session_state.get('current_project', ''),
        help="Give your project a name. This will be used for autosaves.",
        key="project_name_input_manual"
    )
    # Store in session state for autosave to use
    if project_name_input_manual:
        st.session_state['project_name_for_autosave'] = project_name_input_manual
    
    st.header("File Uploads")

    col1, col2 = st.columns(2)

    with col1:
        st.subheader("PDB Structure File")
        pdb_file = st.file_uploader(
            "Upload PDB file",
            type=['pdb'],
            help="Protein structure file in PDB format"
        )
        
        if pdb_file:
            is_valid, error_msg = validate_pdb_file(pdb_file)
            if is_valid:
                st.success(f"✓ Valid PDB file: {pdb_file.name}")
            else:
                st.error(f"✗ {error_msg}")
        
        st.subheader("Multiple Sequence Alignment (MSA)")
        msa_file = st.file_uploader(
            "Upload MSA file",
            type=['txt', 'msa', 'fasta', 'fa', 'aln'],
            help="Multiple sequence alignment file"
        )
        
        if msa_file:
            is_valid, error_msg = validate_msa_file(msa_file)
            if is_valid:
                st.success(f"✓ Valid MSA file: {msa_file.name}")
            else:
                st.error(f"✗ {error_msg}")

    with col2:
        st.subheader("Crossover Points File (Optional)")
        crossover_file = st.file_uploader(
            "Upload crossover points file",
            type=['txt'],
            help="Optional: File containing crossover positions for energy calculation"
        )
        
        if crossover_file:
            is_valid, error_msg = validate_crossover_file(crossover_file)
            if is_valid:
                st.success(f"✓ Valid crossover file: {crossover_file.name}")
            else:
                st.error(f"✗ {error_msg}")
            # Store in session state for use outside this scope
            st.session_state.crossover_file = crossover_file
        else:
            # Clear from session state if not uploaded
            if 'crossover_file' in st.session_state:
                del st.session_state.crossover_file

    # Process chain identifiers
    chain_list = None
    if chains:
        chain_list = [c.strip() for c in chains.split(',')]

    # Calculate contacts (manual mode)
    if st.button("Calculate SCHEMA Contacts", type="primary"):
        if not pdb_file or not msa_file:
            st.error("Please upload both PDB and MSA files")
        else:
            with st.spinner("Calculating SCHEMA contacts..."):
                try:
                    # Save uploaded files temporarily
                    temp_dir = tempfile.mkdtemp()
                    pdb_path = save_uploaded_file(pdb_file, temp_dir)
                    msa_path = save_uploaded_file(msa_file, temp_dir)
                    
                    # Calculate contacts using schemarecomb (automatic PDB alignment)
                    from utils.schemarecomb_bridge import parent_sequences_to_json
                    contacts_result = calculate_contacts(
                        pdb_path,
                        msa_path,
                        chains=chain_list,
                        contact_distance=contact_distance
                    )
                    
                    # Store ParentSequences object as JSON for later use
                    if 'parents_object' in contacts_result:
                        st.session_state.parents_object_json = parent_sequences_to_json(contacts_result['parents_object'])
                    
                    # Store in session state
                    st.session_state.schema_contacts = contacts_result
                    st.session_state.temp_dir = temp_dir
                    st.session_state.pdb_path = pdb_path
                    # Use trimmed MSA file if available (domain boundaries applied)
                    st.session_state.msa_path = contacts_result.get('trimmed_msa_file', msa_path)
                    
                    st.success("✓ SCHEMA contacts calculated successfully!")
                    
                    # Reset downstream RASPP state so new contacts start a fresh design workflow
                    for key in [
                        SESSION_KEYS['raspp_results'],
                        SESSION_KEYS['raspp_raw_results'],
                        SESSION_KEYS['raspp_parents'],
                        SESSION_KEYS['multi_fragment_results'],
                    ]:
                        if key in st.session_state:
                            del st.session_state[key]
                    
                    # Display contact statistics
                    num_contacts = len(contacts_result['contacts'])
                    # Calculate num_residues from residues if available, otherwise from contacts or parents
                    if contacts_result.get('residues') is not None:
                        num_residues = len(contacts_result['residues'])
                    elif contacts_result.get('parents_object') is not None:
                        # Use aligned sequence length from schemarecomb ParentSequences
                        parents_obj = contacts_result['parents_object']
                        if hasattr(parents_obj, 'p0_aligned'):
                            num_residues = len(parents_obj.p0_aligned)
                        elif contacts_result['contacts']:
                            # Fallback: calculate from max contact index
                            max_res = max(max(c[0], c[1]) for c in contacts_result['contacts'])
                            num_residues = max_res + 1
                        else:
                            num_residues = 0
                    elif contacts_result['contacts']:
                        # Fallback: calculate from max contact index
                        max_res = max(max(c[0], c[1]) for c in contacts_result['contacts'])
                        num_residues = max_res + 1
                    else:
                        num_residues = 0
                    num_parents = len(contacts_result['parents'])
                    
                    st.info(f"""
                    **Contact Statistics:**
                    - Number of contacts: {num_contacts}
                    - Number of residues: {num_residues}
                    - Number of parent sequences: {num_parents}
                    """)
                    
                    # Autosave after contacts calculated (manual mode)
                    from utils.session_manager import auto_save
                    auto_save("contacts_calculated")
                    
                    from utils.session_manager import render_save_project_ui
                    render_save_project_ui('save_contacts_btn')
                except Exception as e:
                    st.error(f"Error calculating contacts: {str(e)}")

# Display structure viewer and contact editing if available
contacts_data = st.session_state.get(SESSION_KEYS['schema_contacts'])
if contacts_data and isinstance(contacts_data, dict) and contacts_data.get('contacts') is not None:
    contacts = contacts_data['contacts']
    
    # Check if we have PDB structure and path
    pdb_structure = contacts_data.get('pdb_structure')
    pdb_path = st.session_state.get(SESSION_KEYS['pdb_path'])
    pdb_id = st.session_state.get(SESSION_KEYS['pdb_id'])  # UniProt ID for AlphaFold
    
    if pdb_structure:
        st.header("Structure Viewer & Contact Editor")
        
        try:
            from streamlit_molstar import st_molstar
            from utils.structure_viewer import (
                precompute_nearby_residues, get_contacts_for_residue,
                update_contacts, get_pdb_file_content, create_colored_pdb_content
            )

            # Validate that contact indices match the current structure
            valid_indices = set(aa.index for aa in pdb_structure.amino_acids)
            invalid_contacts = []
            for contact in contacts:
                if len(contact) >= 2:
                    i, j = contact[0], contact[1]
                    if i not in valid_indices or j not in valid_indices:
                        invalid_contacts.append((i, j))
            
            if invalid_contacts:
                st.warning(f"""
                ⚠️ **Contact/Structure Mismatch Detected**
                
                {len(invalid_contacts)} contacts reference residues not in the current structure.
                This usually happens when loading a project saved with different alignment settings.
                
                **Structure residue range:** {min(valid_indices)} - {max(valid_indices)} ({len(valid_indices)} residues)
                **Invalid contact indices include:** {invalid_contacts[:5]}{'...' if len(invalid_contacts) > 5 else ''}
                """)
                
                # Offer to recalculate contacts from the current structure
                if st.button("🔄 Recalculate Contacts from Current Structure", type="primary"):
                    with st.spinner("Recalculating contacts..."):
                        try:
                            from utils.schemarecomb_bridge import contacts_from_pdb_structure
                            new_contacts = contacts_from_pdb_structure(pdb_structure)
                            
                            # Update contacts in contacts_data
                            contacts_data['contacts'] = [(i, j, None, None) for i, j in new_contacts]
                            st.session_state[SESSION_KEYS['schema_contacts']] = contacts_data
                            
                            # Clear edited_contacts to reload from new contacts
                            if 'edited_contacts' in st.session_state:
                                del st.session_state['edited_contacts']
                            
                            st.success(f"✓ Recalculated {len(new_contacts)} contacts from structure")
                            st.rerun()
                        except Exception as e:
                            st.error(f"Error recalculating contacts: {str(e)}")
            
            # Create a unique viewer ID based on structure to avoid key conflicts
            # Use PDB ID or structure hash to make it unique per structure
            import hashlib
            if pdb_id:
                structure_id = pdb_id
            else:
                # Create hash from structure characteristics
                structure_str = f"{len(pdb_structure.amino_acids)}_{pdb_structure.amino_acids[0].index if pdb_structure.amino_acids else 0}"
                structure_id = hashlib.md5(structure_str.encode()).hexdigest()[:8]
            
            # Use a stable key that doesn't change on reruns
            # Include a unique identifier to ensure no conflicts
            # Use structure hash + first few residues to make it unique
            structure_hash = hashlib.md5(
                str([aa.index for aa in pdb_structure.amino_acids[:10]]).encode()
            ).hexdigest()[:8]
            viewer_base_key = f"molstar_{structure_id}_{structure_hash}"
            
            # Create a fingerprint of the current structure's indices to detect changes
            # This ensures we recalculate nearby_residues when the structure changes
            current_structure_fingerprint = str(sorted([aa.index for aa in pdb_structure.amino_acids]))
            
            # Check if we need to recalculate nearby_residues (structure changed or not computed)
            needs_recalc = (
                'nearby_residues' not in st.session_state or
                st.session_state.get('structure_fingerprint') != current_structure_fingerprint
            )
            
            if needs_recalc:
                with st.spinner("Pre-computing distances..."):
                    st.session_state.nearby_residues = precompute_nearby_residues(
                        pdb_structure, distance_threshold=NEARBY_RESIDUE_DISTANCE
                    )
                    st.session_state.structure_fingerprint = current_structure_fingerprint
                    # Also reset edited_contacts since structure changed
                    if 'edited_contacts' in st.session_state:
                        del st.session_state['edited_contacts']
            
            if 'edited_contacts' not in st.session_state:
                # Convert contacts to list of (i, j) tuples for easier manipulation
                # Handle both (i, j) and (i, j, ri, rj) formats
                st.session_state.edited_contacts = []
                for contact in contacts:
                    if len(contact) >= 2:
                        i, j = contact[0], contact[1]
                        # Validate that these indices exist in the current structure
                        valid_indices = set(aa.index for aa in pdb_structure.amino_acids)
                        if i in valid_indices and j in valid_indices:
                            st.session_state.edited_contacts.append((i, j))
                    else:
                        st.session_state.edited_contacts.append(contact)
            
            # Get the sorted list of actual residue indices (alignment column positions)
            # After renumbering, aa.index is the alignment column position, not sequential
            sorted_indices = sorted([aa.index for aa in pdb_structure.amino_acids])
            index_to_aa = {aa.index: aa for aa in pdb_structure.amino_acids}
            num_residues = len(sorted_indices)
            
            # Map from selector position (0, 1, 2...) to actual index
            # This allows the selector to iterate through actual indices
            if 'residue_selector_pos' not in st.session_state:
                st.session_state.residue_selector_pos = 0
            
            # Layout: Structure viewer on left, controls on right
            col1, col2 = st.columns([2, 1])
            
            with col1:
                # Get the actual residue index from the selector position
                selector_pos = st.session_state.residue_selector_pos
                if selector_pos < len(sorted_indices):
                    selected_idx = sorted_indices[selector_pos]
                else:
                    selected_idx = sorted_indices[0] if sorted_indices else 0
                    st.session_state.residue_selector_pos = 0
                nearby = st.session_state.nearby_residues.get(selected_idx, set())
                contacting = get_contacts_for_residue(st.session_state.edited_contacts, selected_idx)
                
                # Display UniProt ID / AlphaFold link above the viewer
                if pdb_id:
                    alphafold_url = f"https://alphafold.ebi.ac.uk/entry/{pdb_id}"
                    st.markdown(f"**Structure:** [{pdb_id}]({alphafold_url}) (AlphaFold)")
                
                # Display structure with Mol*
                # Generate PDB content with color encoding for visualization
                try:
                    pdb_content = create_colored_pdb_content(
                        pdb_structure, selected_idx, nearby, contacting
                    )
                except Exception as e:
                    # If structure conversion fails, try fallback
                    st.warning(f"Could not generate colored PDB from structure object: {str(e)}")
                    # Try to use pdb_path if it exists and is a valid file
                    if pdb_path:
                        try:
                            pdb_content = get_pdb_file_content(pdb_path)
                        except Exception as e2:
                            st.error(f"Could not load PDB content: {str(e2)}")
                            st.info("Please recalculate contacts to regenerate the structure.")
                            pdb_content = None
                    else:
                        st.error("No PDB structure or file available.")
                        pdb_content = None
                
                if not pdb_content:
                    st.error("Could not load PDB structure for visualization.")
                    st.info("Upload a PDB structure or recalculate contacts on this page, then refresh.")
                    st.stop()
                
                # Embedded Mol* viewer (B-factor categories → uncertainty color theme)
                try:
                    from utils.structure_viewer import create_molstar_html_viewer
                    import streamlit.components.v1 as components
                    
                    html_content = create_molstar_html_viewer(
                        pdb_content, pdb_structure, selected_idx, nearby, contacting,
                        viewer_id=f"viewer_{viewer_base_key}"
                    )
                    # Caption bar + Mol* viewport (~620px)
                    components.html(html_content, height=680)
                except Exception as e:
                    # Fall back to basic st_molstar if custom viewer fails
                    st.warning(f"Could not create custom viewer: {str(e)}. Using basic viewer.")
                    import traceback
                    with st.expander("Error details"):
                        st.code(traceback.format_exc())
                    
                    # Use the existing temp_dir from session state if available
                    if 'temp_dir' in st.session_state and st.session_state.temp_dir:
                        temp_dir = st.session_state.temp_dir
                    else:
                        temp_dir = tempfile.mkdtemp()
                        st.session_state.temp_dir = temp_dir
                    
                    # Save PDB content to a file in temp directory (use stable filename)
                    pdb_viewer_path = os.path.join(temp_dir, f"structure_viewer_{structure_id}.pdb")
                    with open(pdb_viewer_path, 'w') as f:
                        f.write(pdb_content)
                    
                    # Create Mol* component with file path
                    st_molstar(
                        pdb_viewer_path,
                        key=viewer_base_key,
                        height=700
                    )
                
                # Instructions
                st.info("💡 Use the controls on the right to select residues and edit contacts.")
            
            with col2:
                st.subheader("Residue Selection")
                
                # Residue selector with dropdown and +/- buttons
                col_minus, col_dropdown, col_plus = st.columns([1, 3, 1])
                
                # Get current position, ensuring it's valid
                current_pos = st.session_state.residue_selector_pos
                if current_pos >= len(sorted_indices):
                    current_pos = 0
                    st.session_state.residue_selector_pos = 0
                
                # Handle - button (previous residue)
                with col_minus:
                    prev_disabled = (current_pos <= 0)
                    if st.button("◀", key="residue_prev", use_container_width=True, 
                                 help="Previous residue", disabled=prev_disabled):
                        st.session_state.residue_selector_pos = current_pos - 1
                        st.rerun()
                
                # Dropdown for direct selection
                with col_dropdown:
                    if sorted_indices:
                        # Format function to show residue number and amino acid
                        def format_residue(idx):
                            if idx in index_to_aa:
                                aa = index_to_aa[idx]
                                return f"{idx} ({aa.letter})"
                            return f"{idx}"
                        
                        # Selectbox without key - controlled entirely by index parameter
                        # This ensures buttons can update the displayed value
                        new_idx = st.selectbox(
                            "Select residue",
                            options=sorted_indices,
                            index=current_pos,
                            label_visibility="collapsed",
                            format_func=format_residue
                        )
                        
                        # If user changed dropdown, update position
                        if new_idx != sorted_indices[current_pos]:
                            st.session_state.residue_selector_pos = sorted_indices.index(new_idx)
                            st.rerun()
                
                # Handle + button (next residue)
                with col_plus:
                    next_disabled = (current_pos >= len(sorted_indices) - 1)
                    if st.button("▶", key="residue_next", use_container_width=True,
                                 help="Next residue", disabled=next_disabled):
                        st.session_state.residue_selector_pos = current_pos + 1
                        st.rerun()
                
                # Show selected residue info (selected_idx was calculated at top of col1)
                if selected_idx in index_to_aa:
                    aa = index_to_aa[selected_idx]
                    st.markdown(f"**Residue {selected_idx}: {aa.name} ({aa.letter})**")
                
                st.markdown("---")
                st.subheader("Contacts")
                
                # Save contacts button - at top for easy access
                if st.button("💾 Save Contacts", type="primary", use_container_width=True):
                    # Convert edited contacts back to original format if needed
                    # Original format might be (i, j, ri, rj) but we're working with (i, j)
                    # For saving, we'll use (i, j, None, None) format
                    formatted_contacts = [(i, j, None, None) for i, j in st.session_state.edited_contacts]
                    
                    # Update session state contacts
                    contacts_data['contacts'] = formatted_contacts
                    st.session_state[SESSION_KEYS['schema_contacts']] = contacts_data
                    
                    # Save to file if temp_dir exists
                    if 'temp_dir' in st.session_state:
                        contact_file_path = os.path.join(st.session_state.temp_dir, "contacts.txt")
                        from utils.schema_wrapper import save_contacts
                        save_contacts(formatted_contacts, contact_file_path)
                        st.success("✓ Contacts saved to file!")
                    else:
                        st.success("✓ Contacts updated in session!")
                    
                    # Autosave after contacts edited
                    from utils.session_manager import auto_save
                    auto_save("contacts_edited")
                
                st.markdown("---")
                
                # Get contacts for selected residue
                contacting_residues = get_contacts_for_residue(
                    st.session_state.edited_contacts, selected_idx
                )
                
                # Get nearby residues (within 8Å)
                nearby_residues = st.session_state.nearby_residues.get(selected_idx, set())
                
                # Show all nearby residues with checkboxes
                # Contacts are checked, others are unchecked
                all_candidates = sorted(list(nearby_residues | contacting_residues))
                
                if all_candidates:
                    for res_j in all_candidates:
                        is_contact = res_j in contacting_residues
                        
                        # Get amino acid info using index_to_aa mapping
                        if res_j in index_to_aa:
                            aa_j = index_to_aa[res_j]
                            label = f"Residue {res_j}: {aa_j.name} ({aa_j.letter})"
                        else:
                            label = f"Residue {res_j}"
                        
                        # Checkbox
                        new_value = st.checkbox(
                            label,
                            value=is_contact,
                            key=f"contact_{selected_idx}_{res_j}"
                        )
                        
                        # Update contacts if changed
                        if new_value != is_contact:
                            st.session_state.edited_contacts = update_contacts(
                                st.session_state.edited_contacts,
                                selected_idx,
                                res_j,
                                new_value
                            )
                            st.rerun()
                else:
                    st.info(f"No nearby residues within {NEARBY_RESIDUE_DISTANCE}Å")
        
        except ImportError:
            st.warning("⚠️ streamlit-molstar not installed. Install with: pip install streamlit-molstar")
            st.info("Structure viewer requires streamlit-molstar. Contact editing will be available after installation.")
        except Exception as e:
            st.error(f"Error displaying structure viewer: {str(e)}")
            import traceback
            with st.expander("Error Details"):
                st.code(traceback.format_exc())

# Display contact map if available
if 'schema_contacts' in st.session_state:
    st.header("Contact Map Visualization")
    
    contacts_data = st.session_state[SESSION_KEYS['schema_contacts']]
    contacts = contacts_data['contacts']
    
    # Calculate num_residues from residues if available, otherwise from contacts or parents
    if contacts_data.get('residues') is not None:
        num_residues = len(contacts_data['residues'])
    elif contacts_data.get('parents_object') is not None:
        # Use aligned sequence length from schemarecomb ParentSequences
        parents_obj = contacts_data['parents_object']
        if hasattr(parents_obj, 'p0_aligned'):
            num_residues = len(parents_obj.p0_aligned)
        elif contacts:
            # Fallback: calculate from max contact index
            max_res = max(max(c[0], c[1]) for c in contacts)
            num_residues = max_res + 1
        else:
            num_residues = None  # Let plot_contact_map infer it
    elif contacts:
        # Fallback: calculate from max contact index
        max_res = max(max(c[0], c[1]) for c in contacts)
        num_residues = max_res + 1
    else:
        num_residues = None  # Let plot_contact_map infer it
    
    fig = plot_contact_map(contacts, num_residues)
    st.plotly_chart(fig, use_container_width=True)
    
    # Save checkpoint option
    with st.expander("💾 Save Project", expanded=False):
        from utils.session_manager import render_save_project_ui
        render_save_project_ui('save_contacts_viz_btn', inside_expander=True)

    # Calculate energies if crossover file is provided
    if 'crossover_file' in st.session_state:
        crossover_file = st.session_state.crossover_file
        if crossover_file:
            st.header("Energy Calculation")
            
            if st.button("Calculate SCHEMA Energies", type="primary"):
                with st.spinner("Calculating SCHEMA energies..."):
                    try:
                        # Save crossover file
                        if 'temp_dir' not in st.session_state:
                            temp_dir = tempfile.mkdtemp()
                            st.session_state.temp_dir = temp_dir
                        else:
                            temp_dir = st.session_state.temp_dir
                        
                        crossover_path = save_uploaded_file(crossover_file, temp_dir)
                        
                        # Save contacts to file temporarily
                        contact_file_path = os.path.join(temp_dir, "contacts.txt")
                        save_contacts(contacts, contact_file_path)
                        
                        # Calculate energies
                        energy_result = calculate_energies(
                            st.session_state.msa_path,
                            contact_file_path,
                            crossover_path,
                            print_E=True,
                            print_m=False
                        )
                        
                        # Store in session state
                        st.session_state.schema_results = energy_result
                        
                        st.success("✓ SCHEMA energies calculated successfully!")
                        
                        # Display results
                        st.subheader("Energy Results")
                        results_df = st.session_state.schema_results['results']
                        
                        if results_df:
                            # Create summary
                            disruptions = [r['disruption'] for r in results_df if r.get('disruption') is not None]
                            if disruptions:
                                st.metric("Average Disruption", f"{sum(disruptions)/len(disruptions):.2f}")
                                st.metric("Min Disruption", f"{min(disruptions):.0f}")
                                st.metric("Max Disruption", f"{max(disruptions):.0f}")
                                
                                # Plot energy distribution
                                fig_dist = plot_energy_distribution(results_df)
                                if fig_dist:
                                    st.plotly_chart(fig_dist, use_container_width=True)
                                
                                # Display table
                                import pandas as pd
                                display_data = []
                                for i, r in enumerate(results_df[:100]):  # Limit to first 100
                                    display_data.append({
                                        'Chimera': r['chimera'],
                                        'Disruption': r['disruption'],
                                        'Sequence': r.get('sequence', '')[:50] + '...' if r.get('sequence') else ''
                                    })
                                df = pd.DataFrame(display_data)
                                st.dataframe(df, use_container_width=True)
                        
                    except Exception as e:
                        st.error(f"Error calculating energies: {str(e)}")

# Download contacts file
if 'schema_contacts' in st.session_state:
    st.header("Export Results")
    
    contacts = st.session_state.schema_contacts['contacts']
    
    # Create download button for contacts
    contact_file_content = ""
    for contact in contacts:
        contact_file_content += f"{contact[0]} {contact[1]} {contact[2]} {contact[3]}\n"
    
    st.download_button(
        label="Download Contacts File",
        data=contact_file_content,
        file_name="schema_contacts.txt",
        mime="text/plain"
    )
