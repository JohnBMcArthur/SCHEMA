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
from utils.session_manager import init_session_state
from utils.workflow_state import clear_downstream_from_contacts
from utils.config import DEFAULTS, SESSION_KEYS, SCHEMA_DEBUG

st.set_page_config(
    page_title="SCHEMA Energy",
    page_icon="⚡",
    layout="wide"
)

# Initialize session state
init_session_state()

st.title("⚡ SCHEMA Energy Calculation")

st.markdown("""
This page identifies **SCHEMA contacts** — pairs of residues whose interaction would be
disrupted if they end up on opposite sides of a crossover in a chimera. Those contacts
power **RASPP Design** (step 2) and the rest of the pipeline.

After contacts are calculated, a **contact map** appears below. You can **download**
the contacts file for your records. In **manual** mode you may also upload a crossover
points file to estimate SCHEMA disruption energies for specific chimeras (optional).

Choose an input method below.
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
    help="Automated: BLAST + MUSCLE + AlphaFold from a single sequence. Manual: bring your own PDB and MSA.",
)

if workflow_mode == "Automated (Enter Sequence)":
    st.markdown("""
**Automated workflow**

Best when you have a **query protein sequence** but not yet an MSA or structure. The app
will:

1. **BLAST** your sequence against the AlphaFold database (EBI) to find diverse homologs
   (identity filters are adjustable in *Automated Workflow Parameters*).
2. **Align** the query plus homologs with EBI MUSCLE to build an MSA.
3. **Download** an AlphaFold PDB for the top BLAST hit (or the best match to your query).
4. **Calculate SCHEMA contacts** by aligning the structure to the MSA and detecting
   residue pairs within the contact distance threshold (default 5.0 Å).

**What you need:** a project name (recommended), your query sequence (≥ 20 amino acids),
and patience — BLAST and alignment can take several minutes.

**If no structure is found:** upload a PDB manually in the section below; the MSA from
step 2 is already in session.

**Next step:** go to **2. RASPP Design** once contacts are calculated successfully.
    """)

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
                        
                        if SCHEMA_DEBUG:
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
                            clear_downstream_from_contacts()
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
                clear_downstream_from_contacts()
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
    st.markdown("""
**Manual workflow**

Best when you already have a **PDB structure** and **multiple sequence alignment (MSA)**
for your parents. The app aligns the structure to the MSA (via schemarecomb), trims to
the shared domain if needed, and computes SCHEMA contacts.

**What you need:**

- **PDB file** — structure of the query (or reference) protein.
- **MSA file** — aligned parent sequences (`.txt`, `.fasta`, `.aln`, etc.).
- **Crossover points file** *(optional)* — whitespace-separated crossover positions if you
  want SCHEMA **disruption energies** for specific chimeras on this page. If you skip this,
  you can still proceed to RASPP and choose crossovers later.

**Steps:**

1. Enter a project name (recommended for autosave).
2. Upload PDB and MSA, validate messages appear, then click **Calculate SCHEMA Contacts**.
3. Review the **contact map** below; download the contacts file if needed.
4. *(Optional)* If you uploaded a crossover file, click **Calculate SCHEMA Energies**.

**Next step:** go to **2. RASPP Design** once contacts are in session.
    """)

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
                    clear_downstream_from_contacts()
                    
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
                except Exception as e:
                    st.error(f"Error calculating contacts: {str(e)}")

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
