"""
Page 5: Project Management
Full interface for managing saved projects/checkpoints.
"""

import streamlit as st
import pandas as pd
from pathlib import Path
from datetime import datetime

# Import utilities
import sys
sys.path.insert(0, str(Path(__file__).parent.parent))

from utils.session_manager import (
    init_session_state, list_checkpoints, load_checkpoint,
    delete_checkpoint, get_checkpoint_metadata, get_checkpoint_summary,
    export_checkpoint, import_checkpoint, validate_checkpoint,
    _sanitize_import_filename,
)
from utils.config import CHECKPOINT_DIR, SESSION_KEYS, is_cloud_hosting
from utils.workflow_state import render_cloud_persistence_notice

st.set_page_config(
    page_title="Projects",
    page_icon="📁",
    layout="wide"
)

# Initialize session state
init_session_state()

st.title("📁 Project Management")

st.markdown("""
**Projects** lets you save, load, export, and import SCHEMA-RASPP work across sessions.
Checkpoints capture session state at key workflow stages (contacts, RASPP results,
applied crossovers, diversity pools, oligopool options, and more).

**When to use Projects**

| Situation | Action |
|-----------|--------|
| Starting a new design | Name your project on **SCHEMA Energy** (autosave uses it) or **Save Current Project** here |
| Resuming later | **Load** a checkpoint from **All Projects** |
| Backing up (especially Streamlit Cloud) | **Export** as JSON or ZIP, then download |
| Sharing with a colleague | Export JSON/ZIP and send the file; they **Import** on their instance |
| Cleaning up | **Validate** before trusting an old file; **Delete** when no longer needed |

**Streamlit Cloud note:** server-side checkpoint folders can be lost on redeploy. After
meaningful progress, **export your project** and keep a local copy.

**Tabs**

1. **All Projects** — browse, search, load, validate, export, or delete saved checkpoints.
2. **Save Current Project** — write the current session to disk with an optional description.
3. **Import Project** — upload a previously exported JSON or ZIP checkpoint.
4. **Project Details** — inspect what data a checkpoint contains before loading.

**Tips**

- Export **JSON** for lightweight backups; **ZIP** includes auxiliary files when present.
- Loading a project replaces relevant session keys — save first if you might want to undo.
- The **Workflow Stage** column shows how far the project progressed (contacts, RASPP, etc.).
""")

if is_cloud_hosting():
    render_cloud_persistence_notice()

# Show current project if loaded
if 'current_project' in st.session_state and st.session_state.get('current_project'):
    st.info(f"**Current Project:** {st.session_state['current_project']}")

# Tabs for different operations
tab1, tab2, tab3, tab4 = st.tabs([
    "📋 All Projects",
    "💾 Save Current Project",
    "📥 Import Project",
    "🔍 Project Details"
])

with tab1:
    st.header("All Saved Projects")
    
    # Refresh button
    if st.button("🔄 Refresh List", key='refresh_projects'):
        st.rerun()
    
    # Get all checkpoints
    checkpoints = list_checkpoints()
    
    if not checkpoints:
        st.info("No saved projects yet. Save a project to get started!")
    else:
        # Create summary table
        summary_data = []
        for project_path, metadata in checkpoints:
            summary = get_checkpoint_summary(Path(project_path).name)
            if summary:
                summary_data.append({
                    'Project Name': summary.get('project_name', Path(project_path).name),
                    'Description': summary.get('description', ''),
                    'Workflow Stage': summary.get('workflow_stage', 'unknown'),
                    'Saved': summary.get('timestamp', '')[:10] if summary.get('timestamp') else '',
                    'Has Contacts': '✓' if summary.get('data_summary', {}).get('has_contacts') else '',
                    'Has Energies': '✓' if summary.get('data_summary', {}).get('has_energies') else '',
                    'Has RASPP': '✓' if summary.get('data_summary', {}).get('has_raspp') else '',
                    'Path': project_path
                })
        
        if summary_data:
            df = pd.DataFrame(summary_data)
            
            # Search and filter
            col1, col2 = st.columns(2)
            with col1:
                search_term = st.text_input("🔍 Search projects", key='search_projects')
            with col2:
                filter_stage = st.selectbox(
                    "Filter by stage",
                    ['All', 'contacts_calculated', 'energies_calculated', 'raspp_complete'],
                    key='filter_stage'
                )
            
            # Apply filters
            filtered_df = df.copy()
            if search_term:
                filtered_df = filtered_df[
                    filtered_df['Project Name'].str.contains(search_term, case=False, na=False) |
                    filtered_df['Description'].str.contains(search_term, case=False, na=False)
                ]
            if filter_stage != 'All':
                filtered_df = filtered_df[filtered_df['Workflow Stage'] == filter_stage]
            
            # Display table
            st.dataframe(
                filtered_df.drop(columns=['Path']),
                use_container_width=True,
                height=400
            )
            
            # Action buttons for selected project
            selected_project = st.selectbox(
                "Select project for actions",
                range(len(filtered_df)),
                format_func=lambda x: f"{filtered_df.iloc[x]['Project Name']} ({filtered_df.iloc[x]['Saved']})",
                key='selected_project_action'
            )
            
            if selected_project is not None:
                project_path = filtered_df.iloc[selected_project]['Path']
                project_name = filtered_df.iloc[selected_project]['Project Name']
                
                st.subheader("Project Actions")
                
                col1, col2 = st.columns(2)
                with col1:
                    if st.button("📂 Load", use_container_width=True, key='load_selected_project'):
                        try:
                            metadata = load_checkpoint(project_path)
                            st.success(f"✓ Loaded: {metadata.get('project_name', project_name)}")
                            st.rerun()
                        except Exception as e:
                            st.error(f"Error loading: {str(e)}")
                
                with col2:
                    if st.button("✅ Validate", use_container_width=True, key='validate_selected_project'):
                        try:
                            validate_checkpoint(project_path)
                            st.success("✓ Checkpoint is valid")
                        except Exception as e:
                            st.error(f"Validation failed: {str(e)}")
                
                st.markdown("**Export**")
                export_format = st.radio(
                    "Format",
                    ['zip', 'json'],
                    key='export_format_selected',
                    horizontal=True,
                    help="JSON is recommended for backups on Streamlit Cloud.",
                )
                if st.button("📦 Prepare export", use_container_width=True, key='export_selected_project'):
                    try:
                        output_path = export_checkpoint(project_path, format=export_format)
                        st.session_state['projects_export_ready'] = {
                            'path': output_path,
                            'name': project_name,
                            'format': export_format,
                        }
                    except Exception as e:
                        st.error(f"Error exporting: {str(e)}")
                
                export_ready = st.session_state.get('projects_export_ready')
                if export_ready and export_ready.get('name') == project_name:
                    try:
                        with open(export_ready['path'], 'rb') as f:
                            file_data = f.read()
                        fmt = export_ready['format']
                        st.download_button(
                            label=f"📥 Download {fmt.upper()}",
                            data=file_data,
                            file_name=f"{project_name}.{fmt}",
                            mime="application/zip" if fmt == 'zip' else "application/json",
                            key='download_exported_selected',
                            use_container_width=True,
                        )
                    except Exception as e:
                        st.error(f"Could not read export file: {str(e)}")
                
                st.markdown("**Delete**")
                pending_delete = st.session_state.get('projects_pending_delete')
                if pending_delete != project_path:
                    if st.button("🗑️ Delete project", use_container_width=True, key='delete_selected_project'):
                        st.session_state['projects_pending_delete'] = project_path
                        st.rerun()
                else:
                    st.warning(f"⚠️ Delete **{project_name}**? This cannot be undone.")
                    c1, c2 = st.columns(2)
                    with c1:
                        if st.button("Confirm delete", type="primary", key='confirm_delete'):
                            try:
                                delete_checkpoint(project_path)
                                st.session_state.pop('projects_pending_delete', None)
                                st.session_state.pop('projects_export_ready', None)
                                st.success(f"✓ Deleted: {project_name}")
                                st.rerun()
                            except Exception as e:
                                st.error(f"Error deleting: {str(e)}")
                    with c2:
                        if st.button("Cancel", key='cancel_delete'):
                            st.session_state.pop('projects_pending_delete', None)
                            st.rerun()

with tab2:
    st.header("Save Current Project")
    
    # Check if there's data to save
    has_data = (
        SESSION_KEYS['schema_contacts'] in st.session_state and st.session_state[SESSION_KEYS['schema_contacts']] or
        SESSION_KEYS['raspp_results'] in st.session_state and st.session_state[SESSION_KEYS['raspp_results']] or
        SESSION_KEYS['schema_results'] in st.session_state and st.session_state[SESSION_KEYS['schema_results']]
    )
    
    if not has_data:
        st.warning("No data to save. Please run calculations on other pages first.")
    else:
        col1, col2 = st.columns(2)
        
        with col1:
            project_name = st.text_input(
                "Project Name",
                value=st.session_state.get('current_project', ''),
                key='save_new_project_name'
            )
        
        with col2:
            description = st.text_area(
                "Description (optional)",
                key='save_new_project_description',
                height=100
            )
        
        # Show what will be saved
        st.subheader("Data to Save")
        col1, col2, col3, col4 = st.columns(4)
        
        with col1:
            has_contacts = SESSION_KEYS['schema_contacts'] in st.session_state and st.session_state[SESSION_KEYS['schema_contacts']]
            st.metric("SCHEMA Contacts", "✓" if has_contacts else "✗")
        
        with col2:
            has_energies = SESSION_KEYS['schema_results'] in st.session_state and st.session_state[SESSION_KEYS['schema_results']]
            st.metric("SCHEMA Energies", "✓" if has_energies else "✗")
        
        with col3:
            has_raspp = SESSION_KEYS['raspp_results'] in st.session_state and st.session_state[SESSION_KEYS['raspp_results']]
            st.metric("RASPP Results", "✓" if has_raspp else "✗")
        
        with col4:
            has_multi = SESSION_KEYS['multi_fragment_results'] in st.session_state and st.session_state[SESSION_KEYS['multi_fragment_results']]
            st.metric("Multi-Fragment", "✓" if has_multi else "✗")

        col5, col6, col7, col8 = st.columns(4)
        with col5:
            has_xo = bool(st.session_state.get("selected_crossover_positions"))
            st.metric("Applied Crossovers", "✓" if has_xo else "✗")
        with col6:
            has_pools = bool(st.session_state.get("diversity_saved_selections"))
            st.metric("Saved Pools", "✓" if has_pools else "✗")
        with col7:
            has_div = bool(st.session_state.get("diversity_analysis_result"))
            st.metric("Diversity Analysis", "✓" if has_div else "✗")
        with col8:
            has_opt = bool(st.session_state.get("library_opt_results"))
            st.metric("Library Optimization", "✓" if has_opt else "✗")
        
        if st.button("💾 Save Project", type="primary", key='save_new_project_btn'):
            from utils.session_manager import save_checkpoint
            try:
                checkpoint_path = save_checkpoint(
                    project_name=project_name if project_name else None,
                    description=description
                )
                st.success(f"✓ Project saved: {Path(checkpoint_path).name}")
                st.session_state['current_project'] = project_name if project_name else Path(checkpoint_path).name
                st.session_state['current_project_path'] = checkpoint_path
                st.rerun()
            except Exception as e:
                st.error(f"Error saving project: {str(e)}")

with tab3:
    st.header("Import Project")
    
    st.markdown("""
    Import a project from a **JSON** or **ZIP** file exported from this app.
    JSON is recommended for backups on Streamlit Cloud.
    """)
    
    imported_file = st.file_uploader(
        "Upload Project File",
        type=['zip', 'json'],
        key='import_project_file_main'
    )
    
    if imported_file:
        st.info(f"File: {imported_file.name} ({imported_file.size / 1024:.1f} KB)")
        
        import_name = st.text_input(
            "Project Name (optional, defaults to filename)",
            value=imported_file.name.rsplit('.', 1)[0],
            key='import_project_name_main'
        )
        
        if st.button("📥 Import Project", type="primary", key='import_project_btn_main'):
            from utils.session_manager import import_checkpoint
            import tempfile
            import shutil
            
            try:
                temp_dir = tempfile.mkdtemp(prefix="schema_import_")
                try:
                    safe_name = _sanitize_import_filename(imported_file.name)
                    temp_file = Path(temp_dir) / safe_name
                    with open(temp_file, 'wb') as f:
                        f.write(imported_file.getbuffer())
                    
                    project_path = import_checkpoint(
                        str(temp_file),
                        project_name=import_name if import_name else None
                    )
                finally:
                    shutil.rmtree(temp_dir, ignore_errors=True)
                
                st.success(f"✓ Project imported: {Path(project_path).name}")
                
                metadata = get_checkpoint_metadata(Path(project_path).name)
                if metadata:
                    st.json(metadata)
                
                if st.button("📂 Load Imported Project", key='load_imported_main'):
                    metadata = load_checkpoint(project_path)
                    st.success(f"✓ Loaded: {metadata.get('project_name', 'Unknown')}")
                    st.rerun()
                    
            except Exception as e:
                st.error(f"Error importing project: {str(e)}")

with tab4:
    st.header("Project Details")
    
    # Select project to view
    checkpoints = list_checkpoints()
    
    if not checkpoints:
        st.info("No saved projects to view.")
    else:
        project_names = [f"{meta.get('project_name', Path(path).name)}" 
                        for path, meta in checkpoints]
        project_paths = [path for path, meta in checkpoints]
        
        selected_idx = st.selectbox(
            "Select project to view details",
            range(len(project_names)),
            format_func=lambda x: project_names[x],
            key='view_project_details'
        )
        
        if selected_idx is not None:
            project_path = project_paths[selected_idx]
            metadata = get_checkpoint_metadata(Path(project_path).name)
            
            if metadata:
                st.subheader(f"Project: {metadata.get('project_name', 'Unknown')}")
                
                # Basic info
                col1, col2 = st.columns(2)
                
                with col1:
                    st.markdown("**Basic Information**")
                    st.write(f"**Description:** {metadata.get('description', 'None')}")
                    st.write(f"**Created:** {metadata.get('timestamp', 'Unknown')}")
                    st.write(f"**Version:** {metadata.get('version', 'Unknown')}")
                    st.write(f"**Workflow Stage:** {metadata.get('workflow_stage', 'Unknown')}")
                
                with col2:
                    st.markdown("**Data Summary**")
                    data_summary = metadata.get('data_summary', {})
                    st.write(f"**Has Contacts:** {'✓' if data_summary.get('has_contacts') else '✗'}")
                    st.write(f"**Has Energies:** {'✓' if data_summary.get('has_energies') else '✗'}")
                    st.write(f"**Has RASPP:** {'✓' if data_summary.get('has_raspp') else '✗'}")
                    st.write(f"**Number of Parents:** {data_summary.get('num_parents', 0)}")
                    st.write(f"**Number of Contacts:** {data_summary.get('num_contacts', 0)}")
                
                # Parameters
                st.subheader("Parameters")
                parameters = metadata.get('parameters', {})
                if parameters:
                    params_df = pd.DataFrame(list(parameters.items()), columns=['Parameter', 'Value'])
                    st.dataframe(params_df, use_container_width=True)
                else:
                    st.info("No parameters saved")
                
                # Files
                st.subheader("Files")
                files = metadata.get('files', {})
                if files:
                    files_df = pd.DataFrame(list(files.items()), columns=['Type', 'Path'])
                    st.dataframe(files_df, use_container_width=True)
                else:
                    st.info("No files information")
                
                # Full metadata JSON
                with st.expander("View Full Metadata (JSON)"):
                    st.json(metadata)
                
                # Actions
                st.subheader("Actions")
                col1, col2, col3 = st.columns(3)
                
                with col1:
                    if st.button("📂 Load This Project", use_container_width=True, key='load_details_project'):
                        try:
                            metadata = load_checkpoint(project_path)
                            st.success(f"✓ Loaded: {metadata.get('project_name', 'Unknown')}")
                            st.rerun()
                        except Exception as e:
                            st.error(f"Error: {str(e)}")
                
                with col2:
                    export_format_details = st.radio(
                        "Export format",
                        ['zip', 'json'],
                        key='export_format_details',
                        horizontal=True,
                    )
                    if st.button("📦 Export This Project", use_container_width=True, key='export_details_project'):
                        try:
                            output_path = export_checkpoint(project_path, format=export_format_details)
                            st.session_state['projects_details_export'] = {
                                'path': output_path,
                                'format': export_format_details,
                                'name': metadata.get('project_name', 'project'),
                            }
                            st.success(f"✓ Exported: {Path(output_path).name}")
                        except Exception as e:
                            st.error(f"Error: {str(e)}")
                    details_export = st.session_state.get('projects_details_export')
                    if details_export and details_export.get('path'):
                        try:
                            with open(details_export['path'], 'rb') as f:
                                data = f.read()
                            fmt = details_export['format']
                            st.download_button(
                                f"📥 Download {fmt.upper()}",
                                data=data,
                                file_name=f"{details_export['name']}.{fmt}",
                                mime="application/zip" if fmt == 'zip' else "application/json",
                                key='download_details_export',
                                use_container_width=True,
                            )
                        except Exception as e:
                            st.error(str(e))
                
                with col3:
                    if st.button("✅ Validate", use_container_width=True, key='validate_details_project'):
                        try:
                            validate_checkpoint(project_path)
                            st.success("✓ Checkpoint is valid")
                        except Exception as e:
                            st.error(f"Validation failed: {str(e)}")
            else:
                st.error("Could not load project metadata")
