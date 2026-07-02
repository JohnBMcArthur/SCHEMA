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
    export_checkpoint, import_checkpoint, validate_checkpoint
)
from utils.config import CHECKPOINT_DIR, SESSION_KEYS

st.set_page_config(
    page_title="Projects",
    page_icon="📁",
    layout="wide"
)

# Initialize session state
init_session_state()

st.title("📁 Project Management")

st.markdown("""
Manage your saved SCHEMA-RASPP projects. Load, delete, export, and import projects.
""")

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
            st.subheader("Project Actions")
            selected_project = st.selectbox(
                "Select project for actions",
                range(len(filtered_df)),
                format_func=lambda x: f"{filtered_df.iloc[x]['Project Name']} ({filtered_df.iloc[x]['Saved']})",
                key='selected_project_action'
            )
            
            if selected_project is not None:
                project_path = filtered_df.iloc[selected_project]['Path']
                project_name = filtered_df.iloc[selected_project]['Project Name']
                
                col1, col2, col3, col4 = st.columns(4)
                
                with col1:
                    if st.button("📂 Load", use_container_width=True, key='load_selected_project'):
                        try:
                            metadata = load_checkpoint(project_path)
                            st.success(f"✓ Loaded: {metadata.get('project_name', project_name)}")
                            st.rerun()
                        except Exception as e:
                            st.error(f"Error loading: {str(e)}")
                
                with col2:
                    if st.button("📦 Export", use_container_width=True, key='export_selected_project'):
                        try:
                            export_format = st.radio(
                                "Format",
                                ['zip', 'json'],
                                key='export_format_selected',
                                horizontal=True
                            )
                            output_path = export_checkpoint(project_path, format=export_format)
                            st.success(f"✓ Exported: {Path(output_path).name}")
                            
                            # Download button
                            with open(output_path, 'rb') as f:
                                file_data = f.read()
                                st.download_button(
                                    label=f"📥 Download {export_format.upper()}",
                                    data=file_data,
                                    file_name=f"{project_name}.{export_format}",
                                    mime="application/zip" if export_format == 'zip' else "application/json",
                                    key='download_exported_selected'
                                )
                        except Exception as e:
                            st.error(f"Error exporting: {str(e)}")
                
                with col3:
                    if st.button("🗑️ Delete", use_container_width=True, key='delete_selected_project'):
                        st.warning(f"⚠️ Are you sure you want to delete '{project_name}'?")
                        if st.button("Confirm Delete", type="primary", key='confirm_delete'):
                            try:
                                delete_checkpoint(project_path)
                                st.success(f"✓ Deleted: {project_name}")
                                st.rerun()
                            except Exception as e:
                                st.error(f"Error deleting: {str(e)}")
                
                with col4:
                    if st.button("✅ Validate", use_container_width=True, key='validate_selected_project'):
                        try:
                            validate_checkpoint(project_path)
                            st.success("✓ Checkpoint is valid")
                        except Exception as e:
                            st.error(f"Validation failed: {str(e)}")

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
        col1, col2, col3 = st.columns(3)
        
        with col1:
            has_contacts = SESSION_KEYS['schema_contacts'] in st.session_state and st.session_state[SESSION_KEYS['schema_contacts']]
            st.metric("SCHEMA Contacts", "✓" if has_contacts else "✗")
        
        with col2:
            has_energies = SESSION_KEYS['schema_results'] in st.session_state and st.session_state[SESSION_KEYS['schema_results']]
            st.metric("SCHEMA Energies", "✓" if has_energies else "✗")
        
        with col3:
            has_raspp = SESSION_KEYS['raspp_results'] in st.session_state and st.session_state[SESSION_KEYS['raspp_results']]
            st.metric("RASPP Results", "✓" if has_raspp else "✗")
        
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
    Import a project from a ZIP or JSON file. ZIP files contain the complete project
    with all files, while JSON files contain metadata and session data (files must be
    provided separately).
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
            
            try:
                # Save uploaded file temporarily
                temp_dir = tempfile.mkdtemp()
                temp_file = Path(temp_dir) / imported_file.name
                with open(temp_file, 'wb') as f:
                    f.write(imported_file.getbuffer())
                
                # Import checkpoint
                project_path = import_checkpoint(
                    str(temp_file),
                    project_name=import_name if import_name else None
                )
                st.success(f"✓ Project imported: {Path(project_path).name}")
                
                # Show import summary
                metadata = get_checkpoint_metadata(Path(project_path).name)
                if metadata:
                    st.json(metadata)
                
                # Option to load immediately
                if st.button("📂 Load Imported Project", key='load_imported_main'):
                    metadata = load_checkpoint(project_path)
                    st.success(f"✓ Loaded: {metadata.get('project_name', 'Unknown')}")
                    st.rerun()
                    
            except Exception as e:
                st.error(f"Error importing project: {str(e)}")
                import traceback
                with st.expander("Error Details"):
                    st.code(traceback.format_exc())

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
                    if st.button("📦 Export This Project", use_container_width=True, key='export_details_project'):
                        try:
                            export_format = st.radio(
                                "Format",
                                ['zip', 'json'],
                                key='export_format_details',
                                horizontal=True
                            )
                            output_path = export_checkpoint(project_path, format=export_format)
                            st.success(f"✓ Exported: {Path(output_path).name}")
                        except Exception as e:
                            st.error(f"Error: {str(e)}")
                
                with col3:
                    if st.button("✅ Validate", use_container_width=True, key='validate_details_project'):
                        try:
                            validate_checkpoint(project_path)
                            st.success("✓ Checkpoint is valid")
                        except Exception as e:
                            st.error(f"Validation failed: {str(e)}")
            else:
                st.error("Could not load project metadata")
