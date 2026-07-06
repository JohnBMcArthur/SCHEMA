"""
SCHEMA-RASPP Streamlit Application
Main entry point — navigation router for the multipage app.
"""

import streamlit as st
from pathlib import Path

st.set_page_config(
    page_title="SCHEMA-RASPP",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

from utils.session_manager import init_session_state, show_autoload_notice
from utils.workflow_state import render_cloud_persistence_notice
from utils.config import is_cloud_hosting

init_session_state()
show_autoload_notice()

_PROJECTS_PAGE = st.Page("pages/5_Projects.py", title="Projects", icon="📁")

_PAGES = [
    st.Page("pages/0_Introduction.py", title="Introduction", icon="🧬"),
    st.Page("pages/1_SCHEMA_Energy.py", title="1. SCHEMA Energy", icon="⚡"),
    st.Page("pages/2_RASPP_Design.py", title="2. RASPP Design", icon="🔬"),
    st.Page("pages/3_Crossover_Analysis.py", title="3. Crossover Analysis", icon="🧩"),
    st.Page("pages/4_Assembly_Analysis.py", title="4. Assembly Analysis", icon="🔧"),
    st.Page("pages/6_Diversity_Analysis.py", title="5. Diversity Analysis", icon="🌿"),
    st.Page("pages/8_Library_Optimization.py", title="6. Library Optimization", icon="📊"),
    st.Page("pages/9_Simulate_with_AI.py", title="7. Simulate with AI", icon="🤖"),
    st.Page("pages/7_Oligopool_Design.py", title="8. Oligopool Design", icon="🧪"),
    _PROJECTS_PAGE,
]

pg = st.navigation(_PAGES)

pg.run()

with st.sidebar:
    if is_cloud_hosting():
        render_cloud_persistence_notice()
    st.markdown("---")
    st.header("📁 Projects")

    if st.session_state.get("current_project"):
        st.info(f"**Current Project:** {st.session_state['current_project']}")

    with st.form("quick_save_project", clear_on_submit=False):
        project_name = st.text_input(
            "Project Name",
            value=st.session_state.get("current_project", ""),
            key="save_project_name_input",
        )

        description = st.text_area(
            "Description (optional)",
            key="save_project_description_input",
            height=80,
        )

        if st.form_submit_button("💾 Save Project", use_container_width=True):
            from utils.session_manager import save_checkpoint

            try:
                checkpoint_path = save_checkpoint(
                    project_name=project_name if project_name else None,
                    description=description,
                )
                st.success(f"✓ Project saved: {Path(checkpoint_path).name}")
                st.session_state["current_project"] = (
                    project_name if project_name else Path(checkpoint_path).name
                )
                st.session_state["current_project_path"] = checkpoint_path
                st.rerun()
            except Exception as e:
                st.error(f"Error saving project: {str(e)}")

    st.markdown("---")
    st.subheader("Load Project")
    from utils.session_manager import list_checkpoints, load_checkpoint, clear_session

    checkpoints = list_checkpoints()

    if checkpoints:
        checkpoint_names = [
            f"{meta.get('project_name', Path(path).name)} - {meta.get('timestamp', '')[:10]}"
            for path, meta in checkpoints[:10]
        ]
        checkpoint_paths = [path for path, meta in checkpoints[:10]]

        selected_idx = st.selectbox(
            "Select project to load",
            range(len(checkpoint_names)),
            format_func=lambda x: checkpoint_names[x],
            key="load_project_select",
        )

        if st.button("📂 Load Project", use_container_width=True):
            try:
                metadata = load_checkpoint(checkpoint_paths[selected_idx])
                st.success(
                    f"✓ Project loaded: {metadata.get('project_name', 'Unknown')}"
                )
                st.rerun()
            except Exception as e:
                st.error(f"Error loading project: {str(e)}")

        if st.button("🧹 Clear Session", use_container_width=True):
            clear_session()
            st.success("Session cleared. Start a new project or load an existing one.")
            st.rerun()
    else:
        st.info("No saved projects yet")

    st.markdown("---")
    st.page_link(_PROJECTS_PAGE, label="📋 Manage All Projects →")
