"""
Session state management utilities for SCHEMA-RASPP application.
Handles initialization, cleanup, and validation of session state.
"""

import streamlit as st
from utils.config import SESSION_KEYS


def init_session_state():
    """
    Initialize all session state variables with default values.
    Should be called at the start of each page.
    """
    defaults = {
        SESSION_KEYS['schema_contacts']: None,
        SESSION_KEYS['schema_results']: None,
        SESSION_KEYS['raspp_results']: None,
        SESSION_KEYS['raspp_raw_results']: None,
        SESSION_KEYS['raspp_parents']: None,
        SESSION_KEYS['multi_fragment_results']: None,
        SESSION_KEYS['msa_path']: None,
        SESSION_KEYS['pdb_path']: None,
        SESSION_KEYS['pdb_id']: None,
        SESSION_KEYS['sequences']: None,
        SESSION_KEYS['temp_dir']: None,
        SESSION_KEYS['crossover_file']: None,
        SESSION_KEYS['uploaded_files']: {},
    }

    for key, value in defaults.items():
        if key not in st.session_state:
            st.session_state[key] = value

    try_autoload_last_project()
    show_autoload_notice()


def has_required_data(workflow_step):
    """
    Check if required data is present for a workflow step.
    
    Args:
        workflow_step: One of 'contacts', 'energies', 'raspp'
        
    Returns:
        tuple: (has_data, error_message)
    """
    if workflow_step == 'contacts':
        if SESSION_KEYS['schema_contacts'] not in st.session_state or \
           st.session_state[SESSION_KEYS['schema_contacts']] is None:
            return False, "No SCHEMA contacts calculated. Please calculate contacts first."
        return True, None
    
    elif workflow_step == 'energies':
        if SESSION_KEYS['schema_contacts'] not in st.session_state or \
           st.session_state[SESSION_KEYS['schema_contacts']] is None:
            return False, "No SCHEMA contacts available. Please calculate contacts first."
        if SESSION_KEYS['crossover_file'] not in st.session_state or \
           st.session_state[SESSION_KEYS['crossover_file']] is None:
            return False, "No crossover file uploaded. Please upload crossover points file."
        return True, None
    
    elif workflow_step == 'raspp':
        if SESSION_KEYS['schema_contacts'] not in st.session_state or \
           st.session_state[SESSION_KEYS['schema_contacts']] is None:
            return False, "No SCHEMA contacts available. Please calculate contacts first."
        return True, None
    
    return False, f"Unknown workflow step: {workflow_step}"


# Checkpoint functionality
import pickle
import shutil
import json
import hashlib
import zipfile
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Tuple, Optional, Any, Union

from utils.config import (
    CHECKPOINT_VERSION,
    CHECKPOINT_DIR as CONFIG_CHECKPOINT_DIR,
    SESSION_KEYS,
    AUTO_SAVE_ENABLED,
    AUTO_LOAD_ON_STARTUP,
    AUTOSAVE_MAX_CHECKPOINTS,
)
from utils.checkpoint_serialization import (
    serialize_raspp_raw_results, deserialize_raspp_raw_results,
    make_json_safe, restore_from_json_safe
)

CHECKPOINT_DIR = Path(CONFIG_CHECKPOINT_DIR)
CHECKPOINT_DIR.mkdir(exist_ok=True)
LAST_PROJECT_FILE = CHECKPOINT_DIR / ".last_project"
AUTOLOAD_SUPPRESSED_FILE = CHECKPOINT_DIR / ".autoload_suppressed"
DIVERSITY_MSA_BUNDLE_NAME = "diversity_msa.fasta"

# Session keys persisted in checkpoint.pkl (domain state, not Streamlit widgets).
CHECKPOINT_SESSION_KEYS: List[str] = list(SESSION_KEYS.values()) + [
    "parents_object_json",
    "query_sequence",
    "contact_distance",
    "chains",
    "num_crossovers",
    "min_fragment_diversity",
    "num_sequences",
    "min_identity",
    "max_identity",
    "min_fragments",
    "max_fragments",
    "restrict_crossovers_at_graph_marks",
    "msa_seg_conservation_min",
    "msa_seg_diversity_max",
    "selected_crossover_positions",
    "crossover_manual_input",
    "_xo_applied_revision",
    "project_name_for_autosave",
    "diversity_analysis_result",
    "diversity_table_rows",
    "diversity_analysis_upload_hash",
    "diversity_analysis_fasta_text",
    "diversity_msa_query_seq_id",
    "diversity_applied_filters",
    "diversity_draft_filters",
    "diversity_ssn_threshold",
    "diversity_network_preview",
    "diversity_excluded_rows",
    "diversity_saved_selections",
    "diversity_msa_checkpoint_path",
    "library_opt_results",
    "oligopool_forward_primer",
    "oligopool_reverse_primer",
    "oligopool_max_length",
]

CHECKPOINT_PATH_KEYS = {
    SESSION_KEYS["msa_path"],
    SESSION_KEYS["pdb_path"],
    SESSION_KEYS["crossover_file"],
    "diversity_msa_checkpoint_path",
}


def _write_last_project(project_dir: Union[str, Path]) -> None:
    try:
        with open(LAST_PROJECT_FILE, "w", encoding="utf-8") as f:
            f.write(str(project_dir))
        if AUTOLOAD_SUPPRESSED_FILE.exists():
            AUTOLOAD_SUPPRESSED_FILE.unlink()
    except Exception:
        pass


def _session_has_workflow_data() -> bool:
    """True if the session already contains restorable workflow data."""
    if st.session_state.get(SESSION_KEYS["schema_contacts"]):
        return True
    if st.session_state.get(SESSION_KEYS["multi_fragment_results"]):
        return True
    if st.session_state.get(SESSION_KEYS["raspp_results"]):
        return True
    if st.session_state.get("selected_crossover_positions"):
        return True
    if st.session_state.get("diversity_analysis_result"):
        return True
    if st.session_state.get("diversity_saved_selections"):
        return True
    if st.session_state.get("library_opt_results"):
        return True
    return False


def _resolve_autoload_path() -> Optional[Path]:
    if AUTOLOAD_SUPPRESSED_FILE.exists():
        return None

    if LAST_PROJECT_FILE.exists():
        try:
            raw = LAST_PROJECT_FILE.read_text(encoding="utf-8").strip()
            if raw:
                candidate = Path(raw)
                if not candidate.is_absolute():
                    candidate = CHECKPOINT_DIR / candidate
                if (candidate / "checkpoint.pkl").exists():
                    return candidate
        except Exception:
            pass

    checkpoints = list_checkpoints()
    if checkpoints:
        return Path(checkpoints[0][0])
    return None


def try_autoload_last_project() -> Optional[Dict[str, Any]]:
    """Restore the last project once per browser session when session is empty."""
    if st.session_state.get("_autoload_attempted"):
        return None
    st.session_state["_autoload_attempted"] = True

    if not AUTO_LOAD_ON_STARTUP:
        return None
    if _session_has_workflow_data():
        return None

    project_path = _resolve_autoload_path()
    if project_path is None:
        return None

    try:
        metadata = load_checkpoint(project_path, validate=False)
        st.session_state["_autoload_message"] = (
            f"Restored project **{metadata.get('project_name', project_path.name)}** "
            f"(saved {str(metadata.get('timestamp', ''))[:19]})."
        )
        return metadata
    except Exception as exc:
        st.session_state["_autoload_error"] = str(exc)
        return None


def show_autoload_notice() -> None:
    message = st.session_state.pop("_autoload_message", None)
    if message:
        st.success(message)
    error = st.session_state.pop("_autoload_error", None)
    if error:
        st.warning(f"Could not restore last project: {error}")


def _shrink_library_opt_results(results: Any) -> Any:
    if not isinstance(results, dict):
        return results
    shrunk = dict(results)
    chimeras = shrunk.get("chimeras")
    if isinstance(chimeras, list):
        shrunk["chimeras"] = [
            {
                "chimera_id": c.get("chimera_id"),
                "esm_score": c.get("esm_score"),
                "sequence_length": len(c.get("sequence") or ""),
                "picks": c.get("picks"),
            }
            for c in chimeras
            if isinstance(c, dict)
        ]
    return shrunk


def _prepare_snapshot_value(key: str, value: Any) -> Any:
    if key == SESSION_KEYS["raspp_raw_results"]:
        return serialize_raspp_raw_results(value)
    if key == "library_opt_results":
        return _shrink_library_opt_results(value)
    return value


def _collect_session_snapshot() -> Dict[str, Any]:
    snapshot: Dict[str, Any] = {}
    for key in CHECKPOINT_SESSION_KEYS:
        if key not in st.session_state:
            continue
        value = st.session_state[key]
        if value is None:
            continue
        snapshot[key] = _prepare_snapshot_value(key, value)
    return snapshot


def _normalize_state_after_load(project_dir: Path) -> None:
    excluded = st.session_state.get("diversity_excluded_rows")
    if isinstance(excluded, dict):
        st.session_state["diversity_excluded_rows"] = {
            int(k) if str(k).isdigit() else k: set(v) if not isinstance(v, set) else v
            for k, v in excluded.items()
        }

    div_path = st.session_state.get("diversity_msa_checkpoint_path")
    if div_path:
        path = Path(div_path)
        if not path.is_absolute():
            path = project_dir / path
        elif not path.exists() and (project_dir / path.name).exists():
            path = project_dir / path.name
        if path.exists():
            st.session_state["diversity_msa_checkpoint_path"] = str(path)
            st.session_state["diversity_analysis_fasta_bytes"] = path.read_bytes()
            if not st.session_state.get("diversity_analysis_fasta_text"):
                st.session_state["diversity_analysis_fasta_text"] = path.read_text(
                    encoding="utf-8", errors="replace"
                )

    if st.session_state.get("selected_crossover_positions"):
        st.session_state["_xo_sync_checkboxes"] = True


def clear_session():
    """
    Clear all application-related session state, including current project.
    Intended for a full reset before starting a new project or loading one.
    """
    # Core app/session keys
    keys_to_clear = list(SESSION_KEYS.values()) + [
        "current_project",
        "current_project_path",
        "project_name_for_autosave",
        "query_sequence",
        "parents_object_json",
        "selected_crossover_positions",
        "crossover_manual_input",
        "_xo_sync_checkboxes",
        "_xo_eligible_positions",
        "_xo_parent_length",
        "_xo_freq_eligible_sig",
        "_xo_chart_pending",
        "_xo_applied_revision",
        "diversity_analysis_result",
        "diversity_table_rows",
        "diversity_analysis_upload_hash",
        "diversity_analysis_fasta_bytes",
        "diversity_analysis_fasta_text",
        "diversity_msa_query_seq_id",
        "diversity_msa_checkpoint_path",
        "diversity_applied_filters",
        "diversity_draft_filters",
        "diversity_graph_visible",
        "diversity_ssn_threshold",
        "diversity_network_preview",
        "diversity_excluded_rows",
        "diversity_saved_selections",
        "diversity_ssn_cache",
        "library_opt_results",
        "_diversity_fragment_indices",
        "_diversity_prev_fragment",
        "_autoload_message",
        "_autoload_error",
    ]
    for key in keys_to_clear:
        if key in st.session_state:
            del st.session_state[key]

    for key in list(st.session_state.keys()):
        if key.startswith("xo_pick_"):
            del st.session_state[key]

    st.session_state["_autoload_attempted"] = True

    try:
        if LAST_PROJECT_FILE.exists():
            LAST_PROJECT_FILE.unlink()
        AUTOLOAD_SUPPRESSED_FILE.write_text("1", encoding="utf-8")
    except Exception:
        pass


def _get_workflow_stage() -> str:
    """Determine current workflow stage from session state."""
    if st.session_state.get("library_opt_results"):
        return "library_opt_complete"
    if st.session_state.get("diversity_saved_selections") or st.session_state.get(
        "diversity_analysis_result"
    ):
        return "diversity_complete"
    if st.session_state.get("selected_crossover_positions"):
        return "crossovers_applied"
    if SESSION_KEYS["raspp_results"] in st.session_state and st.session_state[
        SESSION_KEYS["raspp_results"]
    ]:
        return "raspp_complete"
    elif SESSION_KEYS["schema_results"] in st.session_state and st.session_state[
        SESSION_KEYS["schema_results"]
    ]:
        return "energies_calculated"
    elif SESSION_KEYS["schema_contacts"] in st.session_state and st.session_state[
        SESSION_KEYS["schema_contacts"]
    ]:
        return "contacts_calculated"
    else:
        return "initialized"


def _get_data_summary() -> Dict[str, Any]:
    """Get summary of data in session state."""
    summary = {
        'has_contacts': False,
        'has_energies': False,
        'has_raspp': False,
        'has_crossovers': False,
        'has_diversity': False,
        'has_library_opt': False,
        'num_parents': 0,
        'num_contacts': 0,
    }
    
    if SESSION_KEYS['schema_contacts'] in st.session_state and st.session_state[SESSION_KEYS['schema_contacts']]:
        summary['has_contacts'] = True
        contacts_data = st.session_state[SESSION_KEYS['schema_contacts']]
        if isinstance(contacts_data, dict):
            if 'contacts' in contacts_data:
                summary['num_contacts'] = len(contacts_data['contacts'])
            if 'parents' in contacts_data:
                summary['num_parents'] = len(contacts_data['parents'])
    
    if SESSION_KEYS['schema_results'] in st.session_state and st.session_state[SESSION_KEYS['schema_results']]:
        summary['has_energies'] = True
    
    if SESSION_KEYS['raspp_results'] in st.session_state and st.session_state[SESSION_KEYS['raspp_results']]:
        summary['has_raspp'] = True

    if st.session_state.get('selected_crossover_positions'):
        summary['has_crossovers'] = True

    if st.session_state.get('diversity_analysis_result') or st.session_state.get(
        'diversity_saved_selections'
    ):
        summary['has_diversity'] = True

    if st.session_state.get('library_opt_results'):
        summary['has_library_opt'] = True

    return summary


def _get_parameters() -> Dict[str, Any]:
    """Extract current parameters from session state."""
    params = {}
    
    # Try to get parameters from session state (may be stored separately)
    param_keys = [
        'contact_distance', 'chains', 'num_crossovers', 'min_fragment_diversity',
        'num_sequences', 'min_identity', 'max_identity', 'min_fragments', 'max_fragments'
    ]
    
    for key in param_keys:
        if key in st.session_state:
            params[key] = st.session_state[key]
    
    return params


def save_checkpoint(project_name=None, description="", parameters=None):
    """
    Save current session state to a checkpoint with full project structure.
    Saving is manual; the app may auto-load the current or last project when session data is missing.
    
    Args:
        project_name: Optional project name (default: timestamp-based)
        description: Optional description
        parameters: Optional dict of parameters to save (if None, extracted from session)
        
    Returns:
        str: Path to saved checkpoint directory
    """
    if project_name is None:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        project_name = f"checkpoint_{timestamp}"
    
    # Sanitize project name for directory name
    safe_name = "".join(c if c.isalnum() or c in ('_', '-') else '_' for c in project_name)
    project_dir = CHECKPOINT_DIR / safe_name
    project_dir.mkdir(exist_ok=True)
    
    files_dir = project_dir / "files"
    files_dir.mkdir(exist_ok=True)
    
    # Get parameters
    if parameters is None:
        parameters = _get_parameters()
    
    # Get workflow stage and data summary
    workflow_stage = _get_workflow_stage()
    data_summary = _get_data_summary()
    
    # Save all session state keys
    session_snapshot = _collect_session_snapshot()

    # Copy and update file paths
    file_mapping = {}

    # Diversity homolog MSA (bytes or text in session)
    div_bytes = st.session_state.get("diversity_analysis_fasta_bytes")
    div_text = st.session_state.get("diversity_analysis_fasta_text")
    if div_bytes or div_text:
        dest_div = files_dir / DIVERSITY_MSA_BUNDLE_NAME
        if div_bytes:
            dest_div.write_bytes(div_bytes)
        else:
            dest_div.write_text(div_text, encoding="utf-8")
        file_mapping["diversity_msa"] = f"files/{DIVERSITY_MSA_BUNDLE_NAME}"
        session_snapshot["diversity_msa_checkpoint_path"] = str(dest_div)

    # Applied crossover positions (redundant JSON for inspection)
    xo_positions = st.session_state.get("selected_crossover_positions")
    if xo_positions:
        xo_json_path = files_dir / "crossovers_applied.json"
        xo_json_path.write_text(
            json.dumps(
                {
                    "crossover_positions_1based": list(xo_positions),
                    "crossover_manual_input": st.session_state.get(
                        "crossover_manual_input", ""
                    ),
                },
                indent=2,
            )
            + "\n",
            encoding="utf-8",
        )
        file_mapping["crossovers_applied"] = "files/crossovers_applied.json"
    
    # Handle MSA file
    if SESSION_KEYS['msa_path'] in session_snapshot:
        msa_path = Path(session_snapshot[SESSION_KEYS['msa_path']])
        if msa_path.exists():
            dest_msa = files_dir / f"alignment{msa_path.suffix}"
            shutil.copy2(msa_path, dest_msa)
            file_mapping['msa'] = f"files/alignment{msa_path.suffix}"
            session_snapshot[SESSION_KEYS['msa_path']] = str(dest_msa)
    
    # Handle PDB file
    if SESSION_KEYS['pdb_path'] in session_snapshot:
        pdb_path = Path(session_snapshot[SESSION_KEYS['pdb_path']])
        if pdb_path.exists():
            dest_pdb = files_dir / "structure.pdb"
            shutil.copy2(pdb_path, dest_pdb)
            file_mapping['pdb'] = "files/structure.pdb"
            session_snapshot[SESSION_KEYS['pdb_path']] = str(dest_pdb)
    
    # Handle contact file (if saved separately)
    if SESSION_KEYS['schema_contacts'] in session_snapshot:
        contacts_data = session_snapshot[SESSION_KEYS['schema_contacts']]
        if isinstance(contacts_data, dict) and 'contacts' in contacts_data:
            # Save contacts to file
            contact_file = files_dir / "contacts.txt"
            with open(contact_file, 'w') as f:
                for contact in contacts_data['contacts']:
                    f.write(f"{contact[0]} {contact[1]} {contact[2]} {contact[3]}\n")
            file_mapping['contacts'] = "files/contacts.txt"
    
    # Handle crossover file (if exists)
    if SESSION_KEYS['crossover_file'] in session_snapshot:
        crossover_file = session_snapshot[SESSION_KEYS['crossover_file']]
        # If it's an UploadedFile, we can't save it directly, skip
        # Otherwise, if it's a path, copy it
        if isinstance(crossover_file, str):
            crossover_path = Path(crossover_file)
            if crossover_path.exists():
                dest_xo = files_dir / "crossovers.txt"
                shutil.copy2(crossover_path, dest_xo)
                file_mapping['crossovers'] = "files/crossovers.txt"
                session_snapshot[SESSION_KEYS['crossover_file']] = str(dest_xo)
    
    # Create metadata
    metadata = {
        'version': CHECKPOINT_VERSION,
        'project_name': project_name,
        'description': description,
        'timestamp': datetime.now().isoformat(),
        'app_version': '1.0.0',  # Could be read from __version__ if available
        'workflow_stage': workflow_stage,
        'parameters': parameters,
        'data_summary': data_summary,
        'files': file_mapping
    }
    
    # Save metadata as JSON
    metadata_path = project_dir / "metadata.json"
    with open(metadata_path, 'w') as f:
        json.dump(metadata, f, indent=2, default=str)
    
    # Save version
    version_path = project_dir / "version.txt"
    with open(version_path, 'w') as f:
        f.write(CHECKPOINT_VERSION)
    
    # Save checkpoint data as pickle
    checkpoint_data = {
        'version': CHECKPOINT_VERSION,
        'timestamp': datetime.now().isoformat(),
        'project_name': project_name,
        'description': description,
        'session_data': session_snapshot,
        'metadata': metadata
    }
    
    checkpoint_path = project_dir / "checkpoint.pkl"
    with open(checkpoint_path, 'wb') as f:
        pickle.dump(checkpoint_data, f)
    
    # Calculate checksum
    checksum = _calculate_checksum(checkpoint_path)
    metadata['checksum'] = checksum
    with open(metadata_path, 'w') as f:
        json.dump(metadata, f, indent=2, default=str)
    
    _write_last_project(project_dir)
    return str(project_dir)


def render_save_project_ui(key_prefix: str, inside_expander: bool = False):
    """
    Render the shared Save Project UI block (subheader, project name/description inputs, button).
    Call from any page; uses key_prefix to ensure unique Streamlit keys.
    If inside_expander is True, skips the separator and subheader (expander already has a title).
    """
    if not inside_expander:
        st.markdown("---")
        st.subheader("💾 Save Project")
    col1, col2 = st.columns(2)
    with col1:
        save_project_name = st.text_input(
            "Project Name",
            value=st.session_state.get('current_project', ''),
            key=f'{key_prefix}_project_name'
        )
    with col2:
        save_description = st.text_input(
            "Description (optional)",
            key=f'{key_prefix}_description'
        )
    if st.button("💾 Save Project", type="primary", key=key_prefix):
        try:
            checkpoint_path = save_checkpoint(
                project_name=save_project_name if save_project_name else None,
                description=save_description
            )
            st.success(f"✓ Project saved: {Path(checkpoint_path).name}")
            st.session_state['current_project'] = save_project_name if save_project_name else Path(checkpoint_path).name
            st.session_state['current_project_path'] = checkpoint_path
            _write_last_project(checkpoint_path)
        except Exception as e:
            st.error(f"Error saving project: {str(e)}")


def auto_save(stage: str = ""):
    """
    Automatically save the current session state after a key workflow step.
    
    Uses the user's project name (from 'project_name_for_autosave' in session state)
    combined with a timestamp. Falls back to 'autosave_<stage>' if no name set.
    
    Args:
        stage: Description of the workflow stage (e.g., 'blast_complete', 'contacts_calculated')
    
    Returns:
        str: Path to the autosave checkpoint, or None if autosave is disabled or fails
    """
    if not AUTO_SAVE_ENABLED:
        return None
    
    try:
        # Use user's project name if available, otherwise use stage-based name
        user_project_name = st.session_state.get('project_name_for_autosave', '')
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        if user_project_name:
            # User provided a project name - use it with timestamp
            project_name = f"{user_project_name}_{timestamp}"
            description = f"{user_project_name} - {stage}" if stage else user_project_name
        else:
            # No project name - use autosave with stage
            project_name = f"autosave_{stage}_{timestamp}" if stage else f"autosave_{timestamp}"
            description = f"Autosave after {stage}" if stage else "Autosave"
        
        # Save the checkpoint
        checkpoint_path = save_checkpoint(
            project_name=project_name,
            description=description,
            parameters=None
        )
        
        # Update session state to track this as current project
        st.session_state['current_project'] = project_name
        st.session_state['current_project_path'] = checkpoint_path
        _write_last_project(checkpoint_path)
        
        # Clean up old autosaves (keep only AUTOSAVE_MAX_CHECKPOINTS)
        _cleanup_old_autosaves()
        
        return checkpoint_path
    except Exception as e:
        # Silently fail on autosave - don't interrupt user workflow
        import logging
        logging.warning(f"Autosave failed: {e}")
        return None


def _cleanup_old_autosaves():
    """Remove old autosave checkpoints, keeping only the most recent ones.
    
    Cleans up checkpoints that:
    - Start with 'autosave_' (old-style autosaves)
    - Have a timestamp pattern in the name (new-style: projectname_YYYYMMDD_HHMMSS)
    """
    import re
    timestamp_pattern = re.compile(r'_\d{8}_\d{6}$')  # Matches _YYYYMMDD_HHMMSS at end
    
    try:
        autosave_dirs = []
        for project_dir in CHECKPOINT_DIR.iterdir():
            if not project_dir.is_dir():
                continue
            
            # Include if it's an autosave or has timestamp pattern (auto-generated)
            is_autosave = (
                project_dir.name.startswith('autosave_') or
                timestamp_pattern.search(project_dir.name)
            )
            
            if is_autosave:
                metadata_path = project_dir / "metadata.json"
                if metadata_path.exists():
                    with open(metadata_path) as f:
                        metadata = json.load(f)
                    timestamp = metadata.get('timestamp', '')
                    autosave_dirs.append((project_dir, timestamp))
        
        # Sort by timestamp (newest first)
        autosave_dirs.sort(key=lambda x: x[1], reverse=True)
        
        # Remove old autosaves beyond the limit
        last_path: Optional[Path] = None
        if LAST_PROJECT_FILE.exists():
            try:
                last_path = Path(LAST_PROJECT_FILE.read_text(encoding="utf-8").strip())
            except Exception:
                last_path = None

        for project_dir, _ in autosave_dirs[AUTOSAVE_MAX_CHECKPOINTS:]:
            if last_path and project_dir.resolve() == last_path.resolve():
                continue
            shutil.rmtree(project_dir, ignore_errors=True)
    except Exception:
        pass  # Don't fail if cleanup fails


def _calculate_checksum(file_path: Path) -> str:
    """Calculate SHA256 checksum of a file."""
    sha256_hash = hashlib.sha256()
    with open(file_path, "rb") as f:
        for byte_block in iter(lambda: f.read(4096), b""):
            sha256_hash.update(byte_block)
    return sha256_hash.hexdigest()


def load_checkpoint(checkpoint_path: Union[str, Path], validate=True):
    """
    Load a checkpoint and restore session state.
    Also used for startup autoload when the session is empty.
    
    Args:
        checkpoint_path: Path to checkpoint directory or pickle file
        validate: Whether to validate checkpoint integrity
        
    Returns:
        dict: Checkpoint metadata
    """
    checkpoint_path = Path(checkpoint_path)
    
    # Handle both directory and file paths
    if checkpoint_path.is_file():
        project_dir = checkpoint_path.parent
        checkpoint_file = checkpoint_path
    else:
        project_dir = checkpoint_path
        checkpoint_file = project_dir / "checkpoint.pkl"
    
    if not checkpoint_file.exists():
        raise FileNotFoundError(f"Checkpoint not found: {checkpoint_file}")
    
    # Validate if requested
    if validate:
        validate_checkpoint(checkpoint_file)
    
    # Load checkpoint data
    with open(checkpoint_file, 'rb') as f:
        checkpoint_data = pickle.load(f)
    
    # Check version compatibility
    checkpoint_version = checkpoint_data.get('version', '0.0')
    if checkpoint_version != CHECKPOINT_VERSION:
        # Try to load anyway, but warn
        pass  # Could add version migration logic here
    
    # Restore session state
    session_data = checkpoint_data.get('session_data', {})
    
    for key, value in session_data.items():
        # Handle special deserialization for raspp_raw_results
        if key == SESSION_KEYS['raspp_raw_results']:
            st.session_state[key] = deserialize_raspp_raw_results(value)
        elif key in CHECKPOINT_PATH_KEYS:
            if value and isinstance(value, str):
                path = Path(value)
                if not path.is_absolute():
                    path = project_dir / path
                elif not path.exists() and (project_dir / path.name).exists():
                    path = project_dir / path.name
                st.session_state[key] = str(path) if path.exists() else value
            else:
                st.session_state[key] = value
        else:
            st.session_state[key] = value

    _normalize_state_after_load(project_dir)

    for draft_key in (
        "xo_table_pick_draft",
        "_xo_freq_eligible_sig",
        "xo_freq_table_df",
        "xo_freq_data_editor",
    ):
        st.session_state.pop(draft_key, None)
    for key in list(st.session_state.keys()):
        if key.startswith("xo_pick_"):
            st.session_state.pop(key, None)
    
    # Store project info in session state
    st.session_state['current_project'] = checkpoint_data.get('project_name')
    st.session_state['current_project_path'] = str(project_dir)
    st.session_state['_autoload_attempted'] = True
    
    _write_last_project(project_dir)
    
    return checkpoint_data.get('metadata', checkpoint_data)


def validate_checkpoint(checkpoint_path: Union[str, Path]) -> bool:
    """
    Validate checkpoint integrity.
    
    Args:
        checkpoint_path: Path to checkpoint file or directory
        
    Returns:
        bool: True if valid
        
    Raises:
        ValueError: If checkpoint is invalid
    """
    checkpoint_path = Path(checkpoint_path)
    
    if checkpoint_path.is_file():
        project_dir = checkpoint_path.parent
        checkpoint_file = checkpoint_path
    else:
        project_dir = checkpoint_path
        checkpoint_file = project_dir / "checkpoint.pkl"
    
    if not checkpoint_file.exists():
        raise FileNotFoundError(f"Checkpoint file not found: {checkpoint_file}")
    
    # Check metadata exists
    metadata_path = project_dir / "metadata.json"
    if not metadata_path.exists():
        raise ValueError(f"Metadata file not found: {metadata_path}")
    
    # Load and check metadata
    try:
        with open(metadata_path, 'r') as f:
            metadata = json.load(f)
    except Exception as e:
        raise ValueError(f"Could not read metadata: {str(e)}")
    
    # Verify checksum if present
    if 'checksum' in metadata:
        calculated_checksum = _calculate_checksum(checkpoint_file)
        if calculated_checksum != metadata['checksum']:
            raise ValueError("Checkpoint checksum mismatch - file may be corrupted")
    
    # Check required files exist
    files = metadata.get('files', {})
    for file_type, file_path in files.items():
        file_full_path = project_dir / file_path
        if not file_full_path.exists():
            raise ValueError(f"Required file missing: {file_full_path}")
    
    return True


def delete_checkpoint(project_name: str) -> bool:
    """
    Delete a checkpoint project.
    
    Args:
        project_name: Project name or path
        
    Returns:
        bool: True if deleted successfully
    """
    project_path = Path(project_name)
    
    if not project_path.is_absolute():
        project_path = CHECKPOINT_DIR / project_name
    
    if not project_path.exists():
        return False
    
    # Safety check - must be in checkpoint directory
    try:
        project_path.resolve().relative_to(CHECKPOINT_DIR.resolve())
    except ValueError:
        raise ValueError("Cannot delete checkpoint outside checkpoint directory")
    
    # Remove directory
    shutil.rmtree(project_path)
    return True


def get_checkpoint_metadata(project_name: str) -> Optional[Dict[str, Any]]:
    """
    Get checkpoint metadata without loading full checkpoint.
    
    Args:
        project_name: Project name or path
        
    Returns:
        dict: Metadata or None if not found
    """
    project_path = Path(project_name)
    
    if not project_path.is_absolute():
        project_path = CHECKPOINT_DIR / project_name
    
    metadata_path = project_path / "metadata.json"
    
    if not metadata_path.exists():
        return None
    
    try:
        with open(metadata_path, 'r') as f:
            return json.load(f)
    except Exception:
        return None


def get_checkpoint_summary(project_name: str) -> Optional[Dict[str, Any]]:
    """
    Get quick summary of checkpoint without full load.
    
    Args:
        project_name: Project name or path
        
    Returns:
        dict: Summary with key info
    """
    metadata = get_checkpoint_metadata(project_name)
    
    if not metadata:
        return None
    
    return {
        'project_name': metadata.get('project_name'),
        'timestamp': metadata.get('timestamp'),
        'description': metadata.get('description'),
        'workflow_stage': metadata.get('workflow_stage'),
        'data_summary': metadata.get('data_summary', {}),
    }


def list_checkpoints(sort_by='timestamp', reverse=True) -> List[Tuple[str, Dict[str, Any]]]:
    """
    List all available checkpoint projects.
    
    Args:
        sort_by: Field to sort by ('timestamp', 'project_name')
        reverse: Whether to reverse sort order
        
    Returns:
        list: List of (project_path, metadata) tuples
    """
    checkpoints = []
    
    if not CHECKPOINT_DIR.exists():
        return []
    
    for project_dir in CHECKPOINT_DIR.iterdir():
        if not project_dir.is_dir():
            continue
        
        metadata = get_checkpoint_metadata(project_dir.name)
        if metadata:
            checkpoints.append((str(project_dir), metadata))
    
    # Sort
    if sort_by == 'timestamp':
        checkpoints.sort(key=lambda x: x[1].get('timestamp', ''), reverse=reverse)
    elif sort_by == 'project_name':
        checkpoints.sort(key=lambda x: x[1].get('project_name', ''), reverse=reverse)
    
    return checkpoints


def export_checkpoint(project_name: str, output_path: Optional[str] = None, format='zip') -> str:
    """
    Export checkpoint as zip file or JSON.
    
    Args:
        project_name: Project name or path
        output_path: Output file path (default: project_name.zip/json)
        format: 'zip' or 'json'
        
    Returns:
        str: Path to exported file
    """
    project_path = Path(project_name)
    
    if not project_path.is_absolute():
        project_path = CHECKPOINT_DIR / project_name
    
    if not project_path.exists():
        raise FileNotFoundError(f"Project not found: {project_path}")
    
    if format == 'zip':
        if output_path is None:
            output_path = str(CHECKPOINT_DIR / f"{project_path.name}.zip")
        
        with zipfile.ZipFile(output_path, 'w', zipfile.ZIP_DEFLATED) as zipf:
            for file_path in project_path.rglob('*'):
                if file_path.is_file():
                    arcname = file_path.relative_to(project_path)
                    zipf.write(file_path, arcname)
        
        return output_path
    
    elif format == 'json':
        if output_path is None:
            output_path = str(CHECKPOINT_DIR / f"{project_path.name}.json")
        
        metadata = get_checkpoint_metadata(project_path.name)
        
        # Load checkpoint data
        checkpoint_file = project_path / "checkpoint.pkl"
        with open(checkpoint_file, 'rb') as f:
            checkpoint_data = pickle.load(f)
        
        # Convert to JSON-safe format
        json_data = {
            'metadata': metadata,
            'session_data': {}
        }
        
        # Convert session data to JSON-safe
        for key, value in checkpoint_data.get('session_data', {}).items():
            json_data['session_data'][key] = make_json_safe(value)
        
        with open(output_path, 'w') as f:
            json.dump(json_data, f, indent=2, default=str)
        
        return output_path
    
    else:
        raise ValueError(f"Unknown format: {format}")


def import_checkpoint(file_path: str, project_name: Optional[str] = None) -> str:
    """
    Import checkpoint from zip file or JSON.
    
    Args:
        file_path: Path to zip or JSON file
        project_name: Optional project name (default: from file)
        
    Returns:
        str: Path to imported project directory
    """
    file_path = Path(file_path)
    
    if not file_path.exists():
        raise FileNotFoundError(f"File not found: {file_path}")
    
    if file_path.suffix == '.zip':
        # Extract zip
        if project_name is None:
            project_name = file_path.stem
        
        project_dir = CHECKPOINT_DIR / project_name
        project_dir.mkdir(exist_ok=True)
        
        with zipfile.ZipFile(file_path, 'r') as zipf:
            zipf.extractall(project_dir)
        
        return str(project_dir)
    
    elif file_path.suffix == '.json':
        # Import from JSON (requires manual reconstruction)
        if project_name is None:
            project_name = file_path.stem
        
        project_dir = CHECKPOINT_DIR / project_name
        project_dir.mkdir(exist_ok=True)
        
        with open(file_path, 'r') as f:
            json_data = json.load(f)
        
        # Save metadata
        metadata = json_data.get('metadata', {})
        metadata_path = project_dir / "metadata.json"
        with open(metadata_path, 'w') as f:
            json.dump(metadata, f, indent=2)
        
        # Restore session data (convert back from JSON-safe)
        session_data = {}
        for key, value in json_data.get('session_data', {}).items():
            session_data[key] = restore_from_json_safe(value)
        
        # Save as pickle
        checkpoint_data = {
            'version': metadata.get('version', CHECKPOINT_VERSION),
            'timestamp': metadata.get('timestamp'),
            'project_name': metadata.get('project_name', project_name),
            'description': metadata.get('description', ''),
            'session_data': session_data,
            'metadata': metadata
        }
        
        checkpoint_path = project_dir / "checkpoint.pkl"
        with open(checkpoint_path, 'wb') as f:
            pickle.dump(checkpoint_data, f)
        
        return str(project_dir)
    
    else:
        raise ValueError(f"Unsupported file format: {file_path.suffix}")
