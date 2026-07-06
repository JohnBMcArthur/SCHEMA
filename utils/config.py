"""
Configuration constants for SCHEMA-RASPP application.
Centralizes session state keys and other constants.
"""

# Session state keys
SESSION_KEYS = {
    'schema_contacts': 'schema_contacts',
    'schema_results': 'schema_results',
    'raspp_results': 'raspp_results',
    'raspp_raw_results': 'raspp_raw_results',
    'raspp_parents': 'raspp_parents',
    'multi_fragment_results': 'multi_fragment_results',
    'msa_path': 'msa_path',
    'pdb_path': 'pdb_path',
    'pdb_id': 'pdb_id',
    'sequences': 'sequences',
    'temp_dir': 'temp_dir',
    'crossover_file': 'crossover_file',
    'uploaded_files': 'uploaded_files',
}

# Default values
DEFAULTS = {
    'contact_distance': 5.0,
    'num_sequences': 10,
    'min_identity': 0.5,  # Minimum identity for BLAST hits (50%)
    'max_identity': 0.9,  # Maximum identity for sequence selection (90%)
    'num_crossovers': 3,
    'min_fragment_diversity': 10,
    'min_fragments': 5,
    'max_fragments': 20,
    'restrict_crossovers_at_graph_marks': False,
    'msa_segment_conservation_min': 80.0,
    'msa_segment_diversity_max': 3,
    'chains': ['A'],  # Chain IDs for PDB; use first when single chain
}

# BLAST parameters
BLAST_PARAMS = {
    'max_alignments': 1000,
    'max_scores': 1000,
    'evalue': '1e-5',
    'max_candidates': 20,  # Process up to 20 candidates for sequence selection
}

# AlphaFold parameters
ALPHAFOLD_PARAMS = {
    'max_version': 20,
    'min_version': 1,
    'timeout': 30,
}

# EBI API parameters (set EBI_EMAIL env var or Streamlit secret to override email)
import os as _os

_DEFAULT_EBI_EMAIL = "johnbmcarthur@gmail.com"


def get_ebi_email() -> str:
    email = (_os.environ.get("EBI_EMAIL") or "").strip()
    if email:
        return email
    try:
        import streamlit as st

        if hasattr(st, "secrets") and "EBI_EMAIL" in st.secrets:
            secret_email = str(st.secrets["EBI_EMAIL"]).strip()
            if secret_email:
                return secret_email
    except Exception:
        pass
    return _DEFAULT_EBI_EMAIL


def get_hosting_profile() -> str:
    profile = (_os.environ.get("SCHEMA_HOSTING") or "local").strip().lower()
    try:
        import streamlit as st

        if hasattr(st, "secrets") and "SCHEMA_HOSTING" in st.secrets:
            profile = str(st.secrets["SCHEMA_HOSTING"]).strip().lower()
    except Exception:
        pass
    return profile or "local"


def is_cloud_hosting() -> bool:
    return get_hosting_profile() == "cloud"


EBI_API_PARAMS = {
    "email": _DEFAULT_EBI_EMAIL,
    "max_wait_time": 600,  # 10 minutes for BLAST
    "poll_interval": 5,  # seconds
    "alignment_timeout": 300,  # 5 minutes for alignment
}

# Upload limit (align with .streamlit/config.toml maxUploadSize when set)
MAX_UPLOAD_BYTES = 50 * 1024 * 1024

# Conservative caps when SCHEMA_HOSTING=cloud (Streamlit Community Cloud)
CLOUD_MAX_ESM_SAMPLES = 2000
CLOUD_MAX_OPT_BATCH = 1000
CLOUD_MAX_OPT_ROUNDS = 10

# Set SCHEMA_DEBUG=1 to show BLAST TSV debug UI on SCHEMA Energy page
SCHEMA_DEBUG = _os.environ.get("SCHEMA_DEBUG", "").lower() in ("1", "true", "yes")

# Structure viewer: distance threshold (Å) for "nearby" residues in contact editor
NEARBY_RESIDUE_DISTANCE = 8.0

# Crossover Analysis: alternating cartoon colors per SCHEMA fragment (even / odd)
CROSSOVER_FRAGMENT_COLORS = ("#d32f2f", "#2e7d32")

# Golden Gate Assembly: AA-pair overhang compatibility table
from pathlib import Path as _Path

GGA_COMPATIBILITY_YAML = _Path(__file__).resolve().parent.parent / "data" / "gga_aa_pair_compatibility.yaml"
STUFFER_SEQUENCES_YAML = _Path(__file__).resolve().parent.parent / "data" / "stuffer_sequences.yaml"

# Checkpoint system parameters
CHECKPOINT_VERSION = "2.0"
CHECKPOINT_DIR = "checkpoints"
AUTO_SAVE_ENABLED = True  # Autosave after key workflow steps (event-based only)
AUTO_LOAD_ON_STARTUP = True  # Restore last project when session is empty
AUTOSAVE_MAX_CHECKPOINTS = 10  # Keep only N most recent autosaves
# Autosave is event-triggered only (e.g. after contacts, RASPP run), not interval-based
