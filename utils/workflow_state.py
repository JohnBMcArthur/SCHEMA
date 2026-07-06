"""
Cross-page workflow helpers: downstream invalidation, signatures, and banners.
"""

from __future__ import annotations

import hashlib
import json
from typing import Any, Dict, List, Optional

import streamlit as st

from utils.config import SESSION_KEYS
from utils.gga_overhangs import assembly_gga_options_from_mapping


DOWNSTREAM_FROM_CONTACTS: List[str] = [
    SESSION_KEYS["raspp_results"],
    SESSION_KEYS["raspp_raw_results"],
    SESSION_KEYS["raspp_parents"],
    SESSION_KEYS["multi_fragment_results"],
    "selected_crossover_positions",
    "crossover_manual_input",
    "_xo_applied_revision",
    "diversity_analysis_result",
    "diversity_table_rows",
    "diversity_analysis_upload_hash",
    "diversity_saved_selections",
    "library_opt_results",
    "oligopool_library_result",
    "oligopool_stuffer_result",
]

DOWNSTREAM_FROM_RASPP: List[str] = [
    "selected_crossover_positions",
    "crossover_manual_input",
    "_xo_applied_revision",
    "diversity_analysis_result",
    "diversity_saved_selections",
    "library_opt_results",
    "oligopool_library_result",
    "oligopool_stuffer_result",
]


def clear_session_keys(keys: List[str]) -> None:
    for key in keys:
        st.session_state.pop(key, None)


def clear_downstream_from_contacts() -> None:
    clear_session_keys(DOWNSTREAM_FROM_CONTACTS)


def clear_downstream_from_raspp() -> None:
    clear_session_keys(DOWNSTREAM_FROM_RASPP)


def crossovers_applied() -> bool:
    return bool(st.session_state.get("selected_crossover_positions"))


def pools_saved() -> bool:
    return bool(st.session_state.get("diversity_saved_selections"))


def assembly_gga_signature() -> str:
    return json.dumps(
        assembly_gga_options_from_mapping(st.session_state),
        sort_keys=True,
    )


def upstream_workflow_signature() -> str:
    payload = {
        "crossovers": list(st.session_state.get("selected_crossover_positions") or []),
        "gga": assembly_gga_options_from_mapping(st.session_state),
    }
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, default=str).encode("utf-8")
    ).hexdigest()


def invalidate_diversity_if_upstream_changed() -> bool:
    """
    Clear cached diversity analysis when crossovers or assembly GGA options change.
    Returns True if results were invalidated.
    """
    sig = upstream_workflow_signature()
    prev = st.session_state.get("_diversity_upstream_sig")
    invalidated = False
    if prev is not None and prev != sig and st.session_state.get("diversity_analysis_result"):
        for key in (
            "diversity_analysis_result",
            "diversity_table_rows",
            "diversity_network_preview",
            "diversity_ssn_cache",
        ):
            st.session_state.pop(key, None)
        invalidated = True
    st.session_state["_diversity_upstream_sig"] = sig
    return invalidated


def pools_signature(selections: Optional[Dict[str, Any]]) -> str:
    if not selections:
        return ""
    return hashlib.sha256(
        json.dumps(selections, sort_keys=True, default=str).encode("utf-8")
    ).hexdigest()


def render_cloud_persistence_notice() -> None:
    st.info(
        "On Streamlit Cloud, saved projects on the server are **temporary** "
        "(lost on redeploy). **Export** your project from the **Projects** page "
        "after each session and **re-import** when needed."
    )


def render_crossover_required_banner() -> None:
    if not crossovers_applied():
        st.warning(
            "Crossovers are not applied yet. On **Crossover Analysis**, select positions "
            "and click **Apply crossover selection**."
        )


def render_pools_required_banner() -> None:
    if not pools_saved():
        st.warning(
            "Homolog pools are not saved to session. On **Diversity Analysis** or "
            "**Library Optimization**, click **Save main list to session** (or "
            "**Save filtered list to session**) before using this page."
        )
