"""
Query fragment assembly helpers for the Assembly Analysis page.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from utils.config import SESSION_KEYS


def strip_msa_gaps(sequence: str) -> str:
    """Remove gap characters from an MSA sequence string."""
    return (sequence or "").replace("-", "").replace(".", "")


def split_query_sequence_into_fragments(
    aligned_sequence: str,
    crossovers_1based: Optional[List[int]] = None,
) -> List[Dict[str, Any]]:
    """
    Split the aligned query into fragments at 1-based crossover columns.

    Each crossover position belongs to the **left** fragment (not the right).
    Gap characters are removed from each fragment sequence.
    """
    seq = aligned_sequence or ""
    length = len(seq)
    crossovers = sorted(
        {int(c) for c in (crossovers_1based or []) if 1 <= int(c) <= length}
    )

    cuts = [0] + crossovers + [length]
    fragments: List[Dict[str, Any]] = []
    for i in range(len(cuts) - 1):
        start_0 = cuts[i]
        end_0 = cuts[i + 1]
        aligned_slice = seq[start_0:end_0]
        ungapped = strip_msa_gaps(aligned_slice)
        fragments.append(
            {
                "index": i + 1,
                "aligned_start_1based": start_0 + 1,
                "aligned_end_1based": end_0,
                "aligned_sequence": aligned_slice,
                "sequence": ungapped,
                "length": len(ungapped),
            }
        )
    return fragments


def _sequence_from_parent_entry(entry: Any) -> Optional[str]:
    if isinstance(entry, tuple) and len(entry) >= 2:
        return str(entry[1])
    if isinstance(entry, str):
        return entry
    return None


def _find_query_in_parent_list(parent_list: List[Any]) -> Tuple[Optional[str], str]:
    for entry in parent_list:
        if isinstance(entry, tuple) and len(entry) >= 2:
            name = str(entry[0]).lower()
            if name in ("query", "query_sequence"):
                return str(entry[1]), f"MSA parent '{entry[0]}'"
    seq = _sequence_from_parent_entry(parent_list[0]) if parent_list else None
    if seq is not None:
        return seq, "first parent in alignment"
    return None, "unknown"


def get_aligned_query_sequence() -> Tuple[Optional[str], str]:
    """
    Resolve the aligned query sequence from session state.

    Returns:
        (aligned_sequence or None, source description)
    """
    import streamlit as st

    contacts_data = st.session_state.get(SESSION_KEYS["schema_contacts"])
    if isinstance(contacts_data, dict):
        parents_obj = contacts_data.get("parents_object")
        if parents_obj is not None:
            if hasattr(parents_obj, "p0_aligned"):
                aligned = getattr(parents_obj, "p0_aligned")
                if aligned:
                    return str(aligned), "trimmed query alignment (p0_aligned)"
            if hasattr(parents_obj, "alignment"):
                try:
                    aligned = "".join(
                        aminos[0] for aminos in parents_obj.alignment
                    )
                    if aligned:
                        return aligned, "trimmed query alignment (parents_object)"
                except (TypeError, IndexError):
                    pass

        parents = contacts_data.get("parents")
        if parents:
            seq, source = _find_query_in_parent_list(list(parents))
            if seq:
                return seq, source

    raspp_parents = st.session_state.get(SESSION_KEYS["raspp_parents"])
    if raspp_parents:
        if isinstance(raspp_parents, list):
            if raspp_parents and isinstance(raspp_parents[0], str):
                return str(raspp_parents[0]), "RASPP parents[0]"
            seq, source = _find_query_in_parent_list(raspp_parents)
            if seq:
                return seq, source

    msa_path = st.session_state.get(SESSION_KEYS["msa_path"])
    if msa_path and Path(msa_path).exists():
        try:
            from schema_raspp import schema

            with open(msa_path, "r", encoding="utf-8") as handle:
                alignment_data = schema.readMultipleSequenceAlignmentFile(handle)
            if alignment_data:
                for name, seq in alignment_data:
                    if str(name).lower() in ("query", "query_sequence"):
                        return seq, f"MSA file row '{name}'"
                return alignment_data[0][1], f"MSA file row '{alignment_data[0][0]}'"
        except Exception:
            pass

    ungapped = st.session_state.get("query_sequence")
    if ungapped:
        return str(ungapped), "ungapped query_sequence (no alignment gaps)"

    return None, "not available"
