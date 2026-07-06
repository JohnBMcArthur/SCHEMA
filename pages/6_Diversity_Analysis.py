"""
Page 6: Diversity Analysis
Interactive fragment viewer, filters, and Sequence Similarity Network.
"""

from __future__ import annotations

import hashlib
import json
import sys
import zipfile
from io import BytesIO
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

import pandas as pd
import streamlit as st

sys.path.insert(0, str(Path(__file__).parent.parent))

from utils.assembly_analysis import (
    get_aligned_query_sequence,
    split_query_sequence_into_fragments,
    strip_msa_gaps,
)
from utils.config import CHECKPOINT_DIR, GGA_COMPATIBILITY_YAML
from utils.diversity_analysis import (
    analyze_diversity_msa,
    apply_fragment_filters_and_exclusions,
    build_all_fragment_table_rows,
    fragments_to_compact_json,
    fragments_to_fasta,
    get_saveable_homolog_rows,
)
from utils.diversity_display import metrics_for_row, render_msa_pair_html
from utils.diversity_ssn import (
    compute_distance_matrix,
    create_ssn_figure_from_distances,
    estimate_hamming_slider_max,
    max_from_distance_matrix,
)
from utils.gga_overhangs import (
    assign_golden_gate_overhangs,
    assembly_gga_options_from_mapping,
    build_fragment_terminal_filters,
    load_gga_compatibility,
)
from utils.session_manager import init_session_state
from utils.workflow_state import (
    invalidate_diversity_if_upstream_changed,
    render_crossover_required_banner,
)

st.set_page_config(
    page_title="Diversity Analysis",
    page_icon="🌿",
    layout="wide",
)

init_session_state()

_PROGRESS_VERSION = "2.0"
_FILTER_FIELDS = [
    ("min_length", "max_length", "Length"),
    ("min_pct_length", "max_pct_length", "% of query length"),
    ("min_mutations_non_gap", "max_mutations_non_gap", "Mutations (non-gap)"),
    ("min_mutations_with_gaps", "max_mutations_with_gaps", "Mutations (w/ gaps)"),
    ("min_pct_identity", "max_pct_identity", "% identity to query"),
]


def _init_viewer_state(fragment_indices: List[int]) -> None:
    defaults = {
        "diversity_applied_filters": {},
        "diversity_draft_filters": {},
        "diversity_graph_visible": {},
        "diversity_ssn_threshold": {},
        "diversity_network_preview": {},
        "diversity_excluded_rows": {},
        "diversity_saved_selections": {},
        "diversity_ssn_cache": {},
    }
    for key, value in defaults.items():
        if key not in st.session_state:
            st.session_state[key] = value

    for frag_idx in fragment_indices:
        st.session_state.diversity_applied_filters.setdefault(frag_idx, {})
        st.session_state.diversity_draft_filters.setdefault(frag_idx, {})
        st.session_state.diversity_graph_visible.setdefault(frag_idx, False)
        st.session_state.diversity_ssn_threshold.setdefault(frag_idx, 1)
        st.session_state.diversity_network_preview.setdefault(frag_idx, [])
        st.session_state.diversity_excluded_rows.setdefault(frag_idx, set())


def _reset_viewer_state(fragment_indices: List[int]) -> None:
    st.session_state.diversity_applied_filters = {idx: {} for idx in fragment_indices}
    st.session_state.diversity_draft_filters = {idx: {} for idx in fragment_indices}
    st.session_state.diversity_graph_visible = {idx: False for idx in fragment_indices}
    st.session_state.diversity_ssn_threshold = {idx: 1 for idx in fragment_indices}
    st.session_state.diversity_network_preview = {idx: [] for idx in fragment_indices}
    st.session_state.diversity_excluded_rows = {idx: set() for idx in fragment_indices}
    st.session_state.diversity_saved_selections = {}
    st.session_state.diversity_ssn_cache = {}


def _ssn_row_fingerprint(rows: List[Dict[str, Any]]) -> Tuple[str, ...]:
    return tuple(row["row_id"] for row in rows)


def _invalidate_graph(fragment_idx: int) -> None:
    st.session_state.diversity_graph_visible[fragment_idx] = False
    st.session_state.diversity_ssn_cache.pop(fragment_idx, None)


def _get_or_build_ssn_cache(
    fragment_idx: int,
    nodes: List[Dict[str, Any]],
    *,
    force_rebuild: bool = False,
) -> Optional[Dict[str, Any]]:
    fingerprint = _ssn_row_fingerprint(nodes)
    cache = st.session_state.diversity_ssn_cache.get(fragment_idx)
    if (
        not force_rebuild
        and cache
        and cache.get("fingerprint") == fingerprint
        and cache.get("distances") is not None
    ):
        return cache

    if len(nodes) <= 1:
        distances = [[0]]
        max_hamming = 0
    else:
        distances = compute_distance_matrix(nodes)
        max_hamming = max_from_distance_matrix(distances)

    cache = {
        "fingerprint": fingerprint,
        "distances": distances,
        "max_hamming": max_hamming,
    }
    st.session_state.diversity_ssn_cache[fragment_idx] = cache
    return cache


def _parse_optional_float(text: str) -> Optional[float]:
    text = (text or "").strip()
    if not text:
        return None
    return float(text)


def _parse_optional_int(text: str) -> Optional[int]:
    text = (text or "").strip()
    if not text:
        return None
    return int(text)


def _render_filter_panel(fragment_idx: int) -> None:
    draft = dict(st.session_state.diversity_draft_filters.get(fragment_idx) or {})
    applied = st.session_state.diversity_applied_filters.get(fragment_idx) or {}

    st.markdown("**Filters**")
    st.caption("Leave a bound empty to ignore it. Applies only to this fragment.")

    with st.form(key=f"diversity_filters_form_{fragment_idx}", clear_on_submit=False):
        new_draft: Dict[str, Any] = {}
        for min_key, max_key, label in _FILTER_FIELDS:
            st.markdown(f"**{label}**")
            c_min, c_max = st.columns(2)
            min_val = draft.get(min_key)
            max_val = draft.get(max_key)
            with c_min:
                min_text = st.text_input(
                    "Min",
                    value="" if min_val is None else str(min_val),
                    key=f"div_filt_{fragment_idx}_{min_key}",
                )
            with c_max:
                max_text = st.text_input(
                    "Max",
                    value="" if max_val is None else str(max_val),
                    key=f"div_filt_{fragment_idx}_{max_key}",
                )
            if "pct" in min_key or "identity" in min_key:
                new_draft[min_key] = _parse_optional_float(min_text)
                new_draft[max_key] = _parse_optional_float(max_text)
            else:
                new_draft[min_key] = _parse_optional_int(min_text)
                new_draft[max_key] = _parse_optional_int(max_text)

        btn_update, btn_clear = st.columns(2)
        update_clicked = btn_update.form_submit_button("Update Filters", use_container_width=True)
        clear_clicked = btn_clear.form_submit_button("Clear Filters", use_container_width=True)

    if clear_clicked:
        st.session_state.diversity_draft_filters[fragment_idx] = {}
        st.session_state.diversity_applied_filters[fragment_idx] = {}
        _invalidate_graph(fragment_idx)
        st.rerun()

    if update_clicked:
        cleaned = {k: v for k, v in new_draft.items() if v is not None}
        st.session_state.diversity_draft_filters[fragment_idx] = cleaned
        st.session_state.diversity_applied_filters[fragment_idx] = cleaned
        _invalidate_graph(fragment_idx)
        st.rerun()

    excluded = st.session_state.diversity_excluded_rows.get(fragment_idx) or set()
    if applied:
        st.caption(f"Active numeric filters: {len(applied)} bound(s) set.")
    if excluded:
        st.caption(f"Excluded via network: {len(excluded)} sequence(s).")


def _row_with_fresh_metrics(row: Dict[str, Any]) -> Dict[str, Any]:
    """Return a copy of the row with metrics recomputed from aligned fragments."""
    if row.get("is_query"):
        return dict(row)
    updated = dict(row)
    ng, wg, pct = metrics_for_row(row)
    updated["mutations_non_gap"] = ng
    updated["mutations_with_gaps"] = wg
    updated["pct_identity"] = pct
    return updated


def _rows_to_dataframe(rows: List[Dict[str, Any]]) -> pd.DataFrame:
    table_rows = []
    for row in rows:
        fresh = _row_with_fresh_metrics(row)
        table_rows.append(
            {
                "Sequence ID": row["sequence_id"],
                "Sequence": row["sequence"],
                "Length": row["length_display"],
                "Mutations (non-gap)": fresh["mutations_non_gap"],
                "Mutations (w/ gaps)": fresh["mutations_with_gaps"],
                "% identity": f"{fresh['pct_identity']:.1f}%",
            }
        )
    return pd.DataFrame(table_rows)


def _rows_by_id(rows: List[Dict[str, Any]]) -> Dict[str, Dict[str, Any]]:
    return {row["row_id"]: row for row in rows}


def _preview_rows(
    fragment_idx: int,
    main_rows: List[Dict[str, Any]],
) -> List[Dict[str, Any]]:
    preview_ids = st.session_state.diversity_network_preview.get(fragment_idx) or []
    lookup = _rows_by_id(main_rows)
    return [lookup[row_id] for row_id in preview_ids if row_id in lookup]


def _add_to_network_preview(fragment_idx: int, row_ids: List[str]) -> bool:
    preview: List[str] = list(st.session_state.diversity_network_preview.get(fragment_idx) or [])
    changed = False
    for row_id in row_ids:
        if row_id not in preview:
            preview.append(row_id)
            changed = True
    if changed:
        st.session_state.diversity_network_preview[fragment_idx] = preview
    return changed


def _remove_from_network_preview(fragment_idx: int, row_ids: List[str]) -> None:
    preview: List[str] = list(st.session_state.diversity_network_preview.get(fragment_idx) or [])
    remove = set(row_ids)
    st.session_state.diversity_network_preview[fragment_idx] = [
        row_id for row_id in preview if row_id not in remove
    ]


def _clear_network_preview(fragment_idx: int) -> None:
    st.session_state.diversity_network_preview[fragment_idx] = []


def _exclude_from_main_list(fragment_idx: int, row_ids: List[str]) -> None:
    excluded: Set[str] = set(st.session_state.diversity_excluded_rows.get(fragment_idx) or set())
    for row_id in row_ids:
        if row_id != "__query__":
            excluded.add(row_id)
    st.session_state.diversity_excluded_rows[fragment_idx] = excluded
    _remove_from_network_preview(fragment_idx, row_ids)
    _invalidate_graph(fragment_idx)


def _render_network_preview_panel(fragment_idx: int, main_rows: List[Dict[str, Any]]) -> None:
    preview = [_row_with_fresh_metrics(row) for row in _preview_rows(fragment_idx, main_rows)]
    query_row = next((row for row in main_rows if row.get("is_query")), None)

    st.markdown("---")
    st.subheader("Network preview")
    st.caption(
        "Each homolog is shown in its query alignment frame (Q/H rows). "
        "Yellow = substitution; magenta = gap mutation."
    )

    if not preview:
        st.info("Click or lasso nodes in the network to inspect aligned sequences here.")
        return

    if query_row is None:
        st.warning("Query row missing from main list.")
        return

    homolog_preview = [row for row in preview if not row.get("is_query")]
    if not homolog_preview:
        st.info("Select homolog nodes in the network to inspect alignments here.")
        return

    header_left, header_seq, header_meta = st.columns([0.18, 0.52, 0.30])
    with header_left:
        st.markdown("**ID**")
    with header_seq:
        st.markdown("**Aligned sequence (Q = query frame, H = homolog)**")
    with header_meta:
        st.markdown("**Metrics**")

    remove_targets: List[str] = []
    for row in homolog_preview:
        aligned_q = row.get("aligned_query_fragment") or query_row.get("sequence") or ""
        aligned_h = row.get("aligned_homolog_fragment") or row.get("sequence") or ""
        row_left, row_seq, row_meta = st.columns([0.18, 0.52, 0.30])
        with row_left:
            remove = st.checkbox(
                "Remove",
                key=f"preview_remove_{fragment_idx}_{row['row_id']}",
                label_visibility="collapsed",
            )
            st.caption(row["sequence_id"])
            if remove:
                remove_targets.append(row["row_id"])
        with row_seq:
            st.markdown(render_msa_pair_html(aligned_q, aligned_h), unsafe_allow_html=True)
        with row_meta:
            st.caption(
                f"Length {row['length_display']}\n\n"
                f"Non-gap {row['mutations_non_gap']} | "
                f"W/ gaps {row['mutations_with_gaps']} | "
                f"Identity {row['pct_identity']:.1f}%"
            )

    btn_col1, btn_col2, btn_col3 = st.columns(3)
    with btn_col1:
        if st.button("Remove checked from preview", key=f"preview_remove_btn_{fragment_idx}"):
            if remove_targets:
                _remove_from_network_preview(fragment_idx, remove_targets)
                st.rerun()
    with btn_col2:
        if st.button("Clear preview selection", key=f"preview_clear_btn_{fragment_idx}"):
            _clear_network_preview(fragment_idx)
            st.rerun()
    with btn_col3:
        if st.button(
            "Filter out preview sequences from main list",
            key=f"preview_exclude_btn_{fragment_idx}",
        ):
            homolog_ids = [row["row_id"] for row in preview if not row.get("is_query")]
            if homolog_ids:
                _exclude_from_main_list(fragment_idx, homolog_ids)
                st.rerun()


def _build_progress_payload(
    fragment_indices: List[int],
    table_rows_by_fragment: Dict[int, List[Dict[str, Any]]],
) -> Dict[str, Any]:
    fragment_lists: Dict[str, Any] = {}
    for frag_idx in fragment_indices:
        main_rows = apply_fragment_filters_and_exclusions(
            table_rows_by_fragment.get(frag_idx, []),
            st.session_state.diversity_applied_filters.get(frag_idx),
            st.session_state.diversity_excluded_rows.get(frag_idx),
        )
        query_row = next((row for row in main_rows if row.get("is_query")), None)
        fragment_lists[str(frag_idx)] = {
            "query": query_row,
            "homologs": get_saveable_homolog_rows(main_rows),
        }

    return {
        "version": _PROGRESS_VERSION,
        "fragment_lists": fragment_lists,
        "query_msa_seq_id": st.session_state.get("diversity_msa_query_seq_id"),
        "fragment_filters": {
            str(k): v for k, v in (st.session_state.diversity_applied_filters or {}).items()
        },
        "excluded_row_ids": {
            str(k): sorted(v) for k, v in (st.session_state.diversity_excluded_rows or {}).items()
        },
        "ssn_thresholds": {
            str(k): v for k, v in (st.session_state.diversity_ssn_threshold or {}).items()
        },
        "network_preview_row_ids": {
            str(k): list(v) for k, v in (st.session_state.diversity_network_preview or {}).items()
        },
    }


def _apply_progress_payload(payload: Dict[str, Any], fragment_indices: List[int]) -> None:
    filters = payload.get("fragment_filters") or {}
    thresholds = payload.get("ssn_thresholds") or {}
    excluded = payload.get("excluded_row_ids") or {}
    preview = payload.get("network_preview_row_ids") or payload.get("activated_row_ids") or {}
    lists = payload.get("fragment_lists") or payload.get("fragment_selections") or {}
    query_msa_seq_id = payload.get("query_msa_seq_id")
    if query_msa_seq_id:
        st.session_state.diversity_msa_query_seq_id = query_msa_seq_id

    for frag_idx in fragment_indices:
        key = str(frag_idx)
        st.session_state.diversity_applied_filters[frag_idx] = filters.get(key) or {}
        st.session_state.diversity_draft_filters[frag_idx] = dict(
            st.session_state.diversity_applied_filters[frag_idx]
        )
        st.session_state.diversity_ssn_threshold[frag_idx] = int(thresholds.get(key) or 1)
        st.session_state.diversity_excluded_rows[frag_idx] = set(excluded.get(key) or [])
        st.session_state.diversity_network_preview[frag_idx] = list(preview.get(key) or [])
        st.session_state.diversity_graph_visible[frag_idx] = False
        st.session_state.diversity_ssn_cache.pop(frag_idx, None)

    st.session_state.diversity_saved_selections = lists


def _analysis_needs_refresh(analysis_result: Dict[str, Any]) -> bool:
    """True if cached analysis predates MSA-query-based extraction."""
    return not analysis_result.get("query_msa_seq_id")


def _upload_fasta_bytes(uploaded) -> Optional[bytes]:
    """Return FASTA bytes from the uploader, using a session cache when the buffer is empty."""
    if uploaded is None:
        return st.session_state.get("diversity_analysis_fasta_bytes")
    raw = uploaded.getvalue()
    if raw:
        st.session_state["diversity_analysis_fasta_bytes"] = raw
        return raw
    return st.session_state.get("diversity_analysis_fasta_bytes")


def _upload_fasta_text(uploaded) -> Optional[str]:
    data = _upload_fasta_bytes(uploaded)
    if not data:
        return st.session_state.get("diversity_analysis_fasta_text")
    text = data.decode("utf-8", errors="replace")
    st.session_state["diversity_analysis_fasta_text"] = text
    return text


def _process_new_msa_upload(
    uploaded,
    *,
    fragment_sequences: List[str],
    query_ungapped: str,
    fragment_filters: List[Dict[str, Any]],
    fragment_indices: List[int],
) -> None:
    """Analyze a newly uploaded MSA once; reset viewer state only for a new file hash."""
    if uploaded is None:
        return
    fasta_bytes = _upload_fasta_bytes(uploaded)
    if not fasta_bytes:
        return
    upload_hash = hashlib.sha256(fasta_bytes).hexdigest()
    if st.session_state.get("diversity_analysis_upload_hash") == upload_hash:
        return
    fasta_text = fasta_bytes.decode("utf-8", errors="replace")
    try:
        analysis = analyze_diversity_msa(
            fasta_text,
            fragment_sequences,
            query_ungapped=query_ungapped,
            fragment_filters=fragment_filters,
        )
    except Exception as exc:
        st.error(f"Could not analyze uploaded MSA: {exc}")
        st.stop()
    st.session_state.diversity_analysis_result = analysis
    st.session_state.diversity_analysis_upload_hash = upload_hash
    st.session_state.diversity_analysis_fasta_text = fasta_text
    st.session_state.diversity_msa_query_seq_id = analysis.get("query_msa_seq_id")
    st.session_state.diversity_table_rows = build_all_fragment_table_rows(analysis)
    _reset_viewer_state(fragment_indices)
    from utils.session_manager import auto_save

    auto_save("diversity_msa_analyzed")


def _load_progress_callback() -> None:
    uploaded = st.session_state.get("diversity_upload_progress")
    fragment_indices = st.session_state.get("_diversity_fragment_indices") or []
    if uploaded is None:
        st.session_state["_diversity_progress_error"] = "No progress file uploaded."
        st.session_state.pop("_diversity_progress_success", None)
        return
    try:
        raw = uploaded.getvalue()
        if not raw:
            raise ValueError("Progress upload is empty.")
        payload = json.loads(raw.decode("utf-8"))
        _apply_progress_payload(payload, fragment_indices)
        saved = payload.get("fragment_lists") or payload.get("fragment_selections")
        st.session_state.diversity_saved_selections = saved or {}
        st.session_state.pop("_diversity_progress_error", None)
        st.session_state["_diversity_progress_success"] = "Progress restored."
        from utils.session_manager import auto_save

        auto_save("diversity_progress_loaded")
    except Exception as exc:
        st.session_state.pop("_diversity_progress_success", None)
        st.session_state["_diversity_progress_error"] = str(exc)


if invalidate_diversity_if_upstream_changed():
    st.warning(
        "Crossovers or Assembly Analysis options changed since the last diversity run. "
        "Re-analyze your MSA if results look stale."
    )
render_crossover_required_banner()

st.title("🌿 Diversity Analysis")

st.markdown("""
**Diversity Analysis** mines **natural sequence diversity** for each assembly fragment. You
upload a **large homolog multiple sequence alignment (MSA)**; the app extracts fragment
sequences per homolog, enforces **Golden Gate terminal compatibility** from **Assembly
Analysis**, and lets you filter and explore candidates with a **Sequence Similarity Network
(SSN)** before saving pools for oligo design.

**Two different MSAs — by design**

| | SCHEMA-RASPP MSA (steps 1–3) | Diversity MSA (this page) |
|---|------------------------------|---------------------------|
| **When** | SCHEMA Energy → RASPP Design | Uploaded here |
| **Size** | **Small** — query plus a handful of diverse parents (automated workflow: typically 5–20 BLAST hits; manual: your SCHEMA parent set) | **As large as practical** — hundreds to thousands of homologs when available |
| **Purpose** | Structure–alignment mapping, SCHEMA contacts, crossover recommendations | Variant pools for each fragment in your library |
| **Why small?** | SCHEMA energy and RASPP scans are expensive; a compact, diverse parent set is enough to find good crossover regions | You need breadth of natural variation at each block to build a useful combinatorial library |

The MSA from SCHEMA Energy is **not** reused here on purpose: crossover design benefits from a
focused parent set, while library design benefits from the widest homolog coverage you can
align reliably. Build or download a family-wide alignment (e.g. broader BLAST, UniRef cluster,
or curated family MSA), then upload it below.

**Prerequisites:**

1. **1. SCHEMA Energy** — session **query sequence** (the same protein you will include in
   the diversity MSA).
2. **3. Crossover Analysis** — click **Apply crossover selection**.
3. **4. Assembly Analysis** — review Golden Gate overhangs; terminal residue filters on this
   page are derived from those assignments.

**What this page does:**

1. **Upload homolog MSA (FASTA)** — all rows same alignment length; one row must match the
   session query (IDs containing `query` are preferred).
2. **Extract fragments** — for each homolog, slice MSA columns corresponding to each assembly
   fragment on the query row.
3. **GGA terminal filter** — drop homologs whose fragment ends are incompatible with assigned
   overhangs (see *Assembly fragment filters* expander).
4. **Per-fragment viewer** — main list table, numeric filters, SSN to spot redundancy or
   outliers, and aligned Q/H preview for selected nodes.
5. **Save main list to session** — commits filtered homologs per fragment for **Oligopool
   Design** (and optional **Library Optimization** / **Simulate with AI**).

**MSA requirements**

- **Aligned FASTA** (`.fasta`, `.fa`, `.fas`, `.txt`) with equal-length rows.
- **Include the query** as an aligned row matching the session query exactly (ungapped
  sequence). Rename with `query` in the ID if auto-detection fails.
- **Re-align if needed** — fragment boundaries come from the **query row columns** defined by
  your applied crossovers; homologs must share that alignment frame.
- Prefer **more sequences** over fewer; filter down on this page rather than limiting the
  upload.

**Steps:**

1. Confirm prerequisites (crossovers applied, assembly overhangs look correct).
2. **Upload** your large homolog MSA.
3. For each **fragment**, review the main list; adjust **filters** and use **Build network**
   to explore the SSN. Exclude redundant or unwanted homologs via the network preview.
4. Click **Save main list to session** when pools are ready (required before Oligopool Design).
5. *(Optional)* Download or upload **progress JSON** to resume filtering later.

**Next step:** **7. Oligopool Design** — designs oligos from the saved per-fragment homolog
lists. Run **8. Library Optimization** first if you want to prune the combinatorial space.
""")

crossovers = list(st.session_state.get("selected_crossover_positions") or [])
session_aligned, session_source = get_aligned_query_sequence()

if not crossovers:
    st.warning(
        "No crossover selection applied. Go to **Crossover Analysis**, apply crossovers, "
        "then review overhangs on **Assembly Analysis**."
    )
    st.stop()

if session_aligned is None:
    st.warning("No query sequence in session. Run **SCHEMA Energy** first.")
    st.stop()

if not GGA_COMPATIBILITY_YAML.is_file():
    st.error(f"GGA compatibility file not found: `{GGA_COMPATIBILITY_YAML}`")
    st.stop()

query_ungapped = strip_msa_gaps(session_aligned)
session_fragments = split_query_sequence_into_fragments(session_aligned, crossovers)

try:
    gga_data = load_gga_compatibility(str(GGA_COMPATIBILITY_YAML.resolve()))
    assembly = assign_golden_gate_overhangs(
        session_fragments,
        compatibility=gga_data,
        **assembly_gga_options_from_mapping(st.session_state),
    )
except Exception as exc:
    st.error(f"Could not assign Golden Gate overhangs: {exc}")
    st.stop()

fragment_sequences = [f["sequence"] for f in assembly["fragments"]]

fragment_filters = build_fragment_terminal_filters(assembly["fragments"])
fragment_indices = [int(f["fragment"]) for f in fragment_filters]
st.session_state["_diversity_fragment_indices"] = fragment_indices
_init_viewer_state(fragment_indices)

with st.expander("Assembly fragment filters (from session overhangs)", expanded=False):
    filter_rows = []
    for filt in fragment_filters:
        filter_rows.append(
            {
                "Fragment": filt["fragment"],
                "Query sequence": filt["query_sequence"],
                "5′ overhang": filt.get("overhang_5") or "—",
                "3′ overhang": filt.get("overhang_3") or "—",
                "Allowed N-term": ", ".join(filt["allowed_n_terminal"])
                if filt.get("allowed_n_terminal")
                else "—",
                "Allowed C-term": ", ".join(filt["allowed_c_terminal"])
                if filt.get("allowed_c_terminal")
                else "—",
            }
        )
    st.dataframe(pd.DataFrame(filter_rows), use_container_width=True, hide_index=True)
    st.caption(f"Session query: {session_source} ({len(query_ungapped)} residues, ungapped)")

st.markdown("---")
st.subheader("Upload homolog MSA (FASTA)")
st.caption(
    "Use your **broadest** reliable family alignment here — not the small SCHEMA-RASPP parent "
    "MSA from step 1. The query row must match the session query."
)

uploaded = st.file_uploader(
    "Upload aligned FASTA",
    type=["fasta", "fa", "fas", "txt"],
    key="diversity_msa_upload",
)

if uploaded is None and not st.session_state.get("diversity_analysis_result"):
    st.info(
        "Upload a large homolog FASTA MSA to begin. Include the session query as one aligned "
        "row; fragment extraction and GGA filtering run automatically after upload."
    )
    st.stop()

_process_new_msa_upload(
    uploaded,
    fragment_sequences=fragment_sequences,
    query_ungapped=query_ungapped,
    fragment_filters=fragment_filters,
    fragment_indices=fragment_indices,
)

analysis = st.session_state.get("diversity_analysis_result")
if not analysis:
    st.stop()

if _analysis_needs_refresh(analysis):
    fasta_text = _upload_fasta_text(uploaded)
    if not fasta_text:
        st.warning(
            "Cached analysis is outdated and the FASTA upload is no longer available. "
            "Re-upload the MSA file."
        )
        st.stop()
    try:
        analysis = analyze_diversity_msa(
            fasta_text,
            fragment_sequences,
            query_ungapped=query_ungapped,
            fragment_filters=fragment_filters,
        )
        st.session_state.diversity_analysis_result = analysis
        st.session_state.diversity_msa_query_seq_id = analysis.get("query_msa_seq_id")
        st.session_state.diversity_table_rows = build_all_fragment_table_rows(analysis)
    except Exception as exc:
        st.error(f"Could not refresh analysis: {exc}")
        st.stop()

table_rows_by_fragment = build_all_fragment_table_rows(analysis)
st.session_state.diversity_table_rows = table_rows_by_fragment

query_msa_seq_id = analysis.get("query_msa_seq_id") or st.session_state.get(
    "diversity_msa_query_seq_id"
)

stats = analysis.get("stats") or {}
st.caption(
    f"MSA query row: `{query_msa_seq_id or 'unknown'}` | "
    f"Processed {stats.get('sequences_processed', 0)} homolog(s) "
    f"({stats.get('sequences_failed', 0)} failed) of "
    f"{stats.get('homologs_processed', stats.get('sequences_total', 0) - 1)} homolog rows."
)

st.markdown("---")
st.subheader("Fragment viewer")

fragment_labels = [f"Fragment {idx}" for idx in fragment_indices]
selected_label = st.selectbox(
    "Fragment",
    fragment_labels,
    key="diversity_fragment_dropdown",
)
selected_fragment = fragment_indices[fragment_labels.index(selected_label)]

if st.session_state.get("_diversity_prev_fragment") != selected_fragment:
    st.session_state._diversity_prev_fragment = selected_fragment
    _invalidate_graph(selected_fragment)

all_rows = table_rows_by_fragment.get(selected_fragment, [])
applied_filters = st.session_state.diversity_applied_filters.get(selected_fragment) or {}
excluded_rows = st.session_state.diversity_excluded_rows.get(selected_fragment) or set()
main_rows = apply_fragment_filters_and_exclusions(all_rows, applied_filters, excluded_rows)
homolog_count = sum(1 for row in main_rows if not row.get("is_query"))

table_col, filter_col = st.columns([3, 1])
with table_col:
    st.markdown(
        f"**{selected_label} — main list** — {homolog_count} homolog(s) "
        f"(+ query). This list is saved for **Oligopool Design**."
    )
    st.dataframe(_rows_to_dataframe(main_rows), use_container_width=True, hide_index=True)

with filter_col:
    _render_filter_panel(selected_fragment)

if not st.session_state.diversity_graph_visible.get(selected_fragment):
    if applied_filters or excluded_rows:
        st.info("Filters changed — click **Build network** to refresh the network view.")

st.markdown("---")
st.subheader("Sequence Similarity Network")

ssn_nodes = main_rows
ssn_cache = st.session_state.diversity_ssn_cache.get(selected_fragment)
if ssn_cache and ssn_cache.get("fingerprint") == _ssn_row_fingerprint(ssn_nodes):
    slider_max = max(1, int(ssn_cache.get("max_hamming") or 1))
else:
    slider_max = estimate_hamming_slider_max(ssn_nodes)

default_threshold = st.session_state.diversity_ssn_threshold.get(selected_fragment, 1)
threshold = st.slider(
    "Hamming distance threshold (edge if distance ≤ threshold)",
    min_value=0,
    max_value=slider_max,
    value=min(default_threshold, slider_max),
    key=f"diversity_ssn_slider_{selected_fragment}",
)
st.session_state.diversity_ssn_threshold[selected_fragment] = threshold

build_col, clear_sel_col, _ = st.columns([1, 1, 2])
with build_col:
    build_network = st.button("Build network", key=f"diversity_build_ssn_{selected_fragment}")
with clear_sel_col:
    if st.button("Deselect all", key=f"diversity_clear_preview_{selected_fragment}"):
        _clear_network_preview(selected_fragment)
        st.rerun()

if build_network:
    st.session_state.diversity_graph_visible[selected_fragment] = True
    with st.spinner("Computing pairwise distances…"):
        _get_or_build_ssn_cache(selected_fragment, ssn_nodes, force_rebuild=True)

if st.session_state.diversity_graph_visible.get(selected_fragment):
    cache = _get_or_build_ssn_cache(selected_fragment, ssn_nodes)
    if cache is None:
        st.warning("Could not build network cache.")
    else:
        preview_ids = set(st.session_state.diversity_network_preview.get(selected_fragment) or [])
        fig, _graph = create_ssn_figure_from_distances(
            ssn_nodes,
            cache["distances"],
            threshold,
            preview_ids,
            query_row_id="__query__",
        )
        fig.update_layout(dragmode="lasso")
        selection = st.plotly_chart(
            fig,
            use_container_width=True,
            on_select="rerun",
            key=f"diversity_ssn_chart_{selected_fragment}",
        )
        if selection and "selection" in selection:
            points = selection["selection"].get("points") or []
            selected_indices = sorted(
                {
                    int(point.get("point_index"))
                    for point in points
                    if point.get("curve_number") == 1 and point.get("point_index") is not None
                }
            )
            if selected_indices:
                row_ids = [
                    ssn_nodes[idx]["row_id"]
                    for idx in selected_indices
                    if 0 <= idx < len(ssn_nodes)
                ]
                if _add_to_network_preview(selected_fragment, row_ids):
                    st.rerun()
else:
    st.caption("Network hidden until you click **Build network**.")

_render_network_preview_panel(selected_fragment, main_rows)

st.markdown("---")
st.subheader("Save progress")

progress_payload = _build_progress_payload(fragment_indices, table_rows_by_fragment)
progress_json = json.dumps(progress_payload, indent=2) + "\n"

save_col, reset_col, download_col, upload_col = st.columns(4)
with save_col:
    if st.button("Save main list to session", key="diversity_save_session"):
        st.session_state.diversity_saved_selections = progress_payload["fragment_lists"]
        from utils.session_manager import auto_save

        auto_save("diversity_selections_saved")
        st.success(
            "Main filtered list saved to session. "
            "**Library Optimization**, **Simulate with AI**, and **Oligopool Design** "
            "use this list (Library Optimization can overwrite it after pruning)."
        )
with reset_col:
    if st.button("Reset all filters", key="diversity_reset_all_filters"):
        _reset_viewer_state(fragment_indices)
        st.success(
            "All filters, network exclusions, and the saved fragment list cleared. "
            "The table shows the full analyzed homolog set again."
        )
        st.rerun()
with download_col:
    st.download_button(
        "Download progress JSON",
        data=progress_json,
        file_name="diversity_progress.json",
        mime="application/json",
        key="diversity_download_progress",
    )
with upload_col:
    if st.session_state.get("_diversity_progress_success"):
        st.success(st.session_state.pop("_diversity_progress_success"))
    if st.session_state.get("_diversity_progress_error"):
        st.error(
            f"Could not load progress file: "
            f"{st.session_state.pop('_diversity_progress_error')}"
        )
    st.file_uploader(
        "Upload progress JSON",
        type=["json"],
        key="diversity_upload_progress",
    )
    st.button(
        "Load progress",
        key="diversity_load_progress_btn",
        on_click=_load_progress_callback,
    )

with st.expander("Export filtered fragments (legacy)", expanded=False):
    compact = fragments_to_compact_json(analysis)
    col_json, col_zip, col_disk = st.columns(3)
    with col_json:
        st.download_button(
            "Download compact JSON",
            data=json.dumps(compact, indent=2) + "\n",
            file_name="diversity_filtered.json",
            mime="application/json",
            key="download_diversity_compact_json",
        )
    with col_zip:
        zip_buffer = BytesIO()
        with zipfile.ZipFile(zip_buffer, "w", zipfile.ZIP_DEFLATED) as zf:
            zf.writestr("diversity_filtered.json", json.dumps(compact, indent=2) + "\n")
            for frag in analysis["fragments"]:
                fasta_body = fragments_to_fasta(frag)
                if fasta_body.strip():
                    zf.writestr(f"fragment_{frag['fragment']}.fasta", fasta_body)
        st.download_button(
            "Download ZIP (JSON + FASTA per fragment)",
            data=zip_buffer.getvalue(),
            file_name="diversity_filtered.zip",
            mime="application/zip",
            key="download_diversity_zip",
        )
    with col_disk:
        _repo = Path(__file__).resolve().parent.parent
        _save_dir = _repo / CHECKPOINT_DIR / "diversity_analysis"
        if st.button("Save to checkpoints/diversity_analysis", key="save_diversity_disk"):
            _save_dir.mkdir(parents=True, exist_ok=True)
            json_path = _save_dir / "diversity_filtered.json"
            json_path.write_text(json.dumps(compact, indent=2) + "\n", encoding="utf-8")
            saved_fastas = []
            for frag in analysis["fragments"]:
                fasta_body = fragments_to_fasta(frag)
                if not fasta_body.strip():
                    continue
                fasta_path = _save_dir / f"fragment_{frag['fragment']}.fasta"
                fasta_path.write_text(fasta_body, encoding="utf-8")
                saved_fastas.append(fasta_path.name)
            st.success(
                f"Saved `{json_path.name}`"
                + (
                    f" and {len(saved_fastas)} FASTA file(s) under `{_save_dir}`."
                    if saved_fastas
                    else "."
                )
            )
