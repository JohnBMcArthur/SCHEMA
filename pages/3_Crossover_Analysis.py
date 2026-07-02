"""
Page 3: Crossover analysis
Crossover distribution, MSA markers, and structure colored by recombination fragments.
"""

import json

import streamlit as st
from pathlib import Path

import sys

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent.parent))

from utils.session_manager import init_session_state
from utils.config import SESSION_KEYS, CHECKPOINT_DIR

from utils.visualization import (
    compute_crossover_position_counts,
    plot_crossover_distribution,
    format_msa_with_crossovers_html,
)


def _parse_1based_crossover_text(text):
    if not text or not str(text).strip():
        return []
    out = []
    for part in str(text).split(","):
        part = part.strip()
        if part.isdigit():
            out.append(int(part))
    return sorted(set(out))


def _safe_crossover_save_basename(name: str) -> str:
    base = "".join(
        c if c.isalnum() or c in ("_", "-") else "_"
        for c in (name or "").strip()
    )
    return base or "selected_crossovers"


def _parse_crossover_upload(uploaded_file) -> list:
    """Parse exported TXT or JSON crossover selection file."""
    if hasattr(uploaded_file, "getvalue"):
        raw = uploaded_file.getvalue()
    else:
        raw = uploaded_file.read()
    name = (getattr(uploaded_file, "name", "") or "").lower()
    if name.endswith(".json"):
        data = json.loads(raw.decode("utf-8"))
        if isinstance(data, dict) and "crossover_positions_1based" in data:
            return sorted({int(p) for p in data["crossover_positions_1based"]})
        if isinstance(data, list):
            return sorted({int(p) for p in data})
        raise ValueError("JSON must contain crossover_positions_1based or a list of positions.")
    text = raw.decode("utf-8")
    positions = _parse_1based_crossover_text(text)
    if not positions:
        raise ValueError("No crossover positions found in file.")
    return positions


def _position_checkbox_key(position: int) -> str:
    return f"xo_pick_{int(position)}"


def _sync_position_checkboxes(eligible_positions, picked_positions):
    """Set per-position checkbox session keys (only on commit / load / eligible change)."""
    picked = set(picked_positions or [])
    for pos in eligible_positions:
        st.session_state[_position_checkbox_key(pos)] = pos in picked


def _read_position_checkboxes(eligible_positions):
    return sorted(
        pos
        for pos in eligible_positions
        if st.session_state.get(_position_checkbox_key(pos), False)
    )


def _set_applied_crossover_positions(positions, parent_length=None, eligible_positions=None):
    """Update applied crossover list (MSA, structure, save UI read this). Returns cleaned list."""
    cleaned = sorted({int(p) for p in (positions or []) if int(p) > 0})
    if parent_length is not None:
        cleaned = [p for p in cleaned if p <= int(parent_length)]
    st.session_state["selected_crossover_positions"] = cleaned
    st.session_state["_xo_applied_revision"] = (
        int(st.session_state.get("_xo_applied_revision", 0)) + 1
    )
    st.session_state.pop("_xo_chart_pending", None)

    eligible_set = set(eligible_positions or [])
    if eligible_set:
        extras = sorted(p for p in cleaned if p not in eligible_set)
        st.session_state["crossover_manual_input"] = ", ".join(str(p) for p in extras)
    else:
        st.session_state["crossover_manual_input"] = ", ".join(str(p) for p in cleaned)
    from utils.session_manager import auto_save

    auto_save("crossover_selection")
    return cleaned


def _apply_crossover_selection_callback():
    eligible_positions = st.session_state.get("_xo_eligible_positions") or []
    parent_length = st.session_state.get("_xo_parent_length")
    table_picks = _read_position_checkboxes(eligible_positions)
    manual_1b = _parse_1based_crossover_text(
        st.session_state.get("crossover_manual_input", "")
    )
    chart_1b = list(st.session_state.get("_xo_chart_pending") or [])
    merged = sorted(set(table_picks) | set(manual_1b) | set(chart_1b))
    cleaned = _set_applied_crossover_positions(
        merged, parent_length, eligible_positions
    )
    eligible_set = set(eligible_positions)
    _sync_position_checkboxes(
        eligible_positions,
        [p for p in cleaned if p in eligible_set],
    )


def _clear_crossover_selection_callback():
    eligible_positions = st.session_state.get("_xo_eligible_positions") or []
    parent_length = st.session_state.get("_xo_parent_length")
    _set_applied_crossover_positions([], parent_length, eligible_positions)
    _sync_position_checkboxes(eligible_positions, [])


def _load_crossover_upload_callback():
    uploaded = st.session_state.get("xo_load_crossovers_upload")
    if uploaded is None:
        st.session_state["_xo_load_error"] = "No file uploaded."
        st.session_state.pop("_xo_load_success", None)
        return
    try:
        loaded = _parse_crossover_upload(uploaded)
        parent_length = st.session_state.get("_xo_parent_length")
        eligible_positions = st.session_state.get("_xo_eligible_positions") or []
        _set_applied_crossover_positions(loaded, parent_length, eligible_positions)
        st.session_state["_xo_sync_checkboxes"] = True
        st.session_state.pop("_xo_load_error", None)
        st.session_state["_xo_load_success"] = (
            f"Loaded {len(loaded)} crossover position(s)."
        )
    except Exception as exc:
        st.session_state.pop("_xo_load_success", None)
        st.session_state["_xo_load_error"] = str(exc)


def _load_crossover_saved_callback():
    pick = st.session_state.get("xo_load_crossovers_saved_pick")
    if not pick:
        return
    _repo = Path(__file__).resolve().parent.parent
    _save_dir = _repo / CHECKPOINT_DIR / "crossover_selections"
    try:
        _path = _save_dir / f"{pick}.json"
        with open(_path, "rb") as _f:
            loaded = _parse_crossover_upload(_f)
        parent_length = st.session_state.get("_xo_parent_length")
        eligible_positions = st.session_state.get("_xo_eligible_positions") or []
        _set_applied_crossover_positions(loaded, parent_length, eligible_positions)
        st.session_state["_xo_sync_checkboxes"] = True
        st.session_state.pop("_xo_load_error", None)
        st.session_state["_xo_load_success"] = (
            f"Loaded {len(loaded)} crossover position(s) from `{_path.name}`."
        )
    except Exception as exc:
        st.session_state.pop("_xo_load_success", None)
        st.session_state["_xo_load_error"] = str(exc)


def _sync_draft_checkboxes_if_needed(eligible, applied_positions):
    eligible_sig = tuple(eligible)
    eligible_positions = [p1 for p1, _ in eligible]
    if st.session_state.pop("_xo_sync_checkboxes", False):
        _sync_position_checkboxes(
            eligible_positions,
            [p for p in applied_positions if p in eligible_positions],
        )
    elif eligible and st.session_state.get("_xo_freq_eligible_sig") != eligible_sig:
        st.session_state["_xo_freq_eligible_sig"] = eligible_sig
        _sync_position_checkboxes(
            eligible_positions,
            [p for p in applied_positions if p in eligible_positions],
        )
    return eligible_positions


def _render_crossover_draft_picks(
    eligible,
    table_threshold,
    max_freq,
    applied_positions,
):
    """Draft checkboxes only (fragment: toggles do not rerun MSA / structure)."""
    eligible_positions = _sync_draft_checkboxes_if_needed(eligible, applied_positions)
    freq_by_pos = {p1: fr for p1, fr in eligible}

    if eligible:
        st.caption(
            f"Positions with frequency **>** {table_threshold:.3g} "
            f"(5% of max frequency {max_freq}). "
            "Check positions below, then click **Apply crossover selection**."
        )
        pick_col, freq_col = st.columns([1, 3])
        with freq_col:
            st.dataframe(
                pd.DataFrame(
                    {
                        "Position (1-based)": eligible_positions,
                        "Frequency": [freq_by_pos[p] for p in eligible_positions],
                    }
                ),
                use_container_width=True,
                hide_index=True,
            )
        with pick_col:
            st.markdown("**Select**")
            for pos in eligible_positions:
                st.checkbox(
                    f"{pos} ({freq_by_pos[pos]})",
                    key=_position_checkbox_key(pos),
                    help="Draft — MSA and structure update only after Apply.",
                )
    else:
        st.caption("No positions above 5% of the maximum frequency.")

    st.text_input(
        "Additional crossover positions (comma-separated, 1-based)",
        key="crossover_manual_input",
        help="Optional positions not shown in the frequency table above.",
    )


def _render_crossover_selection_controls(
    parent_length,
    eligible,
):
    """Apply / clear / save / load — outside fragment so MSA and structure always refresh."""
    eligible_positions = [p1 for p1, _ in eligible]
    st.session_state["_xo_eligible_positions"] = eligible_positions
    st.session_state["_xo_parent_length"] = parent_length

    if st.session_state.get("_xo_chart_pending"):
        st.caption(
            "Chart picks pending apply: "
            + ", ".join(str(p) for p in st.session_state["_xo_chart_pending"])
        )

    apply_col, clear_col = st.columns([1, 1])
    with apply_col:
        st.button(
            "Apply crossover selection",
            type="primary",
            key="xo_apply_crossover_selection",
            on_click=_apply_crossover_selection_callback,
        )
    with clear_col:
        st.button(
            "Clear selection",
            key="xo_clear_crossover_selection",
            on_click=_clear_crossover_selection_callback,
        )

    committed = list(st.session_state.get("selected_crossover_positions") or [])
    if committed:
        st.caption("**Applied** positions: " + ", ".join(map(str, committed)))
    else:
        st.caption(
            "No crossover selection applied yet. Check positions or use the text field, "
            "then click **Apply crossover selection**."
        )

    _repo = Path(__file__).resolve().parent.parent
    _save_dir = _repo / CHECKPOINT_DIR / "crossover_selections"
    with st.expander("Save / load selected crossovers", expanded=False):
        if committed:
            _xo_txt = ", ".join(str(p) for p in committed)
            _xo_json = json.dumps(
                {
                    "crossover_positions_1based": committed,
                    "format": "SCHEMA-RASPP / MSA indices (1-based)",
                },
                indent=2,
            )
            d1, d2 = st.columns(2)
            with d1:
                st.download_button(
                    label="Download TXT",
                    data=_xo_txt,
                    file_name="selected_crossovers_1based.txt",
                    mime="text/plain",
                    help="Comma-separated 1-based positions (same as the text field).",
                    key="download_selected_crossovers_txt",
                )
            with d2:
                st.download_button(
                    label="Download JSON",
                    data=_xo_json,
                    file_name="selected_crossovers_1based.json",
                    mime="application/json",
                    key="download_selected_crossovers_json",
                )
            st.caption(
                "Downloads reflect the **applied** selection. "
                "Apply changes above before exporting."
            )
            _basename = st.text_input(
                "File name (no extension)",
                value="selected_crossovers",
                key="xo_save_crossovers_basename",
            )
            if st.button("Save to checkpoints folder", key="xo_save_crossovers_disk"):
                safe = _safe_crossover_save_basename(
                    st.session_state.get("xo_save_crossovers_basename", "")
                )
                _save_dir.mkdir(parents=True, exist_ok=True)
                _txt_path = _save_dir / f"{safe}.txt"
                _json_path = _save_dir / f"{safe}.json"
                _txt_path.write_text(_xo_txt + "\n", encoding="utf-8")
                _json_path.write_text(_xo_json + "\n", encoding="utf-8")
                st.success(
                    f"Saved `{_txt_path.name}` and `{_json_path.name}` under `{_save_dir}`."
                )
        else:
            st.caption("Apply a selection above before downloading or saving.")

        st.markdown("---")
        st.markdown("**Load crossover positions**")
        if st.session_state.get("_xo_load_success"):
            st.success(st.session_state.pop("_xo_load_success"))
        if st.session_state.get("_xo_load_error"):
            st.error(
                f"Could not load crossover file: "
                f"{st.session_state.pop('_xo_load_error')}"
            )
        st.file_uploader(
            "Upload crossover TXT or JSON",
            type=["txt", "json"],
            key="xo_load_crossovers_upload",
        )
        st.button(
            "Load from uploaded file",
            key="xo_load_crossovers_upload_btn",
            on_click=_load_crossover_upload_callback,
        )

        _saved_json = sorted(_save_dir.glob("*.json")) if _save_dir.exists() else []
        if _saved_json:
            st.selectbox(
                "Or load from checkpoints/crossover_selections",
                options=[p.stem for p in _saved_json],
                key="xo_load_crossovers_saved_pick",
            )
            st.button(
                "Load saved crossover file",
                key="xo_load_crossovers_saved_btn",
                on_click=_load_crossover_saved_callback,
            )


if hasattr(st, "fragment"):
    _render_crossover_draft_picks = st.fragment(_render_crossover_draft_picks)

st.set_page_config(
    page_title="Crossover Analysis",
    page_icon="🧩",
    layout="wide",
)

init_session_state()

for _legacy_key in (
    "xo_table_pick_draft",
    "xo_freq_table_df",
    "xo_freq_data_editor",
):
    st.session_state.pop(_legacy_key, None)

st.title("🧩 Crossover analysis")

st.markdown(
    "Inspect crossover placements and the parent structure for designs produced on the "
    "**RASPP Design** page."
)

has_raspp = (
    SESSION_KEYS["raspp_results"] in st.session_state
    and st.session_state[SESSION_KEYS["raspp_results"]] is not None
)
has_contacts = (
    SESSION_KEYS["schema_contacts"] in st.session_state
    and st.session_state[SESSION_KEYS["schema_contacts"]] is not None
)

has_multi_fragment = (
    SESSION_KEYS["multi_fragment_results"] in st.session_state
    and st.session_state[SESSION_KEYS["multi_fragment_results"]] is not None
    and isinstance(st.session_state[SESSION_KEYS["multi_fragment_results"]], dict)
    and len(st.session_state[SESSION_KEYS["multi_fragment_results"]]) > 0
)

if not has_raspp and not has_multi_fragment:
    st.warning("No RASPP results yet. Run multi-fragment analysis on the **RASPP Design** page first.")
else:
    all_designs = []
    filtered_designs = []
    raspp_results = None

    if has_multi_fragment:
        multi_fragment_results = st.session_state[SESSION_KEYS["multi_fragment_results"]]

        for num_fragments, fragment_data in multi_fragment_results.items():
            if "error" not in fragment_data and "results" in fragment_data:
                for design in fragment_data["results"]:
                    design_with_frag_count = design.copy()
                    design_with_frag_count["fragment_count"] = num_fragments
                    design_with_frag_count["num_crossovers"] = len(design.get("crossovers", []))
                    all_designs.append(design_with_frag_count)

        if all_designs:
            all_designs.sort(
                key=lambda x: (
                    x.get("average_energy", float("inf"))
                    if x.get("average_energy") is not None
                    else float("inf")
                )
            )

            fragment_counts = sorted(set(d["fragment_count"] for d in all_designs))
            min_fragments = min(fragment_counts) if fragment_counts else 1
            max_fragments = max(fragment_counts) if fragment_counts else 1

            st.subheader("Filter designs")
            fragment_range = st.slider(
                "Fragment count range",
                min_value=min_fragments,
                max_value=max_fragments,
                value=(min_fragments, max_fragments),
                step=1,
                help="Include only designs whose fragment count falls in this range in the distribution.",
            )

            min_frag, max_frag = fragment_range
            filtered_designs = [
                d for d in all_designs if min_frag <= d["fragment_count"] <= max_frag
            ]

            if len(filtered_designs) < len(all_designs):
                st.info(
                    f"Distribution uses {len(filtered_designs)} of {len(all_designs)} designs "
                    f"({min_frag}–{max_frag} fragments)"
                )
            else:
                st.info(f"Distribution uses all {len(all_designs)} designs ({min_frag}–{max_frag} fragments)")
        else:
            st.warning("No designs found in multi-fragment results.")

        parents = st.session_state.get(SESSION_KEYS["raspp_parents"])
        if parents and len(parents) > 0:
            parent_length = len(parents[0])
        else:
            st.warning("Parent sequences not available.")
            parent_length = None

    elif has_raspp:
        raspp_results = st.session_state[SESSION_KEYS["raspp_results"]]
        all_designs = raspp_results
        filtered_designs = raspp_results
        parents = st.session_state.get(SESSION_KEYS["raspp_parents"])

        if parents and len(parents) > 0:
            parent_length = len(parents[0])
        else:
            st.warning("Parent sequences not available.")
            parent_length = None
    else:
        parent_length = None

    results_for_comparison = []
    if has_multi_fragment and filtered_designs:
        results_for_comparison = filtered_designs
    elif has_multi_fragment and all_designs:
        results_for_comparison = all_designs
    elif raspp_results is not None:
        results_for_comparison = raspp_results

    if parent_length and results_for_comparison:
        st.subheader("Crossover position distribution")
        energy_percent_filter = st.slider(
            "Include only lowest X% energy designs",
            min_value=1,
            max_value=100,
            value=100,
            step=1,
        )

        if energy_percent_filter < 100:
            sorted_by_energy = sorted(
                results_for_comparison,
                key=lambda x: (
                    x.get("average_energy", float("inf"))
                    if x.get("average_energy") is not None
                    else float("inf")
                ),
            )
            num_to_show = max(
                1,
                int(len(sorted_by_energy) * energy_percent_filter / 100),
            )
            filtered_results = sorted_by_energy[:num_to_show]
            st.caption(
                f"Using {len(filtered_results)} of {len(results_for_comparison)} designs."
            )
        else:
            filtered_results = results_for_comparison

        all_crossovers_list = [r.get("crossovers", []) for r in filtered_results]

        _applied = list(st.session_state.get("selected_crossover_positions") or [])
        fig_freq_hist = plot_crossover_distribution(
            all_crossovers_list,
            parent_length,
            selected_positions_1based=_applied,
        )
        if fig_freq_hist:
            clicked_data = st.plotly_chart(
                fig_freq_hist,
                use_container_width=True,
                on_select="rerun",
                key="crossover_dist_chart",
            )

            chart_pending = list(st.session_state.get("_xo_chart_pending") or [])
            if clicked_data and "selection" in clicked_data:
                sel = clicked_data["selection"]
                if "points" in sel:
                    for point in sel["points"]:
                        if "x" in point:
                            chart_pending.append(int(point["x"]) + 1)
                    st.session_state["_xo_chart_pending"] = sorted(set(chart_pending))

            position_counts = compute_crossover_position_counts(
                all_crossovers_list, parent_length
            )
            max_freq = max(position_counts.values()) if position_counts else 0
            table_threshold = 0.05 * max_freq
            eligible = [
                (p0 + 1, position_counts[p0])
                for p0 in sorted(position_counts.keys())
                if position_counts[p0] > table_threshold
            ]

            _render_crossover_draft_picks(
                eligible,
                table_threshold,
                max_freq,
                _applied,
            )
            _render_crossover_selection_controls(parent_length, eligible)

            _xo_applied = tuple(
                st.session_state.get("selected_crossover_positions") or []
            )
            if _xo_applied:
                st.subheader("MSA with crossover markers")
                alignment_data = None
                msa_path = st.session_state.get(SESSION_KEYS["msa_path"])

                if msa_path and Path(msa_path).exists():
                    try:
                        from schema_raspp import schema

                        with open(msa_path, "r") as f:
                            alignment_data = schema.readMultipleSequenceAlignmentFile(f)
                    except Exception as e:
                        st.warning(f"Could not read MSA file: {e}")

                if not alignment_data:
                    contacts_data = st.session_state.get(SESSION_KEYS["schema_contacts"], {})
                    if "parents" in contacts_data:
                        parents_list = contacts_data["parents"]
                        if parents_list:
                            if isinstance(parents_list[0], tuple):
                                alignment_data = parents_list
                            else:
                                alignment_data = [
                                    (f"Seq{i + 1}", seq) for i, seq in enumerate(parents_list)
                                ]

                if alignment_data:
                    msa_html = format_msa_with_crossovers_html(
                        alignment_data,
                        list(_xo_applied),
                        max_sequences=10,
                    )
                    if msa_html:
                        _rev = int(st.session_state.get("_xo_applied_revision", 0))
                        st.markdown(
                            f"<!-- crossover_revision:{_rev} -->\n{msa_html}",
                            unsafe_allow_html=True,
                        )
                else:
                    st.warning("MSA data not available for HTML alignment.")
        else:
            st.info("No crossover positions in the selected designs to plot.")

    _contacts_sd = st.session_state.get(SESSION_KEYS["schema_contacts"])
    _pdb_s = _contacts_sd.get("pdb_structure") if isinstance(_contacts_sd, dict) else None
    if parent_length and _pdb_s:
        st.markdown("---")
        st.subheader("Structure (3Dmol) — red/green by segment")
        _xo_v = list(st.session_state.get("selected_crossover_positions") or [])
        st.caption(
            "3Dmol viewer: continuous cartoon with segments alternating **red** and **green** "
            "between SCHEMA crossovers. Crossovers use the **same alignment-column frame** as "
            "contacts (see `renumber_pdb_to_alignment` in `utils/schemarecomb_bridge.py`); "
            "they are **1-based** there. Residue labels in the viewer follow **1…N display order**, "
            "not raw PDB auth numbers."
        )
        _ref = (
            len(_pdb_s.renumbering_seq)
            if getattr(_pdb_s, "renumbering_seq", None) is not None
            else parent_length
        )
        if _ref != parent_length:
            st.warning(
                f"Structure alignment length ({_ref}) differs from session parent length "
                f"({parent_length}). Coloring uses the structure’s alignment frame."
            )
        if not _xo_v:
            st.info(
                "No crossovers selected — the structure is shown as a single segment. "
                "Select crossover positions above to split red/green coloring."
            )
        st.markdown(
            f"<span style='color:#d32f2f;font-weight:600;'>■</span> Odd segments (1, 3, …) &nbsp; "
            f"<span style='color:#2e7d32;font-weight:600;'>■</span> Even segments (2, 4, …)",
            unsafe_allow_html=True,
        )
        try:
            from utils.structure_viewer import display_crossover_fragment_structure

            _n_frag = display_crossover_fragment_structure(
                _pdb_s,
                _xo_v,
                parent_length,
                view_width=800,
                view_height=650,
                render_revision=int(
                    st.session_state.get("_xo_applied_revision", 0)
                ),
            )
            if _xo_v:
                st.caption(
                    f"{_n_frag} fragment(s); 1-based crossovers: "
                    + ", ".join(str(p) for p in sorted(_xo_v))
                )
        except ImportError as e:
            st.warning(
                "Structure viewer requires **py3Dmol**, **stmol**, and **ipython_genutils** "
                "(stmol dependency). Install with: "
                "`pip install py3Dmol stmol ipython_genutils`"
            )
            st.caption(str(e))
        except Exception as e:
            st.warning(f"Could not render structure viewer: {e}")
    elif (
        parent_length
        and has_contacts
        and isinstance(_contacts_sd, dict)
        and _pdb_s is None
    ):
        st.caption(
            "Structure viewer needs a PDB on **SCHEMA Energy** so contacts include a structure object."
        )

    if parent_length is None and (has_raspp or has_multi_fragment):
        st.info("Need parent sequences to show the crossover distribution.")
    elif not results_for_comparison:
        st.info("No designs available for the crossover distribution.")
