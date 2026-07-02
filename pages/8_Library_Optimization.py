"""
Page 8: Library Optimization
Iterative ESM2 + ridge regression on random chimeras from saved fragment pools.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

import pandas as pd
import plotly.graph_objects as go
import streamlit as st

sys.path.insert(0, str(Path(__file__).parent.parent))

from utils.library_optimization.convergence import (
    create_rank_comparison_figure,
    rank_comparison_figures_from_history,
)
from utils.library_optimization.esm_scoring import esm_dependencies_available
from utils.library_optimization.optimizer import run_optimization_rounds
from utils.library_optimization.pools import (
    blocks_signature,
    blocks_to_saved_selections,
    build_filter_table_rows,
    copy_blocks,
    keep_only_variants_in_block,
    keep_top_homologs_by_coefficient,
    load_block_pools_from_session,
    pools_summary,
    remove_variants_from_block,
)
from utils.session_manager import init_session_state

st.set_page_config(
    page_title="Library Optimization",
    page_icon="📊",
    layout="wide",
)

init_session_state()

_OPT_VERSION = "1.1"


def _init_opt_state() -> None:
    if "library_opt_results" not in st.session_state:
        st.session_state.library_opt_results = None
    if "library_opt_working_blocks" not in st.session_state:
        st.session_state.library_opt_working_blocks = None
    if "library_opt_working_sig" not in st.session_state:
        st.session_state.library_opt_working_sig = None


def _sync_working_blocks(source_blocks: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    sig = blocks_signature(source_blocks)
    if (
        st.session_state.library_opt_working_blocks is None
        or st.session_state.library_opt_working_sig != sig
    ):
        st.session_state.library_opt_working_blocks = copy_blocks(source_blocks)
        st.session_state.library_opt_working_sig = sig
    return st.session_state.library_opt_working_blocks


def _format_mutations(value: Any) -> str:
    return "—" if value is None else str(value)


def _format_pct_identity(value: Any) -> str:
    if isinstance(value, (int, float)):
        return f"{float(value):.1f}%"
    return "—"


def _render_filter_prune_section(
    source_blocks: List[Dict[str, Any]],
    coefficients: Dict[str, List[Dict[str, Any]]],
    selections: Dict[str, Any],
) -> None:
    working_blocks = _sync_working_blocks(source_blocks)
    frag_keys = sorted(
        {b["fragment_key"] for b in working_blocks},
        key=lambda k: int(k),
    )
    if not frag_keys:
        st.warning("No fragments available to filter.")
        return

    st.markdown(
        "Each **coefficient** estimates how a fragment variant shifts chimera **ESM2** scores "
        "relative to the query fragment (query = 0). **Higher is better**: positive values "
        "tend to improve predicted fitness; **negative values** suggest worse contribution "
        "than the query. Variants are sorted by coefficient (highest first); the query stays "
        "at the top. Select homologs below, or keep the top N by coefficient."
    )

    selected_frag = st.selectbox(
        "Fragment",
        options=frag_keys,
        format_func=lambda k: f"Fragment {k}",
        key="library_opt_filter_fragment",
    )

    block = next(b for b in working_blocks if b["fragment_key"] == selected_frag)
    selection_block = selections.get(selected_frag) or selections.get(str(selected_frag))
    table_rows = build_filter_table_rows(block, coefficients, selection_block)
    homolog_rows = [r for r in table_rows if not r["is_query"]]

    st.caption(
        f"**{len(table_rows)}** sequence(s) in this fragment "
        f"({len(homolog_rows)} homolog(s) + query)."
    )

    header = st.columns([0.06, 0.14, 0.28, 0.10, 0.08, 0.10, 0.10, 0.10])
    header[0].markdown("**Select**")
    header[1].markdown("**Sequence ID**")
    header[2].markdown("**Sequence**")
    header[3].markdown("**Coefficient**")
    header[4].markdown("**Length**")
    header[5].markdown("**Mut. (non-gap)**")
    header[6].markdown("**Mut. (w/ gaps)**")
    header[7].markdown("**% identity**")

    selected_ids: List[str] = []
    for row in table_rows:
        cols = st.columns([0.06, 0.14, 0.28, 0.10, 0.08, 0.10, 0.10, 0.10])
        with cols[0]:
            if row["is_query"]:
                st.markdown("—")
            else:
                picked = st.checkbox(
                    "Select",
                    key=f"library_opt_pick_{selected_frag}_{row['row_id']}",
                    label_visibility="collapsed",
                )
                if picked:
                    selected_ids.append(row["row_id"])
        cols[1].markdown(row["sequence_id"])
        cols[2].code(row["sequence"] or "", language=None)
        coef = row["coefficient"]
        cols[3].markdown(f"{coef:.4f}" if row["is_query"] else f"{coef:.4f}")
        cols[4].markdown(str(row["length_display"]))
        cols[5].markdown(_format_mutations(row["mutations_non_gap"]))
        cols[6].markdown(_format_mutations(row["mutations_with_gaps"]))
        cols[7].markdown(_format_pct_identity(row["pct_identity"]))

    btn1, btn2, btn3 = st.columns(3)
    with btn1:
        remove_selected = st.button(
            "Remove selected sequences",
            key="library_opt_remove_selected",
        )
    with btn2:
        keep_selected = st.button(
            "Keep only selected sequences",
            key="library_opt_keep_selected",
        )
    with btn3:
        save_filtered = st.button(
            "Save filtered list to session",
            type="primary",
            key="library_opt_save_filtered",
        )

    top_col, top_btn_col = st.columns([1, 2])
    with top_col:
        top_n = st.number_input(
            "Top homologs to keep (by coefficient)",
            min_value=1,
            max_value=max(1, len(homolog_rows)),
            value=min(5, max(1, len(homolog_rows))),
            step=1,
            key=f"library_opt_top_n_{selected_frag}",
            help="Keeps the query plus the N highest-scoring homologs on this fragment.",
        )
    with top_btn_col:
        st.markdown("")  # align with number input
        keep_top = st.button(
            "Keep top by coefficient",
            key="library_opt_keep_top",
        )

    def _validate_block(block_after: Dict[str, Any]) -> bool:
        if block_after["n_variants"] < 2:
            st.error(
                f"Fragment {selected_frag} must keep at least the query and one homolog."
            )
            return False
        return True

    if remove_selected:
        if not selected_ids:
            st.warning("Select one or more homolog sequences to remove.")
        else:
            updated = remove_variants_from_block(
                working_blocks, selected_frag, selected_ids
            )
            block_after = next(
                b for b in updated if b["fragment_key"] == selected_frag
            )
            if _validate_block(block_after):
                st.session_state.library_opt_working_blocks = updated
                st.success(f"Removed {len(selected_ids)} sequence(s) from fragment {selected_frag}.")
                st.rerun()

    if keep_selected:
        if not selected_ids:
            st.warning("Select homolog sequences to keep.")
        else:
            updated = keep_only_variants_in_block(
                working_blocks, selected_frag, selected_ids
            )
            block_after = next(
                b for b in updated if b["fragment_key"] == selected_frag
            )
            if _validate_block(block_after):
                st.session_state.library_opt_working_blocks = updated
                st.success(
                    f"Kept {len(selected_ids)} selected homolog(s) on fragment {selected_frag}."
                )
                st.rerun()

    if keep_top:
        try:
            updated = keep_top_homologs_by_coefficient(
                working_blocks,
                selected_frag,
                coefficients,
                int(top_n),
            )
        except ValueError as exc:
            st.error(str(exc))
        else:
            block_after = next(
                b for b in updated if b["fragment_key"] == selected_frag
            )
            if _validate_block(block_after):
                st.session_state.library_opt_working_blocks = updated
                kept_homologs = block_after["n_variants"] - 1
                st.success(
                    f"Kept query + top {kept_homologs} homolog(s) by coefficient "
                    f"on fragment {selected_frag}."
                )
                st.rerun()

    if save_filtered:
        filtered_selections = blocks_to_saved_selections(
            working_blocks, coefficients
        )
        st.session_state.diversity_saved_selections = filtered_selections
        st.session_state.library_opt_working_blocks = copy_blocks(working_blocks)
        st.session_state.library_opt_working_sig = blocks_signature(working_blocks)
        from utils.session_manager import auto_save

        auto_save("diversity_selections_saved")
        st.success(
            "Filtered pools saved to session. "
            "**Simulate with AI**, **Oligopool Design**, and later runs on this page "
            "will use this list."
        )
        st.rerun()


def _results_for_export(results: Dict[str, Any]) -> Dict[str, Any]:
    """Trim chimera sequences for lighter JSON if needed — keep full for now."""
    return {
        "version": _OPT_VERSION,
        "params": results.get("params"),
        "rounds": results.get("rounds"),
        "n_chimeras": results.get("n_chimeras"),
        "converged": results.get("converged"),
        "status": results.get("status"),
        "coefficients_by_block": results.get("coefficients_by_block"),
        "blocks_summary": results.get("blocks_summary"),
        "chimera_score_summary": [
            {
                "chimera_id": c.get("chimera_id"),
                "esm_score": c.get("esm_score"),
                "sequence_length": len(c.get("sequence") or ""),
            }
            for c in (results.get("chimeras") or [])
        ],
    }


def build_convergence_figure(
    rounds: List[Dict[str, Any]],
    *,
    threshold: Optional[float] = None,
) -> Optional[go.Figure]:
    """Spearman ρ vs round — one connected line+markers series per fragment."""
    if not rounds:
        return None

    block_keys = sorted(
        {
            frag_key
            for r in rounds
            for frag_key in (r.get("spearman_by_block") or {})
        },
        key=lambda k: int(k),
    )
    if not block_keys:
        return None

    fig = go.Figure()
    for frag_key in block_keys:
        xs, ys = [], []
        for r in rounds:
            rho = (r.get("spearman_by_block") or {}).get(frag_key)
            if rho is not None:
                xs.append(r["round"])
                ys.append(rho)
        if not xs:
            continue
        fig.add_trace(
            go.Scatter(
                x=xs,
                y=ys,
                mode="lines+markers",
                name=f"Fragment {frag_key}",
                line=dict(width=2),
                marker=dict(size=8),
                connectgaps=False,
            )
        )

    if threshold is not None:
        fig.add_hline(
            y=float(threshold),
            line=dict(color="rgba(80,80,80,0.55)", dash="dash", width=1),
            annotation_text=f"Threshold ρ = {threshold:.2f}",
            annotation_position="top right",
        )

    max_round = max(int(r["round"]) for r in rounds)
    fig.update_layout(
        title="Spearman ρ by round (variant rank stability vs previous round)",
        xaxis_title="Round",
        yaxis_title="Spearman ρ",
        yaxis=dict(range=[0, 1.05]),
        height=420,
        legend=dict(
            title="Fragment",
            yanchor="top",
            y=1.0,
            xanchor="left",
            x=1.02,
        ),
        margin=dict(r=180, t=60, b=50),
        xaxis=dict(dtick=1, range=[1.5, max_round + 0.5] if max_round >= 2 else None),
    )
    return fig


def _render_convergence_chart(
    rounds: List[Dict[str, Any]],
    *,
    threshold: Optional[float] = None,
    chart_key: Optional[str] = None,
) -> None:
    fig = build_convergence_figure(rounds, threshold=threshold)
    if fig is None:
        return
    st.markdown("### Spearman ρ by round")
    st.caption(
        "Each line is one fragment: Spearman correlation of variant coefficient ranks "
        "vs the previous round (first point appears at round 2)."
    )
    st.plotly_chart(
        fig,
        use_container_width=True,
        key=chart_key or "library_opt_spearman_convergence",
    )


def _coefficient_history_from_results(results: Dict[str, Any]) -> List[tuple]:
    stored = results.get("coefficient_history")
    if stored:
        return [
            (int(item["round"]), item["coefficients_by_block"])
            for item in stored
            if item.get("coefficients_by_block")
        ]
    rounds = results.get("rounds") or []
    return [
        (int(r["round"]), r["coefficients_by_block"])
        for r in rounds
        if r.get("coefficients_by_block")
    ]


def _render_rank_comparison_plots(results: Dict[str, Any]) -> None:
    history = _coefficient_history_from_results(results)
    figures = rank_comparison_figures_from_history(history)
    if not figures:
        return
    st.markdown("### Rank previous vs rank new")
    st.caption(
        "Within each fragment, variants are ranked by ridge coefficient "
        "(1 = highest). Points on the dashed diagonal kept the same rank."
    )
    for _prev_round, curr_round, fig in figures:
        st.plotly_chart(fig, use_container_width=True, key=f"rank_cmp_round_{curr_round}")


def _coefficients_dataframe(rows: List[Dict[str, Any]]) -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "Sequence ID": r.get("sequence_id"),
                "Coefficient": round(float(r.get("coefficient", 0.0)), 5),
                "Exposure": r.get("exposure_count", 0),
                "Query ref.": "✓" if r.get("is_query") else "",
                "Sequence": r.get("sequence", ""),
            }
            for r in rows
        ]
    )


def _run_optimization_ui(
    blocks: List[Dict[str, Any]],
    *,
    prior: Optional[Dict[str, Any]],
) -> None:
    """Render run/continue controls and execute rounds."""
    continuing = prior is not None
    prior_params = (prior or {}).get("params") or {}
    prior_rounds = (prior or {}).get("rounds") or []

    if continuing:
        st.info(
            f"**Continuing** from {len(prior_rounds)} completed round(s) "
            f"with **{prior.get('n_chimeras', 0)}** scored chimeras. "
            "New rounds add to this pool — nothing already scored is discarded."
        )

    col1, col2, col3 = st.columns(3)
    with col1:
        spearman_threshold = st.slider(
            "Spearman ρ threshold",
            min_value=0.5,
            max_value=0.99,
            value=float(prior_params.get("spearman_threshold", 0.9)),
            step=0.01,
            help="Used when early stopping is enabled.",
            key="library_opt_spearman_threshold",
        )
        stop_on_convergence = st.checkbox(
            "Stop when Spearman threshold is reached",
            value=bool(prior_params.get("stop_on_convergence", True)),
            help=(
                "When checked, the run may end before all requested rounds once every "
                "fragment's rank correlation exceeds the threshold (after minimum rounds). "
                "When unchecked, always runs the full number of rounds requested."
            ),
            key="library_opt_stop_on_convergence",
        )
        batch_size = st.number_input(
            "Chimeras per round",
            min_value=1,
            max_value=10000,
            value=int(prior_params.get("batch_size", 2000)),
            step=1,
            help="Random chimeras sampled and ESM-scored each round.",
            key="library_opt_batch_size",
        )
    with col2:
        if continuing:
            min_rounds = int(prior_params.get("min_rounds", 2))
            rounds_to_run = st.number_input(
                "Additional rounds",
                min_value=1,
                max_value=50,
                value=5,
                step=1,
                help="How many more rounds to run in this step.",
                key="library_opt_additional_rounds",
            )
        else:
            min_rounds = st.number_input(
                "Minimum rounds",
                min_value=1,
                max_value=10,
                value=2,
                key="library_opt_min_rounds",
            )
            rounds_to_run = st.number_input(
                "Rounds to run",
                min_value=1,
                max_value=50,
                value=20,
                step=1,
                help="Maximum rounds for this run (may stop earlier if early stopping is on).",
                key="library_opt_max_rounds",
            )
    with col3:
        ridge_alpha = st.number_input(
            "Ridge α",
            min_value=0.01,
            max_value=100.0,
            value=float(prior_params.get("ridge_alpha", 1.0)),
            step=0.1,
            key="library_opt_ridge_alpha",
        )
        esm_batch_size = st.number_input(
            "ESM batch size",
            min_value=1,
            max_value=32,
            value=int(prior_params.get("esm_batch_size", 8)),
            help="Sequences scored in parallel per GPU/CPU batch.",
            key="library_opt_esm_batch_size",
        )
        pause_after_round = st.checkbox(
            "Pause after each round",
            value=False,
            help=(
                "Run exactly one round per click so you can adjust chimeras per round "
                "or other settings before continuing."
            ),
            key="library_opt_pause_after_round",
        )

    if pause_after_round:
        rounds_to_run = 1

    button_label = "Continue optimization" if continuing else "Start optimization"
    run_clicked = st.button(button_label, type="primary", key="library_opt_run")

    if not run_clicked:
        return

    n_rounds = int(rounds_to_run)
    progress_bar = st.progress(0.0)
    status = st.empty()
    round_status = st.empty()
    spearman_chart_slot = st.empty()
    rank_plots_area = st.container()
    live_rounds: List[Dict[str, Any]] = list(prior_rounds) if continuing else []
    total_rounds_this_run = n_rounds
    start_round_num = (len(prior_rounds) + 1) if continuing else 1

    def on_progress(event: Dict[str, Any]) -> None:
        phase = event.get("phase", "")
        rnd = int(event.get("round") or 0)
        msg = event.get("message") or ""
        round_index_in_run = rnd - start_round_num + 1

        if phase == "esm" and "fraction" in event:
            base = (round_index_in_run - 1) / total_rounds_this_run
            span = 1.0 / total_rounds_this_run
            progress_bar.progress(
                min(1.0, base + span * float(event["fraction"]))
            )
        elif phase == "round_complete":
            progress_bar.progress(
                min(1.0, round_index_in_run / total_rounds_this_run)
            )
            live_rounds.append(
                {
                    "round": rnd,
                    "spearman_by_block": event.get("spearman_by_block") or {},
                }
            )
            partial = event.get("partial_results")
            if partial:
                st.session_state.library_opt_results = partial
                try:
                    from utils.session_manager import auto_save

                    auto_save("library_optimization_round")
                except Exception:
                    pass

            if rnd >= 2:
                spearman_fig = build_convergence_figure(
                    live_rounds,
                    threshold=float(spearman_threshold),
                )
                if spearman_fig is not None:
                    with spearman_chart_slot.container():
                        st.markdown("### Spearman ρ by round")
                        st.caption(
                            "Each line is one fragment: Spearman ρ of coefficient ranks "
                            "vs the previous round."
                        )
                        st.plotly_chart(
                            spearman_fig,
                            use_container_width=True,
                            key=f"spearman_live_round_{rnd}",
                        )
                with rank_plots_area:
                    if rnd == 2 and not continuing:
                        st.markdown("### Rank previous vs rank new")
                        st.caption(
                            "Within each fragment, rank 1 = highest ridge coefficient. "
                            "Diagonal = unchanged rank between rounds."
                        )
            prev_coef = event.get("prev_coefficients_by_block")
            curr_coef = event.get("coefficients_by_block")
            if rnd >= 2 and prev_coef and curr_coef:
                fig = create_rank_comparison_figure(
                    prev_coef,
                    curr_coef,
                    prev_round=rnd - 1,
                    curr_round=rnd,
                )
                with rank_plots_area:
                    st.plotly_chart(
                        fig,
                        use_container_width=True,
                        key=f"rank_cmp_live_round_{rnd}",
                    )
        status.markdown(msg)
        if event.get("spearman_by_block"):
            rho_lines = [
                f"Frag {k}: {v:.3f}" if v is not None else f"Frag {k}: —"
                for k, v in sorted(
                    event["spearman_by_block"].items(),
                    key=lambda item: int(item[0]),
                )
            ]
            round_status.caption("Spearman ρ — " + " | ".join(rho_lines))

    with st.spinner("Running optimization (this may take several minutes)…"):
        try:
            results = run_optimization_rounds(
                blocks,
                batch_size=int(batch_size),
                n_rounds=n_rounds,
                min_rounds=int(min_rounds),
                spearman_threshold=float(spearman_threshold),
                stop_on_convergence=bool(stop_on_convergence),
                ridge_alpha=float(ridge_alpha),
                esm_batch_size=int(esm_batch_size),
                prior=prior,
                progress_callback=on_progress,
            )
            st.session_state.library_opt_results = results
            from utils.session_manager import auto_save

            auto_save("library_optimization_complete")
            progress_bar.progress(1.0)
            n_done = len(results.get("rounds") or [])
            n_chimeras = results.get("n_chimeras", 0)
            if results.get("converged"):
                st.success(
                    f"Spearman threshold reached after {n_done} total round(s) "
                    f"({n_chimeras} chimeras scored). "
                    "Adjust settings below and click **Continue optimization** to score more."
                )
            elif pause_after_round:
                st.success(
                    f"Round {n_done} complete ({n_chimeras} chimeras scored). "
                    "Paused — adjust settings and click **Continue optimization** for the next round."
                )
            else:
                st.info(
                    f"Finished this run: {n_done} total round(s), "
                    f"{n_chimeras} chimeras scored. "
                    "Run more rounds below if needed."
                )
        except Exception as exc:
            st.error(f"Optimization failed: {exc}")


_init_opt_state()

st.title("📊 Library Optimization")

st.markdown(
    "Iteratively sample random chimeras from your **saved fragment pools** "
    "(Diversity Analysis → **Save main list to session**), score with **ESM2-150M**, "
    "and fit **ridge regression** to estimate each block variant's contribution. "
    "Pruned pools can be saved back to session for downstream pages. "
    "You can **pause between rounds**, change chimeras per round, and **continue** "
    "later without losing scored designs."
)

selections = st.session_state.get("diversity_saved_selections") or {}
blocks, load_err = load_block_pools_from_session(selections)

if load_err:
    st.warning(load_err)
    st.stop()

summary = pools_summary(blocks)
st.success(
    f"Loaded **{summary['n_blocks']}** block(s) — "
    f"{summary['variants_per_block']} variant(s) per block "
    f"({summary['total_combinations']:,} possible combinations)."
)

esm_ok, esm_err = esm_dependencies_available()
if not esm_ok:
    st.error(esm_err)
    st.stop()

with st.expander("Input pools (read-only)", expanded=False):
    for block in blocks:
        st.markdown(f"**Fragment {block['fragment_key']}** — {block['n_variants']} variants")
        preview = pd.DataFrame(
            [
                {
                    "Sequence ID": v["sequence_id"],
                    "Length": len(v["sequence"]),
                    "Query": "✓" if v.get("is_query") else "",
                }
                for v in block["variants"][:15]
            ]
        )
        st.dataframe(preview, use_container_width=True, hide_index=True)
        if block["n_variants"] > 15:
            st.caption(f"… and {block['n_variants'] - 15} more")

st.markdown("---")
existing_results: Optional[Dict[str, Any]] = st.session_state.get("library_opt_results")

if existing_results:
    st.subheader("Continue optimization")
else:
    st.subheader("Run optimization")

_run_optimization_ui(blocks, prior=existing_results)

results: Optional[Dict[str, Any]] = st.session_state.get("library_opt_results")

if not results:
    st.info("Configure parameters above and click **Start optimization**.")
    st.stop()

st.markdown("---")
st.subheader("Results")

rounds = results.get("rounds") or []
fit_last = (rounds[-1].get("fit_stats") or {}) if rounds else {}
metric_cols = st.columns(4)
metric_cols[0].metric("Chimeras scored", results.get("n_chimeras", 0))
metric_cols[1].metric("Rounds", len(rounds))
metric_cols[2].metric(
    "Converged",
    "Yes" if results.get("converged") else "No",
)
metric_cols[3].metric(
    "Ridge R²",
    f"{fit_last['r2']:.3f}" if fit_last.get("r2") is not None else "—",
)

_render_convergence_chart(
    rounds,
    threshold=float((results.get("params") or {}).get("spearman_threshold", 0.9)),
)
_render_rank_comparison_plots(results)

coefficients = results.get("coefficients_by_block") or {}

st.markdown("### Coefficients by fragment")
frag_tabs = st.tabs([f"Fragment {k}" for k in sorted(coefficients.keys(), key=int)])
for tab, frag_key in zip(frag_tabs, sorted(coefficients.keys(), key=int)):
    with tab:
        st.dataframe(
            _coefficients_dataframe(coefficients[frag_key]),
            use_container_width=True,
            hide_index=True,
        )

st.markdown("---")
st.subheader("Filter & prune")

_render_filter_prune_section(blocks, coefficients, selections)

st.markdown("---")
st.subheader("Export")

export_payload = _results_for_export(results)
export_json = json.dumps(export_payload, indent=2) + "\n"

dl_col, clear_col = st.columns(2)
with dl_col:
    st.download_button(
        "Download results JSON",
        data=export_json,
        file_name="library_optimization_results.json",
        mime="application/json",
        key="library_opt_download",
    )
with clear_col:
    if st.button("Clear optimization results", key="library_opt_clear"):
        st.session_state.library_opt_results = None
        st.session_state.library_opt_working_blocks = None
        st.session_state.library_opt_working_sig = None
        st.rerun()

with st.expander("Round history", expanded=False):
    for r in rounds:
        spearman = r.get("spearman_by_block") or {}
        rho_str = ", ".join(
            f"{k}:{v:.3f}" if v is not None else f"{k}:—"
            for k, v in sorted(spearman.items(), key=lambda item: int(item[0]))
        )
        st.markdown(
            f"**Round {r['round']}** — {r['n_total_chimeras']} chimeras | "
            f"ρ: {rho_str} | converged: {r.get('converged')}"
        )
