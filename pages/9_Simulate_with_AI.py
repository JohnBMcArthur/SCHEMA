"""
Page 9: Simulate with AI
Random chimera ESM2 simulation vs wildtype (all-query chimera).
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

import pandas as pd
import plotly.graph_objects as go
import statistics
import streamlit as st

sys.path.insert(0, str(Path(__file__).parent.parent))

from utils.library_optimization.esm_scoring import esm_dependencies_available
from utils.library_optimization.pools import load_block_pools_from_session, pools_summary
from utils.library_optimization.simulation import (
    empty_sim_results,
    merge_sim_run,
    pools_signature_for_sim,
    score_wildtype,
    simulate_random_chimeras_esm,
    top_scored_table_rows,
    top_sequences_to_fasta,
)
from utils.session_manager import init_session_state

st.set_page_config(
    page_title="Simulate with AI",
    page_icon="🤖",
    layout="wide",
)

init_session_state()


def _init_sim_state() -> None:
    if "ai_sim_results" not in st.session_state:
        st.session_state.ai_sim_results = None


def _sync_sim_results(blocks: List[Dict[str, Any]]) -> Dict[str, Any]:
    sig = pools_signature_for_sim(blocks)
    results = st.session_state.get("ai_sim_results")
    if results and results.get("pools_signature") != sig:
        st.session_state.ai_sim_results = None
        results = None
    return results or empty_sim_results()


def build_delta_histogram(
    score_deltas: List[float],
    *,
    title: Optional[str] = None,
) -> go.Figure:
    fig = go.Figure()
    if score_deltas:
        fig.add_trace(
            go.Histogram(
                x=score_deltas,
                nbinsx=min(40, max(10, len(score_deltas) // 5)),
                name="Chimeras",
                marker=dict(color="#4C78A8", line=dict(color="#2F4B6B", width=1)),
            )
        )
    fig.add_vline(
        x=0.0,
        line=dict(color="rgba(200,60,60,0.9)", width=2, dash="dash"),
        annotation_text="Wildtype (Δ = 0)",
        annotation_position="top right",
    )
    fig.update_layout(
        title=title or "ESM2 score − wildtype",
        xaxis_title="Δ ESM2 score (chimera − wildtype)",
        yaxis_title="Count",
        height=420,
        bargap=0.05,
        showlegend=False,
        margin=dict(t=60, b=50),
    )
    return fig


def _results_to_display_dataframe(results: Dict[str, Any]) -> pd.DataFrame:
    rows = top_scored_table_rows(results, top_n=20)
    return pd.DataFrame(
        [
            {
                "ID": r["chimera_id"],
                "Sequence": r["sequence"],
                "ESM2 score": round(float(r["esm_score"]), 4),
                "Δ vs wildtype": round(float(r["score_delta"]), 4),
            }
            for r in rows
        ]
    )


def _run_simulation(
    blocks: List[Dict[str, Any]],
    *,
    n_samples: int,
    esm_batch_size: int,
    prior: Optional[Dict[str, Any]],
) -> None:
    sig = pools_signature_for_sim(blocks)
    prior = prior if (prior and prior.get("pools_signature") == sig) else None
    start_index = int((prior or {}).get("n_chimeras") or 0)
    continuing = start_index > 0

    progress = st.progress(0.0)
    status = st.empty()
    hist_slot = st.empty()
    table_slot = st.empty()

    wildtype_score: float = 0.0
    wildtype_sequence: str = ""

    prior_deltas: List[float] = list((prior or {}).get("score_deltas") or [])
    prior_records: List[Dict[str, Any]] = list((prior or {}).get("chimeras") or [])

    def on_batch(
        run_records_so_far: List[Dict[str, Any]],
        run_deltas_so_far: List[float],
        n_done: int,
        n_total: int,
    ) -> None:
        combined_records = prior_records + run_records_so_far
        combined_deltas = prior_deltas + run_deltas_so_far
        progress.progress(n_done / n_total)
        status.markdown(f"ESM scoring {n_done}/{n_total} sequences this run…")
        partial_results = {
            "wildtype_score": wildtype_score,
            "wildtype_sequence": wildtype_sequence,
            "chimeras": combined_records,
            "score_deltas": combined_deltas,
        }
        with hist_slot.container():
            st.plotly_chart(
                build_delta_histogram(
                    combined_deltas,
                    title=f"ESM2 score − wildtype ({len(combined_deltas)} chimera(s))",
                ),
                use_container_width=True,
                key=f"ai_sim_hist_{len(combined_deltas)}",
            )
        with table_slot.container():
            st.markdown("### Top sequences (live)")
            st.dataframe(
                _results_to_display_dataframe(partial_results),
                use_container_width=True,
                hide_index=True,
            )

    status.markdown("Scoring wildtype (all-query chimera)…")
    if prior and prior.get("wildtype_score") is not None:
        wildtype_score = float(prior["wildtype_score"])
        wildtype_sequence = prior.get("wildtype_sequence") or ""
    else:
        wildtype_score, wildtype_chimera = score_wildtype(
            blocks,
            esm_batch_size=esm_batch_size,
        )
        wildtype_sequence = wildtype_chimera.get("sequence") or ""

    progress.progress(0.0)
    status.markdown(f"Sampling {n_samples} random chimera(s)…")

    with st.spinner("Running simulation…"):
        new_records, _raw, new_deltas = simulate_random_chimeras_esm(
            blocks,
            n_samples,
            wildtype_score,
            start_index=start_index,
            esm_batch_size=esm_batch_size,
            on_batch=on_batch,
        )

    merged = merge_sim_run(
        prior,
        pools_signature=sig,
        wildtype_score=wildtype_score,
        wildtype_sequence=wildtype_sequence,
        new_records=new_records,
    )
    st.session_state.ai_sim_results = merged
    progress.progress(1.0)
    status.markdown("Done.")

    if continuing:
        st.success(
            f"Added {len(new_deltas)} chimera(s) "
            f"({merged['n_chimeras']} total in histogram)."
        )
    else:
        st.success(
            f"Scored {merged['n_chimeras']} chimera(s). "
            f"Wildtype ESM2 = {wildtype_score:.4f}."
        )


_init_sim_state()

st.title("🤖 Simulate with AI")

st.markdown(
    "Sample random chimeras from your **saved fragment pools** "
    "(**Diversity Analysis** or **Library Optimization** → save list to session), "
    "score each with **ESM2-150M**, "
    "and plot **Δ score = chimera − wildtype**, where wildtype is the chimera built "
    "from the query sequence at every block. Run additional batches to grow the histogram."
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

sim_results = _sync_sim_results(blocks)
has_results = bool(sim_results.get("score_deltas"))

st.markdown("---")
st.subheader("Run simulation" if not has_results else "Continue simulation")

col1, col2 = st.columns(2)
with col1:
    n_samples = st.number_input(
        "Sequences to generate",
        min_value=1,
        max_value=50000,
        value=200,
        step=1,
        help="Random chimeras to sample and ESM-score in this run.",
        key="ai_sim_n_samples",
    )
with col2:
    esm_batch_size = st.number_input(
        "ESM batch size",
        min_value=1,
        max_value=32,
        value=8,
        help="Sequences scored per GPU/CPU batch (histogram updates after each batch).",
        key="ai_sim_esm_batch_size",
    )

button_label = "Run more sequences" if has_results else "Run simulation"
if st.button(button_label, type="primary", key="ai_sim_run"):
    _run_simulation(
        blocks,
        n_samples=int(n_samples),
        esm_batch_size=int(esm_batch_size),
        prior=sim_results if has_results else None,
    )
    st.rerun()

if has_results:
    st.markdown("---")
    st.subheader("Results")

    deltas = sim_results.get("score_deltas") or []
    wt_score = sim_results.get("wildtype_score")
    metric_cols = st.columns(4)
    metric_cols[0].metric("Chimeras scored", len(deltas))
    metric_cols[1].metric(
        "Wildtype ESM2",
        f"{wt_score:.4f}" if wt_score is not None else "—",
    )
    if deltas:
        metric_cols[2].metric("Mean Δ", f"{statistics.mean(deltas):.4f}")
        metric_cols[3].metric("Max Δ", f"{max(deltas):.4f}")
    else:
        metric_cols[2].metric("Mean Δ", "—")
        metric_cols[3].metric("Max Δ", "—")

    st.plotly_chart(
        build_delta_histogram(
            deltas,
            title=f"ESM2 score − wildtype ({len(deltas)} chimera(s))",
        ),
        use_container_width=True,
        key="ai_sim_hist_final",
    )

    st.caption(
        "Δ > 0: chimera scores higher than wildtype under ESM2; "
        "Δ < 0: lower than wildtype. Higher ESM2 (less negative log-likelihood) is better."
    )

    st.markdown("### Top sequences")
    st.caption("Wildtype first, then the 20 highest Δ vs wildtype across all scored chimeras.")
    st.dataframe(
        _results_to_display_dataframe(sim_results),
        use_container_width=True,
        hide_index=True,
    )

    fasta_body = top_sequences_to_fasta(sim_results, top_n=20)
    st.download_button(
        "Download top sequences (FASTA)",
        data=fasta_body,
        file_name="ai_sim_top_sequences.fasta",
        mime="text/plain",
        key="ai_sim_download_fasta",
        disabled=not fasta_body.strip(),
    )

    if st.button("Clear simulation results", key="ai_sim_clear"):
        st.session_state.ai_sim_results = None
        st.rerun()
else:
    st.info("Configure parameters above and click **Run simulation**.")
