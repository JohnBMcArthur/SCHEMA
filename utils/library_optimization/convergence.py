"""
Convergence metrics for iterative library optimization.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Tuple

import plotly.graph_objects as go
from plotly.colors import qualitative
from scipy.stats import rankdata, spearmanr


def coefficient_ranks_by_block(
    coefficients_by_block: Dict[str, List[Dict[str, Any]]],
) -> Dict[str, Dict[str, int]]:
    """
    Per-fragment variant ranks from ridge coefficients (1 = highest coefficient).
    """
    out: Dict[str, Dict[str, int]] = {}
    for frag_key, rows in coefficients_by_block.items():
        if not rows:
            out[frag_key] = {}
            continue
        coefs = [float(r.get("coefficient", 0.0)) for r in rows]
        ranks = rankdata([-c for c in coefs], method="min")
        out[frag_key] = {
            rows[i]["row_id"]: int(ranks[i]) for i in range(len(rows))
        }
    return out


def create_rank_comparison_figure(
    prev_coefficients: Dict[str, List[Dict[str, Any]]],
    curr_coefficients: Dict[str, List[Dict[str, Any]]],
    *,
    prev_round: int,
    curr_round: int,
) -> go.Figure:
    """Scatter of within-fragment ranks: previous round vs current round."""
    prev_ranks = coefficient_ranks_by_block(prev_coefficients)
    curr_ranks = coefficient_ranks_by_block(curr_coefficients)
    palette = qualitative.Plotly + qualitative.D3

    fig = go.Figure()
    frag_keys = sorted(prev_ranks.keys(), key=lambda k: int(k))
    max_rank = 1

    for idx, frag_key in enumerate(frag_keys):
        xs: List[int] = []
        ys: List[int] = []
        labels: List[str] = []
        for row_id, rank_prev in prev_ranks.get(frag_key, {}).items():
            rank_curr = curr_ranks.get(frag_key, {}).get(row_id)
            if rank_curr is None:
                continue
            xs.append(rank_prev)
            ys.append(rank_curr)
            seq_id = next(
                (
                    str(r.get("sequence_id", row_id))
                    for r in curr_coefficients.get(frag_key, [])
                    if r.get("row_id") == row_id
                ),
                row_id,
            )
            labels.append(seq_id)
            max_rank = max(max_rank, rank_prev, rank_curr)

        if not xs:
            continue

        fig.add_trace(
            go.Scatter(
                x=xs,
                y=ys,
                mode="markers",
                name=f"Fragment {frag_key}",
                marker=dict(size=9, color=palette[idx % len(palette)]),
                text=labels,
                hovertemplate=(
                    "Fragment %{fullData.name}<br>"
                    "Previous rank: %{x}<br>"
                    "New rank: %{y}<br>"
                    "%{text}<extra></extra>"
                ),
            )
        )

    fig.add_trace(
        go.Scatter(
            x=[1, max_rank],
            y=[1, max_rank],
            mode="lines",
            line=dict(color="rgba(120,120,120,0.6)", dash="dash", width=1),
            showlegend=False,
            hoverinfo="skip",
        )
    )

    pad = 0.5
    axis_lo = 1 - pad
    axis_hi = max_rank + pad

    fig.update_layout(
        title=f"Round {curr_round}: rank previous (round {prev_round}) vs rank new",
        xaxis_title=f"Rank after round {prev_round}",
        yaxis_title=f"Rank after round {curr_round}",
        height=720,
        legend=dict(
            title="Fragment",
            yanchor="top",
            y=1.0,
            xanchor="left",
            x=1.02,
        ),
        margin=dict(r=180, t=60, b=60),
    )
    fig.update_xaxes(
        range=[axis_lo, axis_hi],
        constrain="domain",
        automargin=True,
    )
    fig.update_yaxes(
        range=[axis_lo, axis_hi],
        scaleanchor="x",
        scaleratio=1,
        constrain="domain",
        automargin=True,
    )
    if max_rank <= 30:
        fig.update_xaxes(dtick=1)
        fig.update_yaxes(dtick=1)
    return fig


def rank_comparison_figures_from_history(
    coefficient_history: List[Tuple[int, Dict[str, List[Dict[str, Any]]]]],
) -> List[Tuple[int, int, go.Figure]]:
    """Build one rank-comparison figure per round transition (round 2 onward)."""
    figures: List[Tuple[int, int, go.Figure]] = []
    for i in range(1, len(coefficient_history)):
        prev_round, prev_coef = coefficient_history[i - 1]
        curr_round, curr_coef = coefficient_history[i]
        figures.append(
            (
                prev_round,
                curr_round,
                create_rank_comparison_figure(
                    prev_coef,
                    curr_coef,
                    prev_round=prev_round,
                    curr_round=curr_round,
                ),
            )
        )
    return figures


def spearman_rank_correlation_by_block(
    current: Dict[str, List[Dict[str, Any]]],
    previous: Optional[Dict[str, List[Dict[str, Any]]]],
) -> Dict[str, Optional[float]]:
    """
    Spearman ρ between variant coefficient ranks in *current* vs *previous*.

    Returns None for each block on the first comparison round (no previous).
    """
    if previous is None:
        return {frag_key: None for frag_key in current}

    out: Dict[str, Optional[float]] = {}
    for frag_key, curr_rows in current.items():
        prev_rows = previous.get(frag_key)
        if not prev_rows:
            out[frag_key] = None
            continue

        prev_coef = {r["row_id"]: float(r["coefficient"]) for r in prev_rows}
        curr_ids = [r["row_id"] for r in curr_rows]
        prev_vals = [prev_coef.get(rid, 0.0) for rid in curr_ids]
        curr_vals = [float(r["coefficient"]) for r in curr_rows]

        if len(set(curr_vals)) < 2 and len(set(prev_vals)) < 2:
            out[frag_key] = 1.0
            continue

        rho, _ = spearmanr(prev_vals, curr_vals)
        out[frag_key] = float(rho) if rho == rho else None  # NaN check

    return out


def all_blocks_converged(
    spearman_by_block: Dict[str, Optional[float]],
    threshold: float,
) -> bool:
    """True when every block has Spearman ρ >= threshold."""
    values = [v for v in spearman_by_block.values() if v is not None]
    if not values:
        return False
    return all(v >= threshold for v in values)
