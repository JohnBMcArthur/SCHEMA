"""
Sequence Similarity Network (SSN) layout and Plotly rendering for Diversity Analysis.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Set, Tuple

import networkx as nx
import plotly.graph_objects as go
from Bio import pairwise2


def aligned_hamming_distance(seq_a: str, seq_b: str) -> int:
    """Global-align two sequences and count column mismatches (gaps count)."""
    if not seq_a or not seq_b:
        return max(len(seq_a), len(seq_b))
    alignments = pairwise2.align.globalxx(seq_a, seq_b)
    if not alignments:
        return max(len(seq_a), len(seq_b))
    aligned_a, aligned_b = alignments[0][0], alignments[0][1]
    return sum(1 for a, b in zip(aligned_a, aligned_b) if a != b)


def compute_distance_matrix(nodes: List[Dict[str, Any]]) -> List[List[int]]:
    """All-by-all aligned Hamming distances for the given node list."""
    n = len(nodes)
    distances: List[List[int]] = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            dist = aligned_hamming_distance(nodes[i]["sequence"], nodes[j]["sequence"])
            distances[i][j] = dist
            distances[j][i] = dist
    return distances


def max_from_distance_matrix(distances: List[List[int]]) -> int:
    """Maximum off-diagonal value in a symmetric distance matrix."""
    if not distances:
        return 0
    max_dist = 0
    n = len(distances)
    for i in range(n):
        for j in range(i + 1, n):
            max_dist = max(max_dist, distances[i][j])
    return max_dist


def estimate_hamming_slider_max(nodes: List[Dict[str, Any]]) -> int:
    """
    Cheap upper bound for the threshold slider before distances are computed.

    Aligned Hamming is at most the sum of the two sequence lengths.
    """
    if len(nodes) <= 1:
        return 1
    max_len = max(len(node.get("sequence") or "") for node in nodes)
    return max(1, max_len * 2)


def build_ssn_graph_from_distances(
    nodes: List[Dict[str, Any]],
    distances: List[List[int]],
    threshold: int,
) -> nx.Graph:
    """Build graph using a precomputed distance matrix."""
    graph = nx.Graph()
    for index, node in enumerate(nodes):
        graph.add_node(index, **node)

    n = len(nodes)
    for i in range(n):
        for j in range(i + 1, n):
            if distances[i][j] <= threshold:
                graph.add_edge(i, j, weight=distances[i][j])
    return graph


def build_ssn_graph(
    nodes: List[Dict[str, Any]],
    threshold: int,
) -> Tuple[nx.Graph, List[List[int]]]:
    """
    Build an undirected graph: edge when aligned Hamming distance <= threshold.

    Returns the NetworkX graph and a symmetric distance matrix (n x n).
    """
    distances = compute_distance_matrix(nodes)
    graph = build_ssn_graph_from_distances(nodes, distances, threshold)
    return graph, distances


def _spring_positions(graph: nx.Graph) -> Dict[int, Tuple[float, float]]:
    if graph.number_of_nodes() == 0:
        return {}
    if graph.number_of_nodes() == 1:
        return {0: (0.0, 0.0)}
    if graph.number_of_edges() == 0:
        nodes = sorted(graph.nodes())
        return {node: (float(i), 0.0) for i, node in enumerate(nodes)}
    return nx.spring_layout(
        graph,
        seed=42,
        k=1.0 / max(graph.number_of_nodes() ** 0.5, 1),
        iterations=75,
    )


def create_ssn_figure_from_distances(
    nodes: List[Dict[str, Any]],
    distances: List[List[int]],
    threshold: int,
    preview_row_ids: Set[str],
    query_row_id: str = "__query__",
) -> Tuple[go.Figure, nx.Graph]:
    """Render SSN using a precomputed distance matrix (no alignment calls)."""
    graph = build_ssn_graph_from_distances(nodes, distances, threshold)
    positions = _spring_positions(graph)

    edge_x: List[Optional[float]] = []
    edge_y: List[Optional[float]] = []
    for u, v in graph.edges():
        x0, y0 = positions[u]
        x1, y1 = positions[v]
        edge_x.extend([x0, x1, None])
        edge_y.extend([y0, y1, None])

    edge_trace = go.Scatter(
        x=edge_x,
        y=edge_y,
        mode="lines",
        line=dict(width=0.8, color="#888"),
        hoverinfo="none",
        showlegend=False,
    )

    def _node_color(row_id: str) -> str:
        if row_id == query_row_id:
            return "#2563eb"
        if row_id in preview_row_ids:
            return "#dc2626"
        return "#94a3b8"

    node_x = [positions[i][0] for i in range(len(nodes))]
    node_y = [positions[i][1] for i in range(len(nodes))]
    colors = [_node_color(n["row_id"]) for n in nodes]
    hover_text = [
        (
            f"<b>{n['sequence_id']}</b><br>"
            f"Length: {n.get('length_display', n.get('length', ''))}<br>"
            f"Mutations (non-gap): {n.get('mutations_non_gap', '—')}<br>"
            f"Mutations (w/ gaps): {n.get('mutations_with_gaps', '—')}<br>"
            f"% identity: {n.get('pct_identity', '—'):.1f}%"
            if isinstance(n.get("pct_identity"), (int, float))
            else f"<b>{n['sequence_id']}</b>"
        )
        for n in nodes
    ]

    node_trace = go.Scatter(
        x=node_x,
        y=node_y,
        mode="markers",
        marker=dict(
            size=14,
            color=colors,
            line=dict(width=1.5, color="#1e293b"),
        ),
        text=hover_text,
        hoverinfo="text",
        customdata=[n["row_id"] for n in nodes],
        showlegend=False,
    )

    fig = go.Figure(data=[edge_trace, node_trace])
    fig.update_layout(
        title=f"Sequence Similarity Network (Hamming ≤ {threshold})",
        showlegend=False,
        hovermode="closest",
        margin=dict(l=10, r=10, t=40, b=10),
        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        height=520,
        clickmode="event+select",
    )
    return fig, graph


def create_ssn_figure(
    nodes: List[Dict[str, Any]],
    threshold: int,
    preview_row_ids: Set[str],
    query_row_id: str = "__query__",
) -> Tuple[go.Figure, nx.Graph]:
    """Render SSN with force-directed layout; computes distances on the fly."""
    distances = compute_distance_matrix(nodes)
    return create_ssn_figure_from_distances(
        nodes, distances, threshold, preview_row_ids, query_row_id
    )
