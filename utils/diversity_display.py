"""
HTML rendering for aligned diversity fragment sequences.
"""

from __future__ import annotations

from typing import Tuple

_CELL_STYLE = (
    "display:inline-block;width:1ch;text-align:center;"
    "font-family:monospace;font-size:0.85rem;line-height:1.2;"
)
_WRAP_STYLE = "white-space:nowrap;overflow-x:auto;"


def normalize_alignment_chars(sequence: str) -> str:
    """Normalize gap characters and uppercase amino-acid letters for comparison."""
    out: list[str] = []
    for char in sequence or "":
        if char in ("-", "."):
            out.append("-")
        else:
            out.append(char.upper())
    return "".join(out)


def pad_alignment(sequence: str, target_length: int) -> str:
    """Pad an aligned string with trailing gaps to a fixed column width."""
    if target_length <= len(sequence):
        return sequence[:target_length]
    return sequence + ("-" * (target_length - len(sequence)))


def normalize_alignment_pair(aligned_query: str, aligned_target: str) -> Tuple[str, str]:
    """Normalize and equalize length for an aligned query/homolog pair."""
    aligned_query = normalize_alignment_chars(aligned_query)
    aligned_target = normalize_alignment_chars(aligned_target)
    width = max(len(aligned_query), len(aligned_target))
    return pad_alignment(aligned_query, width), pad_alignment(aligned_target, width)


def classify_alignment_column(query_char: str, target_char: str) -> str:
    """Return 'match', 'substitution', 'gap_mutation', or 'gap_match'."""
    q = normalize_alignment_chars(query_char)
    t = normalize_alignment_chars(target_char)
    q_gap = q == "-"
    t_gap = t == "-"
    if not q_gap and not t_gap:
        return "match" if q == t else "substitution"
    if q_gap != t_gap:
        return "gap_mutation"
    return "gap_match"


def _render_char(char: str, background: str | None = None) -> str:
    display = char if char not in ("-", ".") else "-"
    style = _CELL_STYLE
    if background:
        style += f"background-color:{background};"
    return f'<span style="{style}">{display}</span>'


def render_aligned_row_html(
    aligned_query: str,
    aligned_target: str,
    *,
    show_query: bool,
) -> str:
    """Render one aligned row with fixed-width columns."""
    aligned_query, aligned_target = normalize_alignment_pair(aligned_query, aligned_target)
    parts = [f'<div style="{_WRAP_STYLE}">']
    for q, t in zip(aligned_query, aligned_target):
        char = q if show_query else t
        if show_query:
            parts.append(_render_char(char))
            continue
        kind = classify_alignment_column(q, t)
        if kind == "substitution":
            parts.append(_render_char(char, "#fde047"))
        elif kind == "gap_mutation":
            parts.append(_render_char(char, "#f0abfc"))
        else:
            parts.append(_render_char(char))
    parts.append("</div>")
    return "".join(parts)


def render_msa_pair_html(aligned_query: str, aligned_target: str) -> str:
    """Render paired query + homolog rows that share the same alignment columns."""
    return (
        '<div style="line-height:1.35;">'
        '<div><span style="color:#64748b;font-size:0.75rem;font-family:monospace;">Q </span>'
        f"{render_aligned_row_html(aligned_query, aligned_target, show_query=True)}"
        "</div>"
        '<div><span style="color:#64748b;font-size:0.75rem;font-family:monospace;">H </span>'
        f"{render_aligned_row_html(aligned_query, aligned_target, show_query=False)}"
        "</div>"
        "</div>"
    )


def render_query_sequence_html(aligned_query: str, aligned_target: str | None = None) -> str:
    """Single-line query row using the homolog's alignment frame."""
    target = aligned_target if aligned_target is not None else aligned_query
    return render_aligned_row_html(aligned_query, target, show_query=True)


def render_homolog_sequence_html(aligned_query: str, aligned_target: str) -> str:
    """Single-line homolog row paired to the same query alignment frame."""
    return render_aligned_row_html(aligned_query, aligned_target, show_query=False)


def render_aligned_sequence_html(
    aligned_query: str,
    aligned_target: str,
    *,
    show_query: bool = True,
    show_target: bool = True,
) -> str:
    """Render query and/or homolog rows with mutation highlighting."""
    parts = ['<div style="font-family:monospace;font-size:0.85rem;line-height:1.6;">']
    if show_query:
        parts.append('<div><strong>Query:</strong>&nbsp;')
        parts.append(render_aligned_row_html(aligned_query, aligned_target, show_query=True))
        parts.append("</div>")
    if show_target:
        parts.append('<div><strong>Homolog:</strong>&nbsp;')
        parts.append(render_aligned_row_html(aligned_query, aligned_target, show_query=False))
        parts.append("</div>")
    parts.append("</div>")
    return "".join(parts)


def metrics_from_aligned_pair(aligned_query: str, aligned_target: str) -> Tuple[int, int, float]:
    """Return (non-gap mutations, gap-inclusive mutations, pct identity)."""
    aligned_query, aligned_target = normalize_alignment_pair(aligned_query, aligned_target)
    mutations_non_gap = 0
    mutations_with_gaps = 0
    paired = 0
    identical = 0
    for q, t in zip(aligned_query, aligned_target):
        if q != t:
            mutations_with_gaps += 1
        q_gap = q == "-"
        t_gap = t == "-"
        if not q_gap and not t_gap:
            paired += 1
            if q == t:
                identical += 1
            else:
                mutations_non_gap += 1
    pct_identity = (100.0 * identical / paired) if paired else 0.0
    return mutations_non_gap, mutations_with_gaps, pct_identity


def metrics_for_row(row: dict) -> Tuple[int, int, float]:
    """Recompute metrics from stored aligned fragments on a table row."""
    aligned_q = row.get("aligned_query_fragment") or row.get("sequence") or ""
    aligned_h = row.get("aligned_homolog_fragment") or row.get("sequence") or ""
    return metrics_from_aligned_pair(aligned_q, aligned_h)
