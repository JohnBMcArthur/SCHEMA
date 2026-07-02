"""
HTML highlighting for Golden Gate / BsaI DNA sequences in Oligopool tables.
"""

from __future__ import annotations

import html
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

from utils.oligopool_design import BSAI_END, BSAI_START, merge_pieces_at_overhang

GGA_HIGHLIGHT_CSS = """
<style>
.oligopool-highlight-legend span {
    font-family: ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, monospace;
    padding: 1px 4px;
    border-radius: 3px;
    margin-right: 10px;
}
.gga-bsa {
    background-color: #b8d4f0;
    padding: 0 1px;
    border-radius: 2px;
    font-family: ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, monospace;
}
.gga-overhang {
    background-color: #b8e6b8;
    padding: 0 1px;
    border-radius: 2px;
    font-family: ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, monospace;
}
.gga-stuffer {
    background-color: #e4e4e4;
    padding: 0 1px;
    border-radius: 2px;
    font-family: ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, monospace;
}
.gga-primer {
    background-color: #ffd8a8;
    padding: 0 1px;
    border-radius: 2px;
    font-family: ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, monospace;
}
.oligopool-seq-table {
    border-collapse: collapse;
    width: 100%;
    font-size: 0.88rem;
    margin-bottom: 0.5rem;
}
.oligopool-seq-table th {
    background-color: #f5f5f5;
    text-align: left;
    font-weight: 600;
    border: 1px solid #ddd;
    padding: 8px;
    white-space: nowrap;
}
.oligopool-seq-table td {
    border: 1px solid #ddd;
    padding: 8px;
    vertical-align: top;
    word-break: break-word;
}
</style>
"""


def render_highlight_legend() -> str:
    return (
        GGA_HIGHLIGHT_CSS
        + '<div class="oligopool-highlight-legend">'
        '<span class="gga-bsa">BsaI Golden Gate site</span>'
        '<span class="gga-overhang">4 bp overhang</span>'
        "</div>"
    )


def render_full_oligo_legend() -> str:
    return (
        GGA_HIGHLIGHT_CSS
        + '<div class="oligopool-highlight-legend">'
        '<span class="gga-stuffer">Stuffer</span>'
        '<span class="gga-primer">Primer</span>'
        '<span class="gga-bsa">BsaI Golden Gate site</span>'
        '<span class="gga-overhang">4 bp overhang</span>'
        "</div>"
    )


def _unique_four_bp(overhangs: Iterable[str]) -> List[str]:
    seen: set[str] = set()
    out: List[str] = []
    for oh in overhangs:
        key = (oh or "").upper()
        if len(key) == 4 and key not in seen:
            seen.add(key)
            out.append(key)
    return out


def _unique_overhangs(overhangs: Iterable[str]) -> List[str]:
    return _unique_four_bp(overhangs)


def _bsa_site_patterns() -> List[str]:
    patterns = [BSAI_START, BSAI_END, "GGTCTC", "GAGACC"]
    seen: set[str] = set()
    out: List[str] = []
    for pat in sorted(patterns, key=len, reverse=True):
        key = pat.upper()
        if key not in seen:
            seen.add(key)
            out.append(key)
    return out


def _coding_region_in_piece(piece_dna: str) -> Tuple[int, int]:
    """Return (start, end) half-open indices of the coding interior in a BsaI piece."""
    seq = (piece_dna or "").upper()
    start = len(BSAI_START) if seq.startswith(BSAI_START) else 0
    end = len(seq) - len(BSAI_END) if seq.endswith(BSAI_END) else len(seq)
    return start, max(start, end)


def _dedupe_spans(spans: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    seen: set[Tuple[int, int]] = set()
    out: List[Tuple[int, int]] = []
    for span in spans:
        if span not in seen:
            seen.add(span)
            out.append(span)
    return out


def _find_in_region(seq: str, pattern: str, region_start: int, region_end: int) -> Optional[int]:
    """First match of *pattern* within seq[region_start:region_end], or None."""
    pattern = pattern.upper()
    if len(pattern) != 4:
        return None
    region = seq[region_start:region_end]
    pos = region.find(pattern)
    if pos == -1:
        return None
    return region_start + pos


def overhangs_for_fragment(
    frag_index: int,
    enriched_fragments: Sequence[Dict[str, Any]],
) -> List[str]:
    """4 bp overhangs assigned to one fragment's junctions."""
    idx = frag_index - 1
    if idx < 0 or idx >= len(enriched_fragments):
        return []
    frag = enriched_fragments[idx]
    candidates: List[str] = []
    oh5 = frag.get("overhang_5") or {}
    oh3 = frag.get("overhang_3") or {}
    if oh5.get("overhang"):
        candidates.append(str(oh5["overhang"]))
    if oh3.get("overhang"):
        candidates.append(str(oh3["overhang"]))
    return _unique_four_bp(candidates)


def overhangs_for_fragment_range(
    fragment_indices: Sequence[int],
    enriched_fragments: Sequence[Dict[str, Any]],
) -> List[str]:
    """All junction overhangs spanned by a list of fragment indices."""
    candidates: List[str] = []
    for frag_index in fragment_indices:
        candidates.extend(overhangs_for_fragment(frag_index, enriched_fragments))
    return _unique_four_bp(candidates)


def overhang_spans_in_piece(
    piece_dna: str,
    frag_index: int,
    enriched_fragments: Sequence[Dict[str, Any]],
) -> List[Tuple[int, int]]:
    """
    4 bp junction overhangs at the 5′/3′ coding termini of one BsaI piece.

    Only the terminal ~6 bp of coding DNA are searched so interior matches
    (e.g. AACC within a codon) are not highlighted.
    """
    seq = (piece_dna or "").upper()
    if not seq:
        return []

    coding_start, coding_end = _coding_region_in_piece(seq)
    if coding_end <= coding_start:
        return []

    idx = frag_index - 1
    if idx < 0 or idx >= len(enriched_fragments):
        return []

    enriched = enriched_fragments[idx]
    spans: List[Tuple[int, int]] = []
    terminal_window = 6

    oh5 = enriched.get("overhang_5") or {}
    if oh5.get("overhang"):
        pos = _find_in_region(
            seq,
            str(oh5["overhang"]),
            coding_start,
            coding_start + min(terminal_window, coding_end - coding_start),
        )
        if pos is not None:
            spans.append((pos, pos + 4))

    oh3 = enriched.get("overhang_3") or {}
    if oh3.get("overhang"):
        tail_start = max(coding_start, coding_end - terminal_window)
        pos = _find_in_region(seq, str(oh3["overhang"]), tail_start, coding_end)
        if pos is not None:
            spans.append((pos, pos + 4))

    return _dedupe_spans(spans)


def overhang_spans_in_coding(
    coding_dna: str,
    frag_index: int,
    enriched_fragments: Sequence[Dict[str, Any]],
) -> List[Tuple[int, int]]:
    """Junction overhang spans within coding DNA only (no BsaI flanks)."""
    coding = (coding_dna or "").upper()
    if not coding:
        return []
    wrapped = f"{BSAI_START}{coding}{BSAI_END}"
    offset = len(BSAI_START)
    return [
        (start - offset, end - offset)
        for start, end in overhang_spans_in_piece(wrapped, frag_index, enriched_fragments)
        if start >= offset and end <= offset + len(coding)
    ]


def shift_spans(
    spans: Sequence[Tuple[int, int]],
    offset: int,
) -> List[Tuple[int, int]]:
    return [(start + offset, end + offset) for start, end in spans]


def overhang_spans_for_block_packed_oligo(
    oligo: Dict[str, Any],
    enriched_fragments: Sequence[Dict[str, Any]],
    *,
    insert_offset: int = 0,
) -> List[Tuple[int, int]]:
    """Highlight spans for concatenated block inserts inside one packed oligo."""
    spans: List[Tuple[int, int]] = []
    offset = 0
    for block in oligo.get("blocks") or []:
        piece_dna = block.get("piece_dna") or ""
        frag_idx = int(block.get("fragment_index") or 0)
        for start, end in overhang_spans_in_piece(
            piece_dna, frag_idx, enriched_fragments
        ):
            spans.append((insert_offset + offset + start, insert_offset + offset + end))
        offset += len(piece_dna)
    return spans


def overhang_spans_in_merged_insert(
    fragment_indices: Sequence[int],
    piece_by_fragment: Dict[int, str],
    enriched_fragments: Sequence[Dict[str, Any]],
    junction_overhangs: Sequence[str],
) -> List[Tuple[int, int]]:
    """
    Overhang positions inside a merged Stage 1 insert (multiple fragments).

    Simulates piece merging to locate internal junction 4 bp sites, plus
    terminal overhangs on the first (5′) and last (3′) fragment.
    """
    if not fragment_indices:
        return []

    pieces = [piece_by_fragment[i].upper() for i in fragment_indices if i in piece_by_fragment]
    if not pieces:
        return []

    spans: List[Tuple[int, int]] = []
    spans.extend(
        overhang_spans_in_piece(pieces[0], fragment_indices[0], enriched_fragments)
    )

    merged = pieces[0]
    for i in range(1, len(pieces)):
        left_frag = fragment_indices[i - 1]
        oh = ""
        junc_idx = left_frag - 1
        if 0 <= junc_idx < len(junction_overhangs):
            oh = (junction_overhangs[junc_idx] or "").upper()
        if len(oh) == 4 and merged.upper().endswith(oh):
            spans.append((len(merged) - 4, len(merged)))
        merged = merge_pieces_at_overhang(merged, pieces[i], oh)

    if len(fragment_indices) > 1:
        last_frag = fragment_indices[-1]
        idx = last_frag - 1
        enriched = enriched_fragments[idx]
        oh3 = enriched.get("overhang_3") or {}
        if oh3.get("overhang"):
            coding_start, coding_end = _coding_region_in_piece(merged)
            tail_start = max(coding_start, coding_end - 6)
            pos = _find_in_region(merged, str(oh3["overhang"]), tail_start, coding_end)
            if pos is not None:
                spans.append((pos, pos + 4))

    return _dedupe_spans(spans)


def _bsa_spans(seq: str) -> List[Tuple[int, int]]:
    spans: List[Tuple[int, int]] = []
    for site in _bsa_site_patterns():
        start = 0
        while start <= len(seq) - len(site):
            idx = seq.find(site, start)
            if idx == -1:
                break
            spans.append((idx, idx + len(site)))
            start = idx + len(site)
    return spans


def _shift_spans_list(
    spans: Sequence[Tuple[int, int]],
    offset: int,
) -> List[Tuple[int, int]]:
    return [(start + offset, end + offset) for start, end in spans]


def region_spans_for_full_oligo(
    oligo: Dict[str, Any],
    enriched_fragments: Sequence[Dict[str, Any]],
) -> Dict[str, List[Tuple[int, int]]]:
    """
    Classify oligo sequence regions for multi-color highlighting.

    Layout: [5′ stuffer][forward primer][insert][RC(reverse primer)][3′ stuffer]
    """
    seq = (oligo.get("sequence") or "").upper()
    s5 = len(oligo.get("stuffer_5") or "")
    s3 = len(oligo.get("stuffer_3") or "")
    fwd = oligo.get("forward_primer") or ""
    rev_rc = oligo.get("reverse_primer_rc") or ""
    fwd_len = len(fwd)
    rev_len = len(rev_rc)
    insert_len = int(oligo.get("insert_length") or len(oligo.get("merged_insert") or ""))
    insert_start = s5 + fwd_len
    insert_end = insert_start + insert_len

    stuffer_spans: List[Tuple[int, int]] = []
    if s5 > 0:
        stuffer_spans.append((0, s5))
    if s3 > 0:
        stuffer_spans.append((len(seq) - s3, len(seq)))

    primer_spans: List[Tuple[int, int]] = []
    if fwd_len > 0:
        primer_spans.append((s5, s5 + fwd_len))
    if rev_len > 0:
        primer_spans.append((insert_end, insert_end + rev_len))

    insert_seq = seq[insert_start:insert_end]
    bsa_spans = _shift_spans_list(_bsa_spans(insert_seq), insert_start)
    overhang_spans = overhang_spans_for_block_packed_oligo(
        oligo,
        enriched_fragments,
        insert_offset=insert_start,
    )

    return {
        "stuffer": stuffer_spans,
        "primer": primer_spans,
        "bsa": bsa_spans,
        "overhang": overhang_spans,
    }


def highlight_oligo_sequence(
    sequence: str,
    regions: Dict[str, Sequence[Tuple[int, int]]],
) -> str:
    """
    Highlight full oligo with stuffer, primer, BsaI, and overhang regions.

    Later region types take precedence where spans overlap.
    """
    seq = (sequence or "").upper()
    if not seq:
        return ""

    n = len(seq)
    priority = ("stuffer", "primer", "bsa", "overhang")
    css_class = {
        "stuffer": "gga-stuffer",
        "primer": "gga-primer",
        "bsa": "gga-bsa",
        "overhang": "gga-overhang",
    }
    mark: List[Optional[str]] = [None] * n

    for kind in priority:
        for start, end in regions.get(kind) or []:
            for i in range(max(0, start), min(n, end)):
                mark[i] = kind

    parts: List[str] = []
    i = 0
    while i < n:
        tag = mark[i]
        j = i + 1
        while j < n and mark[j] == tag:
            j += 1
        chunk = html.escape(seq[i:j])
        if tag and tag in css_class:
            parts.append(f'<span class="{css_class[tag]}">{chunk}</span>')
        else:
            parts.append(chunk)
        i = j
    return "".join(parts)


def highlight_dna_sequence(
    sequence: str,
    overhang_spans: Optional[Sequence[Tuple[int, int]]] = None,
    *,
    highlight_bsaI: bool = True,
) -> str:
    """
    Return HTML with BsaI sites and 4 bp overhangs highlighted.

    Pass explicit *overhang_spans* (start, end half-open) to avoid false positives.
    Overhang highlights take precedence over BsaI where they overlap.
    """
    seq = (sequence or "").upper()
    if not seq:
        return ""

    n = len(seq)
    mark: List[Optional[str]] = [None] * n

    if highlight_bsaI:
        for start, end in _bsa_spans(seq):
            for i in range(start, min(end, n)):
                if mark[i] is None:
                    mark[i] = "bsa"

    for start, end in overhang_spans or []:
        for i in range(max(0, start), min(n, end)):
            mark[i] = "oh"

    parts: List[str] = []
    i = 0
    while i < n:
        tag = mark[i]
        j = i + 1
        while j < n and mark[j] == tag:
            j += 1
        chunk = html.escape(seq[i:j])
        if tag == "oh":
            parts.append(f'<span class="gga-overhang">{chunk}</span>')
        elif tag == "bsa":
            parts.append(f'<span class="gga-bsa">{chunk}</span>')
        else:
            parts.append(chunk)
        i = j
    return "".join(parts)


def render_html_table(
    rows: Sequence[Dict[str, Any]],
    columns: Sequence[Tuple[str, str]],
    *,
    html_columns: Optional[Dict[str, List[str]]] = None,
) -> str:
    """Render a simple HTML table with optional pre-rendered HTML columns."""
    html_columns = html_columns or {}
    head = "".join(f"<th>{html.escape(header)}</th>" for header, _ in columns)
    body_rows: List[str] = []
    for row_idx, row in enumerate(rows):
        cells: List[str] = []
        for _, key in columns:
            if key in html_columns:
                cells.append(f"<td>{html_columns[key][row_idx]}</td>")
            else:
                cells.append(f"<td>{html.escape(str(row.get(key, '')))}</td>")
        body_rows.append(f"<tr>{''.join(cells)}</tr>")

    return (
        GGA_HIGHLIGHT_CSS
        + '<table class="oligopool-seq-table"><thead><tr>'
        + head
        + "</tr></thead><tbody>"
        + "".join(body_rows)
        + "</tbody></table>"
    )
