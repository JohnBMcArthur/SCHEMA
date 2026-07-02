"""
Diversity Analysis: homolog fragments from uploaded MSA (query row required) + GGA filtering.
"""

from __future__ import annotations

from io import StringIO
from typing import Any, Dict, List, Optional, Set, Tuple

from Bio import SeqIO, pairwise2

from utils.assembly_analysis import strip_msa_gaps
from utils.diversity_display import metrics_for_row, metrics_from_aligned_pair
from utils.gga_overhangs import fragment_passes_terminal_filter


def parse_fasta_msa(fasta_text: str) -> List[Tuple[str, str]]:
    """Parse FASTA content into (id, aligned_sequence) rows."""
    records = list(SeqIO.parse(StringIO(fasta_text), "fasta"))
    if not records:
        raise ValueError("No sequences found in FASTA file.")
    rows = [(rec.id, str(rec.seq)) for rec in records]
    lengths = {len(seq) for _, seq in rows}
    if len(lengths) != 1:
        raise ValueError(
            f"All sequences must have the same alignment length; found {sorted(lengths)}."
        )
    return rows


def _ungapped_normalized(sequence: str) -> str:
    return strip_msa_gaps(sequence).upper()


def find_query_row_in_msa(
    rows: List[Tuple[str, str]],
    query_ungapped: str,
) -> Tuple[str, str]:
    """
    Locate the query sequence row in an uploaded MSA.

    1. Prefer sequence IDs containing 'query' (case-insensitive) that match the session query.
    2. Otherwise scan all rows for an ungapped sequence match.
    3. Raise if the session query is not present in the alignment.
    """
    expected = _ungapped_normalized(query_ungapped)
    if not expected:
        raise ValueError("Session query sequence is empty.")

    query_id_candidates = [
        (seq_id, aligned)
        for seq_id, aligned in rows
        if "query" in seq_id.lower()
    ]
    for seq_id, aligned in query_id_candidates:
        if _ungapped_normalized(aligned) == expected:
            return seq_id, aligned

    for seq_id, aligned in rows:
        if _ungapped_normalized(aligned) == expected:
            return seq_id, aligned

    if query_id_candidates:
        ids = ", ".join(seq_id for seq_id, _ in query_id_candidates)
        raise ValueError(
            "Found sequence ID(s) containing 'query', but none match the session query "
            f"({len(expected)} residues ungapped). Checked: {ids}. "
            "Ensure the MSA query row matches the session query exactly, or rename/remove "
            "misleading IDs."
        )

    raise ValueError(
        "The session query sequence was not found in the uploaded alignment. "
        "Add the query as a row in the FASTA MSA (include 'query' in its sequence ID if helpful). "
        f"The session query has {len(expected)} residues (ungapped)."
    )


def _fragment_query_ranges(
    query_ungapped: str, fragment_sequences: List[str]
) -> List[Tuple[int, int, str, int]]:
    """Return (q_start, q_end, frag_seq, fragment_index) for each fragment."""
    query_ungapped = strip_msa_gaps(query_ungapped)
    ranges: List[Tuple[int, int, str, int]] = []
    search_at = 0
    for index, frag_seq in enumerate(fragment_sequences, start=1):
        if not frag_seq:
            raise ValueError(f"Fragment {index} has an empty sequence.")
        pos = query_ungapped.find(frag_seq, search_at)
        if pos < 0:
            raise ValueError(
                f"Fragment {index} sequence not found in the session query at the expected position."
            )
        end_pos = pos + len(frag_seq) - 1
        ranges.append((pos, end_pos, frag_seq, index))
        search_at = end_pos + 1
    return ranges


def _column_for_ungapped_position(aligned: str, position: int) -> int:
    """Map a 0-based ungapped residue index to its alignment column."""
    residue_index = 0
    for col, char in enumerate(aligned):
        if char in ("-", "."):
            continue
        if residue_index == position:
            return col
        residue_index += 1
    return -1


def msa_fragment_column_ranges(
    msa_query_aligned: str,
    query_ungapped: str,
    fragment_sequences: List[str],
) -> Dict[int, Tuple[int, int, str]]:
    """
    Map each assembly fragment to inclusive MSA column bounds on the query row.

    Returns fragment_index -> (col_start, col_end, ungapped_fragment_sequence).
    """
    ranges: Dict[int, Tuple[int, int, str]] = {}
    for q_start, q_end, frag_seq, frag_index in _fragment_query_ranges(
        query_ungapped, fragment_sequences
    ):
        col_start = _column_for_ungapped_position(msa_query_aligned, q_start)
        col_end = _column_for_ungapped_position(msa_query_aligned, q_end)
        if col_start < 0 or col_end < 0:
            raise ValueError(
                f"Fragment {frag_index} boundaries could not be mapped onto the MSA query row."
            )
        ranges[frag_index] = (col_start, col_end, frag_seq)
    return ranges


def extract_msa_fragment_at_columns(
    msa_query_aligned: str,
    homolog_aligned: str,
    col_start: int,
    col_end: int,
) -> Optional[Dict[str, Any]]:
    """Extract one fragment slice from shared MSA columns."""
    frag_q = msa_query_aligned[col_start : col_end + 1]
    frag_h = homolog_aligned[col_start : col_end + 1]
    sequence = strip_msa_gaps(frag_h)
    if not sequence:
        return None
    mutations_non_gap, mutations_with_gaps, pct_identity = metrics_from_aligned_pair(
        frag_q, frag_h
    )
    return {
        "sequence": sequence,
        "aligned_query_fragment": frag_q,
        "aligned_homolog_fragment": frag_h,
        "mutations_non_gap": mutations_non_gap,
        "mutations_with_gaps": mutations_with_gaps,
        "pct_identity": pct_identity,
    }


def _empty_fragment_buckets(
    fragment_filters: List[Dict[str, Any]],
) -> Dict[int, Dict[str, Any]]:
    buckets: Dict[int, Dict[str, Any]] = {}
    for filt in fragment_filters:
        idx = int(filt["fragment"])
        buckets[idx] = {
            "fragment": idx,
            "query_sequence": filt.get("query_sequence", ""),
            "msa_query_aligned_fragment": None,
            "overhang_5": filt.get("overhang_5"),
            "overhang_3": filt.get("overhang_3"),
            "allowed_n_terminal": filt.get("allowed_n_terminal"),
            "allowed_c_terminal": filt.get("allowed_c_terminal"),
            "entries": [],
            "rejected_unmapped": 0,
            "rejected_terminal_filter": 0,
        }
    return buckets


def analyze_diversity_msa(
    fasta_text: str,
    fragment_sequences: List[str],
    query_ungapped: str,
    fragment_filters: List[Dict[str, Any]],
) -> Dict[str, Any]:
    """
    Analyze homolog rows in an uploaded MSA using shared column structure with the query row.
    """
    rows = parse_fasta_msa(fasta_text)
    query_ungapped = strip_msa_gaps(query_ungapped)
    query_seq_id, msa_query_aligned = find_query_row_in_msa(rows, query_ungapped)

    if _ungapped_normalized(msa_query_aligned) != _ungapped_normalized(query_ungapped):
        raise ValueError("Internal error: MSA query row does not match session query.")

    fragment_columns = msa_fragment_column_ranges(
        msa_query_aligned, query_ungapped, fragment_sequences
    )
    filter_by_index = {int(f["fragment"]): f for f in fragment_filters}
    buckets = _empty_fragment_buckets(fragment_filters)

    for frag_index, (col_start, col_end, _frag_seq) in fragment_columns.items():
        bucket = buckets.get(frag_index)
        if bucket is not None:
            bucket["msa_query_aligned_fragment"] = msa_query_aligned[
                col_start : col_end + 1
            ]

    stats = {
        "sequences_total": len(rows),
        "sequences_processed": 0,
        "sequences_failed": 0,
        "homologs_processed": 0,
    }

    for seq_id, homolog_aligned in rows:
        if seq_id == query_seq_id:
            continue
        stats["homologs_processed"] += 1
        try:
            for frag_index, (col_start, col_end, _frag_seq) in fragment_columns.items():
                bucket = buckets.get(frag_index)
                filt = filter_by_index.get(frag_index)
                if bucket is None or filt is None:
                    continue

                frag_data = extract_msa_fragment_at_columns(
                    msa_query_aligned,
                    homolog_aligned,
                    col_start,
                    col_end,
                )
                if not frag_data:
                    bucket["rejected_unmapped"] += 1
                    continue
                sequence = frag_data["sequence"]
                if not fragment_passes_terminal_filter(sequence, filt):
                    bucket["rejected_terminal_filter"] += 1
                    continue
                bucket["entries"].append(
                    {
                        "sequence_id": seq_id,
                        "sequence": sequence,
                        "length": len(sequence),
                        "aligned_query_fragment": frag_data["aligned_query_fragment"],
                        "aligned_homolog_fragment": frag_data["aligned_homolog_fragment"],
                        "mutations_non_gap": frag_data["mutations_non_gap"],
                        "mutations_with_gaps": frag_data["mutations_with_gaps"],
                        "pct_identity": frag_data["pct_identity"],
                    }
                )
            stats["sequences_processed"] += 1
        except Exception:
            stats["sequences_failed"] += 1
            continue

    fragment_results = [buckets[idx] for idx in sorted(buckets)]
    for bucket in fragment_results:
        bucket["count"] = len(bucket["entries"])

    return {
        "alignment_length": len(rows[0][1]) if rows else 0,
        "num_sequences": len(rows),
        "query_length": len(query_ungapped),
        "query_msa_seq_id": query_seq_id,
        "fragment_filters": fragment_filters,
        "fragments": fragment_results,
        "stats": stats,
    }


def fragments_to_fasta(fragment_result: Dict[str, Any]) -> str:
    """Format filtered entries for one fragment as FASTA."""
    lines: List[str] = []
    frag_idx = fragment_result["fragment"]
    for entry in fragment_result.get("entries") or []:
        lines.append(f">{entry['sequence_id']}_fragment{frag_idx}")
        lines.append(entry["sequence"])
    return "\n".join(lines) + ("\n" if lines else "")


def _query_aligned_fragment_from_bucket(fragment_result: Dict[str, Any]) -> str:
    aligned = fragment_result.get("msa_query_aligned_fragment")
    if aligned:
        return aligned
    for entry in fragment_result.get("entries") or []:
        aligned = entry.get("aligned_query_fragment")
        if aligned:
            return aligned
    return fragment_result.get("query_sequence") or ""


def build_fragment_table_rows(fragment_result: Dict[str, Any]) -> List[Dict[str, Any]]:
    """
    Build deduplicated table rows for one fragment (query row first).

    Identical homolog sequences are merged; sequence IDs are comma-separated.
    """
    query_seq = fragment_result.get("query_sequence") or ""
    query_aligned = _query_aligned_fragment_from_bucket(fragment_result)
    query_len = len(query_seq)
    query_row: Dict[str, Any] = {
        "row_id": "__query__",
        "sequence_id": "query",
        "sequence": query_seq,
        "aligned_query_fragment": query_aligned,
        "aligned_homolog_fragment": query_aligned,
        "length": query_len,
        "length_display": f"{query_len} (100%)",
        "pct_of_query_length": 100.0 if query_len else 0.0,
        "mutations_non_gap": 0,
        "mutations_with_gaps": 0,
        "pct_identity": 100.0,
        "is_query": True,
        "source_count": 1,
    }

    by_sequence: Dict[str, Dict[str, Any]] = {}
    for entry in fragment_result.get("entries") or []:
        seq = entry.get("sequence") or ""
        if not seq:
            continue
        if seq not in by_sequence:
            by_sequence[seq] = {
                "ids": [],
                "aligned_query_fragment": entry.get("aligned_query_fragment") or query_aligned,
                "aligned_homolog_fragment": entry.get("aligned_homolog_fragment") or seq,
            }
        by_sequence[seq]["ids"].append(entry.get("sequence_id") or "unknown")

    homolog_rows: List[Dict[str, Any]] = []
    for seq, info in sorted(by_sequence.items(), key=lambda item: item[0]):
        homolog_len = len(seq)
        pct_len = (100.0 * homolog_len / query_len) if query_len else 0.0
        aligned_q = info["aligned_query_fragment"]
        aligned_h = info["aligned_homolog_fragment"]
        mutations_non_gap, mutations_with_gaps, pct_identity = metrics_from_aligned_pair(
            aligned_q, aligned_h
        )
        homolog_rows.append(
            {
                "row_id": seq,
                "sequence_id": ", ".join(info["ids"]),
                "sequence": seq,
                "aligned_query_fragment": aligned_q,
                "aligned_homolog_fragment": aligned_h,
                "length": homolog_len,
                "length_display": f"{homolog_len} ({pct_len:.0f}%)",
                "pct_of_query_length": pct_len,
                "mutations_non_gap": mutations_non_gap,
                "mutations_with_gaps": mutations_with_gaps,
                "pct_identity": pct_identity,
                "is_query": False,
                "source_count": len(info["ids"]),
            }
        )

    return [query_row] + homolog_rows


def apply_excluded_rows(
    rows: List[Dict[str, Any]],
    excluded_row_ids: Optional[Set[str]],
) -> List[Dict[str, Any]]:
    """Remove homolog rows whose row_id is in excluded_row_ids; query is always kept."""
    if not excluded_row_ids:
        return list(rows)
    return [
        row
        for row in rows
        if row.get("is_query") or row["row_id"] not in excluded_row_ids
    ]


def apply_fragment_filters_and_exclusions(
    rows: List[Dict[str, Any]],
    filters: Optional[Dict[str, Any]],
    excluded_row_ids: Optional[Set[str]] = None,
) -> List[Dict[str, Any]]:
    """Apply numeric UI filters then network-driven exclusions."""
    filtered = apply_fragment_ui_filters(rows, filters)
    return apply_excluded_rows(filtered, excluded_row_ids)


def get_saveable_homolog_rows(rows: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """Homolog rows from the main filtered table (excludes query)."""
    return [row for row in rows if not row.get("is_query")]


def apply_fragment_ui_filters(
    rows: List[Dict[str, Any]],
    filters: Optional[Dict[str, Any]],
) -> List[Dict[str, Any]]:
    """Apply optional min/max UI filters; query row is always retained."""
    if not filters:
        return list(rows)

    query_rows = [row for row in rows if row.get("is_query")]
    homolog_rows = [row for row in rows if not row.get("is_query")]

    def _in_range(value: float, min_key: str, max_key: str) -> bool:
        min_val = filters.get(min_key)
        max_val = filters.get(max_key)
        if min_val is not None and value < min_val:
            return False
        if max_val is not None and value > max_val:
            return False
        return True

    filtered: List[Dict[str, Any]] = []
    for row in homolog_rows:
        if not _in_range(row["length"], "min_length", "max_length"):
            continue
        if not _in_range(row["pct_of_query_length"], "min_pct_length", "max_pct_length"):
            continue
        if not _in_range(row["mutations_non_gap"], "min_mutations_non_gap", "max_mutations_non_gap"):
            continue
        if not _in_range(row["mutations_with_gaps"], "min_mutations_with_gaps", "max_mutations_with_gaps"):
            continue
        if not _in_range(row["pct_identity"], "min_pct_identity", "max_pct_identity"):
            continue
        filtered.append(row)

    return query_rows + filtered


def build_all_fragment_table_rows(analysis: Dict[str, Any]) -> Dict[int, List[Dict[str, Any]]]:
    """Precompute deduplicated table rows for every fragment."""
    return {
        int(frag["fragment"]): build_fragment_table_rows(frag)
        for frag in analysis.get("fragments") or []
    }


def fragments_to_compact_json(analysis: Dict[str, Any]) -> Dict[str, Any]:
    """Compact export: fragment -> list of {sequence_id, sequence}."""
    return {
        "stats": analysis.get("stats"),
        "query_msa_seq_id": analysis.get("query_msa_seq_id"),
        "fragment_filters": analysis.get("fragment_filters"),
        "fragments": [
            {
                "fragment": f["fragment"],
                "count": f["count"],
                "entries": [
                    {"sequence_id": e["sequence_id"], "sequence": e["sequence"]}
                    for e in f.get("entries") or []
                ],
            }
            for f in analysis.get("fragments") or []
        ],
    }


# ---------------------------------------------------------------------------
# Legacy pairwise helpers (kept for SSN Hamming distances in diversity_ssn.py)
# ---------------------------------------------------------------------------


def _pairwise_alignment(query_ungapped: str, target_ungapped: str) -> Tuple[str, str, float]:
    query_ungapped = strip_msa_gaps(query_ungapped)
    target_ungapped = strip_msa_gaps(target_ungapped)
    alignments = pairwise2.align.globalxx(query_ungapped, target_ungapped)
    if not alignments:
        raise ValueError("Pairwise alignment failed.")
    aligned_q, aligned_t, score, _begin, _end = alignments[0]
    return aligned_q, aligned_t, float(score)
