"""
Load block variant pools from Diversity Analysis session selections.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Tuple


def _variant_from_row(row: Dict[str, Any], *, is_query: bool = False) -> Dict[str, Any]:
    return {
        "row_id": row.get("row_id") or ("__query__" if is_query else row.get("sequence", "")),
        "sequence_id": row.get("sequence_id") or ("query" if is_query else "unknown"),
        "sequence": row.get("sequence") or "",
        "is_query": bool(is_query or row.get("is_query")),
    }


def load_block_pools_from_session(
    selections: Optional[Dict[str, Any]],
) -> Tuple[List[Dict[str, Any]], Optional[str]]:
    """
    Parse ``diversity_saved_selections`` into ordered block pools.

    Returns:
        (blocks, error_message) — blocks is empty and error is set on failure.
    """
    if not selections:
        return [], "No saved fragment lists. Save from Diversity Analysis first."

    blocks: List[Dict[str, Any]] = []
    for frag_key in sorted(selections.keys(), key=lambda k: int(k)):
        block = selections[frag_key] or {}
        query = block.get("query")
        homologs = block.get("homologs") or block.get("activated") or []

        variants: List[Dict[str, Any]] = []
        if query and query.get("sequence"):
            variants.append(_variant_from_row(query, is_query=True))
        elif homologs:
            return [], f"Fragment {frag_key}: missing query row in saved selection."

        seen_ids = {v["row_id"] for v in variants}
        for row in homologs:
            if not row.get("sequence"):
                continue
            rid = row.get("row_id") or row["sequence"]
            if rid in seen_ids:
                continue
            variants.append(_variant_from_row(row))
            seen_ids.add(rid)

        if len(variants) < 2:
            return [], (
                f"Fragment {frag_key}: need at least query + 1 homolog "
                f"({len(variants)} variant(s) found)."
            )

        blocks.append(
            {
                "fragment": int(frag_key),
                "fragment_key": str(frag_key),
                "variants": variants,
                "n_variants": len(variants),
            }
        )

    if not blocks:
        return [], "Saved selections contain no fragments."

    return blocks, None


def pools_summary(blocks: List[Dict[str, Any]]) -> Dict[str, Any]:
    """Summary stats for UI display."""
    combos = 1
    for b in blocks:
        combos *= b["n_variants"]
    return {
        "n_blocks": len(blocks),
        "variants_per_block": [b["n_variants"] for b in blocks],
        "total_variants": sum(b["n_variants"] for b in blocks),
        "total_combinations": combos,
    }


def blocks_to_saved_selections(
    blocks: List[Dict[str, Any]],
    coefficients_by_block: Optional[Dict[str, List[Dict[str, Any]]]] = None,
) -> Dict[str, Any]:
    """Convert internal block pools to ``diversity_saved_selections`` shape."""
    out: Dict[str, Any] = {}
    coeff_lookup: Dict[str, Dict[str, float]] = {}
    if coefficients_by_block:
        for frag_key, rows in coefficients_by_block.items():
            coeff_lookup[frag_key] = {
                r["row_id"]: float(r.get("coefficient", 0.0)) for r in rows
            }

    for block in blocks:
        frag_key = block["fragment_key"]
        query_row = None
        homologs = []
        lookup = coeff_lookup.get(frag_key, {})
        for variant in block["variants"]:
            row = {
                "row_id": variant["row_id"],
                "sequence_id": variant["sequence_id"],
                "sequence": variant["sequence"],
                "length": variant.get("length") or len(variant["sequence"]),
                "length_display": variant.get("length_display")
                or str(variant.get("length") or len(variant["sequence"])),
                "is_query": variant.get("is_query", False),
            }
            for field in (
                "mutations_non_gap",
                "mutations_with_gaps",
                "pct_identity",
                "aligned_query_fragment",
                "aligned_homolog_fragment",
            ):
                if field in variant:
                    row[field] = variant[field]
            if coefficients_by_block:
                row["coefficient"] = lookup.get(variant["row_id"], 0.0)
            if variant.get("is_query"):
                query_row = row
            else:
                homologs.append(row)
        out[frag_key] = {"query": query_row, "homologs": homologs}
    return out


def copy_blocks(blocks: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """Deep-enough copy of block pools for in-session pruning."""
    return [
        {
            **block,
            "variants": [dict(v) for v in block["variants"]],
            "n_variants": block["n_variants"],
        }
        for block in blocks
    ]


def blocks_signature(blocks: List[Dict[str, Any]]) -> Tuple[Tuple[str, Tuple[str, ...]], ...]:
    """Stable signature for detecting when source pools changed."""
    return tuple(
        (
            block["fragment_key"],
            tuple(v["row_id"] for v in block["variants"]),
        )
        for block in blocks
    )


def enrich_variants_from_selections(
    block: Dict[str, Any],
    selection_block: Optional[Dict[str, Any]],
) -> List[Dict[str, Any]]:
    """Attach diversity-analysis metadata to block variants by row_id."""
    if not selection_block:
        return [dict(v) for v in block["variants"]]

    lookup: Dict[str, Dict[str, Any]] = {}
    query = selection_block.get("query")
    if query and query.get("row_id"):
        lookup[query["row_id"]] = query
    for row in selection_block.get("homologs") or []:
        if row.get("row_id"):
            lookup[row["row_id"]] = row

    enriched: List[Dict[str, Any]] = []
    for variant in block["variants"]:
        extra = lookup.get(variant["row_id"], {})
        merged = {**extra, **variant}
        if "length" not in merged:
            merged["length"] = len(merged.get("sequence") or "")
        if "length_display" not in merged:
            merged["length_display"] = str(merged["length"])
        enriched.append(merged)
    return enriched


def build_filter_table_rows(
    block: Dict[str, Any],
    coefficients_by_block: Dict[str, List[Dict[str, Any]]],
    selection_block: Optional[Dict[str, Any]],
) -> List[Dict[str, Any]]:
    """
    Rows for the Library Optimization filter table: query first, homologs by coefficient.
    """
    frag_key = block["fragment_key"]
    coef_lookup = {
        r["row_id"]: float(r.get("coefficient", 0.0))
        for r in (coefficients_by_block.get(frag_key) or [])
    }
    variants = enrich_variants_from_selections(block, selection_block)

    rows: List[Dict[str, Any]] = []
    for variant in variants:
        row_id = variant["row_id"]
        is_query = bool(variant.get("is_query"))
        coef = coef_lookup.get(row_id, 0.0 if is_query else 0.0)
        rows.append(
            {
                "row_id": row_id,
                "sequence_id": variant.get("sequence_id") or "unknown",
                "sequence": variant.get("sequence") or "",
                "coefficient": coef,
                "length": variant.get("length") or len(variant.get("sequence") or ""),
                "length_display": variant.get("length_display")
                or str(variant.get("length") or len(variant.get("sequence") or "")),
                "mutations_non_gap": variant.get("mutations_non_gap"),
                "mutations_with_gaps": variant.get("mutations_with_gaps"),
                "pct_identity": variant.get("pct_identity"),
                "is_query": is_query,
            }
        )

    query_rows = [r for r in rows if r["is_query"]]
    homolog_rows = [r for r in rows if not r["is_query"]]
    homolog_rows.sort(key=lambda r: (-r["coefficient"], r["sequence_id"]))
    return query_rows + homolog_rows


def remove_variants_from_block(
    working_blocks: List[Dict[str, Any]],
    fragment_key: str,
    row_ids: List[str],
) -> List[Dict[str, Any]]:
    """Remove homolog variants by row_id; query is never removed."""
    remove = {rid for rid in row_ids if rid != "__query__"}
    updated: List[Dict[str, Any]] = []
    for block in working_blocks:
        if block["fragment_key"] != fragment_key:
            updated.append(block)
            continue
        kept = [
            v
            for v in block["variants"]
            if v.get("is_query") or v["row_id"] not in remove
        ]
        updated.append({**block, "variants": kept, "n_variants": len(kept)})
    return updated


def keep_only_variants_in_block(
    working_blocks: List[Dict[str, Any]],
    fragment_key: str,
    row_ids: List[str],
) -> List[Dict[str, Any]]:
    """Keep only selected homologs plus query."""
    keep = set(row_ids)
    updated: List[Dict[str, Any]] = []
    for block in working_blocks:
        if block["fragment_key"] != fragment_key:
            updated.append(block)
            continue
        kept = [
            v
            for v in block["variants"]
            if v.get("is_query") or v["row_id"] in keep
        ]
        updated.append({**block, "variants": kept, "n_variants": len(kept)})
    return updated


def keep_top_homologs_by_coefficient(
    working_blocks: List[Dict[str, Any]],
    fragment_key: str,
    coefficients_by_block: Dict[str, List[Dict[str, Any]]],
    top_n: int,
) -> List[Dict[str, Any]]:
    """Keep the query and the top *top_n* homologs by ridge coefficient."""
    if top_n < 1:
        raise ValueError("top_n must be at least 1")

    coef_lookup = {
        r["row_id"]: float(r.get("coefficient", 0.0))
        for r in (coefficients_by_block.get(fragment_key) or [])
    }

    updated: List[Dict[str, Any]] = []
    for block in working_blocks:
        if block["fragment_key"] != fragment_key:
            updated.append(block)
            continue

        homologs = [v for v in block["variants"] if not v.get("is_query")]
        homologs.sort(
            key=lambda v: (
                -coef_lookup.get(v["row_id"], 0.0),
                v.get("sequence_id") or "",
            )
        )
        keep_ids = {v["row_id"] for v in homologs[:top_n]}
        kept = [
            v
            for v in block["variants"]
            if v.get("is_query") or v["row_id"] in keep_ids
        ]
        updated.append({**block, "variants": kept, "n_variants": len(kept)})
    return updated


def apply_coefficient_filter(
    blocks: List[Dict[str, Any]],
    coefficients_by_block: Dict[str, List[Dict[str, Any]]],
    min_coefficient: Optional[float],
) -> List[Dict[str, Any]]:
    """Return block pools with homologs below min_coefficient removed (query always kept)."""
    if min_coefficient is None:
        return blocks

    coeff_lookup: Dict[str, Dict[str, float]] = {}
    for frag_key, rows in coefficients_by_block.items():
        coeff_lookup[frag_key] = {
            row["row_id"]: float(row.get("coefficient", 0.0)) for row in rows
        }

    filtered: List[Dict[str, Any]] = []
    for block in blocks:
        frag_key = block["fragment_key"]
        lookup = coeff_lookup.get(frag_key, {})
        kept = []
        for variant in block["variants"]:
            if variant.get("is_query"):
                kept.append(variant)
                continue
            coef = lookup.get(variant["row_id"], 0.0)
            if coef >= min_coefficient:
                kept.append(variant)
        filtered.append({**block, "variants": kept, "n_variants": len(kept)})
    return filtered
