"""
Ridge regression on block-variant indicators (query reference per block).
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Tuple

import numpy as np


def _reference_row_id(block: Dict[str, Any]) -> str:
    for variant in block["variants"]:
        if variant.get("is_query"):
            return variant["row_id"]
    return block["variants"][0]["row_id"]


def build_design_matrix(
    chimeras: List[Dict[str, Any]],
    blocks: List[Dict[str, Any]],
) -> Tuple[np.ndarray, List[Dict[str, Any]], List[str]]:
    """
    Build ridge design matrix with query as reference (coefficient 0) per block.

    Returns:
        X (n_samples x n_features), feature_meta (list of dicts), chimera_ids
    """
    feature_meta: List[Dict[str, Any]] = []
    ref_by_block: Dict[str, str] = {}

    for block in blocks:
        frag_key = block["fragment_key"]
        ref_id = _reference_row_id(block)
        ref_by_block[frag_key] = ref_id
        for variant in block["variants"]:
            if variant["row_id"] == ref_id:
                continue
            feature_meta.append(
                {
                    "fragment": block["fragment"],
                    "fragment_key": frag_key,
                    "row_id": variant["row_id"],
                    "sequence_id": variant["sequence_id"],
                    "is_query": False,
                }
            )

    n = len(chimeras)
    p = len(feature_meta)
    if n == 0 or p == 0:
        return np.zeros((n, max(p, 1))), feature_meta, []

    X = np.zeros((n, p), dtype=np.float64)
    chimera_ids: List[str] = []

    # Map (fragment_key, row_id) -> column index
    col_index: Dict[Tuple[str, str], int] = {}
    for j, meta in enumerate(feature_meta):
        col_index[(meta["fragment_key"], meta["row_id"])] = j

    for i, chimera in enumerate(chimeras):
        chimera_ids.append(chimera.get("chimera_id", f"chimera_{i}"))
        pick_map = {p["fragment_key"]: p["row_id"] for p in chimera.get("picks", [])}
        for block in blocks:
            frag_key = block["fragment_key"]
            row_id = pick_map.get(frag_key)
            if row_id is None or row_id == ref_by_block[frag_key]:
                continue
            key = (frag_key, row_id)
            if key in col_index:
                X[i, col_index[key]] = 1.0

    return X, feature_meta, chimera_ids


def fit_ridge(
    chimeras: List[Dict[str, Any]],
    blocks: List[Dict[str, Any]],
    *,
    alpha: float = 1.0,
) -> Tuple[Dict[str, List[Dict[str, Any]]], Dict[str, Any]]:
    """
    Fit ridge regression: esm_score ~ block variant indicators.

    Returns:
        coefficients_by_block, fit_stats
    """
    y_list = [c.get("esm_score") for c in chimeras]
    valid = [
        (i, c)
        for i, c in enumerate(chimeras)
        if c.get("esm_score") is not None and not np.isnan(c["esm_score"])
    ]
    if len(valid) < 2:
        empty = _coefficients_from_values(blocks, {})
        return empty, {"n_samples": len(valid), "r2": None, "alpha": alpha}

    valid_chimeras = [c for _, c in valid]
    y = np.array([c["esm_score"] for c in valid_chimeras], dtype=np.float64)
    X, feature_meta, _ = build_design_matrix(valid_chimeras, blocks)

    if X.shape[1] == 0:
        coef_by_block = _coefficients_from_values(blocks, {})
        return coef_by_block, {
            "n_samples": len(valid_chimeras),
            "r2": None,
            "alpha": alpha,
            "intercept": float(np.mean(y)),
        }

    # Ridge: (X'X + alpha I)^{-1} X'y
    XtX = X.T @ X
    penalty = alpha * np.eye(X.shape[1])
    beta = np.linalg.solve(XtX + penalty, X.T @ y)
    intercept = float(np.mean(y) - np.mean(X @ beta))

    y_hat = X @ beta + intercept
    ss_res = float(np.sum((y - y_hat) ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else None

    coef_map: Dict[Tuple[str, str], float] = {}
    for j, meta in enumerate(feature_meta):
        coef_map[(meta["fragment_key"], meta["row_id"])] = float(beta[j])

    coef_by_block = _coefficients_from_values(blocks, coef_map, intercept=intercept)

    # Exposure counts: how often each variant appeared in sampled chimeras
    exposure = _variant_exposure(chimeras, blocks)
    for frag_key, rows in coef_by_block.items():
        for row in rows:
            row["exposure_count"] = exposure.get((frag_key, row["row_id"]), 0)

    return coef_by_block, {
        "n_samples": len(valid_chimeras),
        "n_features": X.shape[1],
        "r2": r2,
        "alpha": alpha,
        "intercept": intercept,
    }


def _coefficients_from_values(
    blocks: List[Dict[str, Any]],
    coef_map: Dict[Tuple[str, str], float],
    *,
    intercept: Optional[float] = None,
) -> Dict[str, List[Dict[str, Any]]]:
    out: Dict[str, List[Dict[str, Any]]] = {}
    for block in blocks:
        frag_key = block["fragment_key"]
        ref_id = _reference_row_id(block)
        rows: List[Dict[str, Any]] = []
        for variant in block["variants"]:
            is_query = variant["row_id"] == ref_id
            coef = 0.0 if is_query else coef_map.get((frag_key, variant["row_id"]), 0.0)
            rows.append(
                {
                    "row_id": variant["row_id"],
                    "sequence_id": variant["sequence_id"],
                    "sequence": variant["sequence"],
                    "is_query": is_query,
                    "coefficient": coef,
                }
            )
        rows.sort(key=lambda r: (-r["coefficient"], r["sequence_id"]))
        out[frag_key] = rows
    return out


def _variant_exposure(
    chimeras: List[Dict[str, Any]],
    blocks: List[Dict[str, Any]],
) -> Dict[Tuple[str, str], int]:
    counts: Dict[Tuple[str, str], int] = {}
    for chimera in chimeras:
        for pick in chimera.get("picks", []):
            key = (pick["fragment_key"], pick["row_id"])
            counts[key] = counts.get(key, 0) + 1
    return counts
