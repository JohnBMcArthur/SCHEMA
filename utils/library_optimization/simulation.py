"""
Random chimera ESM2 simulation vs wildtype (all-query chimera).
"""

from __future__ import annotations

from typing import Any, Callable, Dict, List, Optional, Tuple

from utils.library_optimization.esm_scoring import (
    load_esm_model,
    score_chimeras,
    score_sequences_batch,
)
from utils.library_optimization.pools import blocks_signature
from utils.library_optimization.sampling import reindex_chimera_ids, sample_random_chimeras


BatchCallback = Callable[[List[Dict[str, Any]], List[float], int, int], None]


def chimera_score_record(
    chimera: Dict[str, Any],
    esm_score: float,
    wildtype_score: float,
) -> Dict[str, Any]:
    return {
        "chimera_id": chimera.get("chimera_id"),
        "sequence": chimera.get("sequence") or "",
        "esm_score": float(esm_score),
        "score_delta": float(esm_score) - float(wildtype_score),
        "is_wildtype": bool(chimera.get("is_wildtype")),
    }


def wildtype_table_row(
    wildtype_sequence: str,
    wildtype_score: float,
) -> Dict[str, Any]:
    return {
        "chimera_id": "wildtype",
        "sequence": wildtype_sequence,
        "esm_score": float(wildtype_score),
        "score_delta": 0.0,
        "is_wildtype": True,
    }


def top_scored_table_rows(
    results: Dict[str, Any],
    *,
    top_n: int = 20,
) -> List[Dict[str, Any]]:
    """Wildtype first, then top *top_n* chimeras by Δ vs wildtype (highest first)."""
    wt_score = results.get("wildtype_score")
    wt_seq = results.get("wildtype_sequence") or ""
    if wt_score is None:
        return []

    rows: List[Dict[str, Any]] = [wildtype_table_row(wt_seq, float(wt_score))]
    homologs = [
        c
        for c in (results.get("chimeras") or [])
        if not c.get("is_wildtype")
    ]
    homologs.sort(
        key=lambda c: (-float(c.get("score_delta", 0.0)), str(c.get("chimera_id") or "")),
    )
    rows.extend(homologs[:top_n])
    return rows


def top_sequences_to_fasta(
    results: Dict[str, Any],
    *,
    top_n: int = 20,
) -> str:
    """FASTA for wildtype + top *top_n* chimeras (same order as the results table)."""
    lines: List[str] = []
    for row in top_scored_table_rows(results, top_n=top_n):
        seq_id = str(row.get("chimera_id") or "sequence")
        delta = float(row.get("score_delta", 0.0))
        esm = float(row.get("esm_score", 0.0))
        sequence = (row.get("sequence") or "").strip()
        if not sequence:
            continue
        lines.append(f">{seq_id} delta={delta:.4f} esm2={esm:.4f}")
        lines.append(sequence)
    return "\n".join(lines) + ("\n" if lines else "")


def build_wildtype_chimera(blocks: List[Dict[str, Any]]) -> Dict[str, Any]:
    """Chimera using the query (wildtype) variant at every block."""
    parts: List[str] = []
    picks: List[Dict[str, Any]] = []
    for block in blocks:
        query = next(
            (v for v in block["variants"] if v.get("is_query")),
            block["variants"][0],
        )
        parts.append(query["sequence"])
        picks.append(
            {
                "fragment": block["fragment"],
                "fragment_key": block["fragment_key"],
                "row_id": query["row_id"],
                "sequence_id": query["sequence_id"],
                "is_query": True,
            }
        )
    return {
        "chimera_id": "wildtype",
        "sequence": "".join(parts),
        "picks": picks,
        "is_wildtype": True,
    }


def score_wildtype(
    blocks: List[Dict[str, Any]],
    *,
    model: Any = None,
    alphabet: Any = None,
    device: Optional[str] = None,
    esm_batch_size: int = 8,
) -> Tuple[float, Dict[str, Any]]:
    """Return (ESM2 score, wildtype chimera dict)."""
    chimera = build_wildtype_chimera(blocks)
    own_model = model is None
    if own_model:
        model, alphabet, device = load_esm_model()

    try:
        score_chimeras(
            [chimera],
            model=model,
            alphabet=alphabet,
            device=device,
            batch_size=esm_batch_size,
        )
    finally:
        if own_model:
            del model
            try:
                import torch

                if torch.cuda.is_available():
                    torch.cuda.empty_cache()
            except Exception:
                pass

    score = float(chimera.get("esm_score") or 0.0)
    return score, chimera


def simulate_random_chimeras_esm(
    blocks: List[Dict[str, Any]],
    n_samples: int,
    wildtype_score: float,
    *,
    start_index: int = 0,
    seed: Optional[int] = None,
    esm_batch_size: int = 8,
    on_batch: Optional[BatchCallback] = None,
) -> Tuple[List[Dict[str, Any]], List[float], List[float]]:
    """
    Sample and ESM-score random chimeras.

    Returns (scored chimera records, raw scores, score − wildtype).

    *on_batch* is called after each ESM batch with
    (scored_records_this_run, score_deltas_this_run, n_scored_this_run, n_total_this_run).
    """
    if n_samples < 1:
        raise ValueError("n_samples must be at least 1")

    round_seed = (int(seed) + start_index) if seed is not None else None
    chimeras = sample_random_chimeras(blocks, n_samples, seed=round_seed)
    chimeras = reindex_chimera_ids(chimeras, start=start_index)

    model, alphabet, device = load_esm_model()
    raw_scores: List[float] = []
    deltas: List[float] = []
    records: List[Dict[str, Any]] = []

    try:
        sequences = [c["sequence"] for c in chimeras]
        for batch_start in range(0, len(sequences), esm_batch_size):
            batch_chimeras = chimeras[batch_start : batch_start + esm_batch_size]
            batch_seqs = sequences[batch_start : batch_start + esm_batch_size]
            batch_scores = score_sequences_batch(
                model,
                alphabet,
                batch_seqs,
                device,
                batch_size=len(batch_seqs),
            )
            for chimera, score in zip(batch_chimeras, batch_scores):
                raw_scores.append(score)
                deltas.append(score - wildtype_score)
                records.append(chimera_score_record(chimera, score, wildtype_score))
            if on_batch:
                on_batch(records, deltas, len(deltas), len(sequences))
    finally:
        del model
        try:
            import torch

            if torch.cuda.is_available():
                torch.cuda.empty_cache()
        except Exception:
            pass

    return records, raw_scores, deltas


def pools_signature_for_sim(blocks: List[Dict[str, Any]]) -> str:
    """Serializable signature string for session invalidation."""
    return repr(blocks_signature(blocks))


def empty_sim_results() -> Dict[str, Any]:
    return {
        "version": "1.1",
        "pools_signature": None,
        "wildtype_score": None,
        "wildtype_sequence": None,
        "score_deltas": [],
        "chimeras": [],
        "n_chimeras": 0,
    }


def merge_sim_run(
    prior: Optional[Dict[str, Any]],
    *,
    pools_signature: str,
    wildtype_score: float,
    wildtype_sequence: str,
    new_records: List[Dict[str, Any]],
) -> Dict[str, Any]:
    """Append a simulation run to prior session results."""
    base = dict(prior) if prior else empty_sim_results()
    prior_records = list(base.get("chimeras") or [])
    all_records = prior_records + list(new_records)
    all_deltas = [float(r["score_delta"]) for r in all_records]
    return {
        "version": "1.1",
        "pools_signature": pools_signature,
        "wildtype_score": wildtype_score,
        "wildtype_sequence": wildtype_sequence,
        "score_deltas": all_deltas,
        "chimeras": all_records,
        "n_chimeras": len(all_records),
    }
