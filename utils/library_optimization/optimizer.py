"""
Main optimization loop: sample → ESM score → ridge → convergence check.
Supports resuming from prior results and optional early stopping on convergence.
"""

from __future__ import annotations

from typing import Any, Callable, Dict, List, Optional, Tuple

from utils.library_optimization.convergence import (
    all_blocks_converged,
    spearman_rank_correlation_by_block,
)
from utils.library_optimization.esm_scoring import load_esm_model, score_chimeras
from utils.library_optimization.ridge_fit import fit_ridge
from utils.library_optimization.sampling import reindex_chimera_ids, sample_random_chimeras


ProgressCallback = Callable[[Dict[str, Any]], None]


def _coefficient_history_tuples(
    results: Dict[str, Any],
) -> List[Tuple[int, Dict[str, List[Dict[str, Any]]]]]:
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


def optimization_state_from_results(
    results: Dict[str, Any],
) -> Dict[str, Any]:
    """Extract mutable loop state from a stored results dict."""
    rounds = list(results.get("rounds") or [])
    coefficient_history = _coefficient_history_tuples(results)
    prev_coefficients: Optional[Dict[str, List[Dict[str, Any]]]] = None
    if rounds:
        prev_coefficients = rounds[-1].get("coefficients_by_block")
    if prev_coefficients is None:
        prev_coefficients = results.get("coefficients_by_block")

    return {
        "all_chimeras": list(results.get("chimeras") or []),
        "round_history": rounds,
        "coefficient_history": coefficient_history,
        "prev_coefficients": prev_coefficients,
        "converged": bool(results.get("converged")),
        "next_round": len(rounds) + 1,
        "params": dict(results.get("params") or {}),
        "blocks_summary": results.get("blocks_summary"),
    }


def build_results_payload(
    *,
    all_chimeras: List[Dict[str, Any]],
    round_history: List[Dict[str, Any]],
    coefficient_history: List[Tuple[int, Dict[str, List[Dict[str, Any]]]]],
    converged: bool,
    params: Dict[str, Any],
    blocks: List[Dict[str, Any]],
    blocks_summary: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """Assemble the session/export results dict from loop state."""
    final_coefficients = (
        round_history[-1]["coefficients_by_block"] if round_history else {}
    )
    return {
        "params": params,
        "rounds": round_history,
        "coefficient_history": [
            {"round": rnd, "coefficients_by_block": coefs}
            for rnd, coefs in coefficient_history
        ],
        "n_chimeras": len(all_chimeras),
        "chimeras": all_chimeras,
        "coefficients_by_block": final_coefficients,
        "converged": converged,
        "blocks_summary": blocks_summary
        or {
            "n_blocks": len(blocks),
            "variants_per_block": [b["n_variants"] for b in blocks],
        },
        "status": "converged" if converged else "active",
    }


def run_optimization_rounds(
    blocks: List[Dict[str, Any]],
    *,
    batch_size: int,
    n_rounds: int,
    min_rounds: int = 2,
    spearman_threshold: float = 0.9,
    stop_on_convergence: bool = True,
    ridge_alpha: float = 1.0,
    esm_batch_size: int = 8,
    seed: Optional[int] = None,
    prior: Optional[Dict[str, Any]] = None,
    progress_callback: Optional[ProgressCallback] = None,
) -> Dict[str, Any]:
    """
    Run one or more optimization rounds, optionally continuing from *prior* results.

    All chimeras and scores accumulate across rounds. When *prior* is provided,
    sampling and round numbering continue from the last stored round without
    discarding previously scored designs.

    If *stop_on_convergence* is True, the loop may end before *n_rounds* when
    Spearman ρ stabilizes (after *min_rounds*). If False, always runs exactly
    *n_rounds* unless *n_rounds* is 0.
    """
    if int(batch_size) < 1:
        raise ValueError("batch_size must be at least 1")
    if int(n_rounds) < 1:
        raise ValueError("n_rounds must be at least 1")

    if prior:
        state = optimization_state_from_results(prior)
        all_chimeras = state["all_chimeras"]
        round_history = state["round_history"]
        coefficient_history = state["coefficient_history"]
        prev_coefficients = state["prev_coefficients"]
        converged = state["converged"]
        stored_params = state["params"]
        blocks_summary = state["blocks_summary"]
        start_round = state["next_round"]
    else:
        all_chimeras = []
        round_history = []
        coefficient_history = []
        prev_coefficients = None
        converged = False
        stored_params = {}
        blocks_summary = None
        start_round = 1

    min_rounds = int(stored_params.get("min_rounds", min_rounds)) if prior else int(min_rounds)

    params = {
        **stored_params,
        "batch_size": batch_size,
        "min_rounds": min_rounds,
        "spearman_threshold": spearman_threshold,
        "ridge_alpha": ridge_alpha,
        "esm_batch_size": esm_batch_size,
        "seed": seed if seed is not None else stored_params.get("seed"),
        "stop_on_convergence": stop_on_convergence,
        "n_rounds_last_run": n_rounds,
    }

    end_round = start_round + int(n_rounds) - 1

    def _emit(event: Dict[str, Any]) -> None:
        if progress_callback:
            progress_callback(event)

    model, alphabet, device = load_esm_model()

    try:
        for round_num in range(start_round, end_round + 1):
            _emit(
                {
                    "phase": "sampling",
                    "round": round_num,
                    "message": f"Round {round_num}: sampling {batch_size} chimeras…",
                }
            )
            round_seed = None
            base_seed = params.get("seed")
            if base_seed is not None:
                round_seed = int(base_seed) + round_num

            new_batch = sample_random_chimeras(
                blocks,
                batch_size,
                seed=round_seed,
            )
            start_idx = len(all_chimeras)
            new_batch = reindex_chimera_ids(new_batch, start=start_idx)
            all_chimeras.extend(new_batch)

            def esm_progress(frac: float, msg: str) -> None:
                _emit(
                    {
                        "phase": "esm",
                        "round": round_num,
                        "fraction": frac,
                        "message": f"Round {round_num}: {msg}",
                    }
                )

            _emit(
                {
                    "phase": "esm",
                    "round": round_num,
                    "message": f"Round {round_num}: ESM scoring {batch_size} sequences…",
                }
            )
            score_chimeras(
                new_batch,
                model=model,
                alphabet=alphabet,
                device=device,
                batch_size=esm_batch_size,
                progress_callback=esm_progress,
            )

            _emit(
                {
                    "phase": "ridge",
                    "round": round_num,
                    "message": f"Round {round_num}: fitting ridge regression…",
                }
            )
            coefficients, fit_stats = fit_ridge(all_chimeras, blocks, alpha=ridge_alpha)
            spearman = spearman_rank_correlation_by_block(coefficients, prev_coefficients)
            coefficient_history.append((round_num, coefficients))

            round_converged = (
                stop_on_convergence
                and round_num >= min_rounds
                and prev_coefficients is not None
                and all_blocks_converged(spearman, spearman_threshold)
            )

            round_record = {
                "round": round_num,
                "n_new_chimeras": batch_size,
                "n_total_chimeras": len(all_chimeras),
                "spearman_by_block": spearman,
                "fit_stats": fit_stats,
                "converged": round_converged,
                "coefficients_by_block": coefficients,
            }
            round_history.append(round_record)
            prev_snapshot = prev_coefficients
            prev_coefficients = coefficients

            partial = build_results_payload(
                all_chimeras=all_chimeras,
                round_history=round_history,
                coefficient_history=coefficient_history,
                converged=round_converged,
                params=params,
                blocks=blocks,
                blocks_summary=blocks_summary,
            )

            _emit(
                {
                    "phase": "round_complete",
                    "round": round_num,
                    "spearman_by_block": spearman,
                    "fit_stats": fit_stats,
                    "converged": round_converged,
                    "coefficients_by_block": coefficients,
                    "prev_coefficients_by_block": prev_snapshot,
                    "partial_results": partial,
                    "message": f"Round {round_num} complete.",
                }
            )

            if round_converged:
                converged = True
                break
        else:
            converged = False
    finally:
        del model
        try:
            import torch

            if torch.cuda.is_available():
                torch.cuda.empty_cache()
        except Exception:
            pass

    return build_results_payload(
        all_chimeras=all_chimeras,
        round_history=round_history,
        coefficient_history=coefficient_history,
        converged=converged,
        params=params,
        blocks=blocks,
        blocks_summary=blocks_summary,
    )


def run_optimization_loop(
    blocks: List[Dict[str, Any]],
    *,
    batch_size: int = 2000,
    max_rounds: int = 20,
    min_rounds: int = 2,
    spearman_threshold: float = 0.9,
    stop_on_convergence: bool = True,
    ridge_alpha: float = 1.0,
    esm_batch_size: int = 8,
    seed: Optional[int] = None,
    prior: Optional[Dict[str, Any]] = None,
    progress_callback: Optional[ProgressCallback] = None,
) -> Dict[str, Any]:
    """
    Back-compatible wrapper around :func:`run_optimization_rounds`.

    Fresh runs execute up to *max_rounds* rounds. When *prior* is set, runs
    *max_rounds* additional rounds (name kept for callers that pass max_rounds).
    """
    return run_optimization_rounds(
        blocks,
        batch_size=batch_size,
        n_rounds=max_rounds,
        min_rounds=min_rounds,
        spearman_threshold=spearman_threshold,
        stop_on_convergence=stop_on_convergence,
        ridge_alpha=ridge_alpha,
        esm_batch_size=esm_batch_size,
        seed=seed,
        prior=prior,
        progress_callback=progress_callback,
    )
