"""
ESM2 sequence scoring (normalized per-residue log-likelihood).
"""

from __future__ import annotations

from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np

_ESM_MODEL_NAME = "esm2_t30_150M_UR50D"
_MAX_SEQ_LEN = 1022  # ESM2 context limit minus special tokens


def esm_dependencies_available() -> Tuple[bool, Optional[str]]:
    """Return (available, error_message)."""
    try:
        import torch  # noqa: F401
        import esm  # noqa: F401
    except ImportError as exc:
        return False, (
            "ESM scoring requires optional dependencies. Install with:\n"
            "`pip install -r requirements-optimization.txt`"
            f"\n\nMissing: {exc}"
        )
    return True, None


def load_esm_model(
    device: Optional[str] = None,
) -> Tuple[Any, Any, str]:
    """
    Load ESM2-150M model, alphabet, and resolved device string.

    Returns (model, alphabet, device).
    """
    import torch
    import esm

    ok, err = esm_dependencies_available()
    if not ok:
        raise ImportError(err or "ESM dependencies not available.")

    if device is None:
        device = "cuda" if torch.cuda.is_available() else "cpu"

    model, alphabet = esm.pretrained.load_model_and_alphabet(_ESM_MODEL_NAME)
    model = model.to(device)
    model.eval()
    return model, alphabet, device


def _truncate_sequence(sequence: str) -> str:
    seq = (sequence or "").upper().strip()
    if len(seq) > _MAX_SEQ_LEN:
        return seq[:_MAX_SEQ_LEN]
    return seq


def score_sequences_batch(
    model: Any,
    alphabet: Any,
    sequences: Sequence[str],
    device: str,
    *,
    batch_size: int = 8,
    progress_callback: Optional[Callable[[float, str], None]] = None,
) -> List[float]:
    """
    Normalized per-residue log-likelihood for each sequence.

    Uses a single forward pass per sequence (log P(aa_i | full sequence)),
    averaged over residues. Sequences are uppercased and truncated to ESM limits.
    """
    import torch

    batch_converter = alphabet.get_batch_converter()
    cleaned = [_truncate_sequence(s) for s in sequences]
    scores: List[float] = []
    n = len(cleaned)
    if n == 0:
        return scores

    for start in range(0, n, batch_size):
        batch = cleaned[start : start + batch_size]
        labels = [(f"s{i}", seq) for i, seq in enumerate(batch)]
        _, _, tokens = batch_converter(labels)
        tokens = tokens.to(device)

        with torch.no_grad():
            logits = model(tokens)["logits"]

        log_probs = torch.log_softmax(logits, dim=-1)
        for row_idx in range(len(batch)):
            token_row = tokens[row_idx]
            lp_row = log_probs[row_idx]
            # Skip BOS (0) and EOS (last); score interior residues only.
            interior = range(1, token_row.size(0) - 1)
            vals = []
            for pos in interior:
                tok = int(token_row[pos].item())
                vals.append(float(lp_row[pos, tok].item()))
            scores.append(float(np.mean(vals)) if vals else float("nan"))

        if progress_callback:
            frac = min(1.0, (start + len(batch)) / n)
            progress_callback(frac, f"ESM scoring {start + len(batch)}/{n} sequences…")

    return scores


def score_chimeras(
    chimeras: List[Dict[str, Any]],
    *,
    model: Any = None,
    alphabet: Any = None,
    device: Optional[str] = None,
    batch_size: int = 8,
    progress_callback: Optional[Callable[[float, str], None]] = None,
) -> List[Dict[str, Any]]:
    """
    Attach ``esm_score`` to each chimera dict (in place + return list).
    """
    own_model = model is None
    if own_model:
        model, alphabet, device = load_esm_model(device)

    sequences = [c["sequence"] for c in chimeras]
    scores = score_sequences_batch(
        model,
        alphabet,
        sequences,
        device,
        batch_size=batch_size,
        progress_callback=progress_callback,
    )
    for chimera, score in zip(chimeras, scores):
        chimera["esm_score"] = score

    if own_model:
        del model
        try:
            import torch

            if torch.cuda.is_available():
                torch.cuda.empty_cache()
        except Exception:
            pass

    return chimeras
