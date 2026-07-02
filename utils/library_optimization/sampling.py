"""
Random chimera sampling from block variant pools.
"""

from __future__ import annotations

import random
from typing import Any, Dict, List


def sample_random_chimeras(
    blocks: List[Dict[str, Any]],
    n_samples: int,
    *,
    seed: int | None = None,
) -> List[Dict[str, Any]]:
    """
    Uniformly sample *n_samples* chimeras (one variant per block, with replacement).

    Each returned dict has:
        chimera_id, sequence, picks (list of {fragment, row_id, sequence_id, is_query})
    """
    if seed is not None:
        random.seed(seed)

    chimeras: List[Dict[str, Any]] = []
    for i in range(n_samples):
        parts: List[str] = []
        picks: List[Dict[str, Any]] = []
        for block in blocks:
            variant = random.choice(block["variants"])
            parts.append(variant["sequence"])
            picks.append(
                {
                    "fragment": block["fragment"],
                    "fragment_key": block["fragment_key"],
                    "row_id": variant["row_id"],
                    "sequence_id": variant["sequence_id"],
                    "is_query": variant.get("is_query", False),
                }
            )
        chimeras.append(
            {
                "chimera_id": f"chimera_{i}",
                "sequence": "".join(parts),
                "picks": picks,
            }
        )
    return chimeras


def reindex_chimera_ids(chimeras: List[Dict[str, Any]], start: int = 0) -> List[Dict[str, Any]]:
    """Assign stable chimera_id values starting at *start*."""
    out = []
    for offset, chimera in enumerate(chimeras):
        out.append({**chimera, "chimera_id": f"chimera_{start + offset}"})
    return out
