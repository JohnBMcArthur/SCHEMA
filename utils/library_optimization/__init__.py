"""
Iterative library optimization: random chimera sampling, ESM scoring, ridge regression.

Optional dependency — install requirements-optimization.txt for ESM2 scoring.
"""

from utils.library_optimization.pools import load_block_pools_from_session
from utils.library_optimization.optimizer import (
    run_optimization_loop,
    run_optimization_rounds,
)

__all__ = [
    "load_block_pools_from_session",
    "run_optimization_loop",
    "run_optimization_rounds",
]
