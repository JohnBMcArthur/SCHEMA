"""
Design stuffer tails for oligo padding: ATGC sequences replacing N placeholders.

Uses precomputed safe sequences from ``data/stuffer_sequences.yaml``.
"""

from __future__ import annotations

import re
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, List, Tuple

import yaml

from utils.config import STUFFER_SEQUENCES_YAML

_BSAI_MOTIFS = ("GGTCTC", "GAGACC")
_HOMOPOLYMER_RUN = re.compile(r"(.)\1{3,}", re.IGNORECASE)
_COMPLEMENT = {"A": "T", "T": "A", "G": "C", "C": "G"}
_CONTEXT_WINDOW = 8
_BASES = ("A", "T", "G", "C")


def reverse_complement(dna: str) -> str:
    return "".join(_COMPLEMENT.get(b, "N") for b in reversed((dna or "").upper()))


def split_stuffer_lengths(total: int) -> Tuple[int, int]:
    """Split total padding equally between 5′ and 3′ outside tails."""
    if total <= 0:
        return 0, 0
    left = total // 2
    return left, total - left


def assemble_oligo_with_outside_stuffers(
    stuffer_5: str,
    forward_primer: str,
    merged_insert: str,
    reverse_primer_rc: str,
    stuffer_3: str,
) -> str:
    """5′ stuffer + forward primer + insert + RC(reverse primer) + 3′ stuffer."""
    return f"{stuffer_5}{forward_primer}{merged_insert}{reverse_primer_rc}{stuffer_3}"


def build_oligo_stuffer_fields(
    padding: int,
    *,
    stuffer_base: str = "N",
) -> Dict[str, Any]:
    """Return stuffer metadata for outside-split padding."""
    stuffer_5_len, stuffer_3_len = split_stuffer_lengths(padding)
    stuffer_5 = (stuffer_base * stuffer_5_len) if stuffer_5_len > 0 else ""
    stuffer_3 = (stuffer_base * stuffer_3_len) if stuffer_3_len > 0 else ""
    return {
        "stuffer_5": stuffer_5,
        "stuffer_3": stuffer_3,
        "stuffer_5_length": stuffer_5_len,
        "stuffer_3_length": stuffer_3_len,
        "stuffer": stuffer_5 + stuffer_3,
        "stuffer_length": padding,
        "stuffer_designed": stuffer_base.upper() != "N",
    }


def gc_fraction(sequence: str) -> float:
    seq = (sequence or "").upper()
    if not seq:
        return 0.0
    return sum(base in "GC" for base in seq) / len(seq)


def _has_excessive_homopolymer(sequence: str) -> bool:
    return bool(_HOMOPOLYMER_RUN.search(sequence or ""))


def _contains_bsai_motif(sequence: str) -> bool:
    seq = (sequence or "").upper()
    rc = reverse_complement(seq)
    for strand in (seq, rc):
        for motif in _BSAI_MOTIFS:
            if motif in strand:
                return True
    return False


def _has_strong_hairpin(sequence: str, *, min_stem: int = 6, min_loop: int = 3) -> bool:
    """Detect a stem-loop with *min_stem* paired bases and at least *min_loop* loop nt."""
    seq = (sequence or "").upper()
    n = len(seq)
    if n < 2 * min_stem + min_loop:
        return False
    for stem_len in range(min_stem, n // 2 + 1):
        for loop in range(min_loop, min(n, 24)):
            max_start = n - 2 * stem_len - loop
            if max_start < 0:
                continue
            for i in range(max_start + 1):
                stem5 = seq[i : i + stem_len]
                stem3_start = i + stem_len + loop
                stem3 = seq[stem3_start : stem3_start + stem_len]
                if stem5 == reverse_complement(stem3):
                    return True
    return False


def is_internally_safe_stuffer(
    seq: str,
    *,
    gc_min: float = 0.35,
    gc_max: float = 0.65,
    check_hairpin: bool = True,
) -> bool:
    """True when *seq* meets internal stuffer constraints (no primer context)."""
    if not seq:
        return False
    if _has_excessive_homopolymer(seq):
        return False
    if _contains_bsai_motif(seq):
        return False
    if check_hairpin and _has_strong_hairpin(seq):
        return False
    if len(seq) >= 4:
        gc = gc_fraction(seq)
        if not (gc_min <= gc <= gc_max):
            return False
    return True


def _validation_window(
    stuffer: str,
    left_context: str,
    right_context: str,
) -> str:
    left = (left_context or "")[-_CONTEXT_WINDOW:]
    right = (right_context or "")[:_CONTEXT_WINDOW]
    return f"{left}{stuffer}{right}"


def _extends_context_homopolymer(
    stuffer: str,
    left_context: str,
    right_context: str,
) -> bool:
    """Reject stuffers that extend a context homopolymer past 3 nt or run >3 internally."""
    if not stuffer:
        return False
    if _has_excessive_homopolymer(stuffer):
        return True

    left = (left_context or "").upper()
    right = (right_context or "").upper()
    st = stuffer.upper()

    if left:
        base = left[-1]
        leading_st = 0
        for ch in st:
            if ch == base:
                leading_st += 1
            else:
                break
        if leading_st > 0:
            trailing_left = 0
            for ch in reversed(left):
                if ch == base:
                    trailing_left += 1
                else:
                    break
            if trailing_left + leading_st > 3:
                return True

    if right:
        base = right[0]
        trailing_st = 0
        for ch in reversed(st):
            if ch == base:
                trailing_st += 1
            else:
                break
        if trailing_st > 0:
            leading_right = 0
            for ch in right:
                if ch == base:
                    leading_right += 1
                else:
                    break
            if trailing_st + leading_right > 3:
                return True

    return False


def _passes_stuffer_constraints(
    stuffer: str,
    left_context: str,
    right_context: str,
    *,
    check_hairpin: bool = False,
) -> bool:
    """Runtime junction check for a precomputed stuffer (hairpin verified offline)."""
    if not stuffer:
        return True
    window = _validation_window(stuffer, left_context, right_context)
    if _extends_context_homopolymer(stuffer, left_context, right_context):
        return False
    if _contains_bsai_motif(window):
        return False
    if check_hairpin and _has_strong_hairpin(stuffer):
        return False
    return True


@lru_cache(maxsize=1)
def _load_stuffer_library(path: str) -> Tuple[int, Dict[int, List[str]]]:
    """Return (max_length, {length: [candidates]})."""
    yaml_path = Path(path)
    if not yaml_path.is_file():
        return 0, {}
    with open(yaml_path, encoding="utf-8") as handle:
        data = yaml.safe_load(handle) or {}
    max_length = int(data.get("max_length") or 0)
    raw = data.get("sequences") or {}
    library: Dict[int, List[str]] = {}
    for key, candidates in raw.items():
        length = int(key)
        library[length] = [str(c).upper() for c in (candidates or []) if c]
    return max_length, library


def load_stuffer_library() -> Tuple[int, Dict[int, List[str]]]:
    return _load_stuffer_library(str(STUFFER_SEQUENCES_YAML))


def _pick_from_candidates(
    candidates: List[str],
    *,
    left_context: str,
    right_context: str,
    seed: int,
) -> Tuple[str, List[str]]:
    if not candidates:
        return "", ["No precomputed candidates available."]
    start = seed % len(candidates)
    for offset in range(len(candidates)):
        seq = candidates[(start + offset) % len(candidates)]
        if _passes_stuffer_constraints(seq, left_context, right_context):
            return seq, []
    return "", ["No precomputed candidate passed junction checks."]


def design_stuffer_sequence(
    length: int,
    *,
    left_context: str = "",
    right_context: str = "",
    seed: int = 0,
) -> Tuple[str, List[str]]:
    """
    Look up a precomputed safe stuffer sequence of *length*.

    For lengths above the library maximum, concatenates precomputed chunks.
    """
    if length <= 0:
        return "", []

    max_length, library = load_stuffer_library()
    warnings: List[str] = []

    if max_length <= 0:
        warnings.append("Stuffer sequence library not found; leaving N placeholders.")
        return "N" * length, warnings

    if length <= max_length and length in library:
        seq, pick_warnings = _pick_from_candidates(
            library[length],
            left_context=left_context,
            right_context=right_context,
            seed=seed,
        )
        warnings.extend(pick_warnings)
        if seq:
            return seq, warnings
        return "N" * length, warnings + ["Lookup failed; leaving N placeholders."]

    if length > max_length:
        parts: List[str] = []
        remaining = length
        part_idx = 0
        chunk_left = left_context
        chunk_right = ""
        while remaining > 0:
            chunk_len = min(remaining, max_length)
            if remaining <= max_length:
                chunk_right = right_context
            chunk, chunk_warnings = design_stuffer_sequence(
                chunk_len,
                left_context=chunk_left,
                right_context=chunk_right,
                seed=seed + part_idx,
            )
            warnings.extend(chunk_warnings)
            if "N" in chunk:
                return "N" * length, warnings + [
                    f"Could not compose stuffer longer than {max_length} nt."
                ]
            parts.append(chunk)
            remaining -= chunk_len
            part_idx += 1
            chunk_left = ""
        combined = "".join(parts)
        if _passes_stuffer_constraints(combined, left_context, right_context):
            return combined, warnings
        return "N" * length, warnings + ["Concatenated stuffer failed junction checks."]

    warnings.append(f"No precomputed stuffer for length {length}.")
    return "N" * length, warnings


def apply_stuffer_design(library: Dict[str, Any]) -> Dict[str, Any]:
    """Replace N stuffer tails in a packed library with precomputed ATGC sequences."""
    result = {**library, "oligos": [dict(o) for o in (library.get("oligos") or [])]}
    all_warnings: List[Dict[str, Any]] = []
    errors: List[Dict[str, Any]] = []

    for oligo in result["oligos"]:
        idx = int(oligo.get("global_oligo_index") or oligo.get("oligo_index") or 0)
        fwd = oligo.get("forward_primer") or ""
        rev_rc = oligo.get("reverse_primer_rc") or ""
        insert = oligo.get("merged_insert") or ""

        s5_len = int(oligo.get("stuffer_5_length") or 0)
        s3_len = int(oligo.get("stuffer_3_length") or 0)
        if s5_len == 0 and s3_len == 0:
            total = int(oligo.get("stuffer_length") or 0)
            s5_len, s3_len = split_stuffer_lengths(total)

        stuffer_5, w5 = design_stuffer_sequence(
            s5_len,
            left_context="",
            right_context=fwd[:15],
            seed=idx * 1000 + 1,
        )
        stuffer_3, w3 = design_stuffer_sequence(
            s3_len,
            left_context=rev_rc[-15:],
            right_context="",
            seed=idx * 1000 + 2,
        )

        oligo_warnings = w5 + w3
        if "N" in stuffer_5 or "N" in stuffer_3:
            errors.append(
                {
                    "oligo": idx,
                    "error": "Could not fully replace stuffer Ns.",
                }
            )

        oligo["stuffer_5"] = stuffer_5
        oligo["stuffer_3"] = stuffer_3
        oligo["stuffer_5_length"] = len(stuffer_5)
        oligo["stuffer_3_length"] = len(stuffer_3)
        oligo["stuffer"] = stuffer_5 + stuffer_3
        oligo["stuffer_length"] = len(stuffer_5) + len(stuffer_3)
        oligo["stuffer_designed"] = True
        oligo["stuffer_gc_pct"] = gc_fraction(stuffer_5 + stuffer_3)
        oligo["sequence"] = assemble_oligo_with_outside_stuffers(
            stuffer_5, fwd, insert, rev_rc, stuffer_3
        )
        oligo["length"] = len(oligo["sequence"])

        if oligo_warnings:
            all_warnings.append({"oligo": idx, "warnings": oligo_warnings})

    result["stuffer_design_warnings"] = all_warnings
    result["stuffer_design_errors"] = errors
    result["unique_oligo_sequences"] = len({o["sequence"] for o in result["oligos"]})
    return result
