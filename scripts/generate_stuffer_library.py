#!/usr/bin/env python3
"""
Generate precomputed safe stuffer sequences (lengths 1–100) for oligopool padding.

Run from repo root:
    python scripts/generate_stuffer_library.py
"""

from __future__ import annotations

import random
import re
from pathlib import Path

import yaml

MAX_LENGTH = 100
CANDIDATES_PER_LENGTH = 32
GC_MIN = 0.35
GC_MAX = 0.65
OUT_PATH = Path(__file__).resolve().parent.parent / "data" / "stuffer_sequences.yaml"

_BSAI_MOTIFS = ("GGTCTC", "GAGACC")
_HOMOPOLYMER_RUN = re.compile(r"(.)\1{3,}", re.IGNORECASE)
_COMPLEMENT = {"A": "T", "T": "A", "G": "C", "C": "G"}
_BASES = ("A", "T", "G", "C")


def reverse_complement(dna: str) -> str:
    return "".join(_COMPLEMENT.get(b, "N") for b in reversed(dna.upper()))


def gc_fraction(seq: str) -> float:
    return sum(b in "GC" for b in seq.upper()) / len(seq) if seq else 0.0


def has_homopolymer(seq: str) -> bool:
    return bool(_HOMOPOLYMER_RUN.search(seq))


def has_bsai(seq: str) -> bool:
    seq = seq.upper()
    rc = reverse_complement(seq)
    return any(m in s for s in (seq, rc) for m in _BSAI_MOTIFS)


def has_hairpin(seq: str, min_stem: int = 6, min_loop: int = 3) -> bool:
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
                j = i + stem_len + loop
                stem3 = seq[j : j + stem_len]
                if stem5 == reverse_complement(stem3):
                    return True
    return False


def is_safe(seq: str) -> bool:
    if not seq:
        return False
    if has_homopolymer(seq):
        return False
    if has_bsai(seq):
        return False
    if len(seq) >= 15 and has_hairpin(seq):
        return False
    if len(seq) >= 4 and not (GC_MIN <= gc_fraction(seq) <= GC_MAX):
        return False
    return True


def partial_ok(trial: str) -> bool:
    if has_homopolymer(trial):
        return False
    if len(trial) >= 6 and has_bsai(trial[-6:]):
        return False
    return True


def greedy_build(length: int, rng: random.Random, first_base: str) -> str | None:
    target_gc = length // 2
    chars = [first_base]
    while len(chars) < length:
        gc_count = sum(b in "GC" for b in chars)
        pref = ["G", "C", "A", "T"] if gc_count < target_gc else ["A", "T", "G", "C"]
        bases = pref + [b for b in _BASES if b not in pref]
        rng.shuffle(bases)
        placed = False
        for base in bases:
            if len(chars) >= 3 and chars[-1] == chars[-2] == chars[-3] == base:
                continue
            trial = "".join(chars) + base
            if not partial_ok(trial):
                continue
            chars.append(base)
            placed = True
            break
        if not placed:
            return None
    seq = "".join(chars)
    return seq if is_safe(seq) else None


def short_candidates(length: int) -> list[str]:
    if length == 1:
        return list(_BASES)
    out = []
    for a in _BASES:
        for b in _BASES:
            s = a + b
            if is_safe(s) or (len(s) < 4 and not has_homopolymer(s) and not has_bsai(s)):
                out.append(s)
    return list(dict.fromkeys(out))


def generate_for_length(length: int, n: int, rng: random.Random) -> list[str]:
    if length <= 2:
        return short_candidates(length)[:n]
    seen: set[str] = set()
    out: list[str] = []
    for attempt in range(n * 100):
        if len(out) >= n:
            break
        first = _BASES[(attempt + len(out)) % 4]
        seq = greedy_build(length, rng, first)
        if seq and seq not in seen:
            seen.add(seq)
            out.append(seq)
    return out


def main() -> None:
    rng = random.Random(42)
    sequences: dict[str, list[str]] = {}
    shortfalls = []

    for length in range(1, MAX_LENGTH + 1):
        candidates = generate_for_length(length, CANDIDATES_PER_LENGTH, rng)
        sequences[str(length)] = candidates
        if len(candidates) < CANDIDATES_PER_LENGTH:
            shortfalls.append(length)
        print(f"length {length:3d}: {len(candidates)} candidates")

    payload = {
        "version": 1,
        "max_length": MAX_LENGTH,
        "candidates_per_length": CANDIDATES_PER_LENGTH,
        "constraints": {
            "bsai_motifs": list(_BSAI_MOTIFS),
            "max_homopolymer": 3,
            "gc_min": GC_MIN,
            "gc_max": GC_MAX,
        },
        "sequences": sequences,
    }

    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT_PATH, "w", encoding="utf-8") as handle:
        yaml.dump(payload, handle, default_flow_style=False, sort_keys=False, width=120)

    print(f"\nWrote {OUT_PATH}")
    if shortfalls:
        print(f"Warning: {len(shortfalls)} lengths under target count.")
    else:
        print("All lengths filled.")


if __name__ == "__main__":
    main()
