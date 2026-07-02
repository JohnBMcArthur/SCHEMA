"""
Golden Gate Assembly overhang selection from AA-pair compatibility data.
"""

from __future__ import annotations

import re
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, List, Optional

import yaml

from utils.config import GGA_COMPATIBILITY_YAML

# N-terminal Met start overhangs (4 bp); chosen last after junction assignment.
XATG_OVERHANGS = ("AATG", "TATG", "CATG", "GATG")
_XATG_SET = frozenset(XATG_OVERHANGS)

# 3+ consecutive identical bases (e.g. AAA in AAAT, CCC in GCCC) are rejected.
_HOMOPOLYMER_RUN = re.compile(r"(.)\1{2,}", re.IGNORECASE)


def _has_excessive_homopolymer(sequence: str) -> bool:
    return bool(_HOMOPOLYMER_RUN.search(sequence or ""))


def _normalize_candidate(raw: Dict[str, Any]) -> Dict[str, Any]:
    return {
        "overhang": str(raw.get("overhang", "")),
        "reverse_complement": str(raw.get("reverse_complement", "")),
        "frame": raw.get("frame"),
        "example_codon_pair": str(raw.get("example_codon_pair", "")),
        "efficiency": float(raw.get("efficiency", 0.0)),
        "count": int(raw.get("count", 0)),
        "compatible_pairs": list(raw.get("compatible_pairs") or []),
    }


@lru_cache(maxsize=2)
def load_gga_compatibility(yaml_path: str) -> Dict[str, Any]:
    """Load and cache the GGA AA-pair compatibility YAML."""
    path = Path(yaml_path)
    if not path.is_file():
        raise FileNotFoundError(f"GGA compatibility file not found: {path}")
    with open(path, encoding="utf-8") as handle:
        data = yaml.safe_load(handle)
    if not isinstance(data, dict) or "aa_pairs" not in data:
        raise ValueError(f"Invalid GGA compatibility YAML: {path}")
    return data


def get_sorted_candidates(
    compatibility: Dict[str, Any], aa_pair: str
) -> List[Dict[str, Any]]:
    """Return overhang candidates for an AA pair, highest efficiency first."""
    entry = compatibility.get("aa_pairs", {}).get(aa_pair.upper())
    if not entry:
        return []
    candidates = [
        _normalize_candidate(item)
        for item in (entry.get("candidate_overhangs") or [])
    ]
    candidates.sort(key=lambda c: c["efficiency"], reverse=True)
    return [
        c
        for c in candidates
        if not _has_excessive_homopolymer(c["overhang"])
        and not _has_excessive_homopolymer(c["reverse_complement"])
    ]


def select_overhang_for_aa_pair(
    compatibility: Dict[str, Any],
    aa_pair: str,
    used_overhangs: set[str],
) -> Dict[str, Any]:
    """
    Pick the best-efficiency overhang for an AA pair not already used elsewhere.

    Each overhang sequence may be assigned to at most one junction.
    """
    candidates = get_sorted_candidates(compatibility, aa_pair)
    if not candidates:
        return {
            "aa_pair": aa_pair.upper(),
            "overhang": None,
            "error": f"No eligible overhang candidates for AA pair {aa_pair.upper()}",
        }

    chosen = None
    chosen_rank = None
    for rank, candidate in enumerate(candidates):
        overhang_key = candidate["overhang"].upper()
        if overhang_key not in used_overhangs:
            chosen = candidate
            chosen_rank = rank
            break

    if chosen is None:
        return {
            "aa_pair": aa_pair.upper(),
            "overhang": None,
            "error": (
                f"All eligible overhangs for {aa_pair.upper()} are already used "
                f"at other junctions"
            ),
        }

    used_overhangs.add(chosen["overhang"].upper())
    result = {
        "aa_pair": aa_pair.upper(),
        "usage_rank": chosen_rank,
        "overhang": chosen["overhang"],
        "reverse_complement": chosen["reverse_complement"],
        "frame": chosen["frame"],
        "example_codon_pair": chosen["example_codon_pair"],
        "efficiency": chosen["efficiency"],
        "count": chosen["count"],
        "compatible_pairs": chosen["compatible_pairs"],
    }
    if chosen_rank is not None and chosen_rank > 0:
        result["warning"] = (
            f"Higher-efficiency overhang(s) for {aa_pair.upper()} were skipped "
            f"(homopolymer filter or already used)."
        )
    return result


def _reverse_complement(dna: str) -> str:
    comp = {"A": "T", "T": "A", "G": "C", "C": "G"}
    return "".join(comp.get(b, "N") for b in reversed((dna or "").upper()))


def apply_fragment1_prepend_m(fragments: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """Return a shallow copy of *fragments*, optionally prefixing fragment 1 with Met."""
    out = [dict(f) for f in fragments]
    if not out:
        return out
    seq = (out[0].get("sequence") or "").strip()
    if seq and seq[0].upper() != "M":
        new_seq = "M" + seq
        out[0] = {**out[0], "sequence": new_seq, "length": len(new_seq)}
    return out


def find_overhang_candidate(
    compatibility: Dict[str, Any],
    overhang: str,
) -> Optional[Dict[str, Any]]:
    """Look up the highest-efficiency YAML entry for a 4 bp overhang sequence."""
    target = (overhang or "").upper()
    if len(target) != 4:
        return None
    best: Optional[Dict[str, Any]] = None
    best_aa_pair: Optional[str] = None
    for aa_pair, entry in (compatibility.get("aa_pairs") or {}).items():
        for raw in entry.get("candidate_overhangs") or []:
            cand = _normalize_candidate(raw)
            if cand["overhang"].upper() != target:
                continue
            if best is None or cand["efficiency"] > best["efficiency"]:
                best = cand
                best_aa_pair = str(aa_pair).upper()
    if best is None:
        return None
    return {**best, "aa_pair": best_aa_pair}


def collect_xatg_candidates(compatibility: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Best YAML candidate per XATG overhang (synthetic fallback if absent)."""
    by_oh: Dict[str, Dict[str, Any]] = {}
    for aa_pair, entry in (compatibility.get("aa_pairs") or {}).items():
        for raw in entry.get("candidate_overhangs") or []:
            cand = _normalize_candidate(raw)
            oh = cand["overhang"].upper()
            if oh not in _XATG_SET:
                continue
            if _has_excessive_homopolymer(cand["overhang"]) or _has_excessive_homopolymer(
                cand["reverse_complement"]
            ):
                continue
            existing = by_oh.get(oh)
            if existing is None or cand["efficiency"] > existing["efficiency"]:
                by_oh[oh] = {**cand, "aa_pair": str(aa_pair).upper()}

    out: List[Dict[str, Any]] = []
    for oh in XATG_OVERHANGS:
        if oh in by_oh:
            out.append(by_oh[oh])
        else:
            out.append(
                {
                    "overhang": oh,
                    "reverse_complement": _reverse_complement(oh),
                    "frame": 2,
                    "example_codon_pair": "GCA-ATG",
                    "efficiency": 0.0,
                    "count": 0,
                    "compatible_pairs": [
                        "AM", "CM", "DM", "EM", "FM", "GM", "HM", "IM", "KM", "LM",
                        "MM", "NM", "PM", "QM", "RM", "SM", "TM", "VM", "WM", "YM",
                    ],
                    "aa_pair": "XM",
                }
            )
    out.sort(key=lambda c: float(c.get("efficiency") or 0), reverse=True)
    return out


def _selection_from_candidate(
    candidate: Dict[str, Any],
    *,
    used_overhangs: set[str],
    site: str,
    usage_rank: int = 0,
) -> Dict[str, Any]:
    oh = candidate["overhang"].upper()
    if oh in used_overhangs:
        return {
            "aa_pair": candidate.get("aa_pair"),
            "overhang": None,
            "site": site,
            "error": f"Overhang {oh} is already used at another site",
        }
    used_overhangs.add(oh)
    return {
        "aa_pair": candidate.get("aa_pair"),
        "usage_rank": usage_rank,
        "overhang": candidate["overhang"],
        "reverse_complement": candidate.get("reverse_complement"),
        "frame": candidate.get("frame"),
        "example_codon_pair": candidate.get("example_codon_pair"),
        "efficiency": candidate.get("efficiency"),
        "count": candidate.get("count"),
        "compatible_pairs": list(candidate.get("compatible_pairs") or []),
        "site": site,
    }


def select_best_xatg_overhang(
    compatibility: Dict[str, Any],
    used_overhangs: set[str],
) -> Dict[str, Any]:
    """Pick the best unused XATG overhang (assigned after junction overhangs)."""
    candidates = collect_xatg_candidates(compatibility)
    for rank, cand in enumerate(candidates):
        oh = cand["overhang"].upper()
        if oh not in used_overhangs:
            return _selection_from_candidate(
                cand, used_overhangs=used_overhangs, site="n_terminal_xatg", usage_rank=rank
            )
    return {
        "aa_pair": "XM",
        "overhang": None,
        "site": "n_terminal_xatg",
        "error": "All XATG overhangs (AATG, TATG, CATG, GATG) are already used",
    }


def select_fragment1_overhang_5(
    compatibility: Dict[str, Any],
    fragment1_sequence: str,
    used_overhangs: set[str],
    *,
    manual_overhang: Optional[str] = None,
) -> Optional[Dict[str, Any]]:
    """
    Assign fragment 1 5′ overhang after junctions are placed.

    - Starts with M: auto best unused XATG, or *manual_overhang* if provided.
    - Otherwise with *manual_overhang*: vector-matched overhang from the table.
    - Otherwise: no 5′ overhang (same as other fragments without a 5′ junction).
    """
    seq = (fragment1_sequence or "").strip()
    if not seq:
        return None

    if seq[0].upper() == "M":
        if manual_overhang:
            manual = manual_overhang.upper()
            if manual not in _XATG_SET:
                return {
                    "aa_pair": "XM",
                    "overhang": None,
                    "site": "n_terminal_xatg",
                    "error": f"Manual overhang must be one of {', '.join(XATG_OVERHANGS)}",
                }
            cand = find_overhang_candidate(compatibility, manual)
            if cand is None:
                synth = next(c for c in collect_xatg_candidates(compatibility) if c["overhang"] == manual)
                cand = synth
            return _selection_from_candidate(
                cand, used_overhangs=used_overhangs, site="n_terminal_xatg"
            )
        return select_best_xatg_overhang(compatibility, used_overhangs)

    if manual_overhang:
        cand = find_overhang_candidate(compatibility, manual_overhang)
        if cand is None:
            return {
                "aa_pair": None,
                "overhang": None,
                "site": "n_terminal_vector",
                "error": f"Overhang {manual_overhang.upper()} not found in compatibility table",
            }
        return _selection_from_candidate(
            cand, used_overhangs=used_overhangs, site="n_terminal_vector"
        )

    return None


def assembly_gga_options_from_mapping(mapping: Dict[str, Any]) -> Dict[str, Any]:
    """Read fragment-1 GGA options from session state or similar mapping."""
    manual = (mapping.get("assembly_fragment1_manual_overhang") or "").strip().upper()
    return {
        "fragment1_prepend_m": bool(mapping.get("assembly_fragment1_prepend_m")),
        "fragment1_manual_overhang": manual or None,
    }


def build_fragment_terminal_filters(
    assembly_fragments: List[Dict[str, Any]],
) -> List[Dict[str, Any]]:
    """
    Derive allowed N- and C-terminal residues per fragment from assigned overhangs.

    For junction overhang compatible pair ``XY``, fragment left of the junction
    must end with ``X``; fragment right must start with ``Y``.
    """
    filters: List[Dict[str, Any]] = []
    for frag in assembly_fragments:
        oh5 = frag.get("overhang_5")
        oh3 = frag.get("overhang_3")
        allowed_n = None
        allowed_c = None
        if oh5 and oh5.get("compatible_pairs"):
            allowed_n = sorted(
                {str(p)[1] for p in oh5["compatible_pairs"] if len(str(p)) == 2}
            )
        if oh3 and oh3.get("compatible_pairs"):
            allowed_c = sorted(
                {str(p)[0] for p in oh3["compatible_pairs"] if len(str(p)) == 2}
            )
        filters.append(
            {
                "fragment": frag["index"],
                "query_sequence": frag.get("sequence", ""),
                "overhang_5": oh5.get("overhang") if oh5 else None,
                "overhang_3": oh3.get("overhang") if oh3 else None,
                "allowed_n_terminal": allowed_n,
                "allowed_c_terminal": allowed_c,
            }
        )
    return filters


def fragment_passes_terminal_filter(
    sequence: str,
    terminal_filter: Dict[str, Any],
) -> bool:
    """Return True if sequence terminals satisfy GGA compatible-pair constraints."""
    if not sequence:
        return False
    allowed_n = terminal_filter.get("allowed_n_terminal")
    allowed_c = terminal_filter.get("allowed_c_terminal")
    if allowed_n is not None and sequence[0] not in allowed_n:
        return False
    if allowed_c is not None and sequence[-1] not in allowed_c:
        return False
    return True


def assign_golden_gate_overhangs(
    fragments: List[Dict[str, Any]],
    compatibility: Optional[Dict[str, Any]] = None,
    yaml_path: Optional[Path] = None,
    *,
    fragment1_prepend_m: bool = False,
    fragment1_manual_overhang: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Assign 5' and 3' Golden Gate overhangs to each query fragment.

    Rules:
    - Junction overhangs are chosen from the AA pair at each boundary
      (last residue of left fragment + first residue of right fragment).
    - Fragment 1 5′ overhang is assigned **after** junctions:
        - If fragment 1 starts with M (natural or prepended): best unused XATG
          (AATG, TATG, CATG, GATG), or *fragment1_manual_overhang* when set.
        - If not M and *fragment1_manual_overhang* is set: that vector overhang.
        - Otherwise fragment 1 has no 5′ overhang.
    - The last fragment has no 3' overhang.
    - The 3' overhang of fragment N equals the 5' overhang of fragment N+1.
    """
    path = yaml_path or GGA_COMPATIBILITY_YAML
    compat = compatibility or load_gga_compatibility(str(path.resolve()))

    working = [dict(f) for f in fragments]
    if fragment1_prepend_m:
        working = apply_fragment1_prepend_m(working)

    used_overhangs: set[str] = set()
    junctions: List[Dict[str, Any]] = []

    for idx in range(len(working) - 1):
        left = working[idx].get("sequence") or ""
        right = working[idx + 1].get("sequence") or ""
        junction_num = idx + 1

        if not left or not right:
            junctions.append(
                {
                    "junction": junction_num,
                    "left_fragment": working[idx]["index"],
                    "right_fragment": working[idx + 1]["index"],
                    "aa_pair": None,
                    "overhang": None,
                    "error": "Empty fragment sequence at junction",
                }
            )
            continue

        aa_pair = f"{left[-1]}{right[0]}".upper()
        selection = select_overhang_for_aa_pair(compat, aa_pair, used_overhangs)
        junctions.append(
            {
                "junction": junction_num,
                "left_fragment": working[idx]["index"],
                "right_fragment": working[idx + 1]["index"],
                **selection,
            }
        )

    frag1_seq = (working[0].get("sequence") or "") if working else ""
    manual_oh = (fragment1_manual_overhang or "").strip().upper() or None
    frag1_oh5 = select_fragment1_overhang_5(
        compat,
        frag1_seq,
        used_overhangs,
        manual_overhang=manual_oh,
    )

    enriched: List[Dict[str, Any]] = []
    for idx, frag in enumerate(working):
        out = dict(frag)
        if idx == 0:
            out["overhang_5"] = frag1_oh5
        else:
            out["overhang_5"] = junctions[idx - 1]

        if idx == len(working) - 1:
            out["overhang_3"] = None
        else:
            out["overhang_3"] = junctions[idx]

        enriched.append(out)

    return {
        "fragments": enriched,
        "junctions": junctions,
        "fragment1_n_terminal": frag1_oh5,
        "fragment1_prepend_m": fragment1_prepend_m,
        "yaml_path": str(path),
        "metadata": compat.get("metadata", {}),
    }
