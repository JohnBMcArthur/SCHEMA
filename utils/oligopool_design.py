"""
Oligopool design: BsaI-flanked fragment inserts and full oligo assembly.

Stage 1 — per fragment / homolog:
    [BsaI 5′ cassette] + [coding DNA] + [BsaI 3′ cassette]

Stage 2 — pack consecutive Stage 1 pieces into minimal full-length oligos:
    [5′ stuffer] + forward_primer + merged_pieces + reverse_complement(reverse_primer) + [3′ stuffer]
    (each oligo padded to max length when possible; stuffers sit outside primer binding)
"""

from __future__ import annotations

import itertools
import re
from typing import Any, Dict, Iterable, List, Optional, Tuple

from Bio.Data import CodonTable
from Bio.Seq import Seq

from utils.assembly_analysis import (
    get_aligned_query_sequence,
    split_query_sequence_into_fragments,
)
from utils.config import GGA_COMPATIBILITY_YAML
from utils.gga_overhangs import (
    assign_golden_gate_overhangs,
    assembly_gga_options_from_mapping,
    load_gga_compatibility,
)
from utils.stuffer_design import (
    assemble_oligo_with_outside_stuffers,
    build_oligo_stuffer_fields,
)

# NEB BsaI-HFv2 Golden Gate sites (recognition + 1 bp context; no throwaway spacers).
BSAI_START = "GGTCTCC"
BSAI_END = "CGAGACC"

_STANDARD = CodonTable.unambiguous_dna_by_name["Standard"]
_AA_TO_CODONS: Dict[str, List[str]] = {}
for _codon, _aa in _STANDARD.forward_table.items():
    _AA_TO_CODONS.setdefault(_aa, []).append(_codon.upper())

# E. coli high-expression preferred codons (interior residues).
ECOLI_PREFERRED: Dict[str, str] = {
    "A": "GCG",
    "C": "TGC",
    "D": "GAT",
    "E": "GAA",
    "F": "TTC",
    "G": "GGC",
    "H": "CAT",
    "I": "ATT",
    "K": "AAA",
    "L": "CTG",
    "M": "ATG",
    "N": "AAC",
    "P": "CCG",
    "Q": "CAG",
    "R": "CGC",
    "S": "AGC",
    "T": "ACC",
    "V": "GTG",
    "W": "TGG",
    "Y": "TAC",
    "*": "TAA",
}


def _clean_aa(sequence: str) -> str:
    return re.sub(r"[^A-Za-z*]", "", (sequence or "")).upper()


def _junction_info(junc: Optional[Dict[str, Any]]) -> Tuple[Optional[str], int]:
    if not junc or not junc.get("overhang"):
        return None, 0
    frame = junc.get("frame")
    position = int(frame) if frame is not None else 0
    return str(junc["overhang"]).upper(), position


def _pick_codon(aa: str, pattern: str) -> str:
    aa = aa.upper()
    preferred = ECOLI_PREFERRED.get(aa)
    candidates = list(_AA_TO_CODONS.get(aa, ["NNN"]))
    if preferred and preferred in candidates:
        candidates.remove(preferred)
        candidates.insert(0, preferred)
    for codon in candidates:
        if re.fullmatch(pattern.upper(), codon):
            return codon
    return preferred or candidates[0]


def junction_codon_segment(
    junction: Dict[str, Any],
    aa_left: str,
    aa_right: str,
    side: str,
) -> str:
    """
    DNA segment encoding the junction spanning codons for (aa_left, aa_right).

    *side* is ``'front'`` (5′ of fragment) or ``'back'`` (3′ of fragment).
    """
    oh_seq, position = _junction_info(junction)
    if not oh_seq:
        return _pick_codon(aa_left if side == "back" else aa_right, "...")

    overall_pattern = "." * position + oh_seq + "." * (2 - position)
    pattern1 = overall_pattern[:3]
    pattern2 = overall_pattern[3:]
    c1 = _pick_codon(aa_left.upper(), pattern1)
    c2 = _pick_codon(aa_right.upper(), pattern2)
    cdn_seq = c1 + c2

    if side == "front":
        return cdn_seq[position:]
    if position == 2:
        return cdn_seq
    return cdn_seq[: -2 + position] if (-2 + position) != 0 else cdn_seq


def best_ecoli_codon(aa: str) -> str:
    aa = aa.upper()
    return ECOLI_PREFERRED.get(aa, _AA_TO_CODONS.get(aa, ["NNN"])[0])


def _n_terminal_left_aa(oh5: Dict[str, Any], first_aa: str) -> str:
    """Upstream partner AA for fragment 1 5′ junction codon encoding."""
    first_aa = first_aa.upper()
    for pair in oh5.get("compatible_pairs") or []:
        ps = str(pair).upper()
        if len(ps) == 2 and ps[1] == first_aa:
            return ps[0]
    example = str(oh5.get("example_codon_pair") or "")
    if "-" in example:
        left_codon = example.split("-", 1)[0].strip().upper()
        if len(left_codon) == 3:
            from Bio.Data import CodonTable

            aa = CodonTable.unambiguous_dna_by_name["Standard"].forward_table.get(left_codon)
            if aa:
                return aa
    return "G"


def build_fragment_coding_dna(
    aa_sequence: str,
    frag_index: int,
    query_fragments: List[Dict[str, Any]],
    enriched_fragments: List[Dict[str, Any]],
) -> Tuple[str, str]:
    """
    Return (codon display string, coding DNA) for one homolog fragment.

    Interior residues use E. coli–preferred codons; terminal segments follow
    Assembly Analysis junction overhangs and reading frame.
    """
    aa = _clean_aa(aa_sequence)
    if not aa:
        return "", ""

    idx = frag_index - 1
    n = len(aa)
    is_first = idx == 0
    is_last = idx == len(query_fragments) - 1
    oh5_raw = enriched_fragments[idx].get("overhang_5")
    oh5 = oh5_raw if oh5_raw and oh5_raw.get("overhang") else None
    oh3_raw = None if is_last else enriched_fragments[idx].get("overhang_3")
    oh3 = oh3_raw if oh3_raw and oh3_raw.get("overhang") else None

    dna_parts: List[str] = []
    codon_labels: List[str] = []

    if n == 1:
        if oh5 and oh3:
            left = (
                _n_terminal_left_aa(oh5, aa[0])
                if is_first
                else query_fragments[idx - 1]["sequence"][-1]
            )
            right = query_fragments[idx + 1]["sequence"][0]
            seg5 = junction_codon_segment(oh5, left, aa[0], "front")
            seg3 = junction_codon_segment(oh3, aa[0], right, "back")
            dna_parts.extend([seg5, seg3])
            codon_labels.append(f"5′:{seg5} 3′:{seg3}")
        elif oh5:
            left = (
                _n_terminal_left_aa(oh5, aa[0])
                if is_first
                else query_fragments[idx - 1]["sequence"][-1]
            )
            seg = junction_codon_segment(oh5, left, aa[0], "front")
            dna_parts.append(seg)
            codon_labels.append(f"5′:{seg}")
        elif oh3:
            right = query_fragments[idx + 1]["sequence"][0]
            seg = junction_codon_segment(oh3, aa[0], right, "back")
            dna_parts.append(seg)
            codon_labels.append(f"3′:{seg}")
        else:
            c = best_ecoli_codon(aa[0])
            dna_parts.append(c)
            codon_labels.append(c)
    else:
        if oh5:
            if is_first:
                left = _n_terminal_left_aa(oh5, aa[0])
            else:
                left = query_fragments[idx - 1]["sequence"][-1]
            seg = junction_codon_segment(oh5, left, aa[0], "front")
            dna_parts.append(seg)
            codon_labels.append(f"5′:{seg}")
        else:
            c0 = best_ecoli_codon(aa[0])
            dna_parts.append(c0)
            codon_labels.append(c0)

        for i in range(1, n - 1):
            ci = best_ecoli_codon(aa[i])
            dna_parts.append(ci)
            codon_labels.append(ci)

        if oh3:
            right = query_fragments[idx + 1]["sequence"][0]
            seg = junction_codon_segment(oh3, aa[-1], right, "back")
            dna_parts.append(seg)
            codon_labels.append(f"3′:{seg}")
        else:
            cn = best_ecoli_codon(aa[-1])
            dna_parts.append(cn)
            codon_labels.append(cn)

    return " ".join(codon_labels), "".join(dna_parts).upper()


def build_bsaI_flanked_piece(coding_dna: str) -> str:
    """Wrap coding DNA with BsaI Golden Gate cassettes."""
    return f"{BSAI_START}{coding_dna}{BSAI_END}".upper()


def merge_pieces_at_overhang(left: str, right: str, overhang: str) -> str:
    """Concatenate adjacent fragment pieces, merging the shared 4 bp overhang."""
    oh = (overhang or "").upper()
    if not oh:
        return left + right
    lu, ru = left.upper(), right.upper()
    if lu.endswith(oh) and ru.startswith(oh):
        return left + right[len(oh) :]
    return left + right


def junction_overhangs_between_fragments(
    enriched_fragments: List[Dict[str, Any]],
) -> List[str]:
    """Overhang between fragment i and i+1 (length n_fragments - 1)."""
    out: List[str] = []
    for frag in enriched_fragments[:-1]:
        oh3 = frag.get("overhang_3") or {}
        out.append(str(oh3.get("overhang") or ""))
    return out


def reverse_complement(dna: str) -> str:
    """Return the reverse complement of a DNA string."""
    if not dna:
        return ""
    return str(Seq(dna.upper()).reverse_complement())


def merge_piece_range(
    pieces: List[str],
    junction_overhangs: List[str],
    start: int,
    end: int,
) -> str:
    """Merge consecutive pieces from *start* through *end* (inclusive, 0-based)."""
    if start > end or start < 0 or end >= len(pieces):
        raise ValueError("Invalid piece range for merge")
    merged = pieces[start]
    for idx in range(start + 1, end + 1):
        oh = junction_overhangs[idx - 1] if idx - 1 < len(junction_overhangs) else ""
        merged = merge_pieces_at_overhang(merged, pieces[idx], oh)
    return merged


def pack_oligos(
    pieces: List[str],
    junction_overhangs: List[str],
    *,
    forward_primer: str = "",
    reverse_primer: str = "",
    max_oligo_length: int,
    pad_to_max: bool = True,
    stuffer_base: str = "N",
) -> Dict[str, Any]:
    """
    Pack consecutive Stage 1 pieces into the fewest oligos ≤ *max_oligo_length*.

    Each oligo layout:
        [5′ stuffer] + forward_primer + merged_insert + reverse_complement(reverse_primer) + [3′ stuffer]

    When *pad_to_max* is True, outside stuffer tails fill unused space so each oligo
    is exactly *max_oligo_length* (when primers fit within the limit).
    """
    fwd = (forward_primer or "").upper()
    rev_rc = reverse_complement(reverse_primer or "")
    overhead = len(fwd) + len(rev_rc)
    payload_capacity = max_oligo_length - overhead

    if payload_capacity < 0:
        return {
            "oligos": [],
            "n_oligos": 0,
            "payload_capacity": payload_capacity,
            "error": (
                f"Primer sequences ({overhead} nt) exceed maximum oligo length "
                f"({max_oligo_length} nt)."
            ),
        }

    if not pieces:
        padding = payload_capacity if pad_to_max else 0
        stuffer_fields = build_oligo_stuffer_fields(padding, stuffer_base=stuffer_base)
        seq = assemble_oligo_with_outside_stuffers(
            stuffer_fields["stuffer_5"], fwd, "", rev_rc, stuffer_fields["stuffer_3"]
        )
        return {
            "oligos": [
                {
                    "oligo_index": 1,
                    "fragment_indices": [],
                    "fragment_start": None,
                    "fragment_end": None,
                    "merged_insert": "",
                    "insert_length": 0,
                    **stuffer_fields,
                    "forward_primer": fwd,
                    "reverse_primer_rc": rev_rc,
                    "sequence": seq,
                    "length": len(seq),
                }
            ],
            "n_oligos": 1,
            "payload_capacity": payload_capacity,
            "error": None,
        }

    oligos: List[Dict[str, Any]] = []
    idx = 0
    n = len(pieces)

    while idx < n:
        if len(pieces[idx]) > payload_capacity:
            return {
                "oligos": oligos,
                "n_oligos": len(oligos),
                "payload_capacity": payload_capacity,
                "error": (
                    f"Fragment {idx + 1} piece ({len(pieces[idx])} nt) exceeds "
                    f"payload capacity ({payload_capacity} nt) after primers."
                ),
            }

        end = idx
        merged = pieces[idx]
        while end + 1 < n:
            candidate = merge_piece_range(pieces, junction_overhangs, idx, end + 1)
            if len(candidate) <= payload_capacity:
                merged = candidate
                end += 1
            else:
                break

        padding = (payload_capacity - len(merged)) if pad_to_max else 0
        stuffer_fields = build_oligo_stuffer_fields(padding, stuffer_base=stuffer_base)
        sequence = assemble_oligo_with_outside_stuffers(
            stuffer_fields["stuffer_5"],
            fwd,
            merged,
            rev_rc,
            stuffer_fields["stuffer_3"],
        )

        oligos.append(
            {
                "oligo_index": len(oligos) + 1,
                "fragment_indices": list(range(idx + 1, end + 2)),
                "fragment_start": idx + 1,
                "fragment_end": end + 1,
                "merged_insert": merged,
                "insert_length": len(merged),
                **stuffer_fields,
                "forward_primer": fwd,
                "reverse_primer_rc": rev_rc,
                "sequence": sequence,
                "length": len(sequence),
            }
        )
        idx = end + 1

    return {
        "oligos": oligos,
        "n_oligos": len(oligos),
        "payload_capacity": payload_capacity,
        "error": None,
    }


def build_all_unique_stage1_pieces(
    blocks: List[Dict[str, Any]],
    query_fragments: List[Dict[str, Any]],
    enriched_fragments: List[Dict[str, Any]],
) -> List[Dict[str, Any]]:
    """Build Stage 1 inserts for every saved homolog at every fragment site."""
    pieces: List[Dict[str, Any]] = []
    seen: set[tuple[int, str]] = set()
    for block in blocks:
        frag_index = int(block["fragment"])
        for variant in block.get("variants") or []:
            seq_id = str(variant.get("sequence_id") or "unknown")
            key = (frag_index, seq_id)
            if key in seen:
                continue
            seen.add(key)
            insert = build_fragment_insert(
                variant.get("sequence") or "",
                seq_id,
                frag_index,
                query_fragments,
                enriched_fragments,
            )
            pieces.append(
                {
                    **insert,
                    "fragment_index": frag_index,
                    "is_query": bool(variant.get("is_query")),
                }
            )
    return pieces


def pack_block_oligos(
    blocks: List[Dict[str, Any]],
    query_fragments: List[Dict[str, Any]],
    enriched_fragments: List[Dict[str, Any]],
    *,
    forward_primer: str = "",
    reverse_primer: str = "",
    max_oligo_length: int = 300,
    pad_to_max: bool = True,
    stuffer_base: str = "N",
) -> Dict[str, Any]:
    """
    Pack each unique Stage 1 block insert once into the fewest oligos.

    Block pieces are concatenated back-to-back (each retains full BsaI flanks).
    Uses first-fit decreasing bin packing on piece lengths.
    """
    fwd = (forward_primer or "").upper()
    rev_rc = reverse_complement(reverse_primer or "")
    overhead = len(fwd) + len(rev_rc)
    payload_capacity = max_oligo_length - overhead

    stage1_pieces = build_all_unique_stage1_pieces(
        blocks, query_fragments, enriched_fragments
    )

    if payload_capacity < 0:
        return {
            "oligos": [],
            "n_oligos": 0,
            "total_oligos": 0,
            "total_blocks": len(stage1_pieces),
            "unique_oligo_sequences": 0,
            "payload_capacity": payload_capacity,
            "errors": [],
            "stage1_pieces": stage1_pieces,
            "error": (
                f"Primer sequences ({overhead} nt) exceed maximum oligo length "
                f"({max_oligo_length} nt)."
            ),
        }

    if not stage1_pieces:
        return {
            "oligos": [],
            "n_oligos": 0,
            "total_oligos": 0,
            "total_blocks": 0,
            "unique_oligo_sequences": 0,
            "payload_capacity": payload_capacity,
            "errors": [],
            "stage1_pieces": [],
            "error": "No Stage 1 block pieces to pack.",
        }

    sorted_pieces = sorted(
        stage1_pieces,
        key=lambda p: int(p.get("piece_length") or 0),
        reverse=True,
    )
    bins: List[Dict[str, Any]] = []
    errors: List[Dict[str, Any]] = []

    for piece in sorted_pieces:
        piece_len = int(piece.get("piece_length") or len(piece.get("piece_dna") or ""))
        if piece_len > payload_capacity:
            errors.append(
                {
                    "block": f"F{piece['fragment_index']}:{piece['sequence_id']}",
                    "error": (
                        f"Block piece ({piece_len} nt) exceeds payload capacity "
                        f"({payload_capacity} nt) after primers."
                    ),
                }
            )
            continue

        placed = False
        for bin_entry in bins:
            if bin_entry["used"] + piece_len <= payload_capacity:
                bin_entry["pieces"].append(piece)
                bin_entry["used"] += piece_len
                placed = True
                break
        if not placed:
            bins.append({"pieces": [piece], "used": piece_len})

    oligos: List[Dict[str, Any]] = []
    for bin_idx, bin_entry in enumerate(bins, start=1):
        block_pieces = bin_entry["pieces"]
        merged_insert = "".join(p["piece_dna"] for p in block_pieces)
        padding = (payload_capacity - len(merged_insert)) if pad_to_max else 0
        stuffer_fields = build_oligo_stuffer_fields(padding, stuffer_base=stuffer_base)
        sequence = assemble_oligo_with_outside_stuffers(
            stuffer_fields["stuffer_5"],
            fwd,
            merged_insert,
            rev_rc,
            stuffer_fields["stuffer_3"],
        )
        block_labels = ", ".join(
            f"F{p['fragment_index']}:{p['sequence_id']}" for p in block_pieces
        )
        oligos.append(
            {
                "oligo_index": bin_idx,
                "global_oligo_index": bin_idx,
                "blocks": [
                    {
                        "fragment_index": p["fragment_index"],
                        "sequence_id": p["sequence_id"],
                        "piece_dna": p["piece_dna"],
                        "piece_length": p["piece_length"],
                        "is_query": p.get("is_query"),
                    }
                    for p in block_pieces
                ],
                "block_labels": block_labels,
                "n_blocks": len(block_pieces),
                "merged_insert": merged_insert,
                "insert_length": len(merged_insert),
                **stuffer_fields,
                "forward_primer": fwd,
                "reverse_primer_rc": rev_rc,
                "sequence": sequence,
                "length": len(sequence),
            }
        )

    unique_sequences = len({o["sequence"] for o in oligos})
    return {
        "oligos": oligos,
        "n_oligos": len(oligos),
        "total_oligos": len(oligos),
        "total_blocks": len(stage1_pieces),
        "blocks_packed": len(stage1_pieces) - len(errors),
        "unique_oligo_sequences": unique_sequences,
        "payload_capacity": payload_capacity,
        "errors": errors,
        "stage1_pieces": stage1_pieces,
        "error": None,
    }


def _chimera_id_from_variants(variant_by_fragment: Dict[int, Dict[str, Any]]) -> str:
    parts = [
        f"F{frag}-{variant_by_fragment[frag]['sequence_id']}"
        for frag in sorted(variant_by_fragment)
    ]
    return "__".join(parts)


def pack_library_oligos(
    blocks: List[Dict[str, Any]],
    query_fragments: List[Dict[str, Any]],
    enriched_fragments: List[Dict[str, Any]],
    *,
    forward_primer: str = "",
    reverse_primer: str = "",
    max_oligo_length: int = 300,
    pad_to_max: bool = True,
    max_chimeras: Optional[int] = None,
) -> Dict[str, Any]:
    """
    Pack oligos for every chimera (one homolog choice per fragment site).

    Each chimera uses one variant from each block pool; all combinations are
    processed so every saved homolog piece appears in some packed oligo design.
    """
    if not blocks:
        return {
            "oligos": [],
            "chimera_count": 0,
            "total_oligos": 0,
            "unique_oligo_sequences": 0,
            "payload_capacity": 0,
            "errors": [],
            "stage1_pieces": [],
            "error": "No fragment blocks to pack.",
        }

    junction_oh = junction_overhangs_between_fragments(enriched_fragments)
    stage1_pieces = build_all_unique_stage1_pieces(
        blocks, query_fragments, enriched_fragments
    )

    variant_lists = [block.get("variants") or [] for block in blocks]
    if any(not variants for variants in variant_lists):
        return {
            "oligos": [],
            "chimera_count": 0,
            "total_oligos": 0,
            "unique_oligo_sequences": 0,
            "payload_capacity": 0,
            "errors": [],
            "stage1_pieces": stage1_pieces,
            "error": "Every fragment block must have at least one variant.",
        }

    all_oligos: List[Dict[str, Any]] = []
    errors: List[Dict[str, Any]] = []
    payload_capacity: Optional[int] = None
    chimera_count = 0
    global_oligo_idx = 0

    for combo in itertools.product(*variant_lists):
        if max_chimeras is not None and chimera_count >= max_chimeras:
            break

        variant_by_fragment: Dict[int, Dict[str, Any]] = {}
        pieces: List[str] = []
        piece_details: List[Dict[str, Any]] = []
        for block, variant in zip(blocks, combo):
            frag_index = int(block["fragment"])
            variant_by_fragment[frag_index] = variant
            insert = build_fragment_insert(
                variant.get("sequence") or "",
                str(variant.get("sequence_id") or "unknown"),
                frag_index,
                query_fragments,
                enriched_fragments,
            )
            pieces.append(insert["piece_dna"])
            piece_details.append(
                {
                    **insert,
                    "fragment_index": frag_index,
                    "is_query": bool(variant.get("is_query")),
                }
            )

        chimera_id = _chimera_id_from_variants(variant_by_fragment)
        packed = pack_oligos(
            pieces,
            junction_oh,
            forward_primer=forward_primer,
            reverse_primer=reverse_primer,
            max_oligo_length=max_oligo_length,
            pad_to_max=pad_to_max,
        )
        if payload_capacity is None:
            payload_capacity = packed.get("payload_capacity", 0)

        if packed.get("error"):
            errors.append({"chimera_id": chimera_id, "error": packed["error"]})
            chimera_count += 1
            continue

        for oligo in packed.get("oligos") or []:
            global_oligo_idx += 1
            all_oligos.append(
                {
                    **oligo,
                    "global_oligo_index": global_oligo_idx,
                    "chimera_id": chimera_id,
                    "variant_ids": {
                        frag: variant_by_fragment[frag]["sequence_id"]
                        for frag in sorted(variant_by_fragment)
                    },
                    "piece_details": piece_details,
                }
            )
        chimera_count += 1

    unique_sequences = len({o["sequence"] for o in all_oligos})
    return {
        "oligos": all_oligos,
        "chimera_count": chimera_count,
        "total_oligos": len(all_oligos),
        "unique_oligo_sequences": unique_sequences,
        "payload_capacity": payload_capacity or 0,
        "errors": errors,
        "stage1_pieces": stage1_pieces,
        "junction_overhangs": junction_oh,
        "error": None,
    }


def oligos_to_fasta(
    oligos: List[Dict[str, Any]],
    *,
    prefix: str = "oligo",
    include_chimera_id: bool = False,
) -> str:
    """Format packed oligos as a multi-entry FASTA string."""
    lines: List[str] = []
    for oligo in oligos:
        if oligo.get("block_labels"):
            frag_label = f"nblocks_{oligo.get('n_blocks', 1)}"
            block_part = f" blocks={oligo['block_labels']}" if include_chimera_id else ""
        else:
            frag_part = oligo.get("fragment_indices") or []
            frag_label = (
                f"frag_{frag_part[0]}"
                if len(frag_part) == 1
                else f"frag_{frag_part[0]}-{frag_part[-1]}"
                if frag_part
                else "empty"
            )
            block_part = ""
            if include_chimera_id and oligo.get("chimera_id"):
                block_part = f" chimera={oligo['chimera_id']}"
        idx = oligo.get("global_oligo_index") or oligo.get("oligo_index")
        header = (
            f">{prefix}_{idx}_{frag_label}{block_part} "
            f"len={oligo['length']} insert={oligo['insert_length']} "
            f"stuffer={oligo['stuffer_length']}"
        )
        lines.append(header)
        lines.append(oligo["sequence"])
    return "\n".join(lines) + ("\n" if lines else "")


def build_oligo_from_fragment_sequences(
    fragment_aa_sequences: List[str],
    query_fragments: List[Dict[str, Any]],
    enriched_fragments: List[Dict[str, Any]],
    *,
    forward_primer: str = "",
    reverse_primer: str = "",
    max_oligo_length: int = 300,
    pad_to_max: bool = True,
    sequence_ids: Optional[List[str]] = None,
) -> Dict[str, Any]:
    """
    Stage 2 assembly for one chimera (one AA sequence per fragment).

    Builds Stage 1 pieces, then packs them into the minimum number of oligos.
    """
    pieces: List[str] = []
    piece_details: List[Dict[str, Any]] = []
    for idx, frag in enumerate(query_fragments):
        frag_index = int(frag["index"])
        aa = fragment_aa_sequences[idx] if idx < len(fragment_aa_sequences) else ""
        sid = (
            sequence_ids[idx]
            if sequence_ids and idx < len(sequence_ids)
            else f"fragment_{frag_index}"
        )
        insert = build_fragment_insert(
            aa,
            sid,
            frag_index,
            query_fragments,
            enriched_fragments,
        )
        pieces.append(insert["piece_dna"])
        piece_details.append({**insert, "fragment_index": frag_index})

    junction_oh = junction_overhangs_between_fragments(enriched_fragments)
    packed = pack_oligos(
        pieces,
        junction_oh,
        forward_primer=forward_primer,
        reverse_primer=reverse_primer,
        max_oligo_length=max_oligo_length,
        pad_to_max=pad_to_max,
    )

    merged_insert = merge_piece_range(pieces, junction_oh, 0, len(pieces) - 1) if pieces else ""
    full_oligo = assemble_full_oligo(
        pieces,
        junction_oh,
        forward_primer=forward_primer,
        reverse_primer=reverse_primer,
    )
    return {
        "pieces": pieces,
        "piece_details": piece_details,
        "merged_insert": merged_insert,
        "insert_length": len(merged_insert),
        "junction_overhangs": junction_oh,
        "packed": packed,
        "n_oligos": packed["n_oligos"],
        "oligos": packed["oligos"],
        "payload_capacity": packed["payload_capacity"],
        "pack_error": packed.get("error"),
        "full_oligo": full_oligo,
        "full_length": len(full_oligo),
    }


def assemble_full_oligo(
    pieces: List[str],
    junction_overhangs: List[str],
    *,
    forward_primer: str = "",
    reverse_primer: str = "",
) -> str:
    """Single oligo: primer + merged pieces + reverse-complement of reverse primer."""
    if not pieces:
        return f"{forward_primer}{reverse_complement(reverse_primer)}"
    merged = pieces[0]
    for i in range(1, len(pieces)):
        oh = junction_overhangs[i - 1] if i - 1 < len(junction_overhangs) else ""
        merged = merge_pieces_at_overhang(merged, pieces[i], oh)
    return f"{forward_primer}{merged}{reverse_complement(reverse_primer)}"


def load_assembly_context(
    crossovers: List[int],
    *,
    fragment1_prepend_m: bool = False,
    fragment1_manual_overhang: Optional[str] = None,
) -> Tuple[Optional[List[Dict[str, Any]]], Optional[List[Dict[str, Any]]], Optional[str]]:
    """
    Build query fragments and GGA-enriched fragment metadata.

    Returns (query_fragments, enriched_fragments, error_message).
    """
    aligned_query, _source = get_aligned_query_sequence()
    if not aligned_query:
        return None, None, "No query sequence in session. Run SCHEMA Energy first."

    query_fragments = split_query_sequence_into_fragments(aligned_query, crossovers)
    if not query_fragments:
        return None, None, "Could not split query into fragments."

    if not GGA_COMPATIBILITY_YAML.is_file():
        return None, None, f"GGA compatibility file missing: {GGA_COMPATIBILITY_YAML}"

    try:
        gga_data = load_gga_compatibility(str(GGA_COMPATIBILITY_YAML.resolve()))
        assignment = assign_golden_gate_overhangs(
            query_fragments,
            compatibility=gga_data,
            fragment1_prepend_m=fragment1_prepend_m,
            fragment1_manual_overhang=fragment1_manual_overhang,
        )
    except Exception as exc:
        return None, None, f"Golden Gate assignment failed: {exc}"

    enriched = assignment["fragments"]
    return enriched, enriched, None


def build_fragment_insert(
    aa_sequence: str,
    sequence_id: str,
    frag_index: int,
    query_fragments: List[Dict[str, Any]],
    enriched_fragments: List[Dict[str, Any]],
) -> Dict[str, Any]:
    """Stage-1 insert for one homolog at one fragment position."""
    codons_display, coding_dna = build_fragment_coding_dna(
        aa_sequence,
        frag_index,
        query_fragments,
        enriched_fragments,
    )
    piece_dna = build_bsaI_flanked_piece(coding_dna)
    aa = _clean_aa(aa_sequence)
    return {
        "sequence_id": sequence_id,
        "aa_sequence": aa,
        "codons_display": codons_display,
        "coding_dna": coding_dna,
        "piece_dna": piece_dna,
        "piece_length": len(piece_dna),
        "aa_length": len(aa),
    }


def rows_from_selection_block(
    block: Dict[str, Any],
) -> List[Dict[str, Any]]:
    """Flatten query + homolog rows from a diversity_saved_selections block."""
    rows: List[Dict[str, Any]] = []
    query = block.get("query")
    if query and query.get("sequence"):
        rows.append(
            {
                "sequence_id": query.get("sequence_id") or "query",
                "sequence": query.get("sequence") or "",
                "is_query": True,
            }
        )
    for homolog in block.get("homologs") or block.get("activated") or []:
        if homolog.get("sequence"):
            rows.append(
                {
                    "sequence_id": homolog.get("sequence_id") or "unknown",
                    "sequence": homolog.get("sequence") or "",
                    "is_query": False,
                }
            )
    return rows


def query_aa_sequences_all_fragments(
    selections: Dict[str, Any],
    query_fragments: List[Dict[str, Any]],
) -> Tuple[Optional[List[str]], Optional[str]]:
    """
    Query AA sequence for every assembly fragment (1..N) from saved selections.

    Fragment order follows *query_fragments*, not the Step 1 dropdown.
    """
    sequences: List[str] = []
    for frag in query_fragments:
        frag_index = int(frag["index"])
        key = str(frag_index)
        if key not in selections:
            return None, f"Fragment {frag_index} is missing from saved selections."
        block = selections[key] or {}
        query_row = block.get("query") or {}
        seq = query_row.get("sequence") or frag.get("sequence") or ""
        seq = _clean_aa(seq)
        if not seq:
            return None, f"Fragment {frag_index} has no query sequence in saved selections."
        sequences.append(seq)
    return sequences, None


def build_stage1_inserts_all_fragments(
    selections: Dict[str, Any],
    query_fragments: List[Dict[str, Any]],
    enriched_fragments: List[Dict[str, Any]],
    *,
    use_query: bool = True,
    sequence_id: str = "query",
) -> Tuple[List[Dict[str, Any]], Optional[str]]:
    """
    Build Stage 1 BsaI inserts for every fragment position (1..N).

    When *use_query* is True, uses the query row from each saved fragment block.
    """
    inserts: List[Dict[str, Any]] = []
    for frag in query_fragments:
        frag_index = int(frag["index"])
        key = str(frag_index)
        if key not in selections:
            return [], f"Fragment {frag_index} is missing from saved selections."
        block = selections[key] or {}
        if use_query:
            row = block.get("query") or {}
            aa = row.get("sequence") or frag.get("sequence") or ""
            sid = row.get("sequence_id") or sequence_id
        else:
            return [], "Homolog-specific all-fragment inserts are not implemented."
        if not aa:
            return [], f"Fragment {frag_index} has no sequence for Stage 1."
        insert = build_fragment_insert(
            aa,
            sid,
            frag_index,
            query_fragments,
            enriched_fragments,
        )
        inserts.append({**insert, "fragment_index": frag_index})
    return inserts, None
