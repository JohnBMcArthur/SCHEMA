"""
Page 7: Oligopool Design
BsaI Golden Gate oligos from saved fragment pools and Assembly Analysis overhangs.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any, Dict, List

import streamlit as st
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent.parent))

from utils.library_optimization.pools import load_block_pools_from_session, pools_summary
from utils.oligopool_display import (
    GGA_HIGHLIGHT_CSS,
    highlight_dna_sequence,
    highlight_oligo_sequence,
    overhang_spans_for_block_packed_oligo,
    overhang_spans_in_coding,
    overhang_spans_in_piece,
    region_spans_for_full_oligo,
    render_full_oligo_legend,
    render_highlight_legend,
    render_html_table,
)
from utils.gga_overhangs import assembly_gga_options_from_mapping
from utils.oligopool_design import (
    BSAI_END,
    BSAI_START,
    build_fragment_insert,
    load_assembly_context,
    oligos_to_fasta,
    pack_block_oligos,
    reverse_complement,
    rows_from_selection_block,
)
from utils.stuffer_design import apply_stuffer_design
from utils.session_manager import auto_save, init_session_state

st.set_page_config(
    page_title="Oligopool Design",
    page_icon="🧪",
    layout="wide",
)

init_session_state()

st.title("🧪 Oligopool Design")

st.markdown(
    "Design **BsaI Golden Gate** oligos from saved fragment homolog pools. "
    "Terminal codons follow **Assembly Analysis** junction overhangs and reading "
    "frame; interior residues use high-abundance **E. coli** codons."
)

selections = st.session_state.get("diversity_saved_selections") or {}
crossovers = list(st.session_state.get("selected_crossover_positions") or [])

if not selections:
    st.warning(
        "No saved fragment lists found. Save from **Diversity Analysis** "
        "(**Save main list to session**) or **Library Optimization** "
        "(**Save filtered list to session**)."
    )
    st.stop()

blocks, pool_err = load_block_pools_from_session(selections)
if pool_err:
    st.error(pool_err)
    st.stop()

if not crossovers:
    st.warning(
        "No crossover positions applied. Select crossovers on **Crossover Analysis** "
        "and click **Apply crossover selection** before designing oligos."
    )
    st.stop()

query_fragments, enriched_fragments, asm_err = load_assembly_context(
    crossovers,
    **assembly_gga_options_from_mapping(st.session_state),
)
if asm_err or not query_fragments or not enriched_fragments:
    st.error(asm_err or "Could not load assembly context.")
    st.stop()

if len(blocks) != len(query_fragments):
    st.warning(
        f"Saved pools cover **{len(blocks)}** fragment(s) but assembly has "
        f"**{len(query_fragments)}** fragment(s). Ensure saved lists match "
        "applied crossovers."
    )

summary = pools_summary(blocks)
st.caption(
    f"{summary['n_blocks']} fragment pool(s) · "
    f"{summary['total_variants']} total variants · "
    f"up to {summary['total_combinations']:,} chimera combinations"
)

st.markdown("---")
st.subheader("Options")


def _init_oligopool_option_widgets() -> None:
    """Load saved options into widgets on first visit after reload."""
    if st.session_state.pop("_oligopool_sync_widgets_on_load", False):
        st.session_state.oligopool_fwd_widget = (
            st.session_state.get("oligopool_forward_primer") or ""
        )
        st.session_state.oligopool_rev_widget = (
            st.session_state.get("oligopool_reverse_primer") or ""
        )
        st.session_state.oligopool_max_widget = int(
            st.session_state.get("oligopool_max_length") or 300
        )
        return
    st.session_state.setdefault(
        "oligopool_fwd_widget",
        st.session_state.get("oligopool_forward_primer") or "",
    )
    st.session_state.setdefault(
        "oligopool_rev_widget",
        st.session_state.get("oligopool_reverse_primer") or "",
    )
    st.session_state.setdefault(
        "oligopool_max_widget",
        int(st.session_state.get("oligopool_max_length") or 300),
    )


_init_oligopool_option_widgets()

opt_col1, opt_col2, opt_col3 = st.columns(3)
with opt_col1:
    st.text_input(
        "Forward Primer",
        key="oligopool_fwd_widget",
        help="5′ annealing sequence appended to each packed oligo (Stage 2).",
    )
with opt_col2:
    st.text_input(
        "Reverse Primer",
        key="oligopool_rev_widget",
        help="3′ annealing sequence; reverse complement is used in each oligo.",
    )
with opt_col3:
    st.number_input(
        "Maximum oligo length (nt)",
        min_value=60,
        max_value=500,
        step=10,
        key="oligopool_max_widget",
    )

if st.button("Save Options", type="primary"):
    st.session_state.oligopool_forward_primer = (
        str(st.session_state.oligopool_fwd_widget).strip().upper()
    )
    st.session_state.oligopool_reverse_primer = (
        str(st.session_state.oligopool_rev_widget).strip().upper()
    )
    st.session_state.oligopool_max_length = int(st.session_state.oligopool_max_widget)
    auto_save("oligopool_options")
    st.session_state["_oligopool_sync_widgets_on_load"] = True
    st.success("Options saved — they will be restored on the next page load.")
    st.rerun()

forward_primer = str(st.session_state.oligopool_fwd_widget).strip().upper()
reverse_primer = str(st.session_state.oligopool_rev_widget).strip().upper()
max_oligo_length = int(st.session_state.oligopool_max_widget)

st.caption(
    f"BsaI sites (Stage 1): 5′ `{BSAI_START}` · 3′ `{BSAI_END}` "
    "(GGTCTC / GAGACC recognition, no spacer bases)."
)
st.markdown(render_highlight_legend(), unsafe_allow_html=True)

st.markdown("---")
st.subheader("Step 1 — Fragment inserts (Stage 1)")

frag_keys = sorted(selections.keys(), key=lambda k: int(k))
frag_labels = [f"Fragment {k}" for k in frag_keys]
selected_label = st.selectbox(
    "Fragment",
    options=frag_labels,
    index=0,
)
selected_key = frag_keys[frag_labels.index(selected_label)]
frag_index = int(selected_key)

if frag_index > len(query_fragments):
    st.error(
        f"Fragment {frag_index} is not defined in the current assembly "
        f"({len(query_fragments)} fragment(s))."
    )
else:
    query_frag = query_fragments[frag_index - 1]
    enriched_frag = enriched_fragments[frag_index - 1]
    oh5 = enriched_frag.get("overhang_5")
    oh3 = enriched_frag.get("overhang_3")

    info_cols = st.columns(4)
    info_cols[0].metric("Query fragment AA", query_frag.get("length", len(query_frag["sequence"])))
    info_cols[1].metric("Saved homologs", len(rows_from_selection_block(selections[selected_key])) - 1)
    info_cols[2].write(
        f"**5′ junction:** {oh5.get('overhang') if oh5 else '—'}"
        + (f" (frame {oh5.get('frame')})" if oh5 and oh5.get("frame") is not None else "")
    )
    info_cols[3].write(
        f"**3′ junction:** {oh3.get('overhang') if oh3 else '—'}"
        + (f" (frame {oh3.get('frame')})" if oh3 and oh3.get("frame") is not None else "")
    )

    st.markdown(f"**Query reference:** `{query_frag['sequence']}`")

    block = selections[selected_key] or {}
    rows = rows_from_selection_block(block)
    if not rows:
        st.info("No sequences in this fragment pool.")
    else:
        insert_rows: List[Dict[str, Any]] = []
        for row in rows:
            insert = build_fragment_insert(
                row["sequence"],
                row["sequence_id"],
                frag_index,
                query_fragments,
                enriched_fragments,
            )
            insert_rows.append(
                {
                    "Sequence ID": insert["sequence_id"],
                    "Source": "query" if row.get("is_query") else "homolog",
                    "AA sequence": insert["aa_sequence"],
                    "AA length": insert["aa_length"],
                    "Codons": insert["codons_display"],
                    "Coding DNA": insert["coding_dna"],
                    "BsaI piece (Stage 1)": insert["piece_dna"],
                    "Piece length (nt)": insert["piece_length"],
                }
            )

        coding_html: List[str] = []
        piece_html: List[str] = []
        for row in insert_rows:
            coding_html.append(
                highlight_dna_sequence(
                    row["Coding DNA"],
                    overhang_spans_in_coding(
                        row["Coding DNA"], frag_index, enriched_fragments
                    ),
                    highlight_bsaI=False,
                )
            )
            piece_html.append(
                highlight_dna_sequence(
                    row["BsaI piece (Stage 1)"],
                    overhang_spans_in_piece(
                        row["BsaI piece (Stage 1)"], frag_index, enriched_fragments
                    ),
                )
            )

        st.markdown(
            render_html_table(
                insert_rows,
                [
                    ("Sequence ID", "Sequence ID"),
                    ("Source", "Source"),
                    ("AA sequence", "AA sequence"),
                    ("AA length", "AA length"),
                    ("Codons", "Codons"),
                    ("Coding DNA", "Coding DNA"),
                    ("BsaI piece (Stage 1)", "BsaI piece (Stage 1)"),
                    ("Piece length (nt)", "Piece length (nt)"),
                ],
                html_columns={
                    "Coding DNA": coding_html,
                    "BsaI piece (Stage 1)": piece_html,
                },
            ),
            unsafe_allow_html=True,
        )

        with st.expander("Full BsaI-flanked pieces (Stage 1)", expanded=False):
            for row_idx, row in enumerate(insert_rows):
                st.markdown(
                    f"**{row['Sequence ID']}** ({row['Source']}) — "
                    f"{row['Piece length (nt)']} nt"
                )
                st.markdown(
                    GGA_HIGHLIGHT_CSS + f'<p style="word-break:break-all;">{piece_html[row_idx]}</p>',
                    unsafe_allow_html=True,
                )

st.markdown("---")
st.subheader("Step 2 — Oligo packing (Stage 2)")

st.markdown(
    "Pack **each unique Stage 1 block once** into the minimum number of oligos. "
    "Short blocks are combined on the same oligo when they fit under the max length. "
    "Unused length is filled with **N** stuffer tails **outside** the primer binding "
    "regions: "
    "`[5′ N…][forward primer][BsaI blocks…][RC(reverse primer)][N… 3′]`."
)

st.caption(
    f"**{summary['total_variants']}** unique blocks across **{summary['n_blocks']}** fragment sites. "
    "Set primers and max length above, then run packing."
)

pack_options_sig = (forward_primer, reverse_primer, max_oligo_length)

if st.button("Pack oligos", type="primary", key="oligopool_run_step2"):
    with st.spinner("Packing block oligos…"):
        st.session_state.oligopool_library_result = pack_block_oligos(
            blocks,
            query_fragments,
            enriched_fragments,
            forward_primer=forward_primer,
            reverse_primer=reverse_primer,
            max_oligo_length=max_oligo_length,
            pad_to_max=True,
        )
        st.session_state.oligopool_library_options_sig = pack_options_sig
        st.session_state.pop("oligopool_stuffer_result", None)
        st.session_state.pop("oligopool_stuffer_sig", None)

library = st.session_state.get("oligopool_library_result")

if library is None:
    st.info("Configure options above and click **Pack oligos** to run Stage 2.")
elif pack_options_sig != st.session_state.get("oligopool_library_options_sig"):
    st.warning("Options have changed since the last run. Click **Pack oligos** to recalculate.")
    library = None

if library is not None:
    if library.get("error"):
        st.error(library["error"])

if library is not None and not library.get("error"):
    rev_rc = reverse_complement(reverse_primer)
    m1, m2, m3, m4, m5 = st.columns(5)
    m1.metric("Unique blocks", library["total_blocks"])
    m2.metric("Total oligos", library["total_oligos"])
    m3.metric("Unique sequences", library["unique_oligo_sequences"])
    m4.metric("Payload capacity", f"{library['payload_capacity']} nt")
    m5.metric("Primer overhead", f"{len(forward_primer) + len(rev_rc)} nt")

    if library["errors"]:
        st.warning(
            f"{len(library['errors'])} block(s) could not be packed — see errors below."
        )
    else:
        st.success(
            f"Packed **{library['total_blocks']}** unique blocks "
            f"into **{library['total_oligos']:,}** oligos."
        )

    st.markdown("**All Stage 1 pieces (every homolog × fragment site)**")
    stage1_table = []
    stage1_html: List[str] = []
    for piece in library["stage1_pieces"]:
        frag_idx = piece["fragment_index"]
        stage1_table.append(
            {
                "Fragment": frag_idx,
                "Sequence ID": piece["sequence_id"],
                "Source": "query" if piece.get("is_query") else "homolog",
                "AA": piece["aa_sequence"],
                "BsaI piece": piece["piece_dna"],
                "Length (nt)": piece["piece_length"],
            }
        )
        stage1_html.append(
            highlight_dna_sequence(
                piece["piece_dna"],
                overhang_spans_in_piece(piece["piece_dna"], frag_idx, enriched_fragments),
            )
        )
    st.markdown(
        render_html_table(
            stage1_table,
            [
                ("Fragment", "Fragment"),
                ("Sequence ID", "Sequence ID"),
                ("Source", "Source"),
                ("AA", "AA"),
                ("BsaI piece", "BsaI piece"),
                ("Length (nt)", "Length (nt)"),
            ],
            html_columns={"BsaI piece": stage1_html},
        ),
        unsafe_allow_html=True,
    )

    if library["errors"]:
        with st.expander("Packing errors", expanded=True):
            st.dataframe(
                pd.DataFrame(library["errors"]),
                use_container_width=True,
                hide_index=True,
            )

    oligo_rows = []
    oligo_seq_html: List[str] = []

    display_oligos = library["oligos"]
    max_display = 200
    if len(display_oligos) > max_display:
        st.caption(f"Showing first {max_display} of {len(display_oligos):,} oligos in the table.")

    for oligo in display_oligos[:max_display]:
        stuffer_5_len = int(oligo.get("stuffer_5_length") or 0)
        oligo_insert_start = stuffer_5_len + len(forward_primer)
        oh_spans = overhang_spans_for_block_packed_oligo(
            oligo,
            enriched_fragments,
            insert_offset=oligo_insert_start,
        )
        stuffer_status = (
            "designed"
            if oligo.get("stuffer_designed")
            else ("Ns" if int(oligo.get("stuffer_length") or 0) > 0 else "—")
        )
        oligo_rows.append(
            {
                "Oligo": oligo.get("global_oligo_index"),
                "Blocks": oligo.get("block_labels", ""),
                "# blocks": oligo.get("n_blocks", 0),
                "5′ tail": stuffer_5_len,
                "3′ tail": int(oligo.get("stuffer_3_length") or 0),
                "Stuffer": stuffer_status,
                "Sequence": oligo["sequence"],
                "Insert (nt)": oligo["insert_length"],
                "Total (nt)": oligo["length"],
            }
        )
        oligo_seq_html.append(highlight_dna_sequence(oligo["sequence"], oh_spans))

    st.markdown("**Packed oligo library**")
    st.markdown(
        render_html_table(
            oligo_rows,
            [
                ("Oligo", "Oligo"),
                ("Blocks", "Blocks"),
                ("# blocks", "# blocks"),
                ("5′ tail", "5′ tail"),
                ("3′ tail", "3′ tail"),
                ("Stuffer", "Stuffer"),
                ("Sequence", "Sequence"),
                ("Insert (nt)", "Insert (nt)"),
                ("Total (nt)", "Total (nt)"),
            ],
            html_columns={"Sequence": oligo_seq_html},
        ),
        unsafe_allow_html=True,
    )

    if library["oligos"]:
        fasta_body = oligos_to_fasta(
            library["oligos"],
            prefix="oligopool",
            include_chimera_id=True,
        )
        st.download_button(
            "Download packed oligos with N stuffer (FASTA)",
            data=fasta_body,
            file_name="oligopool_blocks_packed_N.fasta",
            mime="text/plain",
            key="oligopool_download_n_fasta",
        )

    st.markdown("---")
    st.subheader("Step 3 — Stuffer sequence design")

    st.markdown(
        "Replace **N** padding tails with precomputed **A/T/G/C** sequences "
        f"(library: lengths 1–{100}, 32 candidates each). "
        "Sequences avoid internal BsaI sites, target ~50% GC, limit homopolymers "
        "to 3 nt, and avoid strong hairpins. Junction checks use your primer context."
    )

    has_stuffer = any(int(o.get("stuffer_length") or 0) > 0 for o in library["oligos"])
    stuffer_sig = pack_options_sig

    if not has_stuffer:
        st.info("No stuffer padding in the current pack (all oligos already at max insert capacity).")
    elif st.button("Replace stuffer Ns", type="primary", key="oligopool_run_step3"):
        with st.spinner("Designing stuffer sequences…"):
            st.session_state.oligopool_stuffer_result = apply_stuffer_design(library)
            st.session_state.oligopool_stuffer_sig = stuffer_sig
        st.rerun()

    stuffer_result = st.session_state.get("oligopool_stuffer_result")
    if stuffer_result and stuffer_sig != st.session_state.get("oligopool_stuffer_sig"):
        st.warning("Packing options changed. Re-run **Pack oligos** before replacing stuffer Ns.")
        stuffer_result = None

    if stuffer_result:
        n_designed = sum(1 for o in stuffer_result["oligos"] if o.get("stuffer_designed"))
        n_errors = len(stuffer_result.get("stuffer_design_errors") or [])
        if n_errors:
            st.warning(
                f"Designed stuffer tails for **{n_designed}** oligo(s); "
                f"**{n_errors}** oligo(s) still contain N stuffer bases."
            )
        else:
            st.success(f"Designed stuffer tails for **{n_designed}** oligo(s).")

        warn_rows = stuffer_result.get("stuffer_design_warnings") or []
        if warn_rows:
            with st.expander("Stuffer design warnings", expanded=False):
                st.dataframe(pd.DataFrame(warn_rows), use_container_width=True, hide_index=True)

        st.markdown(render_full_oligo_legend(), unsafe_allow_html=True)

        designed_rows = []
        designed_html: List[str] = []
        display_designed = stuffer_result["oligos"]
        if len(display_designed) > max_display:
            st.caption(
                f"Showing first {max_display} of {len(display_designed):,} final oligos."
            )

        for oligo in display_designed[:max_display]:
            gc_pct = oligo.get("stuffer_gc_pct")
            regions = region_spans_for_full_oligo(oligo, enriched_fragments)
            designed_rows.append(
                {
                    "Oligo": oligo.get("global_oligo_index"),
                    "Blocks": oligo.get("block_labels", ""),
                    "# blocks": oligo.get("n_blocks", 0),
                    "5′ tail": int(oligo.get("stuffer_5_length") or 0),
                    "3′ tail": int(oligo.get("stuffer_3_length") or 0),
                    "Stuffer GC": f"{100 * gc_pct:.0f}%" if gc_pct is not None else "—",
                    "Sequence": oligo["sequence"],
                    "Insert (nt)": oligo["insert_length"],
                    "Total (nt)": oligo["length"],
                }
            )
            designed_html.append(
                highlight_oligo_sequence(oligo["sequence"], regions)
            )

        st.markdown("**Final oligos (designed stuffers)**")
        st.markdown(
            render_html_table(
                designed_rows,
                [
                    ("Oligo", "Oligo"),
                    ("Blocks", "Blocks"),
                    ("# blocks", "# blocks"),
                    ("5′ tail", "5′ tail"),
                    ("3′ tail", "3′ tail"),
                    ("Stuffer GC", "Stuffer GC"),
                    ("Sequence", "Sequence"),
                    ("Insert (nt)", "Insert (nt)"),
                    ("Total (nt)", "Total (nt)"),
                ],
                html_columns={"Sequence": designed_html},
            ),
            unsafe_allow_html=True,
        )

        fasta_designed = oligos_to_fasta(
            stuffer_result["oligos"],
            prefix="oligopool",
            include_chimera_id=True,
        )
        st.download_button(
            "Download final oligos (FASTA)",
            data=fasta_designed,
            file_name="oligopool_blocks_final.fasta",
            mime="text/plain",
            key="oligopool_download_final_fasta",
        )
