"""
Page 4: Assembly Analysis
Query sequence fragments and Golden Gate overhang assignment.
"""

import json
import sys
from pathlib import Path

import pandas as pd
import streamlit as st

sys.path.insert(0, str(Path(__file__).parent.parent))

from utils.assembly_analysis import (
    get_aligned_query_sequence,
    split_query_sequence_into_fragments,
    strip_msa_gaps,
)
from utils.config import GGA_COMPATIBILITY_YAML
from utils.gga_overhangs import (
    XATG_OVERHANGS,
    assign_golden_gate_overhangs,
    assembly_gga_options_from_mapping,
    load_gga_compatibility,
)
from utils.session_manager import init_session_state

st.set_page_config(
    page_title="Assembly Analysis",
    page_icon="🧪",
    layout="wide",
)

init_session_state()

st.title("🧪 Assembly Analysis")

st.markdown("""
**Assembly Analysis** turns your applied crossover positions into **query protein
fragments** and assigns **Golden Gate (BsaI) overhangs** at each junction. Those
overhangs define how DNA blocks will assemble in one-pot Golden Gate and which codons
are used at fragment boundaries in **Diversity Analysis** and **Oligopool Design**.

**Prerequisites:**

1. **1. SCHEMA Energy** — provides the **query sequence** (ideally from an MSA-aligned
   parent so crossover indices match alignment columns).
2. **3. Crossover Analysis** — you must click **Apply crossover selection**. Draft
   checkboxes alone are not enough; this page reads `selected_crossover_positions`.

**What this page does:**

1. **Splits the query** at each applied crossover (1-based **alignment column**
   indices). Gap characters (`-`) are stripped from displayed fragment sequences.
2. **Assigns junction overhangs** between fragment pairs using the AA-pair compatibility
   table (`data/gga_aa_pair_compatibility.yaml`). Each crossover column belongs to the
   **left** fragment at that junction.
3. **Handles fragment 1 N-terminus** specially (see below) for cloning into an expression
   vector or Met-start constructs.

**Fragment 1 N-terminal options**

| Query starts with M? | Your choices | Result |
|---------------------|--------------|--------|
| **Yes** (or you prepend M) | Auto or pick **AATG / TATG / CATG / GATG** | Universal Met start overhang (`XATG`) chosen last among unused junction overhangs |
| **No** | **Prepend M** checkbox | Adds N-terminal Met, then XATG logic applies |
| **No**, no prepend | Optional **vector 5′ overhang** (4 bp) | Matches your cloning vector; otherwise fragment 1 has no special 5′ site |

Changing these widgets **immediately recalculates** overhangs and updates session state
used by downstream pages. Save your project if you want fragment-1 settings restored later.

**Reading the results**

- **Junction assignments** — one row per internal junction: amino-acid pair at the cut,
  chosen 4 bp overhang, efficiency score, and example codon pair.
- **Per-fragment overhangs** — 5′ and 3′ sticky ends for each fragment. The last fragment
  has no 3′ overhang; fragment 1 may have only a 5′ overhang if configured.
- **Export** — download TXT or JSON from the expander for your records.

**Steps:**

1. Confirm **applied crossovers** and **query fragments** look correct (lengths and
   sequences). If you see a warning about an ungapped query, rerun SCHEMA Energy with an
   MSA so crossover columns align properly.
2. Set **fragment 1 N-terminal** options if needed for your cloning scheme.
3. Review junction and per-fragment overhang tables; resolve any errors (e.g. no
   compatible overhang for an AA pair).
4. *(Optional)* Export fragment/overhang summary.

**Next step:** **5. Diversity Analysis** — upload a homolog MSA, filter variants per
fragment using these overhangs, and **Save main list to session** before Oligopool Design.
""")

crossovers = list(st.session_state.get("selected_crossover_positions") or [])
aligned_query, query_source = get_aligned_query_sequence()

if aligned_query is None:
    st.warning(
        "No query sequence found. Run **SCHEMA Energy** (upload MSA or search sequences) "
        "so the aligned query is available in session."
    )
    st.stop()

st.caption(f"Query source: {query_source} ({len(aligned_query)} alignment columns)")
if "ungapped query_sequence" in query_source and "-" not in aligned_query:
    st.warning(
        "Using the ungapped query from SCHEMA Energy. Crossover positions are "
        "alignment-column indices — run contact calculation with an MSA for correct fragment splits."
    )

if not crossovers:
    st.info(
        "No crossover positions applied yet. Select crossovers on **Crossover Analysis** "
        "and click **Apply crossover selection**."
    )
    full_seq = strip_msa_gaps(aligned_query)
    st.subheader("Query Fragments")
    st.markdown(f"**Fragment 1:** `{full_seq}`")
    st.caption(f"Length: {len(full_seq)} residues (full query, gaps removed)")
    st.stop()

fragments = split_query_sequence_into_fragments(aligned_query, crossovers)

st.subheader("Query Fragments")
st.caption(
    "Applied crossovers (1-based alignment columns): "
    + ", ".join(str(p) for p in crossovers)
)

for frag in fragments:
    st.markdown(f"**Fragment {frag['index']}:** `{frag['sequence']}`")
    st.caption(
        f"Length: {frag['length']} residues | "
        f"Alignment columns {frag['aligned_start_1based']}–{frag['aligned_end_1based']} (1-based)"
    )

st.markdown("---")
st.subheader("Golden Gate Overhangs")

if not GGA_COMPATIBILITY_YAML.is_file():
    st.error(
        f"GGA compatibility file not found at `{GGA_COMPATIBILITY_YAML}`. "
        "Place `gga_aa_pair_compatibility.yaml` in the `data/` folder."
    )
    st.stop()

try:
    gga_data = load_gga_compatibility(str(GGA_COMPATIBILITY_YAML.resolve()))
except Exception as exc:
    st.error(f"Could not load GGA compatibility data: {exc}")
    st.stop()

frag1_raw = fragments[0]["sequence"] if fragments else ""
frag1_starts_m = bool(frag1_raw) and frag1_raw[0].upper() == "M"

st.markdown("**Fragment 1 N-terminal Golden Gate**")
st.caption(
    "When fragment 1 starts with Met, the app assigns the best unused "
    f"XATG overhang ({', '.join(XATG_OVERHANGS)}) after junction overhangs. "
    "You may override with a specific XATG overhang, or prepend Met when absent."
)

prepend_m = False
if frag1_raw and not frag1_starts_m:
    prepend_m = st.checkbox(
        "Prepend M to fragment 1",
        value=bool(st.session_state.get("assembly_fragment1_prepend_m")),
        key="assembly_fragment1_prepend_m_widget",
        help="Adds an N-terminal methionine before overhang assignment.",
    )
elif frag1_starts_m:
    st.caption("Fragment 1 already starts with M.")

effective_seq = f"M{frag1_raw}" if prepend_m and frag1_raw and not frag1_starts_m else frag1_raw
effective_starts_m = bool(effective_seq) and effective_seq[0].upper() == "M"

manual_overhang: str | None = None
if effective_starts_m:
    xatg_choice = st.selectbox(
        "Fragment 1 5′ overhang",
        ["Auto (app chooses best unused XATG)"] + list(XATG_OVERHANGS),
        index=0,
        key="assembly_fragment1_xatg_select",
    )
    if xatg_choice != "Auto (app chooses best unused XATG)":
        manual_overhang = xatg_choice
else:
    use_manual_vector = st.checkbox(
        "Specify 5′ vector overhang manually",
        value=bool(st.session_state.get("assembly_fragment1_use_manual_vector")),
        key="assembly_fragment1_use_manual_vector_widget",
        help="Match a cloning vector overhang from the compatibility table.",
    )
    if use_manual_vector:
        manual_input = st.text_input(
            "Vector 5′ overhang (4 bp)",
            value=st.session_state.get("assembly_fragment1_manual_overhang_input") or "",
            max_chars=4,
            key="assembly_fragment1_manual_overhang_input_widget",
        ).strip().upper()
        if len(manual_input) == 4:
            manual_overhang = manual_input
        elif manual_input:
            st.warning("Enter exactly 4 bases for the vector overhang.")

st.session_state.assembly_fragment1_prepend_m = prepend_m
st.session_state.assembly_fragment1_use_manual_vector = (
    not effective_starts_m and bool(st.session_state.get("assembly_fragment1_use_manual_vector_widget"))
)
st.session_state.assembly_fragment1_manual_overhang = manual_overhang

gga_opts = assembly_gga_options_from_mapping(st.session_state)
assignment = assign_golden_gate_overhangs(fragments, compatibility=gga_data, **gga_opts)
enriched = assignment["fragments"]
junctions = assignment["junctions"]
frag1_nt = assignment.get("fragment1_n_terminal")

if gga_opts.get("fragment1_prepend_m"):
    st.info("Met was prepended to fragment 1 for Golden Gate assignment.")

st.caption(
    f"Overhang table: `{GGA_COMPATIBILITY_YAML.name}` "
    f"(schema {assignment.get('metadata', {}).get('schema_version', '?')})"
)

if frag1_nt and frag1_nt.get("overhang"):
    nt_note = frag1_nt.get("warning") or frag1_nt.get("error") or "—"
    st.markdown(
        f"**Fragment 1 N-terminal site:** `{frag1_nt['overhang']}` "
        f"(efficiency {frag1_nt.get('efficiency', 0):.4f}) — {nt_note}"
    )
elif frag1_nt and frag1_nt.get("error"):
    st.error(f"Fragment 1 N-terminal overhang: {frag1_nt['error']}")

if len(enriched) < 2:
    st.info("Only one fragment — no junction overhangs to assign.")
else:
    junction_rows = []
    for junc in junctions:
        if junc.get("error") and not junc.get("overhang"):
            junction_rows.append(
                {
                    "Junction": f"{junc['left_fragment']} → {junc['right_fragment']}",
                    "AA pair": junc.get("aa_pair") or "—",
                    "Overhang": "—",
                    "Efficiency": "—",
                    "Example codons": "—",
                    "Notes": junc.get("error", ""),
                }
            )
        else:
            warning = junc.get("warning", "")
            junction_rows.append(
                {
                    "Junction": f"{junc['left_fragment']} → {junc['right_fragment']}",
                    "AA pair": junc.get("aa_pair", "—"),
                    "Overhang": junc.get("overhang", "—"),
                    "Rev. complement": junc.get("reverse_complement", "—"),
                    "Efficiency": (
                        f"{junc['efficiency']:.4f}" if junc.get("efficiency") is not None else "—"
                    ),
                    "Example codons": junc.get("example_codon_pair", "—"),
                    "Notes": warning or "—",
                }
            )

    st.markdown("**Junction assignments**")
    st.dataframe(pd.DataFrame(junction_rows), use_container_width=True, hide_index=True)

fragment_rows = []
for frag in enriched:
    oh5 = frag.get("overhang_5")
    oh3 = frag.get("overhang_3")
    fragment_rows.append(
        {
            "Fragment": frag["index"],
            "Sequence": frag["sequence"],
            "Length": frag["length"],
            "5′ overhang": oh5.get("overhang") if oh5 else "—",
            "3′ overhang": oh3.get("overhang") if oh3 else "—",
        }
    )

st.markdown("**Per-fragment overhangs**")
st.caption(
    "Fragment 1 may have a 5′ XATG or vector-matched overhang when configured above. "
    "The 3′ overhang of each fragment is the 5′ overhang of the next fragment. "
    "The last fragment has no 3′ overhang."
)
st.dataframe(pd.DataFrame(fragment_rows), use_container_width=True, hide_index=True)

for frag in enriched:
    oh5 = frag.get("overhang_5")
    oh3 = frag.get("overhang_3")
    parts = [f"**Fragment {frag['index']}:** `{frag['sequence']}`"]
    parts.append(f"5′: `{oh5['overhang']}`" if oh5 and oh5.get("overhang") else "5′: —")
    parts.append(f"3′: `{oh3['overhang']}`" if oh3 and oh3.get("overhang") else "3′: —")
    if oh3 and oh3.get("example_codon_pair"):
        parts.append(f"3′ junction codons: `{oh3['example_codon_pair']}`")
    st.markdown(" | ".join(parts))

with st.expander("Export fragments and overhangs", expanded=False):

    def _junction_export(junc):
        if not junc:
            return None
        return {
            "junction": junc.get("junction"),
            "left_fragment": junc.get("left_fragment"),
            "right_fragment": junc.get("right_fragment"),
            "aa_pair": junc.get("aa_pair"),
            "overhang": junc.get("overhang"),
            "reverse_complement": junc.get("reverse_complement"),
            "efficiency": junc.get("efficiency"),
            "example_codon_pair": junc.get("example_codon_pair"),
            "usage_rank": junc.get("usage_rank"),
        }

    export_fragments = []
    for frag in enriched:
        export_fragments.append(
            {
                "fragment": frag["index"],
                "sequence": frag["sequence"],
                "length": frag["length"],
                "aligned_columns_1based": [
                    frag["aligned_start_1based"],
                    frag["aligned_end_1based"],
                ],
                "overhang_5": _junction_export(frag.get("overhang_5")),
                "overhang_3": _junction_export(frag.get("overhang_3")),
            }
        )

    export_lines = []
    for frag in enriched:
        oh5 = frag.get("overhang_5")
        oh3 = frag.get("overhang_3")
        line = f"Fragment {frag['index']}: {frag['sequence']}"
        if oh5 and oh5.get("overhang"):
            line += f" | 5'={oh5['overhang']}"
        if oh3 and oh3.get("overhang"):
            line += f" | 3'={oh3['overhang']}"
        export_lines.append(line)

    export_text = "\n".join(export_lines)
    export_json = json.dumps(
        {
            "crossover_positions_1based": crossovers,
            "gga_compatibility_file": str(GGA_COMPATIBILITY_YAML),
            "junctions": [_junction_export(j) for j in junctions],
            "fragments": export_fragments,
        },
        indent=2,
    )
    col_txt, col_json = st.columns(2)
    with col_txt:
        st.download_button(
            "Download TXT",
            data=export_text + "\n",
            file_name="query_fragments_overhangs.txt",
            mime="text/plain",
            key="download_query_fragments_txt",
        )
    with col_json:
        st.download_button(
            "Download JSON",
            data=export_json + "\n",
            file_name="query_fragments_overhangs.json",
            mime="application/json",
            key="download_query_fragments_json",
        )
