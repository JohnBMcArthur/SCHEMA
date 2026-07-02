"""
Structure viewer utilities for Mol* integration.
Provides functions for pre-computing distances and managing structure visualization.
"""

import numpy as np
from scipy.spatial import distance
from typing import Any, Dict, List, Optional, Set, Tuple

try:
    from schemarecomb import PDBStructure
    SCHEMARECOMB_AVAILABLE = True
except ImportError:
    SCHEMARECOMB_AVAILABLE = False
    PDBStructure = None

from utils.config import CROSSOVER_FRAGMENT_COLORS

# Pinned Mol* build for embedded ``components.html`` viewers (crossover + SCHEMA Energy).
MOLSTAR_VIEWER_VERSION = "4.6.0"
MOLSTAR_VIEWER_BASE = (
    f"https://cdn.jsdelivr.net/npm/molstar@{MOLSTAR_VIEWER_VERSION}/build/viewer"
)


def precompute_nearby_residues(pdb_structure, distance_threshold=8.0):
    """
    Pre-compute all residues within distance_threshold Å of each residue.
    
    Args:
        pdb_structure: schemarecomb.PDBStructure object
        distance_threshold: Distance threshold in Angstroms (default: 8.0)
        
    Returns:
        dict: Dictionary mapping residue index -> set of nearby residue indices
    """
    if not SCHEMARECOMB_AVAILABLE:
        raise ImportError("schemarecomb is required")
    
    nearby_residues = {}
    aa_coords = {aa.index: aa.coords for aa in pdb_structure.amino_acids}
    
    # Compute distances between all residue pairs
    for i in aa_coords:
        nearby = set()
        coords_i = aa_coords[i]
        
        for j in aa_coords:
            if i == j:
                continue
            coords_j = aa_coords[j]
            # Compute minimum distance between any atoms in residues i and j
            d = distance.cdist(coords_i, coords_j)
            min_dist = np.min(d)
            if min_dist < distance_threshold:
                nearby.add(j)
        
        nearby_residues[i] = nearby
    
    return nearby_residues


def get_contacts_for_residue(contacts: List[Tuple], residue_index: int) -> Set[int]:
    """
    Get all residues that contact a given residue.
    
    Args:
        contacts: List of contact tuples (i, j) or (i, j, ri, rj)
        residue_index: Index of residue to get contacts for
        
    Returns:
        set: Set of residue indices that contact the given residue
    """
    contacting_residues = set()
    for contact in contacts:
        # Handle both (i, j) and (i, j, ri, rj) formats
        i = contact[0]
        j = contact[1]
        if i == residue_index:
            contacting_residues.add(j)
        elif j == residue_index:
            contacting_residues.add(i)
    return contacting_residues


def update_contacts(contacts: List[Tuple], residue_i: int, residue_j: int, 
                   is_contact: bool) -> List[Tuple]:
    """
    Add or remove a contact between two residues.
    
    Args:
        contacts: Current list of contacts (can be (i, j) or (i, j, ri, rj) format)
        residue_i: First residue index
        residue_j: Second residue index
        is_contact: True to add contact, False to remove
        
    Returns:
        list: Updated list of contacts (preserves original format)
    """
    # Determine format
    is_extended_format = len(contacts[0]) == 4 if contacts else False
    
    # Create a set for fast lookup (using just i, j)
    contact_set = set()
    contact_dict = {}  # Store full contact tuples
    
    for contact in contacts:
        i, j = contact[0], contact[1]
        contact_set.add((i, j))
        contact_set.add((j, i))
        if is_extended_format:
            contact_dict[(i, j)] = contact
            contact_dict[(j, i)] = (j, i, contact[3], contact[2])  # Swap ri, rj
    
    if is_contact:
        # Add contact (both directions)
        contact_set.add((residue_i, residue_j))
        contact_set.add((residue_j, residue_i))
        if is_extended_format:
            # Use None for ri, rj if not available
            contact_dict[(residue_i, residue_j)] = (residue_i, residue_j, None, None)
            contact_dict[(residue_j, residue_i)] = (residue_j, residue_i, None, None)
    else:
        # Remove contact (both directions)
        contact_set.discard((residue_i, residue_j))
        contact_set.discard((residue_j, residue_i))
        contact_dict.pop((residue_i, residue_j), None)
        contact_dict.pop((residue_j, residue_i), None)
    
    # Convert back to list, keeping only one direction per pair
    updated_contacts = []
    seen = set()
    for i, j in contact_set:
        if (i, j) not in seen and (j, i) not in seen:
            if is_extended_format:
                updated_contacts.append(contact_dict.get((i, j), (i, j, None, None)))
            else:
                updated_contacts.append((i, j))
            seen.add((i, j))
    
    return updated_contacts


def get_pdb_file_content(pdb_path_or_content) -> str:
    """
    Read PDB file content or return content if already a string.
    
    Args:
        pdb_path_or_content: Path to PDB file or PDB file content string
        
    Returns:
        str: PDB file content
    """
    import os
    
    # Handle None or empty
    if not pdb_path_or_content:
        return ""
    
    # Convert to string if not already
    pdb_str = str(pdb_path_or_content)
    
    # Check if it looks like PDB content (starts with PDB keywords)
    # Check first 200 characters to avoid reading huge files
    preview = pdb_str[:200] if len(pdb_str) > 200 else pdb_str
    if (preview.startswith("HEADER") or 
        preview.startswith("ATOM") or 
        "\nATOM" in preview or 
        "REMARK" in preview or
        preview.startswith("MODEL") or
        "ATOM      1" in preview):
        # It's already PDB content
        return pdb_str
    
    # Check if it's a valid file path (only if it doesn't look like content)
    # Only check if the string is short (likely a path) and doesn't contain newlines
    if '\n' not in pdb_str and len(pdb_str) < 500:
        try:
            if os.path.exists(pdb_str) and os.path.isfile(pdb_str):
                with open(pdb_str, 'r') as f:
                    return f.read()
        except (OSError, IOError):
            # If file access fails, assume it's content
            pass
    
    # Assume it's already content or return as-is
    return pdb_str


def pdb_structure_to_pdb_content(pdb_structure, color_map: Optional[Dict[int, float]] = None) -> str:
    """
    Convert schemarecomb PDBStructure to PDB file content.
    
    IMPORTANT: This function renumbers residues to use aa.index + 1 (1-based)
    to match our application's indexing. The original PDB residue numbers
    are NOT preserved - this ensures consistency with the contact list UI.
    
    Args:
        pdb_structure: schemarecomb.PDBStructure object
        color_map: Optional dictionary mapping residue index -> B-factor value for coloring
        
    Returns:
        str: PDB file content
    """
    if not SCHEMARECOMB_AVAILABLE:
        raise ImportError("schemarecomb is required")
    
    from dataclasses import replace
    
    lines = []
    
    # Add header
    lines.append("HEADER    PROTEIN STRUCTURE")
    lines.append("REMARK   1 GENERATED FROM SCHEMARECOMB")
    lines.append("REMARK   2 RESIDUE NUMBERS ARE APPLICATION INDICES (0-based + 1)")
    
    # Sort amino acids by index
    sorted_aas = sorted(pdb_structure.amino_acids, key=lambda aa: aa.index)
    
    # Write ATOM lines with sequential serial numbers
    # IMPORTANT: Renumber res_index to match aa.index for consistency with UI
    atom_serial = 1
    for aa in sorted_aas:
        # Get color value for this residue (default to original temp factor)
        color_value = color_map.get(aa.index, aa.atoms[0].temp) if color_map else aa.atoms[0].temp
        
        # Use aa.index as the residue number (0-based, to_line() adds 1)
        # This ensures PDB residue numbers match our application's indexing
        consistent_res_index = aa.index
        
        for atom in aa.atoms:
            # Create ATOM line - update serial number, res_index, and optionally B-factor
            atom_with_updates = replace(
                atom, 
                serial_num=atom_serial, 
                res_index=consistent_res_index,
                temp=color_value
            )
            atom_line = atom_with_updates.to_line()
            lines.append(atom_line)
            atom_serial += 1
    
    # Add END
    lines.append("END")
    
    return '\n'.join(lines)


def pdb_structure_to_sequential_cartoon_pdb(
    pdb_structure,
    color_map: Optional[Dict[int, float]],
) -> str:
    """
    Write PDB with **1..N consecutive** ``resSeq`` in alignment order (sorted ``aa.index``).

    Gaps in MSA column indices are not written as gaps in ``resSeq``, so viewers that
    tie the polymer trace to residue numbering (for example older 3Dmol cartoons)
    still draw one continuous chain. Mol* does not require this, but it keeps
    cross-fragment coloring aligned with residue order.

    ``color_map`` keys are **0-based alignment indices** (``aa.index``); values are
    written to the temperature-factor column (e.g. fragment ids for coloring).

    Alignment / SCHEMA positions are not preserved in resSeq; use only for display.
    """
    if not SCHEMARECOMB_AVAILABLE:
        raise ImportError("schemarecomb is required")

    from dataclasses import replace

    sorted_aas = sorted(pdb_structure.amino_acids, key=lambda aa: aa.index)
    lines = [
        "HEADER    SEQUENTIAL DISPLAY FOR CARTOON",
        "REMARK  99 ResSeq 1..N = order along structure in alignment order, not MSA columns.",
        "REMARK  99 B-factor = encoded value for viewer coloring (e.g. fragment scale).",
    ]

    atom_serial = 1
    for display_pos, aa in enumerate(sorted_aas, start=1):
        res_index_0based = display_pos - 1
        if color_map is not None:
            color_val = float(color_map.get(aa.index, 0.0))
        else:
            color_val = float(aa.atoms[0].temp)

        for atom in aa.atoms:
            atom_with_updates = replace(
                atom,
                serial_num=atom_serial,
                res_index=res_index_0based,
                temp=color_val,
            )
            lines.append(atom_with_updates.to_line())
            atom_serial += 1

    lines.append("END")
    return "\n".join(lines)


def create_colored_pdb_content(pdb_structure, selected_idx: int, nearby: Set[int], 
                               contacting: Set[int]) -> str:
    """
    Create PDB content with color encoding in the B-factor column for Mol*.

    Values are chosen for the **uncertainty** color theme (low / mid / high):
    bulk polymer, focal residue plus nearby (8 Å), and SCHEMA contacts (contacts
    override nearby when both apply).

    Args:
        pdb_structure: schemarecomb.PDBStructure object
        selected_idx: Index of selected/focal residue
        nearby: Set of nearby residue indices (within 8Å)
        contacting: Set of contacting residue indices
        
    Returns:
        str: PDB file content with color encoding
    """
    # B-factor categories for Mol* uncertainty coloring (spread for a clear gradient).
    # 0 = bulk polymer, 1 = focal + nearby, 2 = SCHEMA contacts (overrides 1).
    color_map = {}
    for aa in pdb_structure.amino_acids:
        color_map[aa.index] = 12.0
    stick_residues = {selected_idx} | nearby
    for idx in stick_residues:
        color_map[idx] = 55.0
    for idx in contacting:
        color_map[idx] = 98.0
    
    return pdb_structure_to_pdb_content(pdb_structure, color_map=color_map)


def _embedded_molstar_single_pdb_html(
    *,
    viewer_id: str,
    pdb_json: str,
    caption_html: str,
    height_px: int,
    layout_show_sequence: bool,
) -> str:
    """
    Minimal Mol* (molstar ``Viewer``) page: load one PDB string and apply the
    trajectory ``default`` preset with **polymer cartoon** (not Mol*'s size-based
    ``auto`` preset, which uses spacefill for small structures) and **uncertainty**
    coloring from B-factors.
    """
    css_href = f"{MOLSTAR_VIEWER_BASE}/molstar.css"
    js_src = f"{MOLSTAR_VIEWER_BASE}/molstar.js"
    seq_js = "true" if layout_show_sequence else "false"
    return f"""<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8"/>
    <link rel="stylesheet" type="text/css" href="{css_href}" />
    <style>
        body {{ margin: 0; padding: 0; font-family: system-ui, sans-serif; }}
        .mol-caption {{
            font-size: 13px; color: #333; padding: 6px 8px;
            border-bottom: 1px solid #e0e0e0; background: #fafafa;
        }}
        #{viewer_id} {{
            position: relative; width: 100%; height: {height_px}px;
        }}
    </style>
    <script src="{js_src}"></script>
</head>
<body>
    <div class="mol-caption">{caption_html}</div>
    <div id="{viewer_id}"></div>
    <script>
        document.addEventListener('DOMContentLoaded', async function () {{
            const pdbData = {pdb_json};
            try {{
                const viewer = await molstar.Viewer.create('{viewer_id}', {{
                    layoutIsExpanded: true,
                    layoutShowControls: true,
                    layoutShowRemoteState: false,
                    layoutShowSequence: {seq_js},
                    layoutShowLog: false,
                    layoutShowLeftPanel: false,
                    viewportShowExpand: true,
                    viewportShowSelectionMode: false,
                    viewportShowAnimation: false,
                    pdbProvider: 'rcsb',
                    emdbProvider: 'rcsb',
                }});
                const plugin = viewer.plugin;
                const data = await plugin.builders.data.rawData({{ data: pdbData, label: 'structure' }});
                const traj = await plugin.builders.structure.parseTrajectory(data, 'pdb');
                await plugin.builders.structure.hierarchy.applyPreset(traj, 'default', {{
                    representationPreset: 'polymer-cartoon',
                    representationPresetParams: {{
                        theme: {{ globalName: 'uncertainty' }},
                    }},
                }});
            }} catch (err) {{
                console.error('Mol* viewer:', err);
            }}
        }});
    </script>
</body>
</html>
"""


def create_molstar_html_viewer(
    pdb_content: str,
    pdb_structure,
    selected_idx: int,
    nearby: Set[int],
    contacting: Set[int],
    viewer_id: str = "molstar-viewer",
) -> str:
    """
    Embedded **Mol*** viewer (pinned CDN build) for SCHEMA contact exploration.

    ``pdb_content`` should come from :func:`create_colored_pdb_content`, which encodes
    categories in the **B-factor** column. The viewer colors the default cartoon with
    Mol*'s **uncertainty** theme so bulk, focal+nearby, and contact residues read as
    distinct hues on a **polymer cartoon** (not separate stick layers like the old 3Dmol view).
    """
    import html as html_module
    import json

    num_nearby = len(nearby)
    num_contacts = len(contacting)
    caption_html = (
        f"Alignment column index (0-based) <b>{html_module.escape(str(selected_idx))}</b> — "
        f"{num_nearby} nearby (8 Å) — {num_contacts} SCHEMA contacts. "
        "Colors: B-factor categories via Mol* uncertainty theme."
    )

    pdb_json = json.dumps(pdb_content)
    return _embedded_molstar_single_pdb_html(
        viewer_id=viewer_id,
        pdb_json=pdb_json,
        caption_html=caption_html,
        height_px=620,
        layout_show_sequence=True,
    )


def _effective_alignment_length(pdb_structure, session_parent_length: int) -> int:
    """
    Length of the coordinate system used for aa.index after SCHEMA contact alignment.

    When the structure was renumbered to the parent alignment, indices match columns of
    ``renumbering_seq`` (trimmed query + gaps), not raw PDB residue numbers.
    """
    ren = getattr(pdb_structure, "renumbering_seq", None)
    if ren is not None:
        return len(ren)
    return int(session_parent_length)


def _fragment_id_for_residue(
    residue_index_0b: int, crossovers_1based: List[int], alignment_length: int
) -> int:
    """0-based alignment column -> fragment index from 1-based SCHEMA crossover list."""
    import bisect

    xb = sorted(
        {
            int(c)
            for c in (crossovers_1based or [])
            if 1 <= int(c) <= int(alignment_length)
        }
    )
    starts_0b = [c - 1 for c in xb]
    return bisect.bisect_right(starts_0b, residue_index_0b)


def _fragment_color_hex(fragment_id: int) -> str:
    return CROSSOVER_FRAGMENT_COLORS[fragment_id % len(CROSSOVER_FRAGMENT_COLORS)]


def _sequential_fragment_color_ranges(
    pdb_structure, crossovers_1based: List[int], alignment_length: int
) -> Tuple[List[Dict[str, Any]], int]:
    """
    Contiguous display-order spans (resSeq 1..N) per crossover fragment for 3Dmol.

    Each range is ``{"start": i, "end": j, "color": "#..."}`` with inclusive endpoints.
    """
    sorted_aas = sorted(pdb_structure.amino_acids, key=lambda aa: aa.index)
    if not sorted_aas:
        return [], 1

    xb = sorted(
        {
            int(c)
            for c in (crossovers_1based or [])
            if 1 <= int(c) <= int(alignment_length)
        }
    )
    n_frag = len(xb) + 1

    ranges_out: List[Dict[str, Any]] = []
    display_pos = 0
    cur_fid: Optional[int] = None
    range_start = 0

    for aa in sorted_aas:
        display_pos += 1
        fid = _fragment_id_for_residue(aa.index, crossovers_1based, alignment_length)
        if cur_fid is None:
            cur_fid = fid
            range_start = display_pos
        elif fid != cur_fid:
            ranges_out.append(
                {
                    "start": range_start,
                    "end": display_pos - 1,
                    "color": _fragment_color_hex(cur_fid),
                }
            )
            cur_fid = fid
            range_start = display_pos

    assert cur_fid is not None
    ranges_out.append(
        {
            "start": range_start,
            "end": display_pos,
            "color": _fragment_color_hex(cur_fid),
        }
    )
    return ranges_out, n_frag


def create_crossover_rainbow_pdb_content(
    pdb_structure, crossovers_1based: List[int], session_parent_length: int
) -> str:
    """PDB text with B-factor = fragment id (optional fallback for viewers that color by B)."""
    alen = _effective_alignment_length(pdb_structure, session_parent_length)
    cmap: Dict[int, float] = {}
    for aa in pdb_structure.amino_acids:
        fid = _fragment_id_for_residue(aa.index, crossovers_1based, alen)
        cmap[aa.index] = float(fid)
    return pdb_structure_to_pdb_content(pdb_structure, color_map=cmap)


def display_crossover_fragment_structure(
    pdb_structure,
    crossovers_1based: List[int],
    session_parent_length: int,
    view_width: int = 800,
    view_height: int = 650,
    render_revision: int = 0,
) -> int:
    """
    Streamlit + py3Dmol/stmol: full-chain cartoon with red/green alternating fragments.

    Uses sequential 1..N ``resSeq`` in the viewer PDB so the ribbon is continuous.
    Crossovers are **1-based SCHEMA / alignment** indices (same frame as RASPP and MSA).

    Returns:
        Number of fragments (len(crossovers) + 1).
    """
    if not SCHEMARECOMB_AVAILABLE:
        raise ImportError("schemarecomb is required")

    import py3Dmol
    import streamlit.components.v1 as components

    alignment_len = _effective_alignment_length(pdb_structure, session_parent_length)
    ranges, n_frag = _sequential_fragment_color_ranges(
        pdb_structure, crossovers_1based, alignment_len
    )
    pdb_plain = pdb_structure_to_sequential_cartoon_pdb(pdb_structure, color_map=None)

    chain = "A"
    if pdb_structure.amino_acids:
        chain = pdb_structure.amino_acids[0].atoms[0].chain or "A"

    view = py3Dmol.view(width=view_width, height=view_height)
    view.addModel(pdb_plain, "pdb")
    view.setStyle(
        {},
        {
            "line": {"hidden": True},
            "stick": {"hidden": True},
            "cartoon": {"color": "#cccccc", "opacity": 0.85},
        },
    )

    for span in ranges:
        resi = list(range(int(span["start"]), int(span["end"]) + 1))
        if not resi:
            continue
        sel: Dict[str, Any] = {"resi": resi}
        if chain.strip():
            sel["chain"] = chain
        view.setStyle(sel, {"cartoon": {"color": span["color"], "opacity": 0.95}})

    view.zoomTo()
    html = view._make_html()
    if render_revision:
        html = f"<!-- crossover_revision:{render_revision} -->\n{html}"
    components.html(html, height=view_height, width=view_width)
    return n_frag
