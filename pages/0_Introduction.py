"""
Introduction — SCHEMA-RASPP overview and workflow guide.
"""

import streamlit as st
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parent.parent))

from utils.session_manager import init_session_state
from utils.workflow_state import render_cloud_persistence_notice
from utils.config import is_cloud_hosting

st.set_page_config(
    page_title="Introduction",
    page_icon="🧬",
    layout="wide",
)

init_session_state()

st.title("🧬 JBM-SCHEMA-RASPP")
st.markdown("### Structure-Guided Protein Recombination Library Design")

if is_cloud_hosting():
    render_cloud_persistence_notice()

st.markdown("""
This app assists in generating low-cost, high-quality libraries of chimeric proteins.
Crossover sites are chosen with assistance from iterative SCHEMA-RASPP scans. Diversity
within each block is chosen from natural homologs and then filtered for:

1. Compatibility with one-pot, combinatorial Golden Gate Assembly
2. User-defined length and diversity constraints
3. Regression analysis of simulated sequences scored by the Protein Language Model ESM2-150M

Fragments are packaged into oligos, flanked by user-defined amplification tags, which can
be ordered from Twist or other oligopool synthesis companies.

**Workflow (sidebar order):**

| Step | Page | Summary |
|------|------|---------|
| 1 | **SCHEMA Energy** | Small parent MSA + structure → SCHEMA contacts and query sequence |
| 2 | **RASPP Design** | Multi-fragment RASPP scan for low-energy crossover regions |
| 3 | **Crossover Analysis** | Pick crossover positions; **Apply** to commit for downstream |
| 4 | **Assembly Analysis** | Split query into fragments; assign Golden Gate overhangs |
| 5 | **Diversity Analysis** | **Large** homolog MSA → per-fragment pools; **Save main list** |
| 6 | **Library Optimization** | *(Optional)* ESM2 + ridge regression to prune pools; **Save filtered list** |
| 7 | **Simulate with AI** | *(Optional)* Random chimera ESM2 simulation vs wildtype |
| 8 | **Oligopool Design** | BsaI Golden Gate oligos from saved pools → order synthesis |

**Projects** — Save, load, **export**, and import checkpoints (not part of the numbered
workflow). On Streamlit Cloud, export after each session.

**Important steps (easy to miss):**

- On **Crossover Analysis**, selections are not active until you click **Apply crossover selection**.
- The SCHEMA-RASPP MSA (steps 1–3) is intentionally **small**; upload a **much larger**
  homolog MSA on **Diversity Analysis** for library design.
- On **Diversity Analysis** / **Library Optimization**, filtered pools are not used downstream
  until you click **Save main list to session** (or **Save filtered list to session**).
- **Library Optimization** and **Simulate with AI** are optional; you can go straight from
  Diversity to **Oligopool Design** if you do not need ESM-based pruning.
""")
