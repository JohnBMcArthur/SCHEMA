"""
Introduction — SCHEMA-RASPP overview and workflow guide.
"""

import streamlit as st

st.set_page_config(
    page_title="Introduction",
    page_icon="🧬",
    layout="wide",
)

st.title("🧬 JBM-SCHEMA-RASPP")
st.markdown("### Structure-Guided Protein Recombination Library Design")

st.markdown("""
Welcome to the JBM-SCHEMA-RASPP application! This tool helps you design protein
recombination libraries using structure-guided approaches.

**Navigate using the sidebar to:**

1. **SCHEMA Energy** — Search for diversity and define SCHEMA used to calculate recombination "energy" penalties.
2. **RASPP Design** — Design optimal recombination crossovers.
3. **Crossover Analysis** — Crossover distribution, MSA, and structure by fragment
4. **Assembly Analysis** — Query sequence fragments between applied crossovers
5. **Diversity Analysis** — Homolog fragment variants from an uploaded MSA
6. **Library Optimization** — Iterative ESM2 + ridge regression on fragment pools
7. **Simulate with AI** — Random chimera ESM2 simulation vs wildtype
8. **Oligopool Design** — DNA oligos from saved fragment lists
9. **Projects** — Save and load SCHEMA-RASPP projects
""")
