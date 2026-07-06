# SCHEMA-RASPP Streamlit Application

A web-based interface for SCHEMA-RASPP protein recombination library design.

## Overview

This Streamlit app provides an end-to-end workflow:

1. **SCHEMA Energy** — BLAST/MUSCLE automation or manual PDB+MSA; SCHEMA contacts
2. **RASPP Design** — Multi-fragment RASPP crossover optimization
3. **Crossover Analysis** — Distribution, MSA, structure; **Apply crossover selection**
4. **Assembly Analysis** — Query fragments and Golden Gate overhang assignment
5. **Diversity Analysis** — Homolog MSA analysis; **Save main list to session**
6. **Library Optimization** — Optional ESM2 + ridge regression (requires extra deps)
7. **Simulate with AI** — Optional random chimera ESM2 simulation
8. **Oligopool Design** — BsaI Golden Gate oligo design from saved pools
9. **Projects** — Save, load, export, and import projects

## Local installation

```bash
pip install -r requirements.txt
streamlit run app.py
```

Optional ESM2 features:

```bash
pip install -r requirements-optimization.txt
```

## Streamlit Community Cloud deployment

1. Push this repository to GitHub (do not commit `checkpoints/` — it is in `.gitignore`).
2. Create a new app on [Streamlit Community Cloud](https://streamlit.io/cloud) pointing to `app.py`.
3. Add **Secrets** (Settings → Secrets):

```toml
EBI_EMAIL = "johnbmcarthur@gmail.com"
SCHEMA_HOSTING = "cloud"
```

4. Use base `requirements.txt` only unless your tier has enough RAM for PyTorch/ESM2.
5. **Export projects** from the **Projects** page after each session. Server-side checkpoints are **ephemeral** and are lost on redeploy.

## Required workflow clicks

| Step | Page | Action |
|------|------|--------|
| Apply crossovers | Crossover Analysis | **Apply crossover selection** |
| Save homolog pools | Diversity / Library Optimization | **Save main list to session** |
| Backup project | Projects | **Export** (JSON recommended) |

## Data files

Shipped with the repo:

- `data/gga_aa_pair_compatibility.yaml` — Golden Gate overhang compatibility
- `data/stuffer_sequences.yaml` — Precomputed oligopool stuffer sequences

## Example files

See `examples/` for sample PDB, MSA, and crossover inputs.

## References

- Voigt, C. et al., *Nature Structural Biology* 9(7):553-558 (2002)
- Endelman, J. et al., *Protein Engineering, Design & Selection* 17(7):589-594 (2005)

## License

GPL-3.0 (inherited from SCHEMA-RASPP)
