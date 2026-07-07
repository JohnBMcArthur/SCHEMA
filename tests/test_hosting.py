"""Smoke tests for SCHEMA-RASPP hosting readiness."""

from pathlib import Path

import pytest


def test_import_core_modules():
    import utils.config  # noqa: F401
    import utils.workflow_state  # noqa: F401
    import utils.gga_overhangs  # noqa: F401
    import utils.oligopool_design  # noqa: F401
    import utils.stuffer_design  # noqa: F401


def test_sanitize_project_name():
    from utils.session_manager import _sanitize_project_name

    assert _sanitize_project_name("../evil") == "evil"
    assert _sanitize_project_name("My Project") == "My_Project"


def test_gga_fragment1_xatg():
    from utils.gga_overhangs import assign_golden_gate_overhangs

    fragments = [{"index": 1, "sequence": "MKLV"}]
    result = assign_golden_gate_overhangs(fragments, compatibility={"pairs": []})
    assert "fragments" in result


def test_stuffer_library_lookup():
    from utils.stuffer_design import design_stuffer_sequence

    seq, warnings = design_stuffer_sequence(10, seed=1)
    assert len(seq) == 10
    assert set(seq) <= set("ATGC")


def test_data_files_exist():
    root = Path(__file__).resolve().parent.parent
    assert (root / "data" / "gga_aa_pair_compatibility.yaml").is_file()
    assert (root / "data" / "stuffer_sequences.yaml").is_file()
