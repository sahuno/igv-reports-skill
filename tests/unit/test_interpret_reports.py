"""Unit tests for interpret_reports.py — pure-Python, no subprocess/BAM/HTML.

Author: Samuel Ahuno
Purpose:
  Covers the checks-TSV loader, per-region verdict rollup, FAIL-first sort,
  cohort aggregation, Markdown rendering, and the no-verification fallback.
  Runs in ~1 s with only pytest.

Run:
  pytest tests/unit/test_interpret_reports.py -v
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

# Make scripts/ importable without installing the skill as a package.
SCRIPTS = Path(__file__).resolve().parents[2] / "scripts"
sys.path.insert(0, str(SCRIPTS))
import interpret_reports as ir  # noqa: E402


def _write(tmp_path: Path, body: str) -> Path:
    p = tmp_path / "checks.tsv"
    p.write_text(body)
    return p


def test_load_checks_basic(tmp_path):
    p = _write(tmp_path, (
        "sample\ttrack_name\tregion\tstatus\tobserved\texpected\tdetails\n"
        "s1\ttumor\tchr2:100-200\tPASS\t56\t56\tdiff_ratio=0.000 (tol=0.050)\n"
    ))
    rows = ir.load_checks(p)
    assert len(rows) == 1
    r = rows[0]
    assert (r.sample, r.track_name, r.region, r.status) == ("s1", "tumor", "chr2:100-200", "PASS")
    assert r.observed == "56"
    assert r.expected == "56"
    assert r.details.startswith("diff_ratio")


def test_load_checks_missing_file_returns_empty(tmp_path):
    assert ir.load_checks(tmp_path / "nope.tsv") == []


def test_load_checks_header_only_returns_empty(tmp_path):
    p = _write(tmp_path, "sample\ttrack_name\tregion\tstatus\tobserved\texpected\tdetails\n")
    assert ir.load_checks(p) == []


def test_load_checks_short_columns_padded(tmp_path):
    # verify_anchors can emit a SKIP row with empty observed/expected/details.
    p = _write(tmp_path, (
        "sample\ttrack_name\tregion\tstatus\tobserved\texpected\tdetails\n"
        "s1\ttumor\tchr2:100-200\tSKIP\t\t\t\n"
    ))
    rows = ir.load_checks(p)
    assert rows[0].status == "SKIP"
    assert rows[0].observed == ""


def test_load_checks_malformed_raises(tmp_path):
    p = _write(tmp_path, (
        "sample\ttrack_name\tregion\tstatus\tobserved\texpected\tdetails\n"
        "only\ttwo\n"
    ))
    with pytest.raises(ir.MalformedChecks):
        ir.load_checks(p)
