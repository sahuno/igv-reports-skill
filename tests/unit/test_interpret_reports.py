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


def _res(status):
    return ir.AnchorResult("s", "t", "chr1:1-2", status)


def test_region_verdict_all_pass():
    assert ir.region_verdict([_res("PASS"), _res("PASS")]) == "PASS"


def test_region_verdict_all_fail():
    assert ir.region_verdict([_res("FAIL"), _res("FAIL")]) == "FAIL"


def test_region_verdict_mixed_is_review():
    assert ir.region_verdict([_res("PASS"), _res("FAIL")]) == "REVIEW"


def test_region_verdict_all_skip_is_unverified():
    assert ir.region_verdict([_res("SKIP"), _res("SKIP")]) == "UNVERIFIED"


def test_region_verdict_skip_ignored_when_others_present():
    # PASS + SKIP -> PASS (SKIP doesn't downgrade a clean region)
    assert ir.region_verdict([_res("PASS"), _res("SKIP")]) == "PASS"
    # FAIL + SKIP -> FAIL
    assert ir.region_verdict([_res("FAIL"), _res("SKIP")]) == "FAIL"


def test_parse_region_normal():
    assert ir.parse_region("chr7:148884000-148884001") == ("chr7", 148884000)


def test_parse_region_unparseable_falls_back():
    assert ir.parse_region("weird") == ("weird", 0)


def test_sort_key_orders_fail_first_then_genomic():
    entries = [
        ("chr1:500-600", "PASS"),
        ("chr1:100-200", "FAIL"),
        ("chr1:300-400", "FAIL"),
        ("chr1:700-800", "UNVERIFIED"),
        ("chr1:900-999", "REVIEW"),
    ]
    ordered = sorted(entries, key=lambda e: ir.sort_key(e[1], e[0]))
    assert [v for _, v in ordered] == ["FAIL", "FAIL", "REVIEW", "UNVERIFIED", "PASS"]
    # the two FAILs keep genomic order (100 before 300)
    assert ordered[0][0] == "chr1:100-200"
    assert ordered[1][0] == "chr1:300-400"


def test_aggregate_groups_counts_and_sorts():
    rows = [
        ir.AnchorResult("s1", "tumor", "chr1:500-600", "PASS", "10", "10", "ok"),
        ir.AnchorResult("s1", "tumor", "chr1:100-200", "FAIL", "0", "41", "diff"),
        ir.AnchorResult("s1", "meth", "chr1:300-400", "PASS", "9", "9", "ok"),
        ir.AnchorResult("s1", "tumor", "chr1:300-400", "FAIL", "0", "9", "diff"),  # -> REVIEW
        ir.AnchorResult("s2", "tumor", "chr2:1-2", "SKIP", "", "", "not rendered"),  # -> UNVERIFIED
        ir.AnchorResult("*", "*", "*", "SKIP", "", "", "no anchors for sample s3"),  # note
    ]
    sections, notes = ir.aggregate(rows)

    # one note row routed aside, not a region verdict
    assert len(notes) == 1
    assert notes[0].details.startswith("no anchors")

    # two samples, sorted by name
    assert [s.sample for s in sections] == ["s1", "s2"]

    s1 = sections[0]
    # counts: FAIL chr1:100-200, REVIEW chr1:300-400, PASS chr1:500-600
    assert s1.counts == {"FAIL": 1, "REVIEW": 1, "UNVERIFIED": 0, "PASS": 1}
    # FAIL-first ordering
    assert [e.verdict for e in s1.entries] == ["FAIL", "REVIEW", "PASS"]
    assert s1.entries[0].region == "chr1:100-200"

    s2 = sections[1]
    assert s2.counts == {"FAIL": 0, "REVIEW": 0, "UNVERIFIED": 1, "PASS": 0}
    assert s2.entries[0].verdict == "UNVERIFIED"
