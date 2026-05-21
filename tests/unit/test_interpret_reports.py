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


def test_render_markdown_layout():
    rows = [
        ir.AnchorResult("s1", "tumor", "chr1:100-200", "FAIL", "0", "41", "diff_ratio=0.927 (tol=0.050)"),
        ir.AnchorResult("s1", "tumor", "chr1:300-400", "PASS", "88", "90", "diff_ratio=0.022 (tol=0.050)"),
        ir.AnchorResult("s1", "meth", "chr1:300-400", "FAIL", "0", "14", "diff_ratio=1.000 (tol=0.050)"),
        ir.AnchorResult("s1", "tumor", "chr1:500-600", "PASS", "10", "10", "ok"),
    ]
    sections, notes = ir.aggregate(rows)
    md = ir.render_markdown(sections, notes, "demo_cohort", Path("cohort_verify_anchors.tsv"), "2026-05-20")

    # title + provenance line
    assert md.startswith("# Interpretation — demo_cohort")
    assert "cohort_verify_anchors.tsv" in md
    assert "no new analysis" in md

    # summary table present with the Total row
    assert "| Sample | FAIL | REVIEW | UNVERIFIED | PASS | Total |" in md
    assert "| **Total** | 1 | 1 | 0 | 1 | 3 |" in md

    # FAIL section appears before PASS section within s1
    assert md.index("### FAIL") < md.index("### PASS")

    # FAIL region shows the per-track breakdown line
    assert "**chr1:100-200** — FAIL" in md
    assert "`tumor` — observed 0 vs expected 41 — FAIL — diff_ratio=0.927 (tol=0.050)" in md

    # REVIEW region (chr1:300-400) lists both tracks
    assert "**chr1:300-400** — REVIEW" in md

    # PASS region collapsed to one line, no per-track breakdown
    assert "chr1:500-600 — PASS (1/1 anchors)" in md

    # empty buckets omitted: s1 has no UNVERIFIED region
    assert "### UNVERIFIED" not in md


def test_render_markdown_emits_notes_section():
    rows = [ir.AnchorResult("*", "*", "*", "SKIP", "", "", "no anchors for sample s3")]
    sections, notes = ir.aggregate(rows)
    md = ir.render_markdown(sections, notes, "c", Path("checks.tsv"), "2026-05-20")
    assert "## Notes" in md
    assert "no anchors for sample s3" in md


def _run_cli(monkeypatch, argv):
    monkeypatch.setattr(sys, "argv", ["interpret_reports.py", *argv])
    try:
        ir.main()
        return 0
    except SystemExit as e:
        return e.code or 0


def test_cli_writes_interpretation(monkeypatch, tmp_path):
    checks = tmp_path / "checks.tsv"
    checks.write_text(
        "sample\ttrack_name\tregion\tstatus\tobserved\texpected\tdetails\n"
        "s1\ttumor\tchr1:100-200\tPASS\t10\t10\tok\n"
    )
    out = tmp_path / "reports" / "interpretation.md"
    rc = _run_cli(monkeypatch, ["--checks", str(checks), "--out", str(out)])
    assert rc == 0
    assert out.is_file()
    text = out.read_text()
    assert "# Interpretation —" in text
    # cohort name defaults to the --out parent dir name
    assert "Interpretation — reports" in text


def test_cli_missing_checks_writes_fallback(monkeypatch, tmp_path):
    out = tmp_path / "interpretation.md"
    rc = _run_cli(monkeypatch, ["--checks", str(tmp_path / "nope.tsv"), "--out", str(out)])
    assert rc == 0
    assert "No anchor verification results found" in out.read_text()


def test_cli_malformed_checks_exits_2(monkeypatch, tmp_path):
    checks = tmp_path / "checks.tsv"
    checks.write_text(
        "sample\ttrack_name\tregion\tstatus\tobserved\texpected\tdetails\n"
        "only\ttwo\n"
    )
    out = tmp_path / "interpretation.md"
    rc = _run_cli(monkeypatch, ["--checks", str(checks), "--out", str(out)])
    assert rc == 2
    assert not out.exists()


def test_cli_cohort_name_override(monkeypatch, tmp_path):
    checks = tmp_path / "checks.tsv"
    checks.write_text(
        "sample\ttrack_name\tregion\tstatus\tobserved\texpected\tdetails\n"
        "s1\ttumor\tchr1:100-200\tPASS\t10\t10\tok\n"
    )
    out = tmp_path / "interpretation.md"
    rc = _run_cli(monkeypatch, ["--checks", str(checks), "--out", str(out), "--cohort-name", "ATLL run 3"])
    assert rc == 0
    assert "Interpretation — ATLL run 3" in out.read_text()


# ---------------------------------------------------------------------------
# Fix 1 — load_checks header validation
# ---------------------------------------------------------------------------

def test_load_checks_no_header_raises(tmp_path):
    """A file whose first line is a data row (no header) raises MalformedChecks."""
    p = _write(tmp_path, (
        "s1\ttumor\tchr2:100-200\tPASS\t56\t56\tdiff_ratio=0.000\n"
    ))
    with pytest.raises(ir.MalformedChecks, match="unexpected header"):
        ir.load_checks(p)


def test_load_checks_hash_prefixed_header_accepted(tmp_path):
    """A valid header preceded by a '#' is tolerated."""
    p = _write(tmp_path, (
        "#sample\ttrack_name\tregion\tstatus\tobserved\texpected\tdetails\n"
        "s1\ttumor\tchr2:100-200\tPASS\t56\t56\tok\n"
    ))
    rows = ir.load_checks(p)
    assert len(rows) == 1
    assert rows[0].sample == "s1"


def test_load_checks_valid_header_only_returns_empty(tmp_path):
    """A file with only the valid header returns [] (no data rows)."""
    p = _write(tmp_path, "sample\ttrack_name\tregion\tstatus\tobserved\texpected\tdetails\n")
    assert ir.load_checks(p) == []


# ---------------------------------------------------------------------------
# Fix 2 — sort_key raises ValueError on unknown verdict
# ---------------------------------------------------------------------------

def test_sort_key_unknown_verdict_raises_value_error():
    with pytest.raises(ValueError, match="unknown verdict"):
        ir.sort_key("BOGUS", "chr1:1-2")


# ---------------------------------------------------------------------------
# Fix 3 — natural chromosome ordering (chr2 < chr10)
# ---------------------------------------------------------------------------

def test_sort_key_natural_chrom_order():
    """chr2 must sort before chr10 when verdict is the same."""
    entries = [
        ("chr10:100-200", "FAIL"),
        ("chr2:100-200", "FAIL"),
    ]
    ordered = sorted(entries, key=lambda e: ir.sort_key(e[1], e[0]))
    assert ordered[0][0] == "chr2:100-200"
    assert ordered[1][0] == "chr10:100-200"


# ---------------------------------------------------------------------------
# Fix 4 — SKIP rows render without double-space artifact
# ---------------------------------------------------------------------------

def test_render_markdown_skip_row_no_double_space():
    """UNVERIFIED (all-SKIP) region renders clean '— SKIP —' lines."""
    rows = [
        ir.AnchorResult("s1", "tumor", "chr1:100-200", "SKIP", "", "", "not rendered"),
    ]
    sections, notes = ir.aggregate(rows)
    md = ir.render_markdown(sections, notes, "test", Path("checks.tsv"), "2026-05-20")
    # Must contain the clean SKIP line
    assert "- `tumor` — SKIP — not rendered" in md
    # Must NOT contain double-space artifact patterns
    assert "observed  vs" not in md
    assert "expected  —" not in md


# ---------------------------------------------------------------------------
# Fix 6 — happy-path tests
# ---------------------------------------------------------------------------

def test_render_markdown_no_notes_section_when_empty():
    """render_markdown output does NOT contain '## Notes' when notes=[]."""
    rows = [
        ir.AnchorResult("s1", "tumor", "chr1:100-200", "PASS", "10", "10", "ok"),
    ]
    sections, notes = ir.aggregate(rows)
    assert notes == []
    md = ir.render_markdown(sections, notes, "c", Path("checks.tsv"), "2026-05-20")
    assert "## Notes" not in md


def test_render_markdown_all_pass_no_fail_review_unverified_headings():
    """An all-PASS cohort renders no ### FAIL / ### REVIEW / ### UNVERIFIED."""
    rows = [
        ir.AnchorResult("s1", "tumor", "chr1:100-200", "PASS", "10", "10", "ok"),
        ir.AnchorResult("s1", "meth",  "chr1:100-200", "PASS", "5",  "5",  "ok"),
        ir.AnchorResult("s2", "tumor", "chr2:200-300", "PASS", "8",  "8",  "ok"),
    ]
    sections, notes = ir.aggregate(rows)
    md = ir.render_markdown(sections, notes, "c", Path("checks.tsv"), "2026-05-20")
    assert "### FAIL" not in md
    assert "### REVIEW" not in md
    assert "### UNVERIFIED" not in md
