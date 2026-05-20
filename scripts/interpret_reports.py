#!/usr/bin/env python3
"""interpret_reports.py — triage companion for create_report cohorts.

Reads an anchor-verification checks TSV (produced by
`verify_anchors.py verify-cohort` / `verify`) and rolls per-track
PASS/FAIL/SKIP results up into a per-region verdict, written as a single
cohort-wide `interpretation.md`. Lead with verdicts, sort FAIL-first, give
just enough per-track detail to triage without opening the HTML.

This script interprets output of `igv-reports` (`create_report`) maintained
by the IGV team at the Broad Institute (MIT, © 2018-2019 The Broad Institute
and The Regents of the University of California;
<https://github.com/igvteam/igv-reports>). The triage/interpretation layer
in this file is added by this repo, not by upstream. See CREDITS.md.

Author: Samuel Ahuno
Purpose:
  The verifiers (verify_report / verify_cohort / verify_anchors) produce
  machine-first artifacts: one row per (sample, track, region). A researcher
  opening N HTMLs has no single place that answers "which regions need a
  closer look?". This script is the human-first reading surface: one file,
  opened once, regions ordered so the ones needing attention sit on top.

  Verdicts derive SOLELY from the checks TSV's status column (no new
  thresholds, no re-analysis):
    PASS        — every non-SKIP anchor for the region is PASS
    FAIL        — every non-SKIP anchor for the region is FAIL
    REVIEW      — mixed PASS and FAIL
    UNVERIFIED  — all anchors for the region are SKIP (region not rendered,
                  or no tracks matched) — distinct from FAIL: nothing checked.

  Input: the checks TSV written by verify_anchors.write_checks, columns:
    sample  track_name  region  status  observed  expected  details
  (header line is NOT '#'-prefixed; region is 'chrom:start-end').

Typical use:
  python interpret_reports.py \\
      --checks results/<run>/cohort_verify_anchors.tsv \\
      --out    results/<run>/reports/interpretation.md

Skill location:
  /data1/greenbab/users/ahunos/apps/llm_configs/claude/skills/igv-reports/
"""

from __future__ import annotations

import argparse
import dataclasses
import re
import sys
from datetime import date
from pathlib import Path

CHECKS_HEADER = ["sample", "track_name", "region", "status",
                 "observed", "expected", "details"]
EXPECTED_HEADER = "\t".join(CHECKS_HEADER)

# FAIL-first severity ordering for sort and section emission.
VERDICT_ORDER = {"FAIL": 0, "REVIEW": 1, "UNVERIFIED": 2, "PASS": 3}


class MalformedChecks(Exception):
    """A checks-TSV data row could not be parsed (too few columns)."""


# Mirrors verify_anchors.AnchorCheck field-for-field; kept local to avoid
# importing verify_anchors (which pulls in subprocess/samtools resolution).
@dataclasses.dataclass
class AnchorResult:
    sample: str
    track_name: str
    region: str
    status: str
    observed: str = ""
    expected: str = ""
    details: str = ""


def load_checks(path: Path) -> list[AnchorResult]:
    """Parse a verify_anchors checks TSV into AnchorResult rows.

    Returns [] when the file is absent or header-only (the no-verification
    fallback case). Raises MalformedChecks on a data row with < 4 columns,
    consistent with verify_anchors.load_anchors' fail-loud parsing."""
    if not path.exists():
        return []
    results: list[AnchorResult] = []
    with path.open() as fh:
        header_seen = False
        for i, line in enumerate(fh, start=1):
            line = line.rstrip("\n")
            if not line:
                continue
            if not header_seen:
                header_seen = True
                # Validate header content (tolerate a leading '#').
                stripped = line.lstrip("#")
                if stripped != EXPECTED_HEADER:
                    raise MalformedChecks(
                        f"{path}:1: unexpected header: {line!r} "
                        f"(expected {EXPECTED_HEADER!r})"
                    )
                continue
            cols = line.split("\t")
            if len(cols) < 4:
                raise MalformedChecks(
                    f"{path}:{i}: expected >= 4 tab-separated columns "
                    f"({CHECKS_HEADER}), got {len(cols)}: {cols!r}"
                )
            if len(cols) < len(CHECKS_HEADER):
                cols += [""] * (len(CHECKS_HEADER) - len(cols))
            d = dict(zip(CHECKS_HEADER, cols))
            results.append(AnchorResult(
                sample=d["sample"], track_name=d["track_name"],
                region=d["region"], status=d["status"],
                observed=d["observed"], expected=d["expected"],
                details=d["details"],
            ))
    return results


def region_verdict(results: list[AnchorResult]) -> str:
    """Roll a single region's per-track results into one verdict.

    UNVERIFIED if every anchor is SKIP; otherwise consider only non-SKIP
    anchors: PASS if all PASS, FAIL if all FAIL, REVIEW if mixed."""
    non_skip = [r for r in results if r.status != "SKIP"]
    if not non_skip:
        return "UNVERIFIED"
    statuses = {r.status for r in non_skip}
    if statuses == {"PASS"}:
        return "PASS"
    if statuses == {"FAIL"}:
        return "FAIL"
    return "REVIEW"


def parse_region(region: str) -> tuple[str, int]:
    """Split 'chrom:start-end' into (chrom, start) for genomic ordering.
    Falls back to (region, 0) when the string isn't in that shape."""
    try:
        chrom, span = region.rsplit(":", 1)
        start = int(span.split("-", 1)[0])
        return chrom, start
    except (ValueError, IndexError):
        return region, 0


_CHR_NUM_RE = re.compile(r"(\d+)$")


def _chrom_sort_key(chrom: str) -> tuple[int, str]:
    """Natural-ish chrom ordering: trailing integer first (chr2 < chr10), then
    the raw string for non-numeric chroms (chrX, chrM, contigs)."""
    m = _CHR_NUM_RE.search(chrom)
    return (int(m.group(1)) if m else 10**9, chrom)


def sort_key(verdict: str, region: str) -> tuple[int, int, str, int]:
    """FAIL-first severity, then natural genomic order (chrom, start)."""
    order = VERDICT_ORDER.get(verdict)
    if order is None:
        raise ValueError(
            f"sort_key: unknown verdict {verdict!r} "
            f"(known: {list(VERDICT_ORDER)})"
        )
    chrom, start = parse_region(region)
    return (order, *_chrom_sort_key(chrom), start)


@dataclasses.dataclass
class RegionEntry:
    region: str
    verdict: str
    results: list[AnchorResult]


@dataclasses.dataclass
class SampleSection:
    sample: str
    entries: list[RegionEntry]            # sorted FAIL-first then genomic
    counts: dict[str, int]                # verdict -> region count


def aggregate(
    results: list[AnchorResult],
) -> tuple[list[SampleSection], list[AnchorResult]]:
    """Group results into per-sample sections with per-region verdicts and
    summary counts. Returns (sections, notes); `notes` holds the
    verify-cohort housekeeping rows (sample='*' or region='*') that aren't
    region verdicts. Sections are sorted by sample name; entries within a
    section are sorted FAIL-first then genomic."""
    notes = [r for r in results if r.sample == "*" or r.region == "*"]
    real = [r for r in results if r.sample != "*" and r.region != "*"]

    by_sample: dict[str, dict[str, list[AnchorResult]]] = {}
    for r in real:
        by_sample.setdefault(r.sample, {}).setdefault(r.region, []).append(r)

    sections: list[SampleSection] = []
    for sample in sorted(by_sample):
        counts = {"FAIL": 0, "REVIEW": 0, "UNVERIFIED": 0, "PASS": 0}
        entries: list[RegionEntry] = []
        for region, rs in by_sample[sample].items():
            verdict = region_verdict(rs)
            counts[verdict] += 1
            entries.append(RegionEntry(region, verdict, rs))
        entries.sort(key=lambda e: sort_key(e.verdict, e.region))
        sections.append(SampleSection(sample, entries, counts))
    return sections, notes


def render_markdown(
    sections: list[SampleSection],
    notes: list[AnchorResult],
    cohort_name: str,
    source_path: Path,
    today: str,
) -> str:
    """Render the cohort-wide interpretation.md. FAIL/REVIEW/UNVERIFIED
    regions get a per-track breakdown; PASS regions collapse to one line.
    Empty verdict buckets are omitted per sample."""
    out: list[str] = []
    out.append(f"# Interpretation — {cohort_name}")
    out.append("")
    out.append(
        f"_Generated {today} from anchor verification (`{source_path.name}`). "
        "Verdicts derive solely from anchor pass/fail; this file performs no "
        "new analysis. Open the per-sample HTML report for any REVIEW / FAIL "
        "/ UNVERIFIED region._"
    )
    out.append("")

    # --- summary table ---
    out.append("## Summary")
    out.append("")
    out.append("| Sample | FAIL | REVIEW | UNVERIFIED | PASS | Total |")
    out.append("|--------|-----:|-------:|-----------:|-----:|------:|")
    total = {"FAIL": 0, "REVIEW": 0, "UNVERIFIED": 0, "PASS": 0}
    for s in sections:
        c = s.counts
        n = sum(c.values())
        out.append(
            f"| {s.sample} | {c['FAIL']} | {c['REVIEW']} | "
            f"{c['UNVERIFIED']} | {c['PASS']} | {n} |"
        )
        for k in total:
            total[k] += c[k]
    grand = sum(total.values())
    out.append(
        f"| **Total** | {total['FAIL']} | {total['REVIEW']} | "
        f"{total['UNVERIFIED']} | {total['PASS']} | {grand} |"
    )
    out.append("")

    # --- per-sample sections ---
    for s in sections:
        out.append(f"## {s.sample}")
        out.append("")
        for verdict in ("FAIL", "REVIEW", "UNVERIFIED", "PASS"):
            bucket = [e for e in s.entries if e.verdict == verdict]
            if not bucket:
                continue
            out.append(f"### {verdict}")
            for e in bucket:
                if verdict == "PASS":
                    n_ok = sum(1 for r in e.results if r.status != "SKIP")
                    out.append(f"- {e.region} — PASS ({n_ok}/{n_ok} anchors)")
                else:
                    out.append(f"- **{e.region}** — {verdict}")
                    for r in e.results:
                        if r.status == "SKIP":
                            out.append(f"  - `{r.track_name}` — SKIP — {r.details}")
                            continue
                        out.append(
                            f"  - `{r.track_name}` — observed {r.observed} "
                            f"vs expected {r.expected} — {r.status} — {r.details}"
                        )
            out.append("")

    # --- notes (verify-cohort housekeeping rows) ---
    if notes:
        out.append("## Notes")
        out.append("")
        for r in notes:
            out.append(f"- {r.sample}/{r.region}: {r.details or r.status}")
        out.append("")

    return "\n".join(out).rstrip("\n") + "\n"


def render_fallback(cohort_name: str, source_path: Path) -> str:
    """Minimal file written when no anchor verification results exist."""
    return (
        f"# Interpretation — {cohort_name}\n\n"
        f"_No anchor verification results found (`{source_path}`). Verdicts "
        "require `verify_anchors.py verify-cohort` to have run first. Open the "
        "HTML reports directly._\n"
    )


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--checks", required=True,
                    help="verify_anchors verify-cohort/verify checks TSV")
    ap.add_argument("--out", required=True,
                    help="path to write interpretation.md")
    ap.add_argument("--cohort-name", default=None,
                    help="title string (default: --out parent directory name)")
    args = ap.parse_args()

    checks_path = Path(args.checks)
    out_path = Path(args.out)
    cohort_name = args.cohort_name or out_path.resolve().parent.name

    try:
        results = load_checks(checks_path)
    except MalformedChecks as e:
        sys.stderr.write(f"ERROR: {e}\n")
        sys.exit(2)

    out_path.parent.mkdir(parents=True, exist_ok=True)

    if not results:
        out_path.write_text(render_fallback(cohort_name, checks_path))
        sys.stderr.write(
            f"[interpret_reports] no anchor results in {checks_path} — "
            f"wrote fallback {out_path}\n"
        )
        sys.exit(0)

    sections, notes = aggregate(results)
    md = render_markdown(sections, notes, cohort_name, checks_path,
                         date.today().isoformat())
    out_path.write_text(md)
    n_regions = sum(sum(s.counts.values()) for s in sections)
    sys.stderr.write(
        f"[interpret_reports] wrote {out_path} "
        f"({len(sections)} samples, {n_regions} regions)\n"
    )
    sys.exit(0)


if __name__ == "__main__":
    main()
