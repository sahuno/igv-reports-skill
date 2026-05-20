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
import sys
from datetime import date
from pathlib import Path

CHECKS_HEADER = ["sample", "track_name", "region", "status",
                 "observed", "expected", "details"]

# FAIL-first severity ordering for sort and section emission.
VERDICT_ORDER = {"FAIL": 0, "REVIEW": 1, "UNVERIFIED": 2, "PASS": 3}


class MalformedChecks(Exception):
    """A checks-TSV data row could not be parsed (too few columns)."""


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
                header_seen = True  # first non-empty line is the header
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


def sort_key(verdict: str, region: str) -> tuple[int, str, int]:
    """FAIL-first severity, then genomic (chrom string, start int)."""
    chrom, start = parse_region(region)
    return (VERDICT_ORDER[verdict], chrom, start)


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
