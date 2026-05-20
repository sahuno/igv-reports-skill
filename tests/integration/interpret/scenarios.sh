#!/usr/bin/env bash
# scenarios.sh — end-to-end CLI test for scripts/interpret_reports.py.
#
# Author: Samuel Ahuno
# Purpose:
#   interpret_reports.py consumes only a checks TSV (no BAMs, no HTML), so
#   this scenario is fully portable and runs in CI. It writes a synthetic
#   checks TSV exercising every verdict bucket, runs interpret_reports, and
#   asserts: (a) the file is written, (b) the summary Total row equals the
#   region count, (c) FAIL sorts above PASS, (d) the malformed-input path
#   exits 2, (e) the missing-input path writes the fallback and exits 0.
#
# Runtime: < 1 s. No external data dependency — never SKIPs.
set -euo pipefail

EX_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SKILL_DIR="$(cd "${EX_DIR}/../../.." && pwd)"
INTERPRET="${SKILL_DIR}/scripts/interpret_reports.py"

PY="${IGV_REPORTS_PY:-}"
if [[ -z "${PY}" ]]; then
    if [[ -x /home/ahunos/miniforge3/envs/snakemake/bin/python ]]; then
        PY=/home/ahunos/miniforge3/envs/snakemake/bin/python
    else
        PY=$(command -v python3)
    fi
fi

WORK="$(mktemp -d)"
trap 'rc=$?; if [[ $rc -eq 0 ]]; then rm -rf "${WORK}"; else echo "(interpret scenarios.sh exited $rc — artifacts in ${WORK})"; fi' EXIT

CHECKS="${WORK}/cohort_verify_anchors.tsv"
OUT="${WORK}/reports/interpretation.md"

# --- synthetic checks TSV: 3 regions in s1 (FAIL, REVIEW, PASS) + 1 note ---
cat >"${CHECKS}" <<'EOF'
sample	track_name	region	status	observed	expected	details
s1	tumor	chr1:100-200	FAIL	0	41	diff_ratio=1.000 (tol=0.050)
s1	tumor	chr1:300-400	PASS	88	90	diff_ratio=0.022 (tol=0.050)
s1	meth	chr1:300-400	FAIL	0	14	diff_ratio=1.000 (tol=0.050)
s1	tumor	chr1:500-600	PASS	10	10	diff_ratio=0.000 (tol=0.050)
*	*	*	SKIP			no anchors for sample s2
EOF

echo "=== scenario 0: render interpretation.md from a clean checks TSV ==="
"${PY}" "${INTERPRET}" --checks "${CHECKS}" --out "${OUT}"
[[ -f "${OUT}" ]] || { echo "  FAIL: ${OUT} not written"; exit 1; }

# s1 has 3 regions (chr1:300-400 collapses tumor+meth into one REVIEW region)
if ! grep -q "| s1 | 1 | 1 | 0 | 1 | 3 |" "${OUT}"; then
    echo "  FAIL: summary row for s1 wrong"; grep '| s1 ' "${OUT}"; exit 1
fi
echo "  OK   summary: s1 = 1 FAIL, 1 REVIEW, 1 PASS, total 3"

# FAIL section must appear before PASS section
fail_line=$(grep -n '^### FAIL' "${OUT}" | head -1 | cut -d: -f1)
pass_line=$(grep -n '^### PASS' "${OUT}" | head -1 | cut -d: -f1)
if [[ -z "${fail_line}" || -z "${pass_line}" || "${fail_line}" -ge "${pass_line}" ]]; then
    echo "  FAIL: FAIL section not above PASS section (FAIL=${fail_line:-none} PASS=${pass_line:-none})"; exit 1
fi
echo "  OK   FAIL section sorts above PASS section"

# note row routed to Notes
grep -q '^## Notes' "${OUT}" || { echo "  FAIL: Notes section missing"; exit 1; }
grep -q 'no anchors for sample s2' "${OUT}" || { echo "  FAIL: note text missing"; exit 1; }
echo "  OK   verify-cohort housekeeping row routed to Notes"
echo

echo "=== scenario 1: malformed checks TSV exits 2 ==="
BAD="${WORK}/bad.tsv"
printf 'sample\ttrack_name\tregion\tstatus\tobserved\texpected\tdetails\nonly\ttwo\n' >"${BAD}"
set +e
"${PY}" "${INTERPRET}" --checks "${BAD}" --out "${WORK}/bad.md" 2>/dev/null
rc=$?
set -e
[[ "${rc}" -eq 2 ]] || { echo "  FAIL: expected exit 2, got ${rc}"; exit 1; }
[[ ! -f "${WORK}/bad.md" ]] || { echo "  FAIL: output written despite malformed input"; exit 1; }
echo "  OK   malformed input -> exit 2, no output"
echo

echo "=== scenario 2: missing checks TSV writes fallback, exits 0 ==="
"${PY}" "${INTERPRET}" --checks "${WORK}/missing.tsv" --out "${WORK}/fallback.md"
grep -q 'No anchor verification results found' "${WORK}/fallback.md" \
    || { echo "  FAIL: fallback text missing"; exit 1; }
echo "  OK   missing input -> fallback file, exit 0"
echo

echo "=== all interpret scenarios PASSED ==="
