#!/usr/bin/env bash
# cohort_samplesheet.sh — multi-sample run from a TSV samplesheet.
#
# Produces one HTML per row + an index.html landing page.
# Mirrors the layout used in the ATLL HTLV-1 cohort run:
#   results/<run>/
#     ├── inputs/<sample>/sites.<genome>.bed
#     ├── reports/<sample>.<genome>.html
#     ├── reports/index.html
#     └── logs/run_<timestamp>.log

set -euo pipefail

SKILL_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD="${SKILL_DIR}/scripts/build_igvreports.py"

source /home/ahunos/miniforge3/etc/profile.d/conda.sh
conda activate snakemake

# Samplesheet TSV layout (header required, tab-separated):
#   sample <TAB> bam_tumor <TAB> bam_normal <TAB> vcf <TAB> sites_bed [<TAB> extra_tracks]
#
# extra_tracks (optional) is comma-separated and appended after VCF,
# before default annotation tracks.

# --- edit these for your run ---
SAMPLESHEET=samplesheet.tsv
GENOME=hg38
OUTDIR=results/cohort/reports
# --- end edit ---

python "${BUILD}" \
    --samplesheet "${SAMPLESHEET}" \
    --genome      "${GENOME}" \
    --output-dir  "${OUTDIR}" \
    --fail-on-fail   # auto-invokes verify_cohort.py at the end; exits nonzero if any check FAILs.
                     # Add --no-verify to skip verification entirely (not recommended for cohort runs).

echo "Cohort index: ${OUTDIR}/index.html"
echo "Verification: ${OUTDIR}/cohort_verify.tsv + cohort_verify.summary.md"
