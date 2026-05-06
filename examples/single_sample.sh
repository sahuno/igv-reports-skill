#!/usr/bin/env bash
# single_sample.sh — minimal one-HTML build.
#
# Sites BED + tumor BAM + matched-normal BAM + VCF → one self-contained
# HTML report at ./reports/sample.hg38.html, with default tracks
# (CpG islands, gencode, RepeatMasker) auto-resolved from
# databases_config.yaml.

set -euo pipefail

SKILL_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD="${SKILL_DIR}/scripts/build_igvreports.py"

# Activate the snakemake conda env (where create_report lives).
source /home/ahunos/miniforge3/etc/profile.d/conda.sh
conda activate snakemake

# --- edit these for your run ---
SITES=results/inputs/sites.hg38.bed       # plain BED, headerless: chr/start/end/name
TUMOR_BAM=path/to/tumor.bam
NORMAL_BAM=path/to/normal.bam
VCF=path/to/calls.vcf
GENOME=hg38
OUT=results/reports/sample.${GENOME}.html
# --- end edit ---

python "${BUILD}" \
    --sites      "${SITES}" \
    --bam        "${TUMOR_BAM}" "${NORMAL_BAM}" \
    --vcf        "${VCF}" \
    --genome     "${GENOME}" \
    --output     "${OUT}"
