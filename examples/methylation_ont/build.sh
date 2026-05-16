#!/usr/bin/env bash
# build.sh — runnable ONT 5mC/5hmC methylation viewer.
#
# Builds a single self-contained HTML showing per-read base-mod-colored BAM
# tracks plus per-sample 5mC and 5hmC bedGraph tracks (locked y-axis 0..100)
# for a small set of pre-defined promoter windows.
#
# Two paths:
#   A) Recommended — generate the tracks.json from tracks_spec.example.yaml,
#      then call the skill driver with --track-config (driver auto-picks
#      conda or SIF based on SLURM_JOB_ID).
#   B) Direct — call create_report from inside the dedicated igv-reports
#      SIF (/data1/greenbab/users/ahunos/apps/containers/igv-reports_1.16.0.sif)
#      with a pre-built tracks.json. See recipe.md for the verbatim command.

set -euo pipefail

EX_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SKILL_DIR="$(cd "${EX_DIR}/../.." && pwd)"
GEN_TRACKS="${SKILL_DIR}/scripts/generate_tracks_json.py"
BUILD="${SKILL_DIR}/scripts/build_igvreports.py"

# --- edit these for your run --------------------------------------------------
RUN_DIR="${EX_DIR}"                                  # where slices/ live
SITES="${EX_DIR}/sites.hg38.example.bed"
SPEC="${EX_DIR}/tracks_spec.example.yaml"
GENOME=hg38
OUT="${EX_DIR}/methylation_report.${GENOME}.html"
TITLE="COLO829 promoter methylation: DNMT3A_2 + EZH2"
# --- end edit -----------------------------------------------------------------

source /home/ahunos/miniforge3/etc/profile.d/conda.sh
conda activate snakemake

# A) Generate tracks.json from the YAML spec.
TRACKS_JSON="${EX_DIR}/tracks.json"
python "${GEN_TRACKS}" \
    --spec "${SPEC}" \
    --run-dir "${RUN_DIR}" \
    --out "${TRACKS_JSON}"

# B) Build the report.
python "${BUILD}" \
    --sites         "${SITES}" \
    --track-config  "${TRACKS_JSON}" \
    --genome        "${GENOME}" \
    --flanking      0 \
    --type          mutation \
    --info-columns  name \
    --title         "${TITLE}" \
    --output        "${OUT}"

echo "Wrote: ${OUT}"
