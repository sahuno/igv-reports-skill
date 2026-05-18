#!/usr/bin/env bash
# prep_track_demo.sh — convert a plain-gzip GFF3/GTF/BED.gz to bgzip+tabix.
#
# Run once per misconfigured track. The original file is preserved at
# <name>.bak.original_gzip; the canonical filename is replaced with the
# new bgzip + .tbi sidecar so igv-reports / IGV / tabix can index it.
#
# Tells you immediately whether the file was already bgzip (in which case
# only the .tbi index is rebuilt) or needed full conversion.

set -euo pipefail

SKILL_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PREP="${SKILL_DIR}/scripts/prep_track.sh"

# bgzip / tabix come from htslib in the snakemake conda env.
source /home/ahunos/miniforge3/etc/profile.d/conda.sh
conda activate snakemake

# --- edit this for your run ---
TRACK=/data1/greenbab/database/gencode_annotations/hg38/gencode.v47.annotation.gff3.gz
# --- end edit ---

bash "${PREP}" "${TRACK}"
