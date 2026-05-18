# examples/mskcc

Reference invocations that USE MSKCC HPC lab paths (`/data1/greenbab/...`,
`/home/ahunos/miniforge3/...`, the lab `databases_config.yaml`). Useful as
templates for lab members; outsiders should look at
[`../portable/`](../portable/) instead.

| Script / dir | What it does |
|---|---|
| `single_sample.sh` | One HTML, one SNV site, lab paths |
| `cohort_samplesheet.sh` | Cohort build from a TSV samplesheet, lab paths |
| `prep_track_demo.sh` | Convert a plain-gzip GTF to bgzip+tabix |
| `methylation_ont/` | ONT 5mC + 5hmC track-config workflow |
