# igv-reports-skill

Cohort-aware driver + post-render verifiers on top of the IGV team's
[`igv-reports`](https://github.com/igvteam/igv-reports) (`create_report`).
Builds self-contained, offline HTML genomic-region reports — one HTML per
sample for a whole cohort, with embedded BAM/VCF data slices, default
annotation tracks, and a structural + content audit trail.

## What this is (and is not)

| Component | Source | Role |
|---|---|---|
| `create_report` CLI | upstream PyPI package `igv-reports` | does the actual HTML rendering |
| `scripts/build_igvreports.py` | **this repo** | wraps `create_report` with default-track resolution, cohort/samplesheet mode, SIF auto-detect |
| `scripts/verify_{report,cohort,anchors}.py` | **this repo** | post-render structural + content audits (not in upstream) |
| `scripts/prep_track.sh` | **this repo** | bgzip+tabix utility for annotation tracks |

If you only need raw `create_report` (no cohort mode, no verifiers, no
auto-tracks), use upstream directly. This repo is the right pick when you
want a single command to produce per-sample HTMLs across a whole cohort
with verified content.

## Install

```bash
# Upstream engine (provides `create_report` on PATH):
pip install igv-reports

# This driver (clone for the scripts + tests + SKILL.md):
git clone https://github.com/sahuno/igv-reports-skill.git
cd igv-reports-skill
```

Optional editable install (`pip install -e .`) puts the wrapper scripts on
PATH; running `python scripts/build_igvreports.py ...` works either way.

## Quickstart — single sample

```bash
python scripts/build_igvreports.py \
    --genome hg38 \
    --sites sites.hg38.bed \
    --bam tumor.bam normal.bam \
    --fasta /path/to/hg38.fa \
    --no-default-tracks \
    --extra-track /path/to/your_cpg_islands.bed.gz \
    --extra-track /path/to/gencode.v47.annotation.gff3.gz \
    --output report.hg38.html
```

`--no-default-tracks` skips the lab's auto-resolved annotation tracks; pass
your own with `--extra-track`. On MSKCC HPC (or anywhere you ship a
[`databases_config.yaml`](references/databases_config_paths.md)), omit
`--fasta` / `--no-default-tracks` and the driver resolves FASTA + CpG
islands + gencode + RepeatMasker automatically from `--genome`.

## Quickstart — cohort

```bash
python scripts/build_igvreports.py \
    --genome hg38 \
    --samplesheet cohort.tsv \
    --fasta /path/to/hg38.fa \
    --no-default-tracks \
    --output-dir reports/
```

Samplesheet (TSV with header):
```
sample	bam_tumor	bam_normal	vcf	sites_bed
p001	p001.tumor.bam	p001.normal.bam	p001.vcf.gz	p001.sites.hg38.bed
p002	p002.tumor.bam	p002.normal.bam	p002.vcf.gz	p002.sites.hg38.bed
```

Produces `reports/<sample>.hg38.html` for every row plus an `index.html`
linking them all. A structural verifier (`verify_cohort.py`) runs
automatically and writes `cohort_verify.tsv` + `cohort_verify.summary.md`
next to the reports. Add `--anchors-mode generate` (then `verify` on
re-runs) for read-count content verification.

## Environment overrides

| Var | Effect |
|---|---|
| `IGV_REPORTS_DB_CONFIG` | Path to your own databases YAML (same schema as MSKCC lab's; see [`references/databases_config_paths.md`](references/databases_config_paths.md)) |
| `IGV_REPORTS_SIF` | Path to your own `igv-reports` apptainer SIF |
| `SAMTOOLS_SIF_DEFAULT` | Path to your own `samtools` SIF (verifier only) |
| `IGV_REPORTS_BIND` | Colon-separated bind paths for singularity (default `/data1/greenbab`). Empty string disables binding. |

## Tests

```bash
bash tests/run_all.sh                  # unit + smoke + integration (integration skips if no test BAMs)
bash tests/run_all.sh --unit-only      # 63 hermetic parser tests, runs anywhere
bash tests/run_all.sh --no-integration # unit + smoke (smoke needs samtools)
```

CI runs unit + smoke on Python 3.10 / 3.11 / 3.12 on every push to `main`.

## Use as a Claude Code skill

This repo includes a [`SKILL.md`](SKILL.md) for Claude Code users:

```bash
# Symlink the repo root into Claude Code's skills dir:
ln -s ~/code/igv-reports-skill ~/.claude/skills/igv-reports
```

After symlinking, Claude Code discovers `igv-reports` as a callable skill
and reads `SKILL.md` for invocation guidance.

## Repo layout

```
igv-reports-skill/
├── README.md
├── LICENSE                  # MIT
├── SKILL.md                 # Claude Code skill manifest
├── pyproject.toml           # optional `pip install -e .`
├── scripts/                 # build_igvreports, verify_{report,cohort,anchors}, prep_track
├── tests/                   # unit (hermetic), smoke (samtools), integration (lab BAMs, optional)
├── examples/
│   ├── portable/            # off-MSKCC reference invocations
│   └── mskcc/               # MSKCC HPC reference invocations (lab paths)
├── references/              # databases_config schema + per-genome track availability
└── evals/                   # benchmarks / future eval harnesses
```

## License

MIT — see [LICENSE](LICENSE).

## Citation

If you use this for a publication, please cite both the upstream
`igv-reports` package and this driver:

- Upstream: [igvteam/igv-reports](https://github.com/igvteam/igv-reports)
- This driver: `github.com/sahuno/igv-reports-skill`
