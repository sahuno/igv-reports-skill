---
name: igv-reports
description: |
  Build self-contained, offline HTML genomic-region reports with igv-reports (create_report). Each HTML bundles igv.js viewers per region with embedded BAM/VCF data slices and default tracks (CpG islands, gencode, RepeatMasker); a reviewer clicks the variant table to inspect read-level evidence with no internet, no server, no IGV install.

  USE this skill whenever the user wants an HTML, clickable, or browseable viewer of genomic data — phrases like "HTML IGV report", "offline IGV", "self-contained HTML", "clickable viewer", "create_report", "igv-reports", "email this viewer", or any browseable HTML of reads at variants, fusion breakpoints, SV junctions, viral integrations, ChIP peaks, or ROIs. Trigger even when the user doesn't say "igv-reports" — giveaway is HTML/clickable/offline plus genomic regions. Also fire on /igv-reports.

  DO NOT use for static PNG/PDF/SVG IGV screenshots — use the igv-screenshots skill.

  Supports hg38, mm10, mm39, T2T. Defaults: --flanking 300, --standalone, genome-tagged output.
---

# igv-reports

This skill builds **self-contained HTML genomic-region reports** with
[igv-reports](https://github.com/igvteam/igv-reports) (`create_report`).
Each report is a single browseable HTML containing the igv.js viewer plus
embedded data slices for every region. No server, no internet, no IGV
install needed at view time.

The skill has three entry points:
- **build** — one-shot: sites BED + BAM(s) ± VCF → HTML.
- **cohort** — multi-sample driver from a samplesheet → per-sample HTMLs + index.
- **prep-track** — utility: convert plain-gzip GFF/GTF/BED.gz into a
  bgzip + tabix-indexed track that igv-reports can load.

## When to use which entry point

| User request | Entry point |
|---|---|
| "Make an HTML for these 5 SV breakpoints in tumor.bam" | **build** |
| "Give me one HTML per patient for the cohort integration calls" | **cohort** |
| "create_report fails with 'not BGZF' on this gencode" | **prep-track** |

## Defaults (locked in)

- Tracks always loaded, top-to-bottom in the viewer:
  1. CpG islands (BED, plain or bgzipped)
  2. Gencode full annotation (GFF3.gz, **transcripts + exons + CDS + UTRs**, NOT a gene-level-only file)
  3. RepeatMasker (BED.gz, bgzipped + tabix-indexed)
  Plus the user's BAM(s), VCF, and any extra tracks they pass.
- `--flanking 300` bp on either side of each site (good for SV breakpoints
  and point variants alike). Override per call if needed.
- `--standalone` so the HTML is offline-viewable.
- Output filename includes the genome tag — e.g. `cohort.hg38.html` —
  to pass `enforce-genome-tag.sh`.
- Reference FASTA is resolved from `databases_config.yaml`:
  `/data1/greenbab/users/ahunos/apps/llm_configs/claude/profiles/databases/databases_config.yaml`.
  Supported genome IDs: `hg38`, `mm10`, `mm39`, `t2t_CHM13v2_plusY`, `GRCh37`.
- Per-genome default track availability is recorded in
  `references/databases_config_paths.md` — read it before assembling tracks
  so the skill doesn't try to load a track that doesn't exist for the
  selected genome (e.g., mm39 has no rmsk in our database).

## Sites BED format (critical)

igv-reports' BED parser reads fields **by position** and trips on a header
row (`ValueError: invalid literal for int() with base 10: 'start'`). Always
emit a **plain headerless 4-column BED**:

```
chr    start    end    name
chr2   25227855 25342590 DNMT3A_full_gene
```

Tab-separated. The `name` becomes the row label in the report's variant
table — make it specific enough to identify the site after deduping.

The project's `enforce-genome-tag.sh` hook requires a genome tag in the BED
filename: use `sites.hg38.bed`, not `sites.bed`.

## Pitfalls (the skill should encode and/or detect these)

| Symptom | Root cause | Fix |
|---|---|---|
| `ValueError: invalid literal for int()` on first row | Header row in sites BED | Strip header — plain BED |
| `UnicodeDecodeError: byte 0x8b` reading a track | igv-reports reading bgzip as text | Filename must end `.gff3.gz` / `.bed.gz` AND be true bgzip (check with `file <name>` for "extra field") |
| `tabix: not BGZF` | Track was plain-gzipped, not bgzipped | Run **prep-track** entry point |
| `tabix: out of order` while indexing | GFF/GTF/BED records not pos-sorted within chr | **prep-track** does `sort -k1,1 -k4,4n` before bgzip |
| Annotation track empty in viewer | Tabix returns no rows in displayed window — often correct biology (e.g., CGI-distal site). Confirm with `tabix file region` |
| Genome ID lookup fails with `--genome hg38` | igv.js bundled IDs require internet at view + render time. Use `--fasta /path/to/local.fa` instead (always works offline) |

Full pitfalls + create_report flag reference in `references/best_practices.md`.

## How to run — quick recipe

Activate the snakemake conda env first; `create_report` lives there:

```bash
source /home/ahunos/miniforge3/etc/profile.d/conda.sh
conda activate snakemake
```

Then call the bundled driver script:

```bash
python /data1/greenbab/users/ahunos/apps/llm_configs/claude/skills/igv-reports/scripts/build_igvreports.py \
    --sites results/run/inputs/sites.hg38.bed \
    --bam tumor.bam normal.bam \
    --vcf calls.vcf \
    --genome hg38 \
    --output results/run/reports/cohort.hg38.html
```

The driver:
- Resolves the genome's CpG / gencode / rmsk paths from `databases_config.yaml`
  (skipping any that don't exist for the chosen genome).
- Validates the sites BED is headerless and that all rows have `start < end`.
- Calls `create_report` with `--flanking 300 --standalone`.
- Writes a logs/ entry capturing the full command, the flanking value, the
  per-region embedded data sizes, and the resolved track list — required
  by the project's analysis-script audit-trail expectations.

For multi-sample cohorts, use `--samplesheet samplesheet.tsv` instead of
`--bam/--vcf`. Samplesheet format: `sample, bam_tumor, bam_normal, vcf, sites_bed`.
The driver emits one HTML per sample plus a top-level `index.html` that lists
all samples with links. Layout matches the ATLL viral-integration reference
implementation:

```
results/<run>/
├── inputs/<sample>/sites.<genome>.bed
├── reports/<sample>.<genome>.html
├── reports/index.html
└── logs/run_<timestamp>.log
```

## prep-track — fixing a non-bgzip track

If a GFF3/GTF/BED.gz is plain-gzip rather than bgzip, igv-reports fails
silently or with an obscure error. Convert in place with backup:

```bash
bash /data1/greenbab/users/ahunos/apps/llm_configs/claude/skills/igv-reports/scripts/prep_track.sh \
    /path/to/track.gff3.gz
```

The script:
1. Backs up the original to `<name>.bak.original_gzip`.
2. `gunzip -c`s the file.
3. Sorts by `chr` then numeric `pos` (`sort -k1,1 -k4,4n`).
   (Gencode delivers records interleaved by feature type at the same locus —
   tabix requires pos-sorted.)
4. `bgzip`s in place.
5. `tabix -p <gff|gtf|bed>`s.
6. Verifies a sample tabix query returns rows.

Run from the snakemake conda env (bgzip/tabix from htslib).

## When generating an answer.md / run.sh for the user

The driver script (`build_igvreports.py`) deliberately abstracts the
underlying `create_report` flags — it sets `--standalone`, `--fasta`, the
`--flanking 300` default, and the YAML-resolved annotation tracks
internally so the user doesn't have to remember them. That abstraction is
good for ergonomics but bad for auditability: a reviewer reading the
`answer.md` later can't see what flags are actually being invoked without
opening the driver source.

To keep both: when you produce a runnable command for the user, **also
include a code block titled "Equivalent direct create_report invocation"
that shows the fully-expanded command** with all flags and resolved track
paths inline. The user should see the wrapper command they're going to
run AND the underlying command it expands to. Example:

````
## Run

```bash
python build_igvreports.py --genome mm10 --sites peaks.mm10.bed \\
    --bam ./data/ip.bam ./data/input.bam \\
    --output reports/peaks_qc.mm10.html
```

### Equivalent direct create_report invocation

```bash
create_report peaks.mm10.bed \\
    --fasta /data1/greenbab/database/mm10/mm10.fa \\
    --flanking 300 --standalone \\
    --tracks ./data/ip.bam ./data/input.bam \\
        /data1/greenbab/database/mm10/mm10_CpGIslands.bed \\
        /data1/greenbab/database/mm10/annotations/gencode.vM25.annotation.gtf.gz \\
        /data1/greenbab/database/RepeatMaskerDB/.../rmsk_all_repeats_mm10.bed.gz \\
    --title "ChIP-seq peak QC (mm10) — IP vs Input" \\
    --output reports/peaks_qc.mm10.html
```
````

This costs you ~10 lines and gives the reviewer a full audit trail. For
cohort runs, show the expanded form for ONE representative sample only —
the others differ only in BAM/VCF paths.

## Output and workflow logging

Every run logs to `logs/run_<YYYYMMDD_HHMMSS>.log` next to the reports dir.
The log captures:
- Resolved track paths (per genome, after databases_config.yaml lookup).
- The exact `create_report` command.
- The flanking value used (default **300 bp** — this is the value that's
  baked into all the embedded data slices, so audit trails depend on it).
- Per-region embedded data sizes (extracted post-render so the user can
  see which regions inflated the HTML).
- Total HTML size.

This satisfies CLAUDE.md §"Logging and Audit Trail" — every run is
reproducible from the log alone.

## Track choice nuances

For gencode on hg38, the default points at
`gencode.v47.annotation.gff3.gz` (full annotation, bgzip + tabix). This
gives transcript models with exons / CDS / UTRs. The gene-level-only
companion (`gencode.v47.genes.annotation.sorted.gff3.gz`) renders only
solid gene boxes and is fine for high-zoom views, but the full annotation
is the right default for read-level inspection at integration / fusion /
SV junctions.

For mouse genomes, `databases_config.yaml` ships `.gtf.gz` paths instead.
GTFs work in igv-reports if bgzip + tabix-indexed; **prep-track** converts
plain-gzip GTFs the same way it does GFF3s.

For T2T-CHM13, only the FASTA + GTF + CGI are indexed in our DB; rmsk is
absent and is auto-skipped by the driver. The variant table will load
without rmsk; flag this in the run log.

## Common-case examples

The `examples/` directory has runnable templates:

- `single_sample.sh` — one BAM + one VCF + a sites BED → one HTML.
- `cohort_samplesheet.sh` — TSV-driven multi-sample run.
- `prep_track_demo.sh` — convert a plain-gzip gencode to bgzip+tabix.

These are reference implementations; copy and edit them for new runs
rather than starting from scratch.

## See also

- `references/best_practices.md` — full create_report flag reference,
  format gotchas, performance notes. Read this if a run fails in a way
  not listed in the Pitfalls table above.
- `references/databases_config_paths.md` — per-genome track availability
  matrix and exact YAML keys. Read this when adding a new genome or
  diagnosing a missing-track warning.
- `scripts/build_igvreports.py` — the driver. Reads `--samplesheet` or
  `--bam/--vcf` direct-args, resolves tracks, validates the sites BED,
  writes the HTMLs and the run log.
- `scripts/prep_track.sh` — gunzip → sort → bgzip → tabix utility.
- `igv-screenshots` skill — the **static PNG/PDF/SVG** counterpart based
  on igver. Use it instead of this one when the deliverable is a
  publication-quality figure rather than a clickable viewer.
- Reference implementation:
  `/data1/greenbab/projects/ont/Project_17424/results/20260503_hg38plusHTLV1EBV_cohort_integration_igvreports/`
  — 6-patient ATLL cohort viral-integration HTMLs + DNMT3A sanity check;
  this skill was extracted from that work.
