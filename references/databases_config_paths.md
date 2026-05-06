# Per-genome track availability

Quick lookup for which YAML keys in
`/data1/greenbab/users/ahunos/apps/llm_configs/claude/profiles/databases/databases_config.yaml`
hold which tracks for each supported genome. The skill driver reads this YAML
and skips any track that's missing for the chosen genome (with a warning in
the run log).

## Keys

The YAML structure is:

```yaml
reference_genomes:
  local:
    <genome_id>:
      fasta: <path>
      gtf:   <path>     # gencode .gtf.gz (NOT .gff3.gz in current YAML)
      sizes: <path>     # chrom.sizes
      CpGIslands:    <path>  # .bed
      repMaskerBed:  <path>  # .bed.gz (bgzip + tabix)
```

## Genome → track matrix

| Genome key        | fasta | gtf  | CpGIslands | repMaskerBed | Notes |
|-------------------|:-----:|:----:|:----------:|:------------:|-------|
| `hg38`            |  ✅   |  ✅  |     ✅     |      ✅      | Full set. Skill prefers the bgzip+tabix `gencode.v47.annotation.gff3.gz` over the YAML `gtf` (richer transcript model). |
| `mm10`            |  ✅   |  ✅  |     ✅     |      ✅      | Full set. |
| `mm39`            |  ✅   |  ✅  |     ✅     |      ❌      | No rmsk in YAML — skill auto-skips with warning. |
| `t2t_CHM13v2_plusY` | ✅ |  ✅  |     ✅     |      ❌      | No rmsk in YAML. CGI is at `/data1/greenbab/database/CpGIslands/t2t/chm13v2.0_CGI.bed`. |
| `GRCh37`          |  ✅   |  ✅  |     ✅     |      ❌      | Plus tandem-repeats and gaps tracks not loaded by default. |

## hg38 special case — full gencode annotation

The YAML `gtf` key for hg38 points to `gencode.v47.annotation.gtf.gz`
(plain gzip, gene-level via genes-only sibling file). The skill driver
prefers a sibling file `gencode.v47.annotation.gff3.gz` in the same
directory, which IS bgzip + tabix-indexed and has the full transcript /
exon / CDS / UTR detail. If that sibling is missing, fall back to the
YAML `gtf`.

If the user explicitly wants the YAML's `.gtf.gz` regardless, pass
`--gencode-from-yaml`. But the default is the richer GFF3.

## Genome ID aliases

igv-reports `--genome` flag uses igv.js bundled IDs (`hg38`, `mm10`, etc.)
which auto-fetch from CDN. We don't use that flag — the skill always uses
`--fasta <path>` from the YAML so reports work fully offline. The
`--genome` arg in the driver is just a key into the YAML, not the
igv.js bundled-genome name.

| Skill `--genome` value | YAML key            | Common alias |
|-----------------------|---------------------|--------------|
| `hg38`                | `hg38`              | GRCh38       |
| `mm10`                | `mm10`              | GRCm38       |
| `mm39`                | `mm39`              | GRCm39       |
| `t2t` / `chm13`       | `t2t_CHM13v2_plusY` | T2T-CHM13v2  |
| `grch37` / `hg19`     | `GRCh37`            | hg19         |

The driver normalizes the input alias to the canonical YAML key.

## When a track is missing

If a default track isn't available for the chosen genome, the driver:
1. Logs a `WARNING: <track_type> not configured for <genome>` line.
2. Skips that track in the `--tracks` list.
3. Continues — the report builds without that track.

If the user wants to add a missing track:
1. Build / locate the BED or GFF3.
2. Run `prep_track.sh` if it needs bgzip+tabix conversion.
3. Add the path to `databases_config.yaml` under the appropriate key, or
   pass it via `--extra-track <path>` for a one-off run.
