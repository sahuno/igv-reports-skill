# examples/portable

Reference invocations that DO NOT depend on MSKCC lab paths or the lab
`databases_config.yaml`. Each script accepts environment-variable
overrides for input paths, with `${HOME}/data/...` defaults you can edit
in-place or override at call time:

```bash
FASTA=/path/to/hg38.fa TUMOR_BAM=/path/to/tumor.bam \
    bash examples/portable/single_sample.sh
```

| Script | What it does |
|---|---|
| `single_sample.sh` | Builds one HTML for a tumor/normal pair at two SNV sites |
| `cohort_samplesheet.sh` | Builds per-sample HTMLs + index.html from a 2-row samplesheet |

For MSKCC lab equivalents (using `/data1/greenbab/` paths and the lab
databases YAML), see [`../mskcc/`](../mskcc/).
