---
name: methylation_ont
genome: hg38
assay: ONT 5mC + 5hmC (CpG, basemod2)
inputs:
  - per-sample BAM with MM/ML tags (modkit / dorado basecalled)
  - per-sample 5mC and 5hmC bedGraph (already sliced or whole-genome)
  - sites BED defining the windows to show
outputs:
  - one self-contained HTML with N regions x M tracks
flanking: 0   # sites BED already encodes the window
---

# Recipe — ONT 5mC/5hmC promoter methylation viewer

This directory is the canonical worked example for the igv-reports skill's
**ONT methylation** path. Every file here either parameterizes the build
(`*.template.json`, `*.example.*`, `tracks_spec.example.yaml`) or IS the
sliced data the example uses (`annotation_slices/`, `bedgraph_slices/`,
`gencode.v47.promoter_regions.gtf.gz`).

When the agent is asked to build a methylation viewer, the fastest path is:

1. Copy this directory into the user's results dir.
2. Replace `tracks_spec.example.yaml` with the user's sample list (same shape).
3. Run `build.sh`.

## The four things the agent must know

These are the points the generic skill doesn't cover and that the
plain `--tracks` positional API of `build_igvreports.py` cannot express.

### 1. Use `--track-config tracks.json`, not `--tracks`

Methylation viewers need **named** tracks with per-track styling:
- BAM tracks need `colorBy: "basemod2"` for per-read 5mC/5hmC coloring.
- bedGraph tracks need `min: 0, max: 100` (methylation percent, fixed
  across samples — see `rules/igv.md`: IGV's per-track autoscale hides
  real cross-sample differences).
- Annotation tracks need explicit `color` so CGI vs RepeatMasker vs
  EPDnew don't all render identically.

The skill driver (`build_igvreports.py --track-config <json>`) just
passes the JSON through to `create_report`. No track merging, no
default-track injection — the JSON is the source of truth.

### 2. `--flanking 0` when sites encode the desired window

The default flanking of 300 bp is for SV/integration breakpoints where
the BED row is the breakpoint itself. For whole-promoter or whole-gene
methylation views, the BED row already spans the window of interest
(e.g., 13 kb around DNMT3A's promoter). Adding 300 bp adds nothing
useful and shifts the viewer's initial frame.

### 3. bedGraph, not bigwig

Per `rules/igv.md`: igv-reports' Python slicer cannot read `.bw`. The
only working path is to pre-slice each bigwig to bedGraph over the
report regions (`bigWigToBedGraph -chrom=$c -start=$s -end=$e <bw> <bg>`),
write one bedGraph per sample × modification, and reference those from
`tracks.json`. UCSC `bigWigToBedGraph` truncates `/dev/stdout` between
calls — write per-region temp files and `cat >>` them.

### 4. `#chrom\tstart\tend\tname` comment-header sites BED IS accepted

The skill's `references/best_practices.md` historically warned "BED
must be headerless" because of a `ValueError: invalid literal for
int()` from create_report's positional parser. That message was caused
by a non-`#`-prefixed header row. A line starting with `#` is treated
as a comment by create_report and is fine — and matches the lab's
"BED-like outputs must have a `#`-prefixed header" rule in CLAUDE.md.

## Files in this example

| File | Purpose |
|---|---|
| `tracks.template.json` | Hand-edit template. Replace `{{...}}` placeholders inline. |
| `tracks_spec.example.yaml` | YAML spec consumed by `scripts/generate_tracks_json.py`. The recommended way — agent edits the YAML, generator emits the JSON with right colors/min-max/colorBy. |
| `sites.hg38.example.bed` | Two promoter windows: DNMT3A_2 and EZH2. `#`-prefixed header. |
| `build.sh` | Runnable: generates tracks.json from the YAML spec then calls the skill driver with `--track-config`. |
| `annotation_slices/*.bed` | EPDnew / CpGIslands / RepeatMasker pre-sliced to the two example windows. |
| `bedgraph_slices/*.bedgraph` | Real COLO829 5mC + 5hmC tracks for the same windows (4 samples × 2 mods = 8 files). |
| `gencode.v47.promoter_regions.gtf.gz[.tbi]` | bgzip+tabix gencode slice for the two windows. |

## How to adapt for a new run

1. **Edit `tracks_spec.example.yaml`** — change `samples:` to your sample
   list. Pick a `group:` per sample (`normal` / `tumor` / `replicate1` /
   anything that lives under `group_colors:`). Add a `group_colors:`
   entry if you introduce a new group name.

2. **Slice the bedGraphs and annotations** to your regions. Re-using
   the patterns in `rules/igv.md` (UCSC `bigWigToBedGraph` per-region
   temp files, then `cat >>`). Drop them in `bedgraph_slices/` and
   `annotation_slices/` (or change the paths in the YAML).

3. **Edit `sites.<genome>.example.bed`** to your windows. Keep the
   `#chrom\tstart\tend\tname` header.

4. **Run `build.sh`.** It activates the snakemake conda env (where
   `create_report` lives), generates the JSON, then calls
   `build_igvreports.py --track-config <json>`.

## Direct equivalent (audit trail)

`build.sh` calls the driver, which in turn calls `create_report`. The
fully-expanded direct invocation for the example sites is:

```bash
singularity exec --cleanenv --bind /data1/greenbab \
  /data1/greenbab/users/ahunos/apps/containers/igv-reports_1.16.0.sif \
  create_report \
    sites.hg38.example.bed \
    --fasta /data1/greenbab/database/hg38/v0/Homo_sapiens_assembly38.fasta \
    --type mutation \
    --track-config tracks.json \
    --info-columns name \
    --flanking 0 \
    --standalone \
    --title "COLO829 promoter methylation: DNMT3A_2 + EZH2" \
    --output methylation_report.hg38.html
```

The `--cleanenv` flag is mandatory: host RHEL 8 exports `SSL_CERT_FILE`
pointing at a path the SIF doesn't have, and `create_report --standalone`
makes an HTTPS GET (ideogram CDN) that crashes with
`[SSL: CERTIFICATE_VERIFY_FAILED]` if host env leaks in. See
`rules/apptainer_env_leak.md`. The skill driver always passes `--cleanenv`
so users invoking via `build_igvreports.py --apptainer` don't need to
remember it.

The SIF is a dedicated igv-reports 1.16.0 container (~83 MB) pulled
from the Galaxy depot. It's much lighter than the older
`onttools_v3.10.sif` (which bundles dorado + samtools + bedtools + ...
alongside create_report). The scratch reference run used onttools because
that's what was already on disk at the time — both produce the same
HTML, but the dedicated SIF is portable and matches the lab's standard
practice of one-tool-per-SIF (`rules/severus.md` documents pulling
from Galaxy depot to avoid Quay rate limits and Seqera auth).

Both apptainer and conda paths produce the same HTML; the apptainer
path avoids the conda cold-start tax on fresh compute nodes
(`rules/apptainer_vs_conda.md`). The skill driver uses the conda path
by default — pass `--apptainer` to switch.

## Related rules and skills

- `rules/igv.md` — the y-axis-autoscale problem (why min:0,max:100),
  the bigwig-can't-be-sliced problem (why bedGraph), and the UCSC
  `/dev/stdout`-truncation gotcha (why per-region temp files).
- `rules/apptainer_vs_conda.md` — when to prefer the apptainer SIF.
- `CLAUDE.md` §3A — ONT methylation pipeline upstream of this skill
  (pod5 → dorado → modkit pileup → bedGraph + bigwig).
- `references/methylation_ont.md` (sibling reference in this skill) —
  the cheat-sheet version of this recipe.
