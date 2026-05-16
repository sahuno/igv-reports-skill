# cohort_verify_demo — regression test for `verify_cohort.py`

End-to-end check that the cohort verifier catches the four failure modes
it's designed to catch. Self-asserting — exit nonzero on any mismatch.

## What it does

1. Generates a 3-sample samplesheet (TSV) pointing at three real COLO829 ONT BAMs.
2. Calls `build_igvreports.py --samplesheet ... --no-verify` to produce
   `reports/sample_{1,2,3}.hg38.html` + `index.html`.
3. Runs `verify_cohort.py` against the clean cohort (expects all PASS).
4. Runs four corruption scenarios, each asserting the expected check FAILs:

   | Scenario | Corruption | Expected FAILs |
   |---|---|---|
   | A | Delete `sample_3.hg38.html` | `*/cohort_html_coverage`, `sample_3/html_exists` |
   | B | Replace `sample_1.hg38.html` with sample_2's content | `sample_1/sample_tracks_match`, `sample_1/no_cross_sample_contamination`, `sample_1/sample_id_embedded` |
   | C | Drop one `<li>` from `index.html` | `*/index_consistency` |
   | D | Truncate `sample_2.hg38.html` to 1 KB | `sample_2/html_min_size`, `sample_2/region_count` |

5. Cleans up generated `reports/`, samplesheet, sites BED, and logs/ on exit.

## Run

```bash
bash examples/cohort_verify_demo/scenarios.sh
```

Runtime: ~60-90 s on a warm node (3-sample cohort build at 1-bp point-variant
sites + 4 reverify cycles). Per-sample HTML ends up ~3-5 MB. Cold-cache
network reads of the underlying ONT BAMs can extend this to 2-3 min on
first invocation.

Disk: ~15 MB temporary under `reports/`, auto-cleaned via `trap`.

The sites BED uses 1-bp point-variant style coordinates (not 13 kb promoter
windows like the methylation example) so BAM slicing stays fast — we're
testing the verifier, not the renderer. Adapt for other workflows if you
want to exercise wider windows.

## Adapt for other clusters

The BAM paths at the top of `scenarios.sh` (`BAM_S1`/`BAM_S2`/`BAM_S3`) are
hardcoded to COLO829 ONT runs on MSKCC's `/data1/greenbab`. Off-cluster, edit
those paths to point at any three indexed BAMs you have access to — the
verifier doesn't care which BAMs, only that the three rows in the samplesheet
declare *different* BAMs (so scenario B's contamination check has signal).

## Why this lives in `examples/`

This is a test of the verifier, but it depends on real lab BAMs and on
`create_report` actually running, so it's not a unit test — it's an
integration smoke. The skill doesn't yet have a `tests/` dir; when one
exists, move this here.
