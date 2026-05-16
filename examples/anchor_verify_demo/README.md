# anchor_verify_demo — regression test for `verify_anchors.py`

End-to-end check that the anchor-based content verifier catches the four
failure modes it's designed to catch. Self-asserting — exits nonzero on any
mismatch.

## What it does

1. Generates a 2-sample samplesheet (TSV) pointing at two real COLO829 ONT BAMs.
2. Calls `verify_anchors.py generate` to freeze `samtools view -c` counts
   into `anchors.hg38.tsv` (the regression fixture).
3. Calls `build_igvreports.py --samplesheet ... --no-verify` to produce
   `reports/sample_{1,2}.hg38.html` + `index.html`.
4. Runs `verify_anchors.py verify-cohort` against the clean cohort (all PASS).
5. Runs four corruption scenarios, each asserting the expected outcome:

   | Scenario | Corruption | Expected |
   |---|---|---|
   | A | Mutate anchors `expected` to 9999 (real ~56) | `sample_1/chr2/FAIL` (diff_ratio) |
   | B | Set anchor `min=1000` (real ~56) | `sample_1/chr2/FAIL` (min bound) |
   | C | Mangle a session's base64 payload (`H4sI` → `XXXX`) | `sample_1/*/FAIL` (decode), sample_2 PASS |
   | D | Drop an anchor row | row absent from output, others PASS |

6. Cleans up generated `reports/`, samplesheet, sites BED, anchors TSVs, and
   `logs/` on exit (set `KEEP_REPORTS=1` to leave them).

## Run

```bash
bash examples/anchor_verify_demo/scenarios.sh
```

Runtime: **~6-8 min on a warm node** (the cohort build dominates — verify
itself is ~3 s per pass over 2 samples × 2 anchors). Set `REBUILD=1` to force
rebuild of the HTMLs; otherwise existing HTMLs in `reports/` are reused so
iteration on the verifier is seconds.

Disk: ~10 MB temp under `reports/`, auto-cleaned via `trap`.

## Why these scenarios

The four scenarios cover every status the verifier emits:

- **PASS** (scenario 0): observed within tolerance of expected, or within
  `min`/`max` bounds.
- **FAIL — tolerance** (A): observed read count differs from expected beyond
  the per-row tolerance (default 5%). Catches the silent sample-swap case
  where the wrong source BAM was wired into the build pipeline — same track
  name, different read counts.
- **FAIL — bound** (B): `min`/`max` columns let you assert "this integration
  site should have ≥20 reads supporting it" — a stronger claim than
  tolerance, useful for known-positive sites.
- **FAIL — broken decode** (C): the HTML's session entry can't be gunzipped
  or its inner BAM data URL can't be base64-decoded. Catches arbitrary HTML
  tampering or `create_report` version drift that breaks the embedding format.
- **SKIP** (D): an anchor row references a `(sample, region)` pair that the
  HTML doesn't render. Dropped silently because anchor TSVs are intentionally
  re-usable across runs — a region that exists in one cohort's anchors but
  not in another cohort's HTMLs is benign, not a build failure.

## Adapt for other clusters

The BAM paths at the top of `scenarios.sh` (`BAM_S1`/`BAM_S2`) are hardcoded
to COLO829 ONT runs on MSKCC's `/data1/greenbab`. Off-cluster, edit those
paths to point at any two indexed BAMs you have access to. The verifier
doesn't care which BAMs, only that they are different so scenarios A-C have
the contrast they need.

## Why this lives in `examples/`

Same reasoning as the sibling [cohort_verify_demo](../cohort_verify_demo/):
this is a test of the verifier, but it depends on real lab BAMs and on
`create_report` actually running, so it's an integration smoke rather than
a unit test. Move under `tests/` when the skill gets a real test dir.
