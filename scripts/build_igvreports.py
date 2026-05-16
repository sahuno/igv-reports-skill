#!/usr/bin/env python3
"""build_igvreports.py — generic driver for the igv-reports skill.

Author: Samuel Ahuno
Purpose:
  Build self-contained HTML genomic-region reports with create_report
  (igv-reports). Two run modes:

    1. Single — direct CLI: --sites BED + --bam BAM(s) [+--vcf VCF]
       → one HTML at --output.

    2. Cohort — TSV samplesheet: one HTML per row + an index.html.
       Samplesheet columns (tab-separated, with header):
         sample   bam_tumor   bam_normal   vcf   sites_bed
       Optional fifth column: extra_tracks (comma-separated paths).

  Either way, the driver:
    - Resolves CpG islands, gencode, and RepeatMasker paths from
      databases_config.yaml for the chosen genome (skipping any not
      configured for that genome, with a warning).
    - Validates that the sites BED is headerless and well-formed.
    - Calls create_report with --flanking 300 --standalone by default.
    - Writes a logs/ entry capturing the resolved track list, the full
      command, the flanking value, and per-region embedded data sizes.

Usage:
  python build_igvreports.py --sites SITES.hg38.bed \\
      --bam tumor.bam normal.bam --vcf calls.vcf \\
      --genome hg38 --output report.hg38.html

  python build_igvreports.py --samplesheet sheet.tsv \\
      --genome hg38 --output-dir results/cohort/

Skill location:
  /data1/greenbab/users/ahunos/apps/llm_configs/claude/skills/igv-reports/
"""

from __future__ import annotations

import argparse
import logging
import os
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path

try:
    import yaml  # PyYAML
except ImportError:
    print("ERROR: PyYAML not available. Activate the snakemake conda env first:", file=sys.stderr)
    print("  source /home/ahunos/miniforge3/etc/profile.d/conda.sh && conda activate snakemake", file=sys.stderr)
    sys.exit(2)

DEFAULT_DBCONFIG = Path(
    "/data1/greenbab/users/ahunos/apps/llm_configs/claude/profiles/databases/databases_config.yaml"
)
DEFAULT_FLANKING = 300
# Dedicated igv-reports 1.16.0 SIF (~83 MB), pulled from Galaxy depot.
# Cleanest path for HPC: minimal payload vs. the heavier onttools_v3.10.sif
# which incidentally bundles create_report alongside dorado/samtools/etc.
IGVREPORTS_SIF = Path("/data1/greenbab/users/ahunos/apps/containers/igv-reports_1.16.0.sif")

GENOME_ALIASES = {
    "hg38": "hg38",
    "GRCh38": "hg38",
    "mm10": "mm10",
    "GRCm38": "mm10",
    "mm39": "mm39",
    "GRCm39": "mm39",
    "t2t": "t2t_CHM13v2_plusY",
    "chm13": "t2t_CHM13v2_plusY",
    "T2T": "t2t_CHM13v2_plusY",
    "T2T-CHM13": "t2t_CHM13v2_plusY",
    "t2t_CHM13v2_plusY": "t2t_CHM13v2_plusY",
    "GRCh37": "GRCh37",
    "hg19": "GRCh37",
}


def setup_logger(log_path: Path) -> logging.Logger:
    """Dual-handler logger: file + stderr, with timestamp prefix."""
    log_path.parent.mkdir(parents=True, exist_ok=True)
    fmt = logging.Formatter(
        "[%(asctime)s] %(levelname)s: %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
    )
    log = logging.getLogger("igv_reports")
    log.setLevel(logging.INFO)
    log.handlers.clear()
    fh = logging.FileHandler(log_path)
    fh.setFormatter(fmt)
    log.addHandler(fh)
    sh = logging.StreamHandler(sys.stderr)
    sh.setFormatter(fmt)
    log.addHandler(sh)
    return log


def resolve_genome(genome: str) -> str:
    canon = GENOME_ALIASES.get(genome)
    if not canon:
        raise SystemExit(
            f"ERROR: unknown genome '{genome}'. Supported: {sorted(set(GENOME_ALIASES.values()))}"
        )
    return canon


def load_db_config(path: Path) -> dict:
    if not path.exists():
        raise SystemExit(f"ERROR: databases_config.yaml not found at {path}")
    with path.open() as fh:
        cfg = yaml.safe_load(fh)
    return cfg


def resolve_default_tracks(cfg: dict, genome: str, log: logging.Logger) -> list[str]:
    """Return ordered list of default tracks present on disk for this genome.

    Order matters — first entry renders at the bottom of the IGV.js view by
    default? Actually igv-reports renders --tracks in the order passed,
    top-to-bottom. We put annotation tracks LAST so they sit below the
    BAM/VCF data the user is actually inspecting.
    """
    g = cfg.get("reference_genomes", {}).get("local", {}).get(genome, {})
    if not g:
        raise SystemExit(f"ERROR: no entry for genome '{genome}' in databases_config.yaml")

    tracks: list[str] = []

    # CpG islands.
    cgi = g.get("CpGIslands")
    if cgi and Path(cgi).exists():
        tracks.append(cgi)
    else:
        log.warning(f"CpG islands track missing for {genome} (key=CpGIslands, value={cgi})")

    # Gencode. For hg38 prefer the bgzip+tabix .gff3.gz sibling if present.
    gtf = g.get("gtf")
    gencode_track: str | None = None
    if genome == "hg38" and gtf:
        sibling = (Path(gtf).parent / "gencode.v47.annotation.gff3.gz")
        if sibling.exists() and (sibling.parent / (sibling.name + ".tbi")).exists():
            gencode_track = str(sibling)
            log.info(f"  hg38: using full gencode annotation: {sibling}")
    if gencode_track is None and gtf and Path(gtf).exists():
        gencode_track = gtf
    if gencode_track:
        tracks.append(gencode_track)
    else:
        log.warning(f"Gencode track missing for {genome}")

    # RepeatMasker.
    rmsk = g.get("repMaskerBed")
    if rmsk and Path(rmsk).exists():
        tracks.append(rmsk)
    else:
        log.warning(f"RepeatMasker track not configured for {genome}")

    return tracks


def fasta_for(cfg: dict, genome: str) -> str:
    fasta = cfg["reference_genomes"]["local"][genome].get("fasta")
    if not fasta or not Path(fasta).exists():
        raise SystemExit(f"ERROR: FASTA missing for {genome}: {fasta}")
    if not Path(fasta + ".fai").exists():
        raise SystemExit(
            f"ERROR: FASTA index missing for {fasta} — run `samtools faidx {fasta}`"
        )
    return fasta


def validate_sites_bed(bed: Path) -> None:
    """Sanity-check the sites BED before invoking create_report.

    create_report's BED parser is positional. It skips lines starting with
    `#` or `track ` (so the lab's `#chrom\\tstart\\tend\\tname` header is
    fine), but a non-comment header row like `chrom\\tstart\\tend` crashes
    with `ValueError: invalid literal for int()`. We mirror create_report's
    line-skipping logic and emit an informative error if any data row has
    non-numeric start/end."""
    if not bed.exists():
        raise SystemExit(f"ERROR: sites BED not found: {bed}")
    with bed.open() as fh:
        for i, line in enumerate(fh, start=1):
            line = line.rstrip("\n")
            if not line or line.startswith("#") or line.startswith("track "):
                continue
            cols = line.split("\t")
            if len(cols) < 3:
                raise SystemExit(f"ERROR: {bed}:{i}: BED needs >=3 tab-separated columns; got {cols!r}")
            try:
                start = int(cols[1])
                end = int(cols[2])
            except ValueError:
                raise SystemExit(
                    f"ERROR: {bed}:{i}: non-numeric start/end — likely a header row.\n"
                    "       igv-reports' BED parser is positional and chokes on non-comment\n"
                    "       headers. Prefix the header with `#` (skipped by create_report\n"
                    "       and matches the lab's BED-output convention) or strip it."
                )
            if start >= end:
                raise SystemExit(f"ERROR: {bed}:{i}: start ({start}) >= end ({end})")


def find_create_report() -> str:
    cr = shutil.which("create_report")
    if cr:
        return cr
    candidate = Path("/home/ahunos/miniforge3/envs/snakemake/bin/create_report")
    if candidate.exists():
        return str(candidate)
    raise SystemExit(
        "ERROR: create_report not on PATH. Activate snakemake conda env:\n"
        "  source /home/ahunos/miniforge3/etc/profile.d/conda.sh && conda activate snakemake"
    )


def apptainer_create_report_prefix(sif: Path) -> list[str]:
    """Return the `singularity exec --cleanenv --bind /data1/greenbab <sif>
    create_report` prefix. Used when --apptainer is passed; avoids the NFS
    conda cold-start tax (rules/apptainer_vs_conda.md). The SIF is a dedicated
    igv-reports container (igv-reports_1.16.0.sif, ~83 MB) pulled from the
    Galaxy depot.

    --cleanenv: scrubs host env vars so they don't leak into the SIF.
    Specifically: host SSL_CERT_FILE / SSL_CERT_DIR on RHEL 8 point at paths
    that don't exist inside Galaxy-depot SIFs, and create_report's standalone-
    HTML build path performs an HTTPS GET (for the IGV.js ideogram or similar)
    that aborts with `[SSL: CERTIFICATE_VERIFY_FAILED]`. See
    rules/apptainer_env_leak.md. The driver wraps the invocation so users
    don't need to remember the flag."""
    if not sif.exists():
        raise SystemExit(
            f"ERROR: apptainer SIF not found: {sif}\n"
            "       Pull with:\n"
            "         wget -O /data1/greenbab/users/ahunos/apps/containers/igv-reports_1.16.0.sif \\\n"
            "           'https://depot.galaxyproject.org/singularity/igv-reports:1.16.0--pyh7e72e81_0'"
        )
    return ["singularity", "exec", "--cleanenv", "--bind", "/data1/greenbab", str(sif), "create_report"]


def build_one(
    sites: Path,
    bams: list[Path],
    vcf: Path | None,
    extra_tracks: list[Path],
    fasta: str,
    default_tracks: list[str],
    output: Path,
    title: str,
    flanking: int,
    log: logging.Logger,
    track_config: Path | None = None,
    report_type: str | None = None,
    info_columns: list[str] | None = None,
    use_apptainer: bool = False,
) -> Path:
    """Run create_report for one site set and return the HTML path.

    Two track modes:
      * Default — positional `--tracks <path> <path> ...`. Used when
        `track_config` is None. BAM + VCF + extra + default annotations,
        in render order top-to-bottom.
      * track-config — `--track-config <json>`. Used when `track_config`
        is provided. The JSON is the source of truth; default_tracks,
        bams, vcf, extra_tracks are IGNORED (they go in the JSON instead).
        This is the path required for ONT methylation viewers (named
        tracks, per-track color/min/max/colorBy/displayMode).
    """
    validate_sites_bed(sites)
    output.parent.mkdir(parents=True, exist_ok=True)

    create_report_cmd = (
        apptainer_create_report_prefix(IGVREPORTS_SIF) if use_apptainer
        else [find_create_report()]
    )

    cmd: list[str] = list(create_report_cmd) + [
        str(sites),
        "--fasta", fasta,
        "--flanking", str(flanking),
    ]

    if track_config is not None:
        cmd.extend(["--track-config", str(track_config)])
        log.info(f"  track-config: {track_config}  (defaults+bams+vcf bypassed)")
        if bams or vcf or extra_tracks or default_tracks:
            log.warning(
                "--track-config supplied; ignoring --bam/--vcf/--extra-track and "
                "auto-resolved default tracks. Put everything in the JSON instead."
            )
    else:
        # Track ordering: BAMs (data) -> VCF (calls) -> extra -> defaults (annotation, last).
        tracks: list[str] = [str(b) for b in bams]
        if vcf:
            tracks.append(str(vcf))
        tracks.extend(str(t) for t in extra_tracks)
        tracks.extend(default_tracks)
        cmd.extend(["--tracks", *tracks])
        log.info(f"  tracks (in render order):")
        for i, t in enumerate(tracks, start=1):
            log.info(f"    {i:>2}. {t}")

    if report_type:
        cmd.extend(["--type", report_type])
    if info_columns:
        cmd.extend(["--info-columns", *info_columns])

    cmd.extend([
        "--standalone",
        "--title", title,
        "--output", str(output),
    ])

    log.info(f"  cmd: {' '.join(cmd)}")
    log.info(f"  flanking_bp: {flanking}")

    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        log.error(f"create_report FAILED for {sites}")
        log.error(f"stdout: {proc.stdout}")
        log.error(f"stderr: {proc.stderr}")
        raise SystemExit(proc.returncode)

    if output.exists():
        log.info(f"  HTML: {output} ({output.stat().st_size / 1024 / 1024:.2f} MB)")
    return output


def parse_samplesheet(path: Path) -> list[dict]:
    rows: list[dict] = []
    with path.open() as fh:
        header = fh.readline().lstrip("#").rstrip("\n").split("\t")
        for ln in fh:
            cols = ln.rstrip("\n").split("\t")
            if not cols or not cols[0].strip():
                continue
            row = dict(zip(header, cols))
            rows.append(row)
    required = {"sample", "sites_bed"}
    if rows and not required.issubset(rows[0].keys()):
        raise SystemExit(
            f"ERROR: samplesheet must have columns: sample, sites_bed (got {list(rows[0].keys())}).\n"
            "       Optional columns: bam_tumor, bam_normal, vcf, extra_tracks (comma-separated)."
        )
    return rows


def derive_log_path(out_dir: Path, override: Path | None = None) -> Path:
    """Choose a log dir matching the lab's `results/<run>/{reports,logs}/`
    sibling layout when possible. Fall back to `out_dir/logs/` (in-dir) when
    the sibling can't be created — `out_dir.parent` is root, read-only, or
    otherwise unwritable. Honor an explicit `override` unconditionally."""
    if override is not None:
        log_dir = override
    else:
        out_dir = out_dir.resolve()
        sibling = out_dir.parent / "logs"
        try:
            sibling.mkdir(parents=True, exist_ok=True)
            log_dir = sibling
        except (PermissionError, OSError):
            log_dir = out_dir / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    return log_dir / f"run_{datetime.now():%Y%m%d_%H%M%S}.log"


def write_index(report_paths: dict[str, Path], out: Path, title: str) -> Path:
    items = "\n".join(
        f'  <li><a href="{p.name}">{s}</a></li>'
        for s, p in sorted(report_paths.items())
    )
    out.write_text(
        "<!doctype html>\n<html><head><title>"
        + title
        + "</title><style>body{font-family:Arial,sans-serif;margin:2em}li{margin:0.4em 0}</style></head><body>\n"
        f"<h1>{title}</h1>\n<ul>\n{items}\n</ul>\n</body></html>\n"
    )
    return out


def run_cohort_verify(
    samplesheet: Path,
    reports_dir: Path,
    genome: str,
    db_config: Path,
    fail_on_fail: bool,
    log: logging.Logger,
) -> None:
    """Invoke verify_cohort.py at the end of a cohort build. Writes the TSV +
    summary next to the cohort's index.html. Fails the build if
    `fail_on_fail` is set and the verifier exits nonzero."""
    verify_script = Path(__file__).resolve().parent / "verify_cohort.py"
    if not verify_script.exists():
        log.warning(f"verify_cohort: script not found at {verify_script} — skipping")
        return
    tsv_out = reports_dir / "cohort_verify.tsv"
    md_out = reports_dir / "cohort_verify.summary.md"
    cmd = [
        sys.executable, str(verify_script),
        "--samplesheet", str(samplesheet),
        "--reports-dir", str(reports_dir),
        "--genome", genome,
        "--db-config", str(db_config),
        "--out", str(tsv_out),
        "--summary", str(md_out),
    ]
    if fail_on_fail:
        cmd.append("--fail-on-fail")
    log.info(f"verify_cohort: running {' '.join(cmd)}")
    proc = subprocess.run(cmd, capture_output=True, text=True)
    # Mirror the verifier's stdout/stderr into the run log so audit-trail stays single-source.
    for line in (proc.stdout or "").splitlines():
        log.info(f"  verify_cohort > {line}")
    if proc.stderr:
        for line in proc.stderr.splitlines():
            log.warning(f"  verify_cohort (stderr) > {line}")
    log.info(f"verify_cohort: TSV={tsv_out} summary={md_out} exit={proc.returncode}")
    if proc.returncode != 0:
        if fail_on_fail:
            raise SystemExit(proc.returncode)
        log.warning(f"verify_cohort: exited {proc.returncode} but --fail-on-fail not set; continuing")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genome", required=True, help="hg38 | mm10 | mm39 | t2t | GRCh37 (alias-tolerant)")
    ap.add_argument("--db-config", default=str(DEFAULT_DBCONFIG))
    ap.add_argument("--flanking", type=int, default=DEFAULT_FLANKING)
    ap.add_argument("--extra-track", action="append", default=[], help="(repeat) extra track path; rendered above default annotations")

    mode = ap.add_mutually_exclusive_group(required=True)
    mode.add_argument("--samplesheet", help="TSV: sample, [bam_tumor, bam_normal, vcf,] sites_bed[, extra_tracks]")
    mode.add_argument("--sites", help="path to sites BED for single-sample mode")

    ap.add_argument("--bam", nargs="*", default=[], help="BAM/CRAM tracks (single-sample mode)")
    ap.add_argument("--vcf", help="VCF track (single-sample mode)")

    ap.add_argument("--output", help="output HTML path (single-sample mode)")
    ap.add_argument("--output-dir", help="output dir for cohort mode (default: ./reports)")
    ap.add_argument("--title", default=None, help="report title; defaults to sample name + genome")

    ap.add_argument(
        "--track-config",
        help="path to a tracks.json (igv.js track config). When set, the JSON is "
             "passed straight to create_report --track-config and all default "
             "tracks / --bam / --vcf / --extra-track are bypassed. Use this for "
             "ONT methylation viewers — see examples/methylation_ont/.",
    )
    ap.add_argument(
        "--type",
        dest="report_type",
        choices=["mutation", "fusion", "junction"],
        default=None,
        help="create_report --type. Sets viewer behaviour at each site.",
    )
    ap.add_argument(
        "--info-columns",
        nargs="*",
        default=[],
        help="VCF INFO or BED columns to surface in the variant table. "
             "For BED sites, 'name' is the most useful.",
    )
    ap.add_argument(
        "--apptainer",
        action=argparse.BooleanOptionalAction,
        default=None,
        help=f"Run create_report from inside {IGVREPORTS_SIF} (dedicated "
             "igv-reports 1.16.0 SIF, ~83 MB). Skips the NFS conda cold-"
             "start tax (see rules/apptainer_vs_conda.md). Default: "
             "auto-detect — on if SLURM_JOB_ID is set (i.e. running under "
             "SLURM on a compute node, where the cold-start tax bites), "
             "off otherwise. Override either way with --apptainer / --no-apptainer.",
    )
    ap.add_argument(
        "--log-dir",
        help="explicit log directory. Default: sibling 'logs/' of the output "
             "dir (matches results/<run>/{reports,logs}/ lab layout); falls "
             "back to <out_dir>/logs/ when the sibling is unwritable.",
    )
    ap.add_argument(
        "--verify",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Run scripts/verify_cohort.py at the end of cohort builds "
             "(--samplesheet mode). Single-sample (--sites) mode is unaffected "
             "and emits no cohort verify TSV. Default: on. Use --no-verify to "
             "skip. The verifier inherits --fail-on-fail.",
    )
    ap.add_argument(
        "--fail-on-fail",
        action="store_true",
        help="Propagated to verify_cohort.py: exit nonzero if any verifier "
             "check is FAIL. Only meaningful with --verify and --samplesheet.",
    )

    args = ap.parse_args()

    genome = resolve_genome(args.genome)
    cfg = load_db_config(Path(args.db_config))
    fasta = fasta_for(cfg, genome)

    # Logger placed alongside the output. See derive_log_path docstring.
    if args.samplesheet:
        out_dir = Path(args.output_dir or "reports")
    else:
        if not args.output:
            raise SystemExit("ERROR: --output required in single-sample mode")
        out_dir = Path(args.output).parent
    out_dir.mkdir(parents=True, exist_ok=True)
    log_path = derive_log_path(out_dir, Path(args.log_dir) if args.log_dir else None)
    log = setup_logger(log_path)

    log.info(f"=== igv-reports skill, genome={genome} ===")
    log.info(f"db_config: {args.db_config}")
    log.info(f"fasta:     {fasta}")
    log.info(f"flanking:  {args.flanking} bp (default {DEFAULT_FLANKING})")

    # Resolve --apptainer auto-detect. Tri-state:
    #   user said --apptainer        -> True
    #   user said --no-apptainer     -> False
    #   user said nothing (None)     -> True iff SLURM_JOB_ID is in env
    # Rationale: on a fresh SLURM compute node, the NFS conda cold-start tax
    # (~1-2 M page faults, ~2.5 us each) is large; the dedicated SIF skips it.
    # On the login node, conda is usually warm and the simpler path wins.
    # See rules/apptainer_vs_conda.md.
    slurm_job = os.environ.get("SLURM_JOB_ID")
    if args.apptainer is None:
        args.apptainer = bool(slurm_job)
        decision = (
            f"auto-enabled (SLURM_JOB_ID={slurm_job} detected)"
            if args.apptainer else
            "auto-disabled (no SLURM_JOB_ID; conda env path)"
        )
        log.info(f"apptainer: {decision}")
    else:
        log.info(f"apptainer: {args.apptainer} (explicit)")

    default_tracks = resolve_default_tracks(cfg, genome, log)
    log.info(f"default tracks resolved: {len(default_tracks)}")
    for t in default_tracks:
        log.info(f"  - {t}")

    extra_tracks = [Path(p) for p in args.extra_track]

    track_config = Path(args.track_config) if args.track_config else None
    if track_config is not None and not track_config.exists():
        raise SystemExit(f"ERROR: --track-config file not found: {track_config}")

    if args.sites:
        title = args.title or f"{Path(args.sites).stem} ({genome})"
        build_one(
            sites=Path(args.sites),
            bams=[Path(b) for b in args.bam],
            vcf=Path(args.vcf) if args.vcf else None,
            extra_tracks=extra_tracks,
            fasta=fasta,
            default_tracks=default_tracks,
            output=Path(args.output),
            title=title,
            flanking=args.flanking,
            log=log,
            track_config=track_config,
            report_type=args.report_type,
            info_columns=args.info_columns,
            use_apptainer=args.apptainer,
        )
    else:
        rows = parse_samplesheet(Path(args.samplesheet))
        log.info(f"cohort: {len(rows)} samples from {args.samplesheet}")
        report_paths: dict[str, Path] = {}
        for row in rows:
            sample = row["sample"]
            log.info(f"=== {sample} ===")
            sites = Path(row["sites_bed"])
            bams = [Path(row[k]) for k in ("bam_tumor", "bam_normal") if row.get(k)]
            vcf = Path(row["vcf"]) if row.get("vcf") else None
            sample_extras = list(extra_tracks)
            if row.get("extra_tracks"):
                sample_extras += [Path(p.strip()) for p in row["extra_tracks"].split(",") if p.strip()]
            out_html = out_dir / f"{sample}.{genome}.html"
            title = args.title or f"{sample} ({genome})"
            build_one(
                sites=sites, bams=bams, vcf=vcf, extra_tracks=sample_extras,
                fasta=fasta, default_tracks=default_tracks,
                output=out_html, title=title, flanking=args.flanking, log=log,
                track_config=track_config,
                report_type=args.report_type,
                info_columns=args.info_columns,
                use_apptainer=args.apptainer,
            )
            report_paths[sample] = out_html
        idx = write_index(report_paths, out_dir / "index.html", f"igv-reports cohort ({genome})")
        log.info(f"Wrote cohort index: {idx}")

        if args.verify:
            run_cohort_verify(
                samplesheet=Path(args.samplesheet),
                reports_dir=out_dir,
                genome=genome,
                db_config=Path(args.db_config),
                fail_on_fail=args.fail_on_fail,
                log=log,
            )
        else:
            log.info("verify_cohort: skipped (--no-verify)")

    log.info(f"=== DONE: build_igvreports.py completed successfully ===")


if __name__ == "__main__":
    main()
