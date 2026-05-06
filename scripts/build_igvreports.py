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
    """create_report's BED parser is positional; a header row → ValueError.
    Catch this before invoking create_report so the error is informative."""
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
                    "       igv-reports' BED parser is positional and chokes on headers.\n"
                    "       Strip the header and re-run."
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
) -> Path:
    """Run create_report for one site set and return the HTML path."""
    validate_sites_bed(sites)
    output.parent.mkdir(parents=True, exist_ok=True)

    # Track ordering: BAMs (data) → VCF (calls) → extra → defaults (annotation, last).
    tracks: list[str] = [str(b) for b in bams]
    if vcf:
        tracks.append(str(vcf))
    tracks.extend(str(t) for t in extra_tracks)
    tracks.extend(default_tracks)

    cmd: list[str] = [
        find_create_report(),
        str(sites),
        "--fasta", fasta,
        "--flanking", str(flanking),
        "--tracks", *tracks,
        "--standalone",
        "--title", title,
        "--output", str(output),
    ]
    log.info(f"  cmd: {' '.join(cmd)}")
    log.info(f"  flanking_bp: {flanking}")
    log.info(f"  tracks (in render order):")
    for i, t in enumerate(tracks, start=1):
        log.info(f"    {i:>2}. {t}")
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

    args = ap.parse_args()

    genome = resolve_genome(args.genome)
    cfg = load_db_config(Path(args.db_config))
    fasta = fasta_for(cfg, genome)

    # Logger placed alongside the output.
    if args.samplesheet:
        out_dir = Path(args.output_dir or "reports")
        out_dir.mkdir(parents=True, exist_ok=True)
        log_path = out_dir.parent / "logs" / f"run_{datetime.now():%Y%m%d_%H%M%S}.log"
    else:
        if not args.output:
            raise SystemExit("ERROR: --output required in single-sample mode")
        out_dir = Path(args.output).parent
        out_dir.mkdir(parents=True, exist_ok=True)
        log_path = out_dir.parent / "logs" / f"run_{datetime.now():%Y%m%d_%H%M%S}.log"
    log = setup_logger(log_path)

    log.info(f"=== igv-reports skill, genome={genome} ===")
    log.info(f"db_config: {args.db_config}")
    log.info(f"fasta:     {fasta}")
    log.info(f"flanking:  {args.flanking} bp (default {DEFAULT_FLANKING})")

    default_tracks = resolve_default_tracks(cfg, genome, log)
    log.info(f"default tracks resolved: {len(default_tracks)}")
    for t in default_tracks:
        log.info(f"  - {t}")

    extra_tracks = [Path(p) for p in args.extra_track]

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
            )
            report_paths[sample] = out_html
        idx = write_index(report_paths, out_dir / "index.html", f"igv-reports cohort ({genome})")
        log.info(f"Wrote cohort index: {idx}")

    log.info(f"=== DONE: build_igvreports.py completed successfully ===")


if __name__ == "__main__":
    main()
