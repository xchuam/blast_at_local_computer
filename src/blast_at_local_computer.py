"""Command line interface focused on genome download preparation tasks."""
from __future__ import annotations

import argparse
import os
import re
import sys
from pathlib import Path
from typing import Iterable, List, Sequence

from blast_at_local_tools import (
    ftp_download,
    ftp_re_download,
    ftp_to_md5,
    ftp_to_rsync,
    g_unzip,
    genome_download,
    genome_re_download,
    md5_check,
    md5_download,
    md5_re_download,
    metadata_enrich,
)

DEFAULT_EMAIL = "a@email.address"
CLI_RELEASE_LABEL = "CLI-0.1"
SUPPORTED_FILE_TYPES = ("genome", "gff", "gtf", "protein")
ACCESSION_PATTERN = re.compile(r"\bGC[AF]_\d+\.\d+\b", re.IGNORECASE)


def _normalize_dir_path(path: Path, create: bool = True) -> str:
    """Return a string path that always ends with a path separator."""
    expanded = path.expanduser()
    if create:
        expanded.mkdir(parents=True, exist_ok=True)
    elif not expanded.exists():
        raise FileNotFoundError(f"Required directory not found: {expanded}")
    path_str = str(expanded)
    if not path_str.endswith(os.sep):
        path_str += os.sep
    return path_str


def _read_items_from_file(path: Path) -> List[str]:
    if not path.exists():
        raise FileNotFoundError(f"Input file not found: {path}")
    with path.open("r", encoding="utf-8") as handle:
        return [line.strip() for line in handle if line.strip()]


def _dedupe_preserve(values: Sequence[str]) -> List[str]:
    seen = set()
    ordered: List[str] = []
    for value in values:
        item = value.strip()
        if not item or item in seen:
            continue
        seen.add(item)
        ordered.append(item)
    return ordered


def _extract_accessions(text: str) -> List[str]:
    return [match.upper() for match in ACCESSION_PATTERN.findall(text or "")]


def _read_accessions_from_file(path: Path) -> List[str]:
    if not path.exists():
        raise FileNotFoundError(f"Input file not found: {path}")

    extracted: List[str] = []
    raw_lines: List[str] = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped:
                continue
            raw_lines.append(stripped)
            extracted.extend(_extract_accessions(stripped))

    if extracted:
        return _dedupe_preserve(extracted)
    return _dedupe_preserve(raw_lines)


def _parse_file_types(raw_file_types: str | None) -> List[str]:
    if raw_file_types is None:
        return ["genome"]

    parsed = [item.strip().lower() for item in raw_file_types.split(",") if item.strip()]
    if not parsed:
        return ["genome"]

    invalid = [item for item in parsed if item not in SUPPORTED_FILE_TYPES]
    if invalid:
        supported = ",".join(SUPPORTED_FILE_TYPES)
        raise SystemExit(f"Unsupported --file-types value(s): {','.join(invalid)}. Supported: {supported}")
    return _dedupe_preserve(parsed)


def _infer_run_root(path: Path) -> Path:
    expanded = path.expanduser()
    if expanded.name.lower() in {"download", "download_genome", "download_md5"}:
        return expanded.parent
    return expanded


def _write_resolved_accessions(accessions: Sequence[str], run_root: Path) -> Path:
    temp_dir = run_root / "temp"
    temp_dir.mkdir(parents=True, exist_ok=True)
    output = temp_dir / "resolved_accessions.txt"
    with output.open("w", encoding="utf-8") as handle:
        if accessions:
            handle.write("\n".join(accessions) + "\n")
    return output


def _collect_gca_ids(args: argparse.Namespace) -> Sequence[str]:
    items: List[str] = []
    if args.gca_list:
        items.extend(_read_accessions_from_file(Path(args.gca_list)))
    if getattr(args, "assembly_table", None):
        items.extend(_read_accessions_from_file(Path(args.assembly_table)))
    if args.gca:
        inline: List[str] = []
        for raw_value in args.gca:
            matches = _extract_accessions(raw_value)
            if matches:
                inline.extend(matches)
            elif raw_value.strip():
                inline.append(raw_value.strip())
        items.extend(inline)
    items = _dedupe_preserve(items)
    if not items:
        raise ValueError("No GCA identifiers were provided")
    return items


def _determine_email(provided: str | None) -> str:
    if provided:
        return provided
    env_email = os.environ.get("NCBI_EMAIL")
    if env_email:
        return env_email
    return DEFAULT_EMAIL


def _warn_high_workers(command: str, workers: int) -> None:
    if command in {"genome-download", "md5-download"} and workers > 5:
        print(
            (
                "\n"
                "================================ WARNING ================================\n"
                "NCBI endpoints may throttle heavy parallel traffic; we highly suggest\n"
                "using fewer than 5 workers for stable throughput and fair usage.\n"
                f"Current {command} workers: {workers}\n"
                "=========================================================================\n"
            ),
            file=sys.stderr,
        )
    if command == "gunzip" and workers > 15:
        print(
            (
                "\n"
                "================================ WARNING ================================\n"
                "HDD may not be able to handle more than 15 simultaneous decompressions.\n"
                f"Current gunzip workers: {workers}\n"
                "=========================================================================\n"
            ),
            file=sys.stderr,
        )


def cmd_metadata_download(args: argparse.Namespace) -> None:
    try:
        gca_ids = _collect_gca_ids(args)
    except ValueError as exc:
        raise SystemExit(str(exc))
    download_path = _normalize_dir_path(Path(args.download_path))
    email = _determine_email(args.email)
    run_root = Path(download_path)
    resolved_file = _write_resolved_accessions(gca_ids, run_root)
    error_log_path = run_root / "logs" / "link_download_error.txt"
    ftp_download(
        gca_ids,
        worker=args.workers,
        download_path=download_path,
        email=email,
        error_log_path=str(error_log_path),
    )
    print(f"Resolved accession list written to: {resolved_file}")


def cmd_metadata_retry(args: argparse.Namespace) -> None:
    download_path = _normalize_dir_path(Path(args.download_path))
    email = _determine_email(args.email)
    run_root = Path(download_path)
    error_file = (
        Path(args.error_file).expanduser()
        if args.error_file
        else run_root / "logs" / "link_download_error.txt"
    )
    error_log_path = run_root / "logs" / "link_download_error.txt"
    ftp_re_download(
        error_file=str(error_file),
        worker=args.workers,
        download_path=download_path,
        email=email,
        error_log_path=str(error_log_path),
    )


def cmd_to_rsync(args: argparse.Namespace) -> None:
    ftp_path = _normalize_dir_path(Path(args.ftp_path), create=False)
    rsync_path = _normalize_dir_path(Path(args.rsync_path))
    file_types = _parse_file_types(args.file_types)
    ftp_to_rsync(
        ftp_path=ftp_path,
        rsync_path=rsync_path,
        process_num=args.processes,
        file_types=file_types,
    )


def cmd_genome_download(args: argparse.Namespace) -> None:
    _warn_high_workers("genome-download", args.workers)
    genome_path = _normalize_dir_path(Path(args.genome_path))
    rsync_path = _normalize_dir_path(Path(args.rsync_path), create=False)
    file_types = _parse_file_types(args.file_types)
    run_root = _infer_run_root(Path(genome_path))
    logs_path = _normalize_dir_path(run_root / "logs")
    genome_download(
        genome_path=genome_path,
        rsync_path=rsync_path,
        worker=args.workers,
        file_types=file_types,
        logs_path=logs_path,
    )


def cmd_genome_retry(args: argparse.Namespace) -> None:
    genome_path = _normalize_dir_path(Path(args.genome_path), create=False)
    rsync_path = _normalize_dir_path(Path(args.rsync_path), create=False)
    file_types = _parse_file_types(args.file_types)
    run_root = _infer_run_root(Path(genome_path))
    logs_path = _normalize_dir_path(run_root / "logs")
    genome_re_download(
        genome_path=genome_path,
        rsync_path=rsync_path,
        worker=args.workers,
        file_types=file_types,
        logs_path=logs_path,
    )


def cmd_md5_address(args: argparse.Namespace) -> None:
    ftp_path = _normalize_dir_path(Path(args.ftp_path), create=False)
    md5_address_path = _normalize_dir_path(Path(args.md5_address_path))
    ftp_to_md5(
        ftp_path=ftp_path,
        md5_address_path=md5_address_path,
        process_num=args.processes,
    )


def cmd_md5_download(args: argparse.Namespace) -> None:
    _warn_high_workers("md5-download", args.workers)
    md5_address_path = _normalize_dir_path(Path(args.md5_address_path), create=False)
    md5_download_path = _normalize_dir_path(Path(args.md5_download_path))
    run_root = _infer_run_root(Path(md5_download_path))
    error_log_path = run_root / "logs" / "md5_download_error.txt"
    md5_download(
        md5_download_path=md5_download_path,
        md5_address_path=md5_address_path,
        worker=args.workers,
        error_log_path=str(error_log_path),
    )


def cmd_md5_retry(args: argparse.Namespace) -> None:
    md5_address_path = _normalize_dir_path(Path(args.md5_address_path), create=False)
    md5_download_path = _normalize_dir_path(Path(args.md5_download_path), create=False)
    run_root = _infer_run_root(Path(md5_download_path))
    error_log_path = run_root / "logs" / "md5_download_error.txt"
    md5_re_download(
        md5_download_path=md5_download_path,
        md5_address_path=md5_address_path,
        worker=args.workers,
        error_log_path=str(error_log_path),
    )


def cmd_md5_check(args: argparse.Namespace) -> None:
    generated_path = _normalize_dir_path(Path(args.generated_path), create=False)
    download_path = _normalize_dir_path(Path(args.download_path), create=False)
    run_root = _infer_run_root(Path(download_path))
    output_path = (
        Path(args.not_match_output).expanduser()
        if args.not_match_output
        else run_root / "logs" / "md5_not_match.txt"
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    md5_check(
        md5_generate_path=generated_path,
        md5_download_path=download_path,
        process_num=args.processes,
        not_match_output=str(output_path),
    )


def cmd_gunzip(args: argparse.Namespace) -> None:
    _warn_high_workers("gunzip", args.workers)
    genome_path = _normalize_dir_path(Path(args.genome_path), create=False)
    g_unzip(genome_path=genome_path, worker=args.workers)


def cmd_metadata_enrich(args: argparse.Namespace) -> None:
    json_path = _normalize_dir_path(Path(args.json_path), create=False)
    output_tsv = (
        Path(args.output_tsv).expanduser()
        if args.output_tsv
        else Path.cwd() / "metadata_enriched.tsv"
    )
    error_file = (
        Path(args.error_file).expanduser()
        if args.error_file
        else Path.cwd() / "metadata_enrich_error.txt"
    )
    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    error_file.parent.mkdir(parents=True, exist_ok=True)
    email = _determine_email(args.email)
    metadata_enrich(
        json_path=json_path,
        output_tsv=str(output_tsv),
        worker=args.workers,
        email=email,
        error_file=str(error_file),
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Utilities for downloading and verifying genome data from NCBI "
            f"({CLI_RELEASE_LABEL})"
        ),
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    metadata = subparsers.add_parser(
        "metadata-download", help="Download assembly metadata and FTP links"
    )
    metadata.add_argument("--gca", nargs="*", help="One or more GCA identifiers")
    metadata.add_argument("--gca-list", help="File containing GCA identifiers")
    metadata.add_argument(
        "--assembly-table",
        help="NCBI assembly TSV/table path; GCA/GCF accessions will be extracted automatically",
    )
    metadata.add_argument(
        "--download-path",
        default="Data/",
        help="Directory used to store metadata outputs",
    )
    metadata.add_argument("--workers", type=int, default=2, help="Thread count")
    metadata.add_argument("--email", help="Email address for Entrez access")
    metadata.set_defaults(func=cmd_metadata_download)

    metadata_retry = subparsers.add_parser(
        "metadata-retry", help="Retry failed metadata downloads"
    )
    metadata_retry.add_argument(
        "--error-file",
        help="File containing failed metadata downloads (default: <download-path>/logs/link_download_error.txt)",
    )
    metadata_retry.add_argument(
        "--download-path",
        default="Data/",
        help="Directory used to store metadata outputs",
    )
    metadata_retry.add_argument("--workers", type=int, default=2, help="Thread count")
    metadata_retry.add_argument("--email", help="Email address for Entrez access")
    metadata_retry.set_defaults(func=cmd_metadata_retry)

    metadata_enrich_parser = subparsers.add_parser(
        "metadata-enrich",
        help="Enrich metadata JSON files with BioSample attributes",
    )
    metadata_enrich_parser.add_argument(
        "--json-path",
        default="Data/jsons/",
        help="Directory containing metadata JSON files",
    )
    metadata_enrich_parser.add_argument(
        "--output-tsv",
        help="Output TSV path (default: $PWD/metadata_enriched.tsv)",
    )
    metadata_enrich_parser.add_argument(
        "--workers",
        type=int,
        default=2,
        help="Thread count for BioSample enrichment",
    )
    metadata_enrich_parser.add_argument("--email", help="Email address for Entrez access")
    metadata_enrich_parser.add_argument(
        "--error-file",
        help="Error log path (default: $PWD/metadata_enrich_error.txt)",
    )
    metadata_enrich_parser.set_defaults(func=cmd_metadata_enrich)

    rsync_parser = subparsers.add_parser(
        "make-rsync", help="Convert FTP links to rsync-style address files"
    )
    rsync_parser.add_argument(
        "--ftp-path", default="Data/ftp/", help="Directory containing FTP link files"
    )
    rsync_parser.add_argument(
        "--rsync-path",
        default="Data/rsync/",
        help="Directory to store address files for downstream downloads",
    )
    rsync_parser.add_argument(
        "--file-types",
        default="genome",
        help="Comma-separated asset types to include: genome,gff,gtf,protein",
    )
    rsync_parser.add_argument("--processes", type=int, default=1, help="Process count")
    rsync_parser.set_defaults(func=cmd_to_rsync)

    genome = subparsers.add_parser(
        "genome-download",
        help="Download genomes (auto-converts rsync/ftp links to HTTPS)",
    )
    genome.add_argument(
        "--genome-path", default="Data/download_genome/", help="Output directory"
    )
    genome.add_argument(
        "--rsync-path",
        default="Data/rsync/",
        help="Directory containing genome URL manifests (rsync/ftp/https)",
    )
    genome.add_argument(
        "--file-types",
        default="genome",
        help="Comma-separated asset types to download: genome,gff,gtf,protein",
    )
    genome.add_argument("--workers", type=int, default=2, help="Thread count")
    genome.set_defaults(func=cmd_genome_download)

    genome_retry = subparsers.add_parser(
        "genome-retry", help="Retry failed genome downloads from URL manifests"
    )
    genome_retry.add_argument(
        "--genome-path", default="Data/download_genome/", help="Output directory"
    )
    genome_retry.add_argument(
        "--rsync-path",
        default="Data/rsync/",
        help="Directory containing genome URL manifests (rsync/ftp/https)",
    )
    genome_retry.add_argument(
        "--file-types",
        default="genome",
        help="Comma-separated asset types to retry: genome,gff,gtf,protein",
    )
    genome_retry.add_argument("--workers", type=int, default=2, help="Thread count")
    genome_retry.set_defaults(func=cmd_genome_retry)

    md5_address = subparsers.add_parser(
        "md5-address", help="Generate rsync addresses for MD5 files"
    )
    md5_address.add_argument(
        "--ftp-path", default="Data/ftp/", help="Directory containing FTP link files"
    )
    md5_address.add_argument(
        "--md5-address-path",
        default="Data/md5_address/",
        help="Directory to store MD5 manifest URL files",
    )
    md5_address.add_argument("--processes", type=int, default=1, help="Process count")
    md5_address.set_defaults(func=cmd_md5_address)

    md5_dl = subparsers.add_parser(
        "md5-download",
        help="Download MD5 checksum files (auto-converts rsync/ftp links to HTTPS)",
    )
    md5_dl.add_argument(
        "--md5-address-path",
        default="Data/md5_address/",
        help="Directory containing MD5 URL manifests (rsync/ftp/https)",
    )
    md5_dl.add_argument(
        "--md5-download-path",
        default="Data/download_md5/",
        help="Directory to store downloaded MD5 files",
    )
    md5_dl.add_argument("--workers", type=int, default=2, help="Thread count")
    md5_dl.set_defaults(func=cmd_md5_download)

    md5_retry = subparsers.add_parser("md5-retry", help="Retry failed MD5 downloads")
    md5_retry.add_argument(
        "--md5-address-path",
        default="Data/md5_address/",
        help="Directory containing MD5 URL manifests (rsync/ftp/https)",
    )
    md5_retry.add_argument(
        "--md5-download-path",
        default="Data/download_md5/",
        help="Directory to store downloaded MD5 files",
    )
    md5_retry.add_argument("--workers", type=int, default=2, help="Thread count")
    md5_retry.set_defaults(func=cmd_md5_retry)

    md5_checker = subparsers.add_parser(
        "md5-check", help="Validate downloaded genomes using MD5 sums"
    )
    md5_checker.add_argument(
        "--generated-path",
        default="Data/generated_md5/",
        help="Directory containing generated MD5 sums",
    )
    md5_checker.add_argument(
        "--download-path",
        default="Data/download_md5/",
        help="Directory with downloaded MD5 sums",
    )
    md5_checker.add_argument("--processes", type=int, default=1, help="Process count")
    md5_checker.add_argument(
        "--not-match-output",
        help="Mismatch output file path (default: <download-path-parent>/logs/md5_not_match.txt)",
    )
    md5_checker.set_defaults(func=cmd_md5_check)

    gunzip_parser = subparsers.add_parser(
        "gunzip", help="Decompress downloaded genome archives"
    )
    gunzip_parser.add_argument(
        "--genome-path", default="Data/download_genome/", help="Directory with genomes"
    )
    gunzip_parser.add_argument(
        "--workers",
        type=int,
        default=4,
        help="Thread count for parallel decompression",
    )
    gunzip_parser.set_defaults(func=cmd_gunzip)

    return parser


def main(argv: Iterable[str] | None = None) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
