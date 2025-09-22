"""Command line interface focused on genome download preparation tasks."""
from __future__ import annotations

import argparse
import os
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
)

DEFAULT_EMAIL = "a@email.address"


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


def _collect_gca_ids(args: argparse.Namespace) -> Sequence[str]:
    items: List[str] = []
    if args.gca_list:
        items.extend(_read_items_from_file(Path(args.gca_list)))
    if args.gca:
        items.extend(args.gca)
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


def cmd_metadata_download(args: argparse.Namespace) -> None:
    try:
        gca_ids = _collect_gca_ids(args)
    except ValueError as exc:
        raise SystemExit(str(exc))
    download_path = _normalize_dir_path(Path(args.download_path))
    email = _determine_email(args.email)
    ftp_download(gca_ids, worker=args.workers, download_path=download_path, email=email)


def cmd_metadata_retry(args: argparse.Namespace) -> None:
    download_path = _normalize_dir_path(Path(args.download_path))
    email = _determine_email(args.email)
    ftp_re_download(
        error_file=args.error_file,
        worker=args.workers,
        download_path=download_path,
        email=email,
    )


def cmd_to_rsync(args: argparse.Namespace) -> None:
    ftp_path = _normalize_dir_path(Path(args.ftp_path), create=False)
    rsync_path = _normalize_dir_path(Path(args.rsync_path))
    ftp_to_rsync(ftp_path=ftp_path, rsync_path=rsync_path, process_num=args.processes)


def cmd_genome_download(args: argparse.Namespace) -> None:
    genome_path = _normalize_dir_path(Path(args.genome_path))
    rsync_path = _normalize_dir_path(Path(args.rsync_path), create=False)
    genome_download(genome_path=genome_path, rsync_path=rsync_path, worker=args.workers)


def cmd_genome_retry(args: argparse.Namespace) -> None:
    genome_path = _normalize_dir_path(Path(args.genome_path), create=False)
    rsync_path = _normalize_dir_path(Path(args.rsync_path), create=False)
    genome_re_download(genome_path=genome_path, rsync_path=rsync_path, worker=args.workers)


def cmd_md5_address(args: argparse.Namespace) -> None:
    ftp_path = _normalize_dir_path(Path(args.ftp_path), create=False)
    md5_address_path = _normalize_dir_path(Path(args.md5_address_path))
    ftp_to_md5(
        ftp_path=ftp_path,
        md5_address_path=md5_address_path,
        process_num=args.processes,
    )


def cmd_md5_download(args: argparse.Namespace) -> None:
    md5_address_path = _normalize_dir_path(Path(args.md5_address_path), create=False)
    md5_download_path = _normalize_dir_path(Path(args.md5_download_path))
    md5_download(
        md5_download_path=md5_download_path,
        md5_address_path=md5_address_path,
        worker=args.workers,
    )


def cmd_md5_retry(args: argparse.Namespace) -> None:
    md5_address_path = _normalize_dir_path(Path(args.md5_address_path), create=False)
    md5_download_path = _normalize_dir_path(Path(args.md5_download_path), create=False)
    md5_re_download(
        md5_download_path=md5_download_path,
        md5_address_path=md5_address_path,
        worker=args.workers,
    )


def cmd_md5_check(args: argparse.Namespace) -> None:
    generated_path = _normalize_dir_path(Path(args.generated_path), create=False)
    download_path = _normalize_dir_path(Path(args.download_path), create=False)
    md5_check(
        md5_generate_path=generated_path,
        md5_download_path=download_path,
        process_num=args.processes,
        show_not_match=args.show_not_match,
    )


def cmd_gunzip(args: argparse.Namespace) -> None:
    genome_path = _normalize_dir_path(Path(args.genome_path), create=False)
    g_unzip(genome_path=genome_path)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Utilities for downloading and verifying genome data from NCBI",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    metadata = subparsers.add_parser(
        "metadata-download", help="Download assembly metadata and FTP links"
    )
    metadata.add_argument("--gca", nargs="*", help="One or more GCA identifiers")
    metadata.add_argument("--gca-list", help="File containing GCA identifiers")
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
        default="link_download_error.txt",
        help="File containing failed metadata downloads",
    )
    metadata_retry.add_argument(
        "--download-path",
        default="Data/",
        help="Directory used to store metadata outputs",
    )
    metadata_retry.add_argument("--workers", type=int, default=2, help="Thread count")
    metadata_retry.add_argument("--email", help="Email address for Entrez access")
    metadata_retry.set_defaults(func=cmd_metadata_retry)

    rsync_parser = subparsers.add_parser(
        "make-rsync", help="Convert FTP links to rsync addresses"
    )
    rsync_parser.add_argument(
        "--ftp-path", default="Data/ftp/", help="Directory containing FTP link files"
    )
    rsync_parser.add_argument(
        "--rsync-path", default="Data/rsync/", help="Directory to store rsync addresses"
    )
    rsync_parser.add_argument("--processes", type=int, default=1, help="Process count")
    rsync_parser.set_defaults(func=cmd_to_rsync)

    genome = subparsers.add_parser(
        "genome-download", help="Download genomes using rsync"
    )
    genome.add_argument(
        "--genome-path", default="Data/download_genome/", help="Output directory"
    )
    genome.add_argument(
        "--rsync-path", default="Data/rsync/", help="Directory with rsync addresses"
    )
    genome.add_argument("--workers", type=int, default=45, help="Thread count")
    genome.set_defaults(func=cmd_genome_download)

    genome_retry = subparsers.add_parser(
        "genome-retry", help="Retry failed genome downloads"
    )
    genome_retry.add_argument(
        "--genome-path", default="Data/download_genome/", help="Output directory"
    )
    genome_retry.add_argument(
        "--rsync-path", default="Data/rsync/", help="Directory with rsync addresses"
    )
    genome_retry.add_argument("--workers", type=int, default=45, help="Thread count")
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
        help="Directory to store MD5 rsync addresses",
    )
    md5_address.add_argument("--processes", type=int, default=1, help="Process count")
    md5_address.set_defaults(func=cmd_md5_address)

    md5_dl = subparsers.add_parser("md5-download", help="Download MD5 checksum files")
    md5_dl.add_argument(
        "--md5-address-path",
        default="Data/md5_address/",
        help="Directory containing MD5 rsync addresses",
    )
    md5_dl.add_argument(
        "--md5-download-path",
        default="Data/download_md5/",
        help="Directory to store downloaded MD5 files",
    )
    md5_dl.add_argument("--workers", type=int, default=45, help="Thread count")
    md5_dl.set_defaults(func=cmd_md5_download)

    md5_retry = subparsers.add_parser("md5-retry", help="Retry failed MD5 downloads")
    md5_retry.add_argument(
        "--md5-address-path",
        default="Data/md5_address/",
        help="Directory containing MD5 rsync addresses",
    )
    md5_retry.add_argument(
        "--md5-download-path",
        default="Data/download_md5/",
        help="Directory to store downloaded MD5 files",
    )
    md5_retry.add_argument("--workers", type=int, default=45, help="Thread count")
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
        "--show-not-match",
        action="store_true",
        help="Display unmatched genome identifiers",
    )
    md5_checker.set_defaults(func=cmd_md5_check)

    gunzip_parser = subparsers.add_parser(
        "gunzip", help="Decompress downloaded genome archives"
    )
    gunzip_parser.add_argument(
        "--genome-path", default="Data/download_genome/", help="Directory with genomes"
    )
    gunzip_parser.set_defaults(func=cmd_gunzip)

    return parser


def main(argv: Iterable[str] | None = None) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
