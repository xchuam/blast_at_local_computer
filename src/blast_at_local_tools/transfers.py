"""Helpers for converting FTP manifests to rsync inputs and downloading data."""

from __future__ import annotations

import os
import re
import shutil
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from multiprocessing import Process
from pathlib import Path
from typing import Iterable, List, Sequence
from urllib.error import HTTPError
from urllib.parse import urlparse
from urllib.request import Request, urlopen

import numpy as np

ASSET_SUFFIX_BY_TYPE = {
    "genome": "_genomic.fna.gz",
    "gff": "_genomic.gff.gz",
    "gtf": "_genomic.gtf.gz",
    "protein": "_protein.faa.gz",
}
ASSET_TYPE_BY_SUFFIX = {suffix: asset for asset, suffix in ASSET_SUFFIX_BY_TYPE.items()}
ACCESSION_PATTERN = re.compile(r"(GC[AF]_\d+\.\d+)")


def _write_tsv(path: Path, header: Sequence[str], rows: Sequence[Sequence[str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        handle.write("\t".join(header) + "\n")
        for row in rows:
            handle.write("\t".join(row) + "\n")


def _to_rsync_url(address: str) -> str:
    clean = address.strip().replace("\\", "/")
    if clean.startswith("rsync://"):
        return clean
    if clean.startswith("ftp://"):
        return "rsync://" + clean[len("ftp://") :]
    if clean.startswith("ftp.ncbi.nlm.nih.gov/"):
        return f"rsync://{clean}"
    return clean


def _expand_asset_addresses(address: str, file_types: Sequence[str]) -> List[str]:
    clean = _to_rsync_url(address)
    if not clean:
        return []

    parsed = urlparse(clean)
    path = parsed.path.rstrip("/")
    name = os.path.basename(path)
    suffixes = tuple(ASSET_TYPE_BY_SUFFIX)

    if name.endswith(suffixes):
        base = clean[: -len(next(suffix for suffix in suffixes if name.endswith(suffix)))]
    else:
        root = clean.rstrip("/")
        label = os.path.basename(urlparse(root).path.rstrip("/"))
        base = f"{root}/{label}"

    return [f"{base}{ASSET_SUFFIX_BY_TYPE[file_type]}" for file_type in file_types]


def _address_type(address: str) -> str:
    filename = _address_filename(address)
    for suffix, file_type in ASSET_TYPE_BY_SUFFIX.items():
        if filename.endswith(suffix):
            return file_type
    return "unknown"


def _address_accession(address: str) -> str:
    filename = _address_filename(address)
    match = ACCESSION_PATTERN.search(filename)
    if match:
        return match.group(1)
    return ""


def ftp_modify(
    readins: Iterable[str],
    batch: int,
    ftp_path: str = "Data/ftp/",
    rsync_path: str = "Data/rsync/",
    file_types: Sequence[str] = ("genome",),
) -> None:
    """Convert FTP URLs stored in ``ftp_path`` files to rsync addresses."""

    os.makedirs(rsync_path, exist_ok=True)
    output_path = os.path.join(rsync_path, f"address{batch}.txt")
    with open(output_path, "a", encoding="utf-8") as destination:
        for filename in readins:
            ftp_file = os.path.join(ftp_path, filename)
            with open(ftp_file, "r", encoding="utf-8") as source:
                address = source.read().strip()
            for expanded in _expand_asset_addresses(address, file_types):
                destination.write(f"{expanded}\n")


def ftp_to_rsync(
    ftp_path: str = "Data/ftp/",
    rsync_path: str = "Data/rsync/",
    process_num: int = 1,
    file_types: Sequence[str] = ("genome",),
) -> None:
    """Transform FTP manifests into rsync address files using ``process_num`` workers."""

    file_names = os.listdir(ftp_path)
    missions = np.array_split(file_names, process_num)

    start_time = time.time()
    print(time.asctime())

    jobs: List[Process] = []
    for index, chunk in enumerate(missions):
        p = Process(
            target=ftp_modify,
            args=(list(chunk), index, ftp_path, rsync_path, list(file_types)),
        )
        p.start()
        jobs.append(p)

    for job in jobs:
        job.join()

    end_time = time.time()
    print(f"final time usage {end_time - start_time}")
    print(time.asctime())


def _load_rsync_addresses(rsync_path: str) -> List[str]:
    addresses: List[str] = []
    for filename in os.listdir(rsync_path):
        file_path = os.path.join(rsync_path, filename)
        with open(file_path, "r", encoding="utf-8") as handle:
            addresses.extend(line.strip() for line in handle if line.strip())
    return addresses


def _to_https_url(address: str) -> str:
    clean = address.strip().replace("\\", "/")
    if clean.startswith("https://"):
        return clean
    if clean.startswith("http://"):
        return "https://" + clean[len("http://") :]
    if clean.startswith("rsync://"):
        return "https://" + clean[len("rsync://") :]
    if clean.startswith("ftp://"):
        return "https://" + clean[len("ftp://") :]
    if clean.startswith("ftp.ncbi.nlm.nih.gov/"):
        return f"https://{clean}"
    return clean


def _address_filename(address: str) -> str:
    parsed = urlparse(_to_https_url(address))
    filename = os.path.basename(parsed.path.rstrip("/"))
    if not filename:
        raise ValueError(f"Cannot determine filename from address: {address}")
    return filename


def _download_one(url: str, output_path: str | Path) -> tuple[bool, str, str]:
    req = Request(url, headers={"User-Agent": "blast-at-local-computer/CLI-0.1"})
    try:
        with urlopen(req, timeout=60) as response, open(output_path, "wb") as handle:
            shutil.copyfileobj(response, handle)
        return True, "", ""
    except HTTPError as exc:  # pragma: no cover - network and remote IO dependent
        if os.path.exists(output_path):
            os.remove(output_path)
        if exc.code == 404:
            return False, "missing", f"HTTPError 404: {exc.reason}"
        return False, "error", f"HTTPError {exc.code}: {exc.reason}"
    except Exception as exc:  # pragma: no cover - network and remote IO dependent
        if os.path.exists(output_path):
            os.remove(output_path)
        return False, "error", f"{type(exc).__name__}: {exc}"


def _existing_relative_files(download_root: Path) -> set[str]:
    existing: set[str] = set()
    for file_path in download_root.rglob("*"):
        if file_path.is_file():
            existing.add(str(file_path.relative_to(download_root)))
    return existing


def _known_missing_relative_files(log_path: Path) -> set[str]:
    known: set[str] = set()
    if not log_path.exists():
        return known
    with log_path.open("r", encoding="utf-8") as handle:
        for index, line in enumerate(handle):
            if index == 0:
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            file_type = fields[1].strip()
            filename = fields[2].strip()
            if file_type and filename:
                known.add(str(Path(file_type) / filename))
    return known


def genome_download(
    genome_path: str = "Data/download_genome/",
    rsync_path: str = "Data/rsync/",
    worker: int = 2,
    file_types: Sequence[str] = ("genome",),
    logs_path: str | None = None,
) -> None:
    """Download selected asset files, converting rsync/ftp links to HTTPS automatically."""

    download_root = Path(genome_path)
    download_root.mkdir(parents=True, exist_ok=True)
    logs_dir = Path(logs_path) if logs_path else download_root.parent / "logs"
    logs_dir.mkdir(parents=True, exist_ok=True)

    selected = {item.lower() for item in file_types}
    addresses = [
        address
        for address in _load_rsync_addresses(rsync_path)
        if _address_type(address) in selected
    ]
    if not addresses:
        print(f"No address records found in {rsync_path}")
        return

    start_time = time.time()
    print(time.asctime())

    failures: List[tuple[str, str, str]] = []
    missing_assets: List[tuple[str, str, str, str, str]] = []
    succeeded = 0
    with ThreadPoolExecutor(max_workers=max(1, worker)) as executor:
        futures = {}
        for address in addresses:
            url = _to_https_url(address)
            file_type = _address_type(address)
            output_dir = download_root / file_type
            output_dir.mkdir(parents=True, exist_ok=True)
            output = str(output_dir / _address_filename(address))
            future = executor.submit(_download_one, url, output)
            futures[future] = (address, url)

        for future in as_completed(futures):
            address, url = futures[future]
            ok, status, error = future.result()
            if ok:
                succeeded += 1
            elif status == "missing":
                missing_assets.append(
                    (
                        _address_accession(address),
                        _address_type(address),
                        _address_filename(address),
                        url,
                        error,
                    )
                )
            else:
                failures.append((address, url, error))

    _write_tsv(
        logs_dir / "missing_assets.tsv",
        ("accession", "file_type", "filename", "url", "reason"),
        missing_assets,
    )
    _write_tsv(
        logs_dir / "download_failures.tsv",
        ("address", "url", "reason"),
        failures,
    )

    print(f"Total selected files: {len(addresses)}")
    print(f"Downloaded successfully: {succeeded}")
    print(f"Missing assets (warn only): {len(missing_assets)}")
    print(f"Download failures: {len(failures)}")
    if failures:
        for address, url, error in failures:
            print(f"[genome-download-failed] {address} -> {url} :: {error}")
    print(f"Missing-asset log: {logs_dir / 'missing_assets.tsv'}")
    print(f"Failure log: {logs_dir / 'download_failures.tsv'}")

    end_time = time.time()
    print(f"final time usage {end_time - start_time}")
    print(time.asctime())


def genome_re_download(
    genome_path: str = "Data/download_genome/",
    rsync_path: str = "Data/rsync/",
    worker: int = 2,
    file_types: Sequence[str] = ("genome",),
    logs_path: str | None = None,
) -> None:
    """Retry missing selected asset files, converting rsync/ftp links to HTTPS automatically."""

    download_root = Path(genome_path)
    download_root.mkdir(parents=True, exist_ok=True)
    logs_dir = Path(logs_path) if logs_path else download_root.parent / "logs"
    logs_dir.mkdir(parents=True, exist_ok=True)

    selected = {item.lower() for item in file_types}
    addresses = [
        address
        for address in _load_rsync_addresses(rsync_path)
        if _address_type(address) in selected
    ]
    existing = _existing_relative_files(download_root)

    address_map = {}
    for address in addresses:
        file_type = _address_type(address)
        name = _address_filename(address)
        rel_path = str(Path(file_type) / name)
        address_map[rel_path] = address

    known_missing = _known_missing_relative_files(logs_dir / "missing_assets.tsv")
    pending = sorted(set(address_map).difference(existing).difference(known_missing))
    print(f"{len(pending)} selected files still need to be downloaded.")
    if not pending:
        return

    start_time = time.time()
    print(time.asctime())

    failures: List[tuple[str, str, str]] = []
    missing_assets: List[tuple[str, str, str, str, str]] = []
    succeeded = 0
    with ThreadPoolExecutor(max_workers=max(1, worker)) as executor:
        futures = {}
        for rel_path in pending:
            address = address_map[rel_path]
            url = _to_https_url(address)
            output = download_root / rel_path
            output.parent.mkdir(parents=True, exist_ok=True)
            future = executor.submit(_download_one, url, output)
            futures[future] = (address, url)

        for future in as_completed(futures):
            address, url = futures[future]
            ok, status, error = future.result()
            if ok:
                succeeded += 1
            elif status == "missing":
                missing_assets.append(
                    (
                        _address_accession(address),
                        _address_type(address),
                        _address_filename(address),
                        url,
                        error,
                    )
                )
            else:
                failures.append((address, url, error))

    _write_tsv(
        logs_dir / "missing_assets.tsv",
        ("accession", "file_type", "filename", "url", "reason"),
        missing_assets,
    )
    _write_tsv(
        logs_dir / "download_failures.tsv",
        ("address", "url", "reason"),
        failures,
    )

    print(f"Retried files: {len(pending)}")
    print(f"Downloaded successfully: {succeeded}")
    print(f"Missing assets (warn only): {len(missing_assets)}")
    print(f"Download failures: {len(failures)}")
    if failures:
        for address, url, error in failures:
            print(f"[genome-retry-failed] {address} -> {url} :: {error}")
    print(f"Missing-asset log: {logs_dir / 'missing_assets.tsv'}")
    print(f"Failure log: {logs_dir / 'download_failures.tsv'}")

    end_time = time.time()
    print(f"final time usage {end_time - start_time}")
    print(time.asctime())


def _gunzip_one(gz_file: Path) -> tuple[str, bool, str]:
    result = subprocess.run(
        ["gunzip", str(gz_file)],
        capture_output=True,
        text=True,
        check=False,
    )
    if result.returncode == 0:
        return str(gz_file), True, ""
    return str(gz_file), False, result.stderr.strip()


def g_unzip(genome_path: str = "Data/download_genome/", worker: int = 4) -> None:
    """Gunzip all ``*.gz`` archives in ``genome_path`` using parallel workers."""

    gz_files = sorted(Path(genome_path).rglob("*.gz"))
    if not gz_files:
        print(f"No .gz archives found in {genome_path}")
        return

    start_time = time.time()
    print(time.asctime())

    failures: List[tuple[str, str]] = []
    succeeded = 0
    with ThreadPoolExecutor(max_workers=max(1, worker)) as executor:
        futures = [executor.submit(_gunzip_one, gz_file) for gz_file in gz_files]
        for future in as_completed(futures):
            file_path, ok, error = future.result()
            if ok:
                succeeded += 1
            else:
                failures.append((file_path, error))

    print(f"Total archives: {len(gz_files)}")
    print(f"Successfully decompressed: {succeeded}")
    print(f"Failed decompressions: {len(failures)}")
    if failures:
        for file_path, error in failures:
            print(f"[gunzip-failed] {file_path}: {error}")

    end_time = time.time()
    print(f"final time usage {end_time - start_time}")
    print(time.asctime())
