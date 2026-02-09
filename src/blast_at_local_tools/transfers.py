"""Helpers for converting FTP manifests to rsync inputs and downloading data."""

from __future__ import annotations

import os
import shutil
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from multiprocessing import Process
from pathlib import Path
from typing import Iterable, List
from urllib.parse import urlparse
from urllib.request import Request, urlopen

import numpy as np

def ftp_modify(
    readins: Iterable[str],
    batch: int,
    ftp_path: str = "Data/ftp/",
    rsync_path: str = "Data/rsync/",
) -> None:
    """Convert FTP URLs stored in ``ftp_path`` files to rsync addresses."""

    os.makedirs(rsync_path, exist_ok=True)
    output_path = os.path.join(rsync_path, f"address{batch}.txt")
    with open(output_path, "a", encoding="utf-8") as destination:
        for filename in readins:
            ftp_file = os.path.join(ftp_path, filename)
            with open(ftp_file, "r", encoding="utf-8") as source:
                address = (
                    source.read()
                    .replace("ftp", "rsync", 1)
                    .replace("\\", "/")
                    .replace("_genomic.fna.gz", "_genomic.fna.gz\n")
                )
            destination.write(address)


def ftp_to_rsync(
    ftp_path: str = "Data/ftp/",
    rsync_path: str = "Data/rsync/",
    process_num: int = 1,
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
            args=(list(chunk), index, ftp_path, rsync_path),
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


def _download_one(url: str, output_path: str) -> tuple[bool, str]:
    req = Request(url, headers={"User-Agent": "blast-at-local-computer/CLI-0.1"})
    try:
        with urlopen(req, timeout=300) as response, open(output_path, "wb") as handle:
            shutil.copyfileobj(response, handle)
        return True, ""
    except Exception as exc:  # pragma: no cover - network and remote IO dependent
        if os.path.exists(output_path):
            os.remove(output_path)
        return False, f"{type(exc).__name__}: {exc}"


def genome_download(
    genome_path: str = "Data/download_genome/",
    rsync_path: str = "Data/rsync/",
    worker: int = 2,
) -> None:
    """Download genome FASTA files, converting rsync/ftp links to HTTPS automatically."""

    os.makedirs(genome_path, exist_ok=True)
    addresses = _load_rsync_addresses(rsync_path)
    if not addresses:
        print(f"No address records found in {rsync_path}")
        return

    start_time = time.time()
    print(time.asctime())

    failures: List[tuple[str, str]] = []
    succeeded = 0
    with ThreadPoolExecutor(max_workers=max(1, worker)) as executor:
        futures = {}
        for address in addresses:
            url = _to_https_url(address)
            output = os.path.join(genome_path, _address_filename(address))
            future = executor.submit(_download_one, url, output)
            futures[future] = (address, url)

        for future in as_completed(futures):
            address, url = futures[future]
            ok, error = future.result()
            if ok:
                succeeded += 1
            else:
                failures.append((address, f"{url} :: {error}"))

    print(f"Total genome files: {len(addresses)}")
    print(f"Downloaded successfully: {succeeded}")
    print(f"Download failures: {len(failures)}")
    if failures:
        for address, error in failures:
            print(f"[genome-download-failed] {address} -> {error}")

    end_time = time.time()
    print(f"final time usage {end_time - start_time}")
    print(time.asctime())


def genome_re_download(
    genome_path: str = "Data/download_genome/",
    rsync_path: str = "Data/rsync/",
    worker: int = 2,
) -> None:
    """Retry missing genome files, converting rsync/ftp links to HTTPS automatically."""

    os.makedirs(genome_path, exist_ok=True)
    addresses = _load_rsync_addresses(rsync_path)
    existing = {name for name in os.listdir(genome_path)}

    address_map = {_address_filename(addr): addr for addr in addresses}
    pending = set(address_map).difference(existing)
    print(f"{len(pending)} genome data still need to be downloaded.")
    if not pending:
        return

    start_time = time.time()
    print(time.asctime())

    failures: List[tuple[str, str]] = []
    succeeded = 0
    with ThreadPoolExecutor(max_workers=max(1, worker)) as executor:
        futures = {}
        for name in sorted(pending):
            address = address_map[name]
            url = _to_https_url(address)
            output = os.path.join(genome_path, name)
            future = executor.submit(_download_one, url, output)
            futures[future] = (address, url)

        for future in as_completed(futures):
            address, url = futures[future]
            ok, error = future.result()
            if ok:
                succeeded += 1
            else:
                failures.append((address, f"{url} :: {error}"))

    print(f"Retried files: {len(pending)}")
    print(f"Downloaded successfully: {succeeded}")
    print(f"Download failures: {len(failures)}")
    if failures:
        for address, error in failures:
            print(f"[genome-retry-failed] {address} -> {error}")

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

    gz_files = sorted(Path(genome_path).glob("*.gz"))
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
