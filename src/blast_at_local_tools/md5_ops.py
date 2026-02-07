"""Utilities for generating, downloading, and validating MD5 checksum files."""

from __future__ import annotations

import os
import shutil
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from multiprocessing import Manager, Process
from pathlib import Path
from typing import Iterable, List, Sequence
from urllib.parse import urlparse
from urllib.request import Request, urlopen

import numpy as np
import pandas as pd

def md5_address(
    readins: Iterable[str],
    batch: int,
    ftp_path: str = "Data/ftp/",
    md5_address_path: str = "Data/md5_address/",
) -> None:
    """Convert FTP manifests into MD5 checksum file addresses."""

    os.makedirs(md5_address_path, exist_ok=True)
    output_path = os.path.join(md5_address_path, f"{batch}.txt")
    with open(output_path, "a", encoding="utf-8") as destination:
        for filename in readins:
            ftp_file = os.path.join(ftp_path, filename)
            with open(ftp_file, "r", encoding="utf-8") as source:
                parts = source.read().replace("ftp", "rsync", 1).replace("\\", "/").split("/")
            address_parts = parts[:-1] + ["md5checksums.txt\n"]
            destination.write("/".join(address_parts))


def ftp_to_md5(
    ftp_path: str = "Data/ftp/",
    md5_address_path: str = "Data/md5_address/",
    process_num: int = 1,
) -> None:
    """Generate MD5 checksum manifest files using ``process_num`` workers."""

    file_names = os.listdir(ftp_path)
    missions = np.array_split(file_names, process_num)

    start_time = time.time()
    print(time.asctime())

    jobs: List[Process] = []
    for index, chunk in enumerate(missions):
        p = Process(
            target=md5_address,
            args=(list(chunk), index, ftp_path, md5_address_path),
        )
        p.start()
        jobs.append(p)

    for job in jobs:
        job.join()

    end_time = time.time()
    print(f"final time usage {end_time - start_time}")
    print(time.asctime())


def _load_md5_addresses(md5_address_path: str) -> List[str]:
    addresses: List[str] = []
    for filename in os.listdir(md5_address_path):
        file_path = os.path.join(md5_address_path, filename)
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


def _md5_name(address: str) -> str:
    parsed = urlparse(_to_https_url(address))
    name = Path(parsed.path).parent.name
    if not name:
        raise ValueError(f"Cannot determine md5 filename stem from address: {address}")
    return name


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


def md5_download(
    md5_download_path: str = "Data/download_md5/",
    md5_address_path: str = "Data/md5_address/",
    worker: int = 2,
) -> None:
    """Download MD5 checksum files, converting rsync/ftp links to HTTPS automatically."""

    os.makedirs(md5_download_path, exist_ok=True)
    addresses = _load_md5_addresses(md5_address_path)
    if not addresses:
        print(f"No MD5 address records found in {md5_address_path}")
        return

    start_time = time.time()
    print(time.asctime())

    failures: List[tuple[str, str, str]] = []
    succeeded = 0
    with ThreadPoolExecutor(max_workers=max(1, worker)) as executor:
        futures = {}
        for address in addresses:
            url = _to_https_url(address)
            name = _md5_name(address)
            output = os.path.join(md5_download_path, f"{name}_md5.txt")
            future = executor.submit(_download_one, url, output)
            futures[future] = (address, url)

        for future in as_completed(futures):
            address, url = futures[future]
            ok, error = future.result()
            if ok:
                succeeded += 1
            else:
                failures.append((address, url, error))

    with open("md5_download_error.txt", "w", encoding="utf-8") as handle:
        for address, url, error in failures:
            handle.write(f"{address}\t{url}\t{error}\n")

    print(f"Total MD5 files: {len(addresses)}")
    print(f"Downloaded successfully: {succeeded}")
    print(f"Download failures: {len(failures)}")
    if failures:
        print("Failure details were written to md5_download_error.txt")

    end_time = time.time()
    print(f"final time usage {end_time - start_time}")
    print(time.asctime())


def md5_re_download(
    md5_download_path: str = "Data/download_md5/",
    md5_address_path: str = "Data/md5_address/",
    worker: int = 2,
) -> None:
    """Retry missing MD5 files, converting rsync/ftp links to HTTPS automatically."""

    os.makedirs(md5_download_path, exist_ok=True)
    downloaded = {name.replace("_md5.txt", "") for name in os.listdir(md5_download_path)}

    addresses = _load_md5_addresses(md5_address_path)
    address_map = {_md5_name(addr): addr for addr in addresses}
    pending = set(address_map).difference(downloaded)
    print(f"{len(pending)} md5 files still need to be downloaded.")
    if not pending:
        return

    start_time = time.time()
    print(time.asctime())

    failures: List[tuple[str, str, str]] = []
    succeeded = 0
    with ThreadPoolExecutor(max_workers=max(1, worker)) as executor:
        futures = {}
        for name in sorted(pending):
            address = address_map[name]
            url = _to_https_url(address)
            output = os.path.join(md5_download_path, f"{name}_md5.txt")
            future = executor.submit(_download_one, url, output)
            futures[future] = (address, url)

        for future in as_completed(futures):
            address, url = futures[future]
            ok, error = future.result()
            if ok:
                succeeded += 1
            else:
                failures.append((address, url, error))

    with open("md5_download_error.txt", "w", encoding="utf-8") as handle:
        for address, url, error in failures:
            handle.write(f"{address}\t{url}\t{error}\n")

    print(f"Retried MD5 files: {len(pending)}")
    print(f"Downloaded successfully: {succeeded}")
    print(f"Download failures: {len(failures)}")
    if failures:
        print("Failure details were written to md5_download_error.txt")

    end_time = time.time()
    print(f"final time usage {end_time - start_time}")
    print(time.asctime())


def md5_sum(
    addresses: Sequence[str],
    batch: int,
    genome_path: str = "Data/download_genome/",
    md5_generate_path: str = "Data/generate_md5/",
) -> None:
    """Compute MD5 sums for the genome files listed in ``addresses``."""

    for address in addresses:
        try:
            file_path = os.path.join(genome_path, address)
            result = subprocess.run(
                ["md5sum", file_path], capture_output=True, text=True, check=False
            )
            output_file = os.path.join(md5_generate_path, f"md5_{batch}.txt")
            with open(output_file, "a", encoding="utf-8") as handle:
                handle.write(result.stdout)
        except Exception as ex:  # pragma: no cover - dependent on environment
            error_report = f"{file_path} <> {type(ex)} {ex.args}\n"
            with open("md5_generated_error.txt", "a", encoding="utf-8") as handle:
                handle.write(error_report)


def md5_generate(
    md5_generate_path: str = "Data/generate_md5/",
    genome_path: str = "Data/download_genome/",
    process_num: int = 1,
) -> None:
    """Generate MD5 checksum files for genomes using ``process_num`` workers."""

    os.makedirs(md5_generate_path, exist_ok=True)
    genome_names = os.listdir(genome_path)
    missions = np.array_split(genome_names, process_num)

    start_time = time.time()
    print(time.asctime())

    jobs: List[Process] = []
    for index, chunk in enumerate(missions):
        p = Process(
            target=md5_sum,
            args=(list(chunk), index, genome_path, md5_generate_path),
        )
        p.start()
        jobs.append(p)

    for job in jobs:
        job.join()

    end_time = time.time()
    print(f"final time usage {end_time - start_time}")
    print(time.asctime())


def md5sum_check(
    tuples: Sequence[tuple[str, str, str]],
    shared_list,
    md5_download_path: str = "Data/download_md5/",
) -> None:
    """Check whether the downloaded MD5 matches the locally generated value."""

    for entry in tuples:
        genome_id = entry[0].split("/")[-1].replace("_genomic.fna.gz", "")
        download_file = os.path.join(md5_download_path, f"{genome_id}_md5.txt")
        download_md5 = pd.read_csv(
            download_file,
            sep=" ",
            header=None,
            encoding="utf-8",
        ).set_index([2])

        if download_md5.loc[f"./{genome_id}_genomic.fna.gz", 0] != entry[1]:
            shared_list.append(genome_id)


def md5_check(
    md5_generate_path: str = "Data/generate_md5/",
    md5_download_path: str = "Data/download_md5/",
    process_num: int = 1,
    not_match_output: str = "md5_not_match.txt",
) -> None:
    """Compare generated MD5 checksums with the downloaded checksum files."""

    generated_md5 = pd.DataFrame()
    for filename in os.listdir(md5_generate_path):
        file_path = os.path.join(md5_generate_path, filename)
        frame = pd.read_csv(file_path, sep=" ", header=None, encoding="utf-8").set_index([2])
        generated_md5 = pd.concat([generated_md5, frame])

    missions = np.array_split(list(generated_md5.itertuples(index=True, name=None)), process_num)

    with Manager() as manager:
        shared = manager.list()
        jobs: List[Process] = []
        for chunk in missions:
            p = Process(
                target=md5sum_check,
                args=(chunk, shared, md5_download_path),
            )
            p.start()
            jobs.append(p)

        for job in jobs:
            job.join()

        mismatches = sorted(set(shared))

    print(f"{len(mismatches)} MD5 values cannot match.")
    output_path = Path(not_match_output).expanduser()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as handle:
        if mismatches:
            handle.write("\n".join(mismatches) + "\n")
    print(f"Mismatch IDs written to: {output_path}")
