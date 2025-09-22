"""Helpers for converting FTP manifests to rsync inputs and downloading data."""

from __future__ import annotations

import os
import subprocess
import time
from multiprocessing import Process
from typing import Iterable, List

import numpy as np

from tasker import RsyncTasker


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


def genome_download(
    genome_path: str = "Data/download_genome/",
    rsync_path: str = "Data/rsync/",
    worker: int = 45,
) -> None:
    """Download genome FASTA files referenced by the rsync manifests."""

    os.makedirs(genome_path, exist_ok=True)
    addresses = _load_rsync_addresses(rsync_path)

    start_time = time.time()
    print(time.asctime())

    tasker = RsyncTasker(genome_path)
    tasks = [((address,), 1) for address in addresses]
    tasker.excute(tasks, max_workers=worker)

    end_time = time.time()
    print(f"final time usage {end_time - start_time}")
    print(time.asctime())


def genome_re_download(
    genome_path: str = "Data/download_genome/",
    rsync_path: str = "Data/rsync/",
    worker: int = 45,
) -> None:
    """Retry genome downloads that are missing from ``genome_path``."""

    os.makedirs(genome_path, exist_ok=True)
    addresses = _load_rsync_addresses(rsync_path)
    existing = {name for name in os.listdir(genome_path)}

    address_map = {addr.split("/")[-1]: addr for addr in addresses}
    pending = set(address_map).difference(existing)
    print(f"{len(pending)} genome data still need to be downloaded.")

    start_time = time.time()
    print(time.asctime())

    tasker = RsyncTasker(genome_path)
    tasks = [((address_map[name],), 1) for name in pending]
    tasker.excute(tasks, max_workers=worker)

    end_time = time.time()
    print(f"final time usage {end_time - start_time}")
    print(time.asctime())


def g_unzip(genome_path: str = "Data/download_genome/") -> None:
    """Gunzip all ``*.gz`` archives in ``genome_path``."""

    gunzip_command = f"gunzip {os.path.join(genome_path, '*.gz')}"
    subprocess.run(gunzip_command, shell=True, check=False)
