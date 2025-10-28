"""Utilities for generating, downloading, and validating MD5 checksum files."""

from __future__ import annotations

import os
import subprocess
import time
from multiprocessing import Manager, Process
from typing import Iterable, List, Sequence

import numpy as np
import pandas as pd

from tasker import RsyncTasker


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


def md5_download(
    md5_download_path: str = "Data/download_md5/",
    md5_address_path: str = "Data/md5_address/",
    worker: int = 45,
) -> None:
    """Download MD5 checksum files listed in the manifests."""

    os.makedirs(md5_download_path, exist_ok=True)
    with open("md5_download_error.txt", "w", encoding="utf-8") as handle:
        handle.write("")

    addresses = _load_md5_addresses(md5_address_path)

    start_time = time.time()
    print(time.asctime())

    tasker = RsyncTasker(md5_download_path)
    tasks = []
    for address in addresses:
        name = address.split("/")[-2] + "_md5.txt"
        tasks.append(((address, os.path.join(md5_download_path, name)), 1))
    tasker.excute(tasks, max_workers=worker)

    end_time = time.time()
    print(f"final time usage {end_time - start_time}")
    print(time.asctime())


def md5_re_download(
    md5_download_path: str = "Data/download_md5/",
    md5_address_path: str = "Data/md5_address/",
    worker: int = 45,
) -> None:
    """Retry MD5 downloads that are missing from ``md5_download_path``."""

    os.makedirs(md5_download_path, exist_ok=True)
    downloaded = {name.replace("_md5.txt", "") for name in os.listdir(md5_download_path)}

    addresses = _load_md5_addresses(md5_address_path)
    address_map = {addr.split("/")[-2]: addr for addr in addresses}
    pending = set(address_map).difference(downloaded)
    print(f"{len(pending)} md5 files still need to be downloaded.")

    start_time = time.time()
    print(time.asctime())

    tasker = RsyncTasker(md5_download_path)
    tasks = []
    for name in pending:
        address = address_map[name]
        dest = os.path.join(md5_download_path, f"{name}_md5.txt")
        tasks.append(((address, dest), 1))
    tasker.excute(tasks, max_workers=worker)

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
    show_not_match: bool = False,
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

        mismatches = list(shared)

    print(f"{len(mismatches)} MD5 values cannot match.")
    if show_not_match:
        print(mismatches)
