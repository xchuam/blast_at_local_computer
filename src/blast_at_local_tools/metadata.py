"""Metadata download helpers built on top of Bio.Entrez."""

from __future__ import annotations

import json
import os
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from typing import Iterable, Sequence

import pandas as pd
from Bio import Entrez


def _ensure_directories(base_path: str) -> tuple[str, str, str]:
    """Create the json/ftp/genome subdirectories under ``base_path``."""

    json_dir = os.path.join(base_path, "jsons")
    ftp_dir = os.path.join(base_path, "ftp")
    genome_dir = os.path.join(base_path, "download_genome")
    for directory in (json_dir, ftp_dir, genome_dir):
        os.makedirs(directory, exist_ok=True)
    return json_dir, ftp_dir, genome_dir


def get_assembly_summary(entrez_id: Sequence[str]) -> dict:
    """Return the Entrez assembly summary for the supplied identifier list."""

    with Entrez.esummary(db="assembly", id=entrez_id, report="full") as handle:
        return Entrez.read(handle)


def get_assemblies(
    gca_id: str,
    download: bool = False,
    path: str = "Data/",
    email: str = "a@email.address",
) -> None:
    """Download assembly metadata and optionally the genome FASTA for ``gca_id``."""

    try:
        Entrez.email = email
        with Entrez.esearch(db="assembly", term=gca_id, retmax=1) as handle:
            record = Entrez.read(handle)
        ids = record["IdList"]

        summary = get_assembly_summary(ids)
        json_dir, ftp_dir, genome_dir = _ensure_directories(path)

        json_path = os.path.join(json_dir, f"{gca_id}.json")
        with open(json_path, "w", encoding="utf-8") as outfile:
            json.dump(summary, outfile)

        url = summary["DocumentSummarySet"]["DocumentSummary"][0]["FtpPath_GenBank"]
        label = os.path.basename(url)
        fasta_url = os.path.join(url, f"{label}_genomic.fna.gz")

        ftp_path = os.path.join(ftp_dir, f"{gca_id}.txt")
        with open(ftp_path, "w", encoding="utf-8") as handle:
            handle.write(fasta_url)

        if download:
            dest = os.path.join(genome_dir, f"{label}_genomic.fna.gz")
            urllib.request.urlretrieve(fasta_url, dest)
    except Exception as ex:  # pragma: no cover - network failures are contextual
        error_report = f"{gca_id}<{type(ex)} {ex.args}\n"
        with open("link_download_error.txt", "a", encoding="utf-8") as handle:
            handle.write(error_report)


def _reset_error_file(path: str) -> None:
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("")


def ftp_download(
    gca_list: Iterable[str],
    worker: int = 2,
    download_path: str = "Data/",
    email: str = "a@email.address",
) -> None:
    """Download assembly metadata and FTP links for ``gca_list``."""

    os.makedirs(download_path, exist_ok=True)
    _reset_error_file("link_download_error.txt")

    start_time = time.time()
    print(time.asctime())
    with ThreadPoolExecutor(max_workers=worker) as executor:
        for gca in gca_list:
            executor.submit(get_assemblies, gca, False, download_path, email)
    end_time = time.time()
    print(f"final time usage {end_time - start_time}")
    print(time.asctime())


def ftp_re_download(
    error_file: str = "link_download_error.txt",
    worker: int = 2,
    download_path: str = "Data/",
    email: str = "a@email.address",
) -> None:
    """Retry metadata downloads listed in ``error_file``."""

    error_df = pd.read_csv(error_file, sep="<", names=["gca", "rest"], encoding="utf-8")
    print(f"{len(error_df)} genome metadata still need to be downloaded.")

    _reset_error_file("link_download_error.txt")

    start_time = time.time()
    print(time.asctime())
    with ThreadPoolExecutor(max_workers=worker) as executor:
        for gca in error_df["gca"]:
            executor.submit(get_assemblies, gca, False, download_path, email)
    end_time = time.time()
    print(f"final time usage {end_time - start_time}")
    print(time.asctime())
