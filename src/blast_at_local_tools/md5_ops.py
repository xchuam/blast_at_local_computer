"""Utilities for generating, downloading, and validating MD5 checksum files."""

from __future__ import annotations

import os
import re
import shutil
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from multiprocessing import Manager, Process
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple
from urllib.parse import urlparse
from urllib.request import Request, urlopen

import numpy as np

KNOWN_ASSET_SUFFIXES = (
    "_genomic.fna.gz",
    "_genomic.gff.gz",
    "_genomic.gtf.gz",
    "_protein.faa.gz",
)
ACCESSION_PATTERN = re.compile(r"(GC[AF]_\d+\.\d+)")

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


def _infer_run_root(path: Path) -> Path:
    if path.name.lower() in {"download", "download_md5", "download_genome"}:
        return path.parent
    return path


def _default_md5_error_log(md5_download_path: str) -> Path:
    run_root = _infer_run_root(Path(md5_download_path).expanduser())
    return run_root / "logs" / "md5_download_error.txt"


def _parse_md5_manifest(md5_file: Path) -> Dict[str, str]:
    result: Dict[str, str] = {}
    if not md5_file.exists():
        return result

    with md5_file.open("r", encoding="utf-8") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped:
                continue
            parts = stripped.split(maxsplit=1)
            if len(parts) < 2:
                continue
            checksum, raw_path = parts
            file_name = os.path.basename(raw_path.lstrip("./"))
            if file_name and file_name not in result:
                result[file_name] = checksum
    return result


def _extract_md5_stem(filename: str) -> str:
    for suffix in KNOWN_ASSET_SUFFIXES:
        if filename.endswith(suffix):
            return filename[: -len(suffix)]
    return ""


def md5_download(
    md5_download_path: str = "Data/download_md5/",
    md5_address_path: str = "Data/md5_address/",
    worker: int = 2,
    error_log_path: str | None = None,
) -> None:
    """Download MD5 checksum files, converting rsync/ftp links to HTTPS automatically."""

    os.makedirs(md5_download_path, exist_ok=True)
    log_path = Path(error_log_path).expanduser() if error_log_path else _default_md5_error_log(md5_download_path)
    log_path.parent.mkdir(parents=True, exist_ok=True)
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

    with log_path.open("w", encoding="utf-8") as handle:
        for address, url, error in failures:
            handle.write(f"{address}\t{url}\t{error}\n")

    print(f"Total MD5 files: {len(addresses)}")
    print(f"Downloaded successfully: {succeeded}")
    print(f"Download failures: {len(failures)}")
    if failures:
        print(f"Failure details were written to {log_path}")

    end_time = time.time()
    print(f"final time usage {end_time - start_time}")
    print(time.asctime())


def md5_re_download(
    md5_download_path: str = "Data/download_md5/",
    md5_address_path: str = "Data/md5_address/",
    worker: int = 2,
    error_log_path: str | None = None,
) -> None:
    """Retry missing MD5 files, converting rsync/ftp links to HTTPS automatically."""

    os.makedirs(md5_download_path, exist_ok=True)
    log_path = Path(error_log_path).expanduser() if error_log_path else _default_md5_error_log(md5_download_path)
    log_path.parent.mkdir(parents=True, exist_ok=True)
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

    with log_path.open("w", encoding="utf-8") as handle:
        for address, url, error in failures:
            handle.write(f"{address}\t{url}\t{error}\n")

    print(f"Retried MD5 files: {len(pending)}")
    print(f"Downloaded successfully: {succeeded}")
    print(f"Download failures: {len(failures)}")
    if failures:
        print(f"Failure details were written to {log_path}")

    end_time = time.time()
    print(f"final time usage {end_time - start_time}")
    print(time.asctime())


def md5_sum(
    file_paths: Sequence[str],
    batch: int,
    md5_generate_path: str = "Data/generate_md5/",
    error_log_path: str = "md5_generated_error.txt",
) -> None:
    """Compute MD5 sums for downloaded archives listed in ``file_paths``."""

    output_file = Path(md5_generate_path) / f"md5_{batch}.txt"
    for file_path in file_paths:
        try:
            result = subprocess.run(
                ["md5sum", file_path], capture_output=True, text=True, check=False
            )
            if result.returncode == 0:
                with output_file.open("a", encoding="utf-8") as handle:
                    handle.write(result.stdout)
            else:
                error_report = f"{file_path} <> md5sum_return_code={result.returncode}\n"
                with Path(error_log_path).open("a", encoding="utf-8") as handle:
                    handle.write(error_report)
        except Exception as ex:  # pragma: no cover - dependent on environment
            error_report = f"{file_path} <> {type(ex)} {ex.args}\n"
            with Path(error_log_path).open("a", encoding="utf-8") as handle:
                handle.write(error_report)


def md5_generate(
    md5_generate_path: str = "Data/generate_md5/",
    genome_path: str = "Data/download_genome/",
    process_num: int = 1,
    error_log_path: str | None = None,
) -> None:
    """Generate MD5 checksum files for downloaded ``*.gz`` assets using ``process_num`` workers."""

    os.makedirs(md5_generate_path, exist_ok=True)
    run_root = _infer_run_root(Path(genome_path).expanduser())
    resolved_error_log = Path(error_log_path).expanduser() if error_log_path else run_root / "logs" / "md5_generated_error.txt"
    resolved_error_log.parent.mkdir(parents=True, exist_ok=True)
    resolved_error_log.write_text("", encoding="utf-8")

    archive_files = sorted(str(path) for path in Path(genome_path).rglob("*.gz") if path.is_file())
    if not archive_files:
        print(f"No .gz files found under {genome_path}")
        return

    missions = np.array_split(archive_files, process_num)

    start_time = time.time()
    print(time.asctime())

    jobs: List[Process] = []
    for index, chunk in enumerate(missions):
        p = Process(
            target=md5_sum,
            args=(list(chunk), index, md5_generate_path, str(resolved_error_log)),
        )
        p.start()
        jobs.append(p)

    for job in jobs:
        job.join()

    end_time = time.time()
    print(f"final time usage {end_time - start_time}")
    print(time.asctime())
    print(f"Generation error log: {resolved_error_log}")


def md5sum_check(
    tuples: Sequence[Tuple[str, str]],
    shared_list,
    md5_download_path: str = "Data/download_md5/",
) -> None:
    """Check generated hashes against downloaded md5checksums by exact file name."""

    manifest_cache: Dict[str, Dict[str, str]] = {}

    for file_path, generated_hash in tuples:
        filename = os.path.basename(str(file_path))
        stem = _extract_md5_stem(filename)
        if not stem:
            match = ACCESSION_PATTERN.search(filename)
            if not match:
                shared_list.append(filename)
                continue
            stem = filename.split(match.group(1))[0] + match.group(1)

        md5_file = os.path.join(md5_download_path, f"{stem}_md5.txt")
        if md5_file not in manifest_cache:
            manifest_cache[md5_file] = _parse_md5_manifest(Path(md5_file))

        expected_hash = manifest_cache[md5_file].get(filename)
        if not expected_hash or expected_hash != generated_hash:
            shared_list.append(filename)


def _load_generated_md5_records(md5_generate_path: str) -> List[Tuple[str, str]]:
    records: List[Tuple[str, str]] = []
    for filename in sorted(os.listdir(md5_generate_path)):
        file_path = Path(md5_generate_path) / filename
        if not file_path.is_file():
            continue
        with file_path.open("r", encoding="utf-8") as handle:
            for line in handle:
                stripped = line.strip()
                if not stripped:
                    continue
                parts = stripped.split(maxsplit=1)
                if len(parts) < 2:
                    continue
                records.append((parts[1], parts[0]))
    return records


def md5_check(
    md5_generate_path: str = "Data/generate_md5/",
    md5_download_path: str = "Data/download_md5/",
    process_num: int = 1,
    not_match_output: str = "md5_not_match.txt",
) -> None:
    """Compare generated MD5 checksums with the downloaded checksum files."""

    generated_records = _load_generated_md5_records(md5_generate_path)
    if not generated_records:
        print(f"No generated MD5 records found in {md5_generate_path}")
        output_path = Path(not_match_output).expanduser()
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text("", encoding="utf-8")
        print(f"Mismatch IDs written to: {output_path}")
        return

    missions = np.array_split(generated_records, process_num)

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
