"""NCBI Datasets CLI batching helpers."""

from __future__ import annotations

import csv
import os
import re
import shlex
import shutil
import subprocess
import time
import zipfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List, Sequence

ACCESSION_PATTERN = re.compile(r"\b(?:GC[AF]_\d+\.\d+|PRJ[A-Z]{2}\d+)\b", re.IGNORECASE)
DEFAULT_DATASETS_BINARY = Path("/data/share_data/Softwares/NCBI_Datasets/datasets")
MAX_DATASETS_ACCESSIONS_PER_BATCH = 1000

SUPPORTED_INCLUDE_VALUES = (
    "genome",
    "rna",
    "protein",
    "cds",
    "gff3",
    "gtf",
    "gbff",
    "seq-report",
    "all",
    "none",
)
INCLUDE_ALIASES = {
    "fasta": "genome",
    "fna": "genome",
    "genomic": "genome",
    "genomic-fasta": "genome",
    "genomic_fasta": "genome",
    "gff": "gff3",
    "seq_report": "seq-report",
    "sequence-report": "seq-report",
    "sequence_report": "seq-report",
}


@dataclass(frozen=True)
class DatasetsBatch:
    """A single NCBI Datasets download batch."""

    index: int
    accessions: tuple[str, ...]
    input_file: Path
    zip_file: Path
    log_file: Path
    extract_dir: Path


@dataclass(frozen=True)
class BatchResult:
    """Execution result for one batch."""

    batch: DatasetsBatch
    status: str
    returncode: int
    message: str


@dataclass(frozen=True)
class DatasetsRunSummary:
    """Summary of prepared or executed NCBI Datasets batches."""

    total_accessions: int
    total_batches: int
    completed_batches: int
    skipped_batches: int
    failed_batches: int
    dry_run: bool
    output_dir: Path
    resolved_accessions_file: Path
    batch_manifest_file: Path
    command_manifest_file: Path
    failures_file: Path


def resolve_datasets_binary(value: str | None = None) -> str:
    """Return the datasets executable path to use."""

    if value:
        return value

    env_value = os.environ.get("NCBI_DATASETS_CLI")
    if env_value:
        return env_value

    if DEFAULT_DATASETS_BINARY.exists():
        return str(DEFAULT_DATASETS_BINARY)

    discovered = shutil.which("datasets")
    if discovered:
        return discovered

    return "datasets"


def dedupe_preserve_order(values: Iterable[str]) -> List[str]:
    """Deduplicate non-empty strings while preserving their first-seen order."""

    seen: set[str] = set()
    ordered: List[str] = []
    for value in values:
        item = value.strip().upper()
        if not item or item in seen:
            continue
        seen.add(item)
        ordered.append(item)
    return ordered


def _append_unique(
    destination: List[str],
    seen: set[str],
    values: Iterable[str],
    limit: int | None = None,
) -> bool:
    """Append unique accessions and return True when ``limit`` is reached."""

    for value in values:
        item = value.strip().upper()
        if not item or item in seen:
            continue
        seen.add(item)
        destination.append(item)
        if limit is not None and len(destination) >= limit:
            return True
    return False


def _dedupe_preserve_case(values: Iterable[str]) -> List[str]:
    seen: set[str] = set()
    ordered: List[str] = []
    for value in values:
        item = value.strip()
        key = item.lower()
        if not item or key in seen:
            continue
        seen.add(key)
        ordered.append(item)
    return ordered


def extract_accessions_from_text(text: str) -> List[str]:
    """Extract Assembly or BioProject accessions from free text."""

    return [match.upper() for match in ACCESSION_PATTERN.findall(text or "")]


def normalize_include_values(raw_values: str | Sequence[str] | None) -> List[str]:
    """Normalize user include values to the names accepted by ``datasets``."""

    if raw_values is None:
        return ["genome"]

    if isinstance(raw_values, str):
        pieces = [item.strip() for item in raw_values.split(",")]
    else:
        pieces = []
        for value in raw_values:
            pieces.extend(item.strip() for item in str(value).split(","))

    normalized: List[str] = []
    for piece in pieces:
        if not piece:
            continue
        lowered = piece.lower()
        value = INCLUDE_ALIASES.get(lowered, lowered)
        if value not in SUPPORTED_INCLUDE_VALUES:
            supported = ",".join(SUPPORTED_INCLUDE_VALUES)
            raise ValueError(f"Unsupported include value '{piece}'. Supported: {supported}")
        normalized.append(value)

    normalized = _dedupe_preserve_case(normalized)
    if not normalized:
        return ["genome"]
    if ("all" in normalized or "none" in normalized) and len(normalized) > 1:
        raise ValueError("'all' and 'none' cannot be combined with other include values")
    return normalized


def read_accessions_from_list(path: Path, limit: int | None = None) -> List[str]:
    """Read one accession per line, with forgiving extraction from each line."""

    if not path.exists():
        raise FileNotFoundError(f"Input file not found: {path}")

    found: List[str] = []
    seen: set[str] = set()
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped:
                continue
            matches = extract_accessions_from_text(stripped)
            if matches:
                reached_limit = _append_unique(found, seen, matches, limit)
            else:
                reached_limit = _append_unique(found, seen, (stripped,), limit)
            if reached_limit:
                break
    return found


def _normalized_header_map(fieldnames: Sequence[str] | None) -> dict[str, str]:
    mapping: dict[str, str] = {}
    for field in fieldnames or ():
        key = " ".join(field.strip().lower().replace("_", " ").split())
        mapping[key] = field
    return mapping


def _first_accession(value: str | None) -> str:
    matches = extract_accessions_from_text(value or "")
    return matches[0] if matches else ""


def _pick_table_accessions(
    primary: str,
    paired: str,
    source_preference: str,
) -> List[str]:
    if source_preference == "both":
        return [item for item in (primary, paired) if item]
    if source_preference == "as-is":
        return [primary] if primary else []
    if source_preference == "refseq":
        if primary.startswith("GCF_"):
            return [primary]
        if paired.startswith("GCF_"):
            return [paired]
        return [primary] if primary else []
    if source_preference == "genbank":
        if primary.startswith("GCA_"):
            return [primary]
        if paired.startswith("GCA_"):
            return [paired]
        return [primary] if primary else []
    raise ValueError(f"Unsupported source_preference: {source_preference}")


def read_accessions_from_assembly_table(
    path: Path,
    accession_column: str | None = None,
    paired_accession_column: str | None = None,
    source_preference: str = "refseq",
    limit: int | None = None,
) -> List[str]:
    """Read accessions from an NCBI assembly TSV table.

    The default parser uses the ``Assembly Accession`` column and, when present,
    the paired accession column to avoid downloading both GenBank and RefSeq rows
    for the same assembly unless ``source_preference='both'`` is requested.
    """

    if not path.exists():
        raise FileNotFoundError(f"Input file not found: {path}")
    if source_preference not in {"refseq", "genbank", "as-is", "both"}:
        raise ValueError("source_preference must be one of: refseq, genbank, as-is, both")

    with path.open("r", encoding="utf-8", newline="") as handle:
        sample = handle.read(4096)
        handle.seek(0)
        delimiter = "\t" if "\t" in sample else ","
        reader = csv.DictReader(handle, delimiter=delimiter)
        header_map = _normalized_header_map(reader.fieldnames)
        if not reader.fieldnames:
            return read_accessions_from_list(path, limit=limit)

        accession_field = (
            accession_column
            or header_map.get("assembly accession")
            or header_map.get("accession")
            or header_map.get("assembly")
        )
        paired_field = (
            paired_accession_column
            or header_map.get("assembly paired assembly accession")
            or header_map.get("paired assembly accession")
            or header_map.get("paired accession")
        )
        if not accession_field:
            handle.seek(0)
            return read_accessions_from_list(path, limit=limit)

        accessions: List[str] = []
        seen: set[str] = set()
        for row in reader:
            primary = _first_accession(row.get(accession_field))
            paired = _first_accession(row.get(paired_field)) if paired_field else ""
            reached_limit = _append_unique(
                accessions,
                seen,
                _pick_table_accessions(primary, paired, source_preference),
                limit,
            )
            if reached_limit:
                break
    return accessions


def split_batches(accessions: Sequence[str], batch_size: int) -> List[List[str]]:
    """Split accessions into batches accepted by the datasets CLI."""

    if batch_size < 1:
        raise ValueError("batch_size must be >= 1")
    if batch_size > MAX_DATASETS_ACCESSIONS_PER_BATCH:
        raise ValueError(
            f"batch_size must be <= {MAX_DATASETS_ACCESSIONS_PER_BATCH} for NCBI Datasets"
        )
    return [
        list(accessions[index : index + batch_size])
        for index in range(0, len(accessions), batch_size)
    ]


def prepare_datasets_batches(
    accessions: Sequence[str],
    output_dir: Path,
    batch_size: int = MAX_DATASETS_ACCESSIONS_PER_BATCH,
    prefix: str = "ncbi_dataset",
) -> tuple[List[DatasetsBatch], Path, Path]:
    """Write resolved accession and batch input files."""

    cleaned = dedupe_preserve_order(accessions)
    if not cleaned:
        raise ValueError("No NCBI Assembly or BioProject accessions were provided")

    batches_dir = output_dir / "batches"
    packages_dir = output_dir / "packages"
    logs_dir = output_dir / "logs"
    extracted_dir = output_dir / "extracted"
    for directory in (batches_dir, packages_dir, logs_dir, extracted_dir):
        directory.mkdir(parents=True, exist_ok=True)

    resolved_file = output_dir / "resolved_accessions.txt"
    resolved_file.write_text("\n".join(cleaned) + "\n", encoding="utf-8")

    batches: List[DatasetsBatch] = []
    for index, batch_accessions in enumerate(split_batches(cleaned, batch_size), start=1):
        batch_name = f"{prefix}_batch_{index:05d}"
        input_file = batches_dir / f"{batch_name}.txt"
        input_file.write_text("\n".join(batch_accessions) + "\n", encoding="utf-8")
        batches.append(
            DatasetsBatch(
                index=index,
                accessions=tuple(batch_accessions),
                input_file=input_file,
                zip_file=packages_dir / f"{batch_name}.zip",
                log_file=logs_dir / f"{batch_name}.log",
                extract_dir=extracted_dir / batch_name,
            )
        )

    manifest_file = output_dir / "batches_manifest.tsv"
    _write_batch_manifest(manifest_file, batches, {})
    return batches, resolved_file, manifest_file


def build_datasets_download_command(
    datasets_bin: str,
    batch: DatasetsBatch,
    include_values: Sequence[str],
    dehydrated: bool = False,
    no_progressbar: bool = True,
    fast_zip_validation: bool = False,
    api_key: str | None = None,
    extra_args: Sequence[str] | None = None,
) -> List[str]:
    """Build one ``datasets download genome accession`` command."""

    command = [
        datasets_bin,
        "download",
        "genome",
        "accession",
        "--inputfile",
        str(batch.input_file),
        "--include",
        ",".join(include_values),
        "--filename",
        str(batch.zip_file),
    ]
    if dehydrated:
        command.append("--dehydrated")
    if no_progressbar:
        command.append("--no-progressbar")
    if fast_zip_validation:
        command.append("--fast-zip-validation")
    if api_key:
        command.extend(["--api-key", api_key])
    if extra_args:
        command.extend(extra_args)
    return command


def _write_batch_manifest(
    path: Path,
    batches: Sequence[DatasetsBatch],
    statuses: dict[int, str],
) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(("batch", "accession_count", "input_file", "zip_file", "log_file", "status"))
        for batch in batches:
            writer.writerow(
                (
                    batch.index,
                    len(batch.accessions),
                    batch.input_file,
                    batch.zip_file,
                    batch.log_file,
                    statuses.get(batch.index, "prepared"),
                )
            )


def write_command_manifest(
    path: Path,
    batches: Sequence[DatasetsBatch],
    commands: Sequence[Sequence[str]],
) -> None:
    """Write shell-quoted commands for audit and manual execution."""

    with path.open("w", encoding="utf-8") as handle:
        handle.write("#!/usr/bin/env bash\n")
        handle.write("set -euo pipefail\n\n")
        for batch, command in zip(batches, commands):
            handle.write(f"# batch {batch.index}: {len(batch.accessions)} accessions\n")
            handle.write(shlex.join([str(part) for part in command]) + "\n")


def _extract_zip(batch: DatasetsBatch) -> None:
    batch.extract_dir.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(batch.zip_file) as archive:
        archive.extractall(batch.extract_dir)


def _extract_and_maybe_rehydrate(
    batch: DatasetsBatch,
    datasets_bin: str,
    rehydrate: bool,
    log,
) -> BatchResult | None:
    log.write("\nExtracting zip archive...\n")
    try:
        _extract_zip(batch)
    except Exception as exc:
        log.write(f"Extraction failed: {type(exc).__name__}: {exc}\n")
        return BatchResult(batch, "failed", 1, "zip extraction failed")

    if not rehydrate:
        return None

    rehydrate_command = [
        datasets_bin,
        "rehydrate",
        "--directory",
        str(batch.extract_dir),
    ]
    log.write("\nRehydrate command:\n")
    log.write(shlex.join(rehydrate_command) + "\n\n")
    rehydrate_result = subprocess.run(
        rehydrate_command,
        capture_output=True,
        text=True,
        check=False,
    )
    if rehydrate_result.stdout:
        log.write("REHYDRATE STDOUT:\n")
        log.write(rehydrate_result.stdout)
        log.write("\n")
    if rehydrate_result.stderr:
        log.write("REHYDRATE STDERR:\n")
        log.write(rehydrate_result.stderr)
        log.write("\n")
    if rehydrate_result.returncode != 0:
        log.write(f"Rehydrate return code: {rehydrate_result.returncode}\n")
        return BatchResult(batch, "failed", rehydrate_result.returncode, "datasets rehydrate failed")
    return None


def _run_one_batch(
    command: Sequence[str],
    batch: DatasetsBatch,
    datasets_bin: str,
    skip_existing: bool,
    extract: bool,
    rehydrate: bool,
) -> BatchResult:
    start_time = time.time()
    with batch.log_file.open("w", encoding="utf-8") as log:
        log.write("Command:\n")
        log.write(shlex.join([str(part) for part in command]) + "\n\n")
        if skip_existing and batch.zip_file.exists() and batch.zip_file.stat().st_size > 0:
            log.write("Zip already exists; download skipped.\n")
            if extract or rehydrate:
                failure = _extract_and_maybe_rehydrate(batch, datasets_bin, rehydrate, log)
                if failure is not None:
                    return failure
                elapsed = time.time() - start_time
                log.write(f"\nCompleted from existing zip in {elapsed:.2f} seconds\n")
                return BatchResult(batch, "completed", 0, "zip already exists; extracted")
            return BatchResult(batch, "skipped", 0, "zip already exists")

        result = subprocess.run(command, capture_output=True, text=True, check=False)
        if result.stdout:
            log.write("STDOUT:\n")
            log.write(result.stdout)
            log.write("\n")
        if result.stderr:
            log.write("STDERR:\n")
            log.write(result.stderr)
            log.write("\n")
        if result.returncode != 0:
            log.write(f"Return code: {result.returncode}\n")
            return BatchResult(batch, "failed", result.returncode, "datasets download failed")

        if extract or rehydrate:
            failure = _extract_and_maybe_rehydrate(batch, datasets_bin, rehydrate, log)
            if failure is not None:
                return failure

        elapsed = time.time() - start_time
        log.write(f"\nCompleted in {elapsed:.2f} seconds\n")
    return BatchResult(batch, "completed", 0, "")


def run_datasets_batches(
    accessions: Sequence[str],
    output_dir: Path,
    include_values: Sequence[str],
    datasets_bin: str | None = None,
    batch_size: int = MAX_DATASETS_ACCESSIONS_PER_BATCH,
    prefix: str = "ncbi_dataset",
    workers: int = 1,
    dehydrated: bool = False,
    dry_run: bool = False,
    skip_existing: bool = True,
    extract: bool = False,
    rehydrate: bool = False,
    no_progressbar: bool = True,
    fast_zip_validation: bool = False,
    api_key: str | None = None,
    extra_args: Sequence[str] | None = None,
) -> DatasetsRunSummary:
    """Prepare and optionally execute NCBI Datasets genome download batches."""

    output_dir.mkdir(parents=True, exist_ok=True)
    include_values = normalize_include_values(include_values)
    datasets_executable = resolve_datasets_binary(datasets_bin)
    batches, resolved_file, manifest_file = prepare_datasets_batches(
        accessions=accessions,
        output_dir=output_dir,
        batch_size=batch_size,
        prefix=prefix,
    )
    if rehydrate:
        extract = True

    commands = [
        build_datasets_download_command(
            datasets_executable,
            batch,
            include_values,
            dehydrated=dehydrated,
            no_progressbar=no_progressbar,
            fast_zip_validation=fast_zip_validation,
            api_key=api_key,
            extra_args=extra_args,
        )
        for batch in batches
    ]
    command_manifest = output_dir / "download_commands.sh"
    write_command_manifest(command_manifest, batches, commands)

    failures_file = output_dir / "logs" / "failed_batches.tsv"
    statuses: dict[int, str] = {}
    results: List[BatchResult] = []

    if dry_run:
        statuses = {batch.index: "dry-run" for batch in batches}
    else:
        with ThreadPoolExecutor(max_workers=max(1, workers)) as executor:
            futures = [
                executor.submit(
                    _run_one_batch,
                    command,
                    batch,
                    datasets_executable,
                    skip_existing,
                    extract,
                    rehydrate,
                )
                for batch, command in zip(batches, commands)
            ]
            for future in as_completed(futures):
                result = future.result()
                results.append(result)
                statuses[result.batch.index] = result.status

    _write_batch_manifest(manifest_file, batches, statuses)

    failures = [result for result in results if result.status == "failed"]
    failures_file.parent.mkdir(parents=True, exist_ok=True)
    with failures_file.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(("batch", "input_file", "zip_file", "log_file", "returncode", "message"))
        for result in failures:
            writer.writerow(
                (
                    result.batch.index,
                    result.batch.input_file,
                    result.batch.zip_file,
                    result.batch.log_file,
                    result.returncode,
                    result.message,
                )
            )

    completed = sum(1 for result in results if result.status == "completed")
    skipped = sum(1 for result in results if result.status == "skipped")
    failed = len(failures)
    return DatasetsRunSummary(
        total_accessions=sum(len(batch.accessions) for batch in batches),
        total_batches=len(batches),
        completed_batches=completed,
        skipped_batches=skipped,
        failed_batches=failed,
        dry_run=dry_run,
        output_dir=output_dir,
        resolved_accessions_file=resolved_file,
        batch_manifest_file=manifest_file,
        command_manifest_file=command_manifest,
        failures_file=failures_file,
    )
