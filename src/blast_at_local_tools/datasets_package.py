"""Helpers for working with extracted NCBI Datasets genome packages."""

from __future__ import annotations

import csv
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List, Sequence

from .blast_db import make_database_from_files

ACCESSION_DIR_PATTERN = re.compile(r"GC[AF]_\d+\.\d+", re.IGNORECASE)

FILE_TYPE_BY_NAME = {
    "genomic.gff": "gff3",
    "genomic.gtf": "gtf",
    "genomic.gbff": "gbff",
    "protein.faa": "protein",
    "rna.fna": "rna",
    "cds_from_genomic.fna": "cds",
    "sequence_report.jsonl": "seq-report",
}
FILE_TYPE_SUFFIXES = (
    ("_genomic.fna", "genome"),
    ("_genomic.gff", "gff3"),
    ("_genomic.gtf", "gtf"),
    ("_genomic.gbff", "gbff"),
    ("_protein.faa", "protein"),
)
BLAST_FASTA_TYPES = {"genome", "rna", "cds", "protein"}


@dataclass(frozen=True)
class DatasetFile:
    """One file found inside an extracted NCBI Datasets package."""

    accession: str
    file_type: str
    path: Path


def _infer_accession(path: Path) -> str:
    for parent in (path.parent, *path.parents):
        match = ACCESSION_DIR_PATTERN.fullmatch(parent.name)
        if match:
            return parent.name.upper()
    match = ACCESSION_DIR_PATTERN.search(path.name)
    return match.group(0).upper() if match else ""


def _infer_file_type(path: Path) -> str:
    name = path.name
    lowered = name.lower()
    if lowered in FILE_TYPE_BY_NAME:
        return FILE_TYPE_BY_NAME[lowered]
    for suffix, file_type in FILE_TYPE_SUFFIXES:
        if lowered.endswith(suffix):
            return file_type
    return ""


def discover_datasets_files(
    dataset_root: str | Path,
    file_types: Sequence[str] | None = None,
) -> List[DatasetFile]:
    """Discover genome package files under an extracted NCBI Datasets directory."""

    root = Path(dataset_root)
    if not root.exists():
        raise FileNotFoundError(f"Dataset root not found: {root}")

    selected = {item.lower() for item in file_types} if file_types else None
    discovered: List[DatasetFile] = []
    for path in sorted(root.rglob("*")):
        if not path.is_file():
            continue
        file_type = _infer_file_type(path)
        if not file_type:
            continue
        if selected is not None and file_type not in selected:
            continue
        discovered.append(
            DatasetFile(
                accession=_infer_accession(path),
                file_type=file_type,
                path=path,
            )
        )
    return discovered


def write_datasets_file_manifest(
    dataset_root: str | Path,
    output_tsv: str | Path,
    file_types: Sequence[str] | None = None,
) -> Path:
    """Write a TSV manifest of files discovered in extracted datasets packages."""

    output = Path(output_tsv)
    output.parent.mkdir(parents=True, exist_ok=True)
    files = discover_datasets_files(dataset_root, file_types=file_types)
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(("accession", "file_type", "path"))
        for item in files:
            writer.writerow((item.accession, item.file_type, item.path))
    return output


def datasets_fasta_files_for_blast(
    dataset_root: str | Path,
    file_type: str = "genome",
) -> List[str]:
    """Return FASTA files from extracted datasets packages for BLAST DB creation."""

    normalized = file_type.lower()
    if normalized not in BLAST_FASTA_TYPES:
        supported = ",".join(sorted(BLAST_FASTA_TYPES))
        raise ValueError(f"Unsupported BLAST FASTA type '{file_type}'. Supported: {supported}")
    return [
        str(item.path)
        for item in discover_datasets_files(dataset_root, file_types=(normalized,))
    ]


def make_blast_databases_from_datasets(
    dataset_root: str | Path,
    blastdb_path: str = "Example/output/blast_db/",
    file_type: str = "genome",
    dbtype: str = "nucl",
    makeblastdb_bin: str = "makeblastdb",
    process_num: int = 1,
) -> List[str]:
    """Create BLAST databases from FASTA files inside extracted datasets packages."""

    fasta_files = datasets_fasta_files_for_blast(dataset_root, file_type=file_type)
    make_database_from_files(
        fasta_files,
        blastdb_path=blastdb_path,
        dbtype=dbtype,
        makeblastdb_bin=makeblastdb_bin,
        process_num=process_num,
    )
    return fasta_files
