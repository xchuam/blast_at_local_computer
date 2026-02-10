"""Metadata download helpers built on top of Bio.Entrez."""

from __future__ import annotations

import json
import os
import re
import time
import urllib.request
import xml.etree.ElementTree as ET
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Dict, Iterable, List, Sequence

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
    error_log_path: str = "link_download_error.txt",
) -> None:
    """Download assembly metadata and optionally the genome FASTA for ``gca_id``."""

    try:
        Entrez.email = email
        with Entrez.esearch(db="assembly", term=gca_id, retmax=1) as handle:
            record = Entrez.read(handle)
        ids = record["IdList"]
        if not ids:
            raise ValueError("No assembly entry found in Entrez")

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
        error_path = Path(error_log_path).expanduser()
        error_path.parent.mkdir(parents=True, exist_ok=True)
        with error_path.open("a", encoding="utf-8") as handle:
            handle.write(error_report)


def _reset_error_file(path: str) -> None:
    target = Path(path).expanduser()
    target.parent.mkdir(parents=True, exist_ok=True)
    with target.open("w", encoding="utf-8") as handle:
        handle.write("")


def ftp_download(
    gca_list: Iterable[str],
    worker: int = 2,
    download_path: str = "Data/",
    email: str = "a@email.address",
    error_log_path: str = "link_download_error.txt",
) -> None:
    """Download assembly metadata and FTP links for ``gca_list``."""

    os.makedirs(download_path, exist_ok=True)
    _reset_error_file(error_log_path)

    start_time = time.time()
    print(time.asctime())
    with ThreadPoolExecutor(max_workers=worker) as executor:
        for gca in gca_list:
            executor.submit(get_assemblies, gca, False, download_path, email, error_log_path)
    end_time = time.time()
    print(f"final time usage {end_time - start_time}")
    print(time.asctime())


def ftp_re_download(
    error_file: str = "link_download_error.txt",
    worker: int = 2,
    download_path: str = "Data/",
    email: str = "a@email.address",
    error_log_path: str = "link_download_error.txt",
) -> None:
    """Retry metadata downloads listed in ``error_file``."""

    error_path = Path(error_file).expanduser()
    if not error_path.exists():
        print(f"Error file not found: {error_path}")
        return
    if error_path.stat().st_size == 0:
        print("0 genome metadata still need to be downloaded.")
        return

    error_df = pd.read_csv(error_path, sep="<", names=["gca", "rest"], encoding="utf-8")
    error_df = error_df[error_df["gca"].notna()]
    print(f"{len(error_df)} genome metadata still need to be downloaded.")
    if error_df.empty:
        return

    _reset_error_file(error_log_path)

    start_time = time.time()
    print(time.asctime())
    with ThreadPoolExecutor(max_workers=worker) as executor:
        for gca in error_df["gca"]:
            executor.submit(get_assemblies, str(gca).strip(), False, download_path, email, error_log_path)
    end_time = time.time()
    print(f"final time usage {end_time - start_time}")
    print(time.asctime())


def _extract_summary_record(raw_json: Dict[str, Any], fallback_id: str) -> Dict[str, str]:
    records = raw_json.get("DocumentSummarySet", {}).get("DocumentSummary", [])
    if not records:
        return {
            "assembly_id": fallback_id,
            "assembly_accession": fallback_id,
            "biosample_id": "",
            "biosample_accession": "",
        }

    summary = records[0]
    return {
        "assembly_id": fallback_id,
        "assembly_accession": str(summary.get("AssemblyAccession") or fallback_id),
        "biosample_id": str(summary.get("BioSampleId") or "").strip(),
        "biosample_accession": str(summary.get("BioSampleAccn") or "").strip(),
    }


def _normalize_attr_key(raw_key: str) -> str:
    cleaned = raw_key.strip().lower().replace("-", "_").replace(" ", "_")
    return "_".join(part for part in cleaned.split("_") if part)


def _extract_biosample_attributes(xml_payload: str) -> Dict[str, str]:
    attributes: Dict[str, str] = {}
    root = ET.fromstring(xml_payload)
    for element in root.findall(".//Attribute"):
        raw_key = (
            element.attrib.get("harmonized_name")
            or element.attrib.get("attribute_name")
            or element.attrib.get("display_name")
            or ""
        )
        raw_value = (element.text or "").strip()
        if not raw_key or not raw_value:
            continue
        key = _normalize_attr_key(raw_key)
        if key and key not in attributes:
            attributes[key] = raw_value
    return attributes


def _pick_first(attributes: Dict[str, str], keys: Sequence[str]) -> str:
    for key in keys:
        value = attributes.get(key, "").strip()
        if value:
            return value
    return ""


def _extract_year(collection_date_raw: str) -> str:
    match = re.search(r"(19|20)\d{2}", collection_date_raw)
    return match.group(0) if match else ""


def _read_json_record(json_file: Path) -> Dict[str, str]:
    with json_file.open("r", encoding="utf-8") as handle:
        raw = json.load(handle)
    return _extract_summary_record(raw, fallback_id=json_file.stem)


def _enrich_one_json(json_file: Path, email: str) -> tuple[Dict[str, str], str]:
    record = _read_json_record(json_file)
    base_row = {
        "assembly_id": record["assembly_id"],
        "assembly_accession": record["assembly_accession"],
        "biosample_id": record["biosample_id"],
        "biosample_accession": record["biosample_accession"],
        "country_location": "",
        "isolate": "",
        "host": "",
        "isolation_source": "",
        "collection_date_raw": "",
        "year": "",
    }

    biosample_id = record["biosample_id"]
    if not biosample_id:
        return base_row, f"{record['assembly_id']}\tmissing BioSampleId"

    Entrez.email = email
    with Entrez.efetch(db="biosample", id=biosample_id, retmode="xml") as handle:
        xml_payload = handle.read()
    if isinstance(xml_payload, bytes):
        xml_payload = xml_payload.decode("utf-8", errors="replace")

    attributes = _extract_biosample_attributes(xml_payload)
    collection_date_raw = _pick_first(
        attributes,
        (
            "collection_date",
            "collectiondate",
            "sampling_date",
            "date_of_collection",
        ),
    )

    base_row.update(
        {
            "country_location": _pick_first(
                attributes,
                (
                    "geo_loc_name",
                    "geographic_location",
                    "geographic_location_country_and_or_sea",
                    "country",
                    "location",
                ),
            ),
            "isolate": _pick_first(attributes, ("isolate",)),
            "host": _pick_first(attributes, ("host", "host_name")),
            "isolation_source": _pick_first(
                attributes,
                ("isolation_source", "isolation_site", "source"),
            ),
            "collection_date_raw": collection_date_raw,
            "year": _extract_year(collection_date_raw),
        }
    )
    return base_row, ""


def metadata_enrich(
    json_path: str = "Data/jsons/",
    output_tsv: str = "metadata_enriched.tsv",
    worker: int = 2,
    email: str = "a@email.address",
    error_file: str = "metadata_enrich_error.txt",
) -> None:
    """Fetch BioSample XML metadata and create an enriched TSV from assembly JSON files."""

    json_dir = Path(json_path)
    output_path = Path(output_tsv)
    error_path = Path(error_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    error_path.parent.mkdir(parents=True, exist_ok=True)

    json_files = sorted(json_dir.glob("*.json"))
    columns = [
        "assembly_id",
        "assembly_accession",
        "biosample_id",
        "biosample_accession",
        "country_location",
        "isolate",
        "host",
        "isolation_source",
        "collection_date_raw",
        "year",
    ]

    if not json_files:
        pd.DataFrame(columns=columns).to_csv(output_path, sep="\t", index=False)
        error_path.write_text("", encoding="utf-8")
        print(f"No metadata JSON files found in {json_dir}")
        print(f"Created empty TSV: {output_path}")
        print(f"Error log written to: {error_path}")
        return

    rows: List[Dict[str, str]] = []
    errors: List[str] = []

    start_time = time.time()
    print(time.asctime())

    with ThreadPoolExecutor(max_workers=max(1, worker)) as executor:
        futures = {executor.submit(_enrich_one_json, json_file, email): json_file for json_file in json_files}
        for future in as_completed(futures):
            json_file = futures[future]
            try:
                row, error_msg = future.result()
            except Exception as ex:  # pragma: no cover - network failures are contextual
                row = {
                    "assembly_id": json_file.stem,
                    "assembly_accession": json_file.stem,
                    "biosample_id": "",
                    "biosample_accession": "",
                    "country_location": "",
                    "isolate": "",
                    "host": "",
                    "isolation_source": "",
                    "collection_date_raw": "",
                    "year": "",
                }
                error_msg = f"{json_file.stem}\t{type(ex).__name__}: {ex}"
            rows.append(row)
            if error_msg:
                errors.append(error_msg)

    data_frame = pd.DataFrame(rows, columns=columns).sort_values(
        by=["assembly_accession", "assembly_id"],
        kind="mergesort",
    )
    data_frame.to_csv(output_path, sep="\t", index=False)
    with error_path.open("w", encoding="utf-8") as handle:
        if errors:
            handle.write("\n".join(errors) + "\n")

    end_time = time.time()
    print(f"Processed JSON files: {len(json_files)}")
    print(f"Rows written: {len(data_frame)}")
    print(f"Records with enrichment errors: {len(errors)}")
    print(f"Enriched metadata TSV: {output_path}")
    print(f"Error log written to: {error_path}")
    print(f"final time usage {end_time - start_time}")
    print(time.asctime())
