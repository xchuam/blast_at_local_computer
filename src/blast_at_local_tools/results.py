"""Helpers for extracting sequences and tabular data from BLAST results."""

from __future__ import annotations

import os
import subprocess
import time
from multiprocessing import Process
from typing import List, Sequence

import numpy as np
import pandas as pd
from Bio import SeqIO


def extract_seq(
    df: str,
    process_id: int,
    blast_op_path: str,
    blast_db_path: str,
    extract_seq_path: str,
    dtypes,
    blastdbcmd_bin: str,
) -> None:
    """Extract FASTA sequences for hits in a single BLAST output file."""

    df_path = os.path.join(blast_op_path, df)
    db_path = os.path.join(blast_db_path, df[:-4], df[:-4])

    try:
        blast_df = pd.read_csv(df_path, sep="\t", header=None, dtype=dtypes)
        for idx in blast_df.index:
            s_start = blast_df.loc[idx, 8]
            s_end = blast_df.loc[idx, 9]
            s_total = blast_df.loc[idx, 12]
            if {1, s_total}.isdisjoint({s_start, s_end}):
                strand = "plus"
                out_range = f"{s_start}-{s_end}"
                if s_start > s_end:
                    strand = "minus"
                    out_range = f"{s_end}-{s_start}"

                result = subprocess.run(
                    [
                        blastdbcmd_bin,
                        "-db",
                        db_path,
                        "-entry",
                        blast_df.loc[idx, 1],
                        "-line_length",
                        "60",
                        "-range",
                        out_range,
                        "-strand",
                        strand,
                    ],
                    capture_output=True,
                    text=True,
                    check=False,
                )

                output_file = os.path.join(extract_seq_path, f"{process_id}_{blast_df.loc[idx, 0]}.fas")
                with open(output_file, "a", encoding="utf-8") as handle:
                    handle.write(result.stdout)
    except Exception as ex:  # pragma: no cover - depends on BLAST output
        error_report = f"{df_path} <> {type(ex)} {ex.args}\n"
        with open("extract_seq_error.txt", "a", encoding="utf-8") as handle:
            handle.write(error_report)


def extract_seq_list(
    df_list: Sequence[str],
    process_id: int,
    blast_op_path: str,
    blast_db_path: str,
    extract_seq_path: str,
    dtypes,
    blastdbcmd_bin: str,
) -> None:
    """Extract sequences for a list of BLAST output files."""

    for df in df_list:
        extract_seq(df, process_id, blast_op_path, blast_db_path, extract_seq_path, dtypes, blastdbcmd_bin)


def merge_seq(query: str, extract_seq_path: str) -> None:
    """Merge per-thread FASTA files for ``query`` into a single file."""

    seq_file_pattern = os.path.join(extract_seq_path, f"*_{query}.fas")
    total_extract_seq = os.path.join(extract_seq_path, f"{query}_total.fas")
    subprocess.run(f"cat {seq_file_pattern} > {total_extract_seq}", shell=True, check=False)
    subprocess.run(f"rm {seq_file_pattern}", shell=True, check=False)


def blast_result_seq(
    blastdb_path: str = "Data/blast_db/",
    blast_output_path: str = "Data/blast_output/",
    query_path: str = "Data/example_query.fas",
    extract_seq_path: str = "Data/extract_seq/",
    dtypes = None,
    blastdbcmd_bin: str = "blastdbcmd",
    process_num: int = 1,
) -> None:
    """Extract BLAST hit sequences in parallel."""

    if dtypes is None:
        dtypes = {0: str, 1: str, 2: float, 3: int, 4: int, 5: int, 6: int, 7: int, 8: int, 9: int, 10: float, 11: float, 12: int}

    os.makedirs(extract_seq_path, exist_ok=True)
    with open("extract_seq_error.txt", "w", encoding="utf-8") as handle:
        handle.write("")

    file_names = os.listdir(blast_output_path)
    missions = np.array_split(file_names, process_num)
    query_list = list(SeqIO.to_dict(SeqIO.parse(query_path, "fasta")).keys())

    start_time = time.time()
    print(time.asctime())

    jobs: List[Process] = []
    for index, chunk in enumerate(missions):
        p = Process(
            target=extract_seq_list,
            args=(list(chunk), index, blast_output_path, blastdb_path, extract_seq_path, dtypes, blastdbcmd_bin),
        )
        p.start()
        jobs.append(p)

    for job in jobs:
        job.join()

    for query in query_list:
        merge_seq(query, extract_seq_path)

    end_time = time.time()
    print(f"final time usage {end_time - start_time}")
    print(time.asctime())


def extract_tab(
    df: str,
    process_id: int,
    blast_op_path: str,
    extract_tab_path: str,
    dtypes,
) -> None:
    """Extract tabular BLAST results for a single output file."""

    df_path = os.path.join(blast_op_path, df)
    try:
        blast_df = pd.read_csv(df_path, sep="\t", header=None, dtype=dtypes)
        for idx in set(blast_df.index):
            query = blast_df.loc[idx, 0]
            save_path = os.path.join(extract_tab_path, f"{process_id}_{query}.tab")
            blast_df.loc[[idx]].to_csv(save_path, header=False, sep="\t", mode="a", index=False)
    except Exception as ex:  # pragma: no cover - depends on BLAST output
        error_report = f"{df_path} <> {type(ex)} {ex.args}\n"
        with open("extract_tab_error.txt", "a", encoding="utf-8") as handle:
            handle.write(error_report)


def extract_tab_list(
    df_list: Sequence[str],
    process_id: int,
    blast_op_path: str,
    extract_tab_path: str,
    dtypes,
) -> None:
    """Extract tabular BLAST results for a list of output files."""

    for df in df_list:
        extract_tab(df, process_id, blast_op_path, extract_tab_path, dtypes)


def merge_tab(query: str, extract_tab_path: str, columns: Sequence[str]) -> None:
    """Merge per-thread tabular outputs for ``query`` into a single file."""

    tab_pattern = os.path.join(extract_tab_path, f"*_{query}.tab")
    total_extract_tab = os.path.join(extract_tab_path, f"{query}_total.tab")
    header = "\t".join(columns)
    subprocess.run(f"cat {tab_pattern} > {total_extract_tab}", shell=True, check=False)
    subprocess.run(f"rm {tab_pattern}", shell=True, check=False)
    subprocess.run(f"sed -i '1i\\{header}' {total_extract_tab}", shell=True, check=False)


def blast_result_df(
    blast_output_path: str = "Data/blast_output/",
    query_path: str = "Data/example_query.fas",
    extract_tab_path: str = "Data/extract_tab/",
    dtypes = None,
    columns: Sequence[str] | None = None,
    process_num: int = 1,
) -> None:
    """Extract BLAST hits into per-query tabular files."""

    if dtypes is None:
        dtypes = {0: str, 1: str, 2: float, 3: int, 4: int, 5: int, 6: int, 7: int, 8: int, 9: int, 10: float, 11: float, 12: int}
    if columns is None:
        columns = [
            "query_id",
            "subject_id",
            "pct_identity",
            "alignment_length",
            "mismatches",
            "gap_opens",
            "q_start",
            "q_end",
            "s_start",
            "s_end",
            "evalue",
            "bit_score",
            "s_total_length",
        ]

    os.makedirs(extract_tab_path, exist_ok=True)
    with open("extract_tab_error.txt", "w", encoding="utf-8") as handle:
        handle.write("")

    file_names = os.listdir(blast_output_path)
    missions = np.array_split(file_names, process_num)
    query_list = list(SeqIO.to_dict(SeqIO.parse(query_path, "fasta")).keys())

    start_time = time.time()
    print(time.asctime())

    jobs: List[Process] = []
    for index, chunk in enumerate(missions):
        p = Process(
            target=extract_tab_list,
            args=(list(chunk), index, blast_output_path, extract_tab_path, dtypes),
        )
        p.start()
        jobs.append(p)

    for job in jobs:
        job.join()

    for query in query_list:
        merge_tab(query, extract_tab_path, columns)

    end_time = time.time()
    print(f"final time usage {end_time - start_time}")
    print(time.asctime())
