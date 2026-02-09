"""BLAST database management helpers."""

from __future__ import annotations

import os
import shutil
import subprocess
import time
from multiprocessing import Process
from typing import Iterable, List

import numpy as np


def make_a_db(
    in_abs_pth: str,
    db_dir_abs_pth: str,
    dbtype: str,
    makeblastdb_bin: str,
) -> None:
    """Create a BLAST database for a single genome file."""

    strain_name = os.path.splitext(os.path.basename(in_abs_pth))[0]
    save_dir = os.path.join(db_dir_abs_pth, strain_name)
    os.makedirs(save_dir, exist_ok=True)
    save_name = os.path.join(save_dir, strain_name)

    subprocess.run(
        [
            makeblastdb_bin,
            "-dbtype",
            dbtype,
            "-in",
            in_abs_pth,
            "-out",
            save_name,
            "-parse_seqids",
            "-logfile",
            "/dev/null",
        ],
        check=False,
    )


def make_db_by_ls(
    in_abs_pth_lst: Iterable[str],
    db_dir_abs_pth_: str,
    dbtype_: str,
    makeblastdb_bin_: str,
) -> None:
    """Create BLAST databases for each file in ``in_abs_pth_lst``."""

    for genome in in_abs_pth_lst:
        make_a_db(genome, db_dir_abs_pth_, dbtype_, makeblastdb_bin_)


def make_database(
    genome_path: str = "Data/download_genome/",
    blastdb_path: str = "Data/blast_db/",
    dbtype: str = "nucl",
    makeblastdb_bin: str = "makeblastdb",
    process_num: int = 1,
) -> None:
    """Create BLAST databases from genomes using ``process_num`` workers."""

    start_time = time.time()
    print(time.asctime())

    os.makedirs(blastdb_path, exist_ok=True)
    missions = [
        os.path.join(genome_path, name)
        for name in os.listdir(genome_path)
    ]
    chunks = np.array_split(missions, process_num)

    jobs: List[Process] = []
    for chunk in chunks:
        p = Process(
            target=make_db_by_ls,
            args=(list(chunk), blastdb_path, dbtype, makeblastdb_bin),
        )
        p.start()
        jobs.append(p)

    for job in jobs:
        job.join()

    end_time = time.time()
    print(f"final time usage {end_time - start_time}")
    print(time.asctime())


def database_remove_old(
    GCA_list_remove: Iterable[str],
    ftp_path: str = "Data/ftp/",
    blastdb_path: str = "Data/blast_db/",
    achive_removed_blastdb_path: str = "Data/blast_db_removed/",
) -> None:
    """Archive BLAST databases for assemblies that should be removed."""

    os.makedirs(achive_removed_blastdb_path, exist_ok=True)

    for gca in GCA_list_remove:
        ftp_file = os.path.join(ftp_path, f"{gca}.txt")
        with open(ftp_file, "r", encoding="utf-8") as handle:
            db_name = handle.read().replace(".fna.gz", "").split("/")[-1]

        src = os.path.join(blastdb_path, db_name)
        dst = os.path.join(achive_removed_blastdb_path, db_name)
        shutil.move(src, dst)
