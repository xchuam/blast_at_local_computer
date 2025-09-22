"""Run BLAST searches in parallel across locally prepared databases."""

from __future__ import annotations

import os
import subprocess
import time
from multiprocessing import Lock, Manager, Process
from typing import List, Sequence

import numpy as np


def sequential_blast_high(
    catted_fasta_pth: str,
    db_pth: str,
    blast_res_dir: str,
    blast_bin: str,
    e_value: float,
    threads: int,
) -> None:
    """Run a BLAST search for ``catted_fasta_pth`` against a single database."""

    strain_name = os.path.basename(db_pth)
    db_name = os.path.join(db_pth, strain_name)
    save_dest = os.path.join(blast_res_dir, f"{strain_name}.tab")

    subprocess.run(
        [
            blast_bin,
            "-query",
            catted_fasta_pth,
            "-db",
            db_name,
            "-out",
            save_dest,
            "-evalue",
            str(e_value),
            "-max_hsps",
            "1",
            "-dust",
            "no",
            "-outfmt",
            "6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore slen",
            "-mt_mode",
            "1",
            "-num_threads",
            str(threads),
        ],
        check=False,
    )


def sequential_blast_high_s(
    catted_fasta_pth_: str,
    db_pth_lst: Sequence[str],
    blast_res_dir_: str,
    blast_bin_: str,
    e_value_: float,
    pro_id_: int,
    smart_threads_dict,
    st_lock_,
    available_threads: int,
) -> None:
    """Run BLAST for ``catted_fasta_pth_`` across ``db_pth_lst`` databases."""

    for db in db_pth_lst:
        with st_lock_:
            if sum(smart_threads_dict.values()) >= available_threads:
                usable_threads = smart_threads_dict[pro_id_]
            else:
                smart_threads_dict[pro_id_] += available_threads - sum(smart_threads_dict.values())
                usable_threads = smart_threads_dict[pro_id_]

        sequential_blast_high(
            catted_fasta_pth_,
            db,
            blast_res_dir_,
            blast_bin_,
            e_value_,
            usable_threads,
        )

    with st_lock_:
        smart_threads_dict[pro_id_] = 0


def blast(
    blastdb_path: str = "Data/blast_db/",
    query_path: str = "Data/example_query.fas",
    blast_output_path: str = "Data/blast_output/",
    blast_bin: str = "blastn",
    E_value: float = 1e-5,
    process_num: int = 1,
) -> None:
    """Run BLAST searches across all databases in ``blastdb_path``."""

    start_time = time.time()
    print(time.asctime())

    os.makedirs(blast_output_path, exist_ok=True)
    missions = [os.path.join(blastdb_path, name) for name in os.listdir(blastdb_path)]
    chunks = np.array_split(missions, process_num)

    with Manager() as manager:
        smart_threads = manager.dict()
        lock = Lock()
        for index in range(process_num):
            smart_threads[index] = 1

        jobs: List[Process] = []
        for pro_id, db_list in enumerate(chunks):
            p = Process(
                target=sequential_blast_high_s,
                args=(
                    query_path,
                    list(db_list),
                    blast_output_path,
                    blast_bin,
                    E_value,
                    pro_id,
                    smart_threads,
                    lock,
                    process_num,
                ),
            )
            p.start()
            jobs.append(p)

        for job in jobs:
            job.join()

    end_time = time.time()
    print(f"final time usage {end_time - start_time}")
    print(time.asctime())
