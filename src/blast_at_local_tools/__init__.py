"""High-level helpers for downloading genomes, managing BLAST databases, and
post-processing results.

This package exposes the public functions that used to live in the monolithic
``blast_at_local_tools.py`` module.  The implementation has been split into
focused submodules but the import surface remains unchanged so existing scripts
can continue to ``import blast_at_local_tools``.
"""

from .metadata import (
    ftp_download,
    ftp_re_download,
    get_assemblies,
    get_assembly_summary,
    metadata_enrich,
)
from .datasets_cli import (
    MAX_DATASETS_ACCESSIONS_PER_BATCH,
    normalize_include_values,
    read_accessions_from_assembly_table,
    read_accessions_from_list,
    resolve_datasets_binary,
    run_datasets_batches,
)
from .datasets_package import (
    datasets_fasta_files_for_blast,
    discover_datasets_files,
    make_blast_databases_from_datasets,
    write_datasets_file_manifest,
)
from .transfers import (
    ftp_modify,
    ftp_to_rsync,
    genome_download,
    genome_re_download,
    g_unzip,
)
from .md5_ops import (
    md5_address,
    ftp_to_md5,
    md5_download,
    md5_re_download,
    md5_check,
    md5_generate,
    md5_sum,
    md5sum_check,
)
from .blast_db import (
    database_remove_old,
    find_fasta_files,
    make_a_db,
    make_database,
    make_database_from_files,
    make_db_by_ls,
)
from .blast_pipeline import blast, sequential_blast_high, sequential_blast_high_s
from .results import (
    blast_result_df,
    blast_result_seq,
    extract_seq,
    extract_seq_list,
    extract_tab,
    extract_tab_list,
    merge_seq,
    merge_tab,
)

__all__ = [
    "blast",
    "blast_result_df",
    "blast_result_seq",
    "database_remove_old",
    "datasets_fasta_files_for_blast",
    "discover_datasets_files",
    "extract_seq",
    "extract_seq_list",
    "extract_tab",
    "extract_tab_list",
    "ftp_download",
    "ftp_modify",
    "ftp_re_download",
    "ftp_to_md5",
    "ftp_to_rsync",
    "g_unzip",
    "genome_download",
    "genome_re_download",
    "get_assemblies",
    "get_assembly_summary",
    "find_fasta_files",
    "MAX_DATASETS_ACCESSIONS_PER_BATCH",
    "make_db_by_ls",
    "make_a_db",
    "make_blast_databases_from_datasets",
    "make_database",
    "make_database_from_files",
    "metadata_enrich",
    "md5_address",
    "md5_check",
    "md5_download",
    "md5_generate",
    "md5_sum",
    "md5_re_download",
    "md5sum_check",
    "merge_seq",
    "merge_tab",
    "normalize_include_values",
    "read_accessions_from_assembly_table",
    "read_accessions_from_list",
    "resolve_datasets_binary",
    "run_datasets_batches",
    "sequential_blast_high",
    "sequential_blast_high_s",
    "write_datasets_file_manifest",
]
