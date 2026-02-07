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
from .blast_db import make_a_db, make_db_by_ls, make_database, database_remove_old
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
    "make_db_by_ls",
    "make_a_db",
    "make_database",
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
    "sequential_blast_high",
    "sequential_blast_high_s",
]
