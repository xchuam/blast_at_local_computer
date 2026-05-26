#!/usr/bin/env bash
set -euo pipefail

# batch 1: 3 accessions
/data/share_data/Softwares/NCBI_Datasets/datasets download genome accession --inputfile tmp/ecoli_three_input_download/batches/ncbi_dataset_batch_00001.txt --include genome,gff3,gtf --filename tmp/ecoli_three_input_download/packages/ncbi_dataset_batch_00001.zip --no-progressbar
