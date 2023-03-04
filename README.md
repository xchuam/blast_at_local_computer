# BLAST at local computer
Here we provide a tool **Blast_at_local_computer.ipynb** for the following works:
* download genome data from NCBI
* check downloaded genome data files by MD5 value
* make a local custom genome database
* update the local custom genome database
* run BLAST search against the local custom genome database
* extract BLAST hit sequences

This Jupyter notebook **Blast_at_local_computer.ipynb** should be ran in a linux system, in which the following packages should be installed:
* `parallel`
* `rsync` 
* `ncbi-blast+` 

There are several advantages of this tool:
* The genome downloading function can be used to easily download **a large set of genome data** from NCBI.
* Compared to the online NCBI BLAST, the local custom genome database can also include non-NCBI genome data.
* Compared to the online NCBI BLAST, there is **no maximum limitation** of the hit number. 
* The local custom genome database can be **easily updated**. 
* Compared to other local BLAST tools, this tool can provide not only BLAST result table but also **BLAST resulting hit sequences** for further analysis, such as phylogenetic tree construction. 
