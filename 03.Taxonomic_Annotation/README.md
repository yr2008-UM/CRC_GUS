# Pipeline used for creating the kraken database.</br>

## Repo Contents
* [krakenDB.pl](krakenDB.pl): A script for the analysis process.

## System Requirements
### Hardware Requirements
The perl scripts requires only a standard computer with enough RAM to support the in-memory operations.

### OS Requirements
The scripts have been tested on Linux system, but macOS is theoretically feasible as well.

### Software Requirements
* [perl 5](https://www.perl.org), Version 26, subversion 2 (v5.26.2)
* [Kraken2](https://ccb.jhu.edu/software/kraken2/index.shtml), Version 2.0.7-beta

## Installation
No additional installation is required. Simply download the script for use.

## DEMO
### Step 1: Download necessary files from NCBI.
* [nt.gz](https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nt.gz), Date: 17/10/2023, using Aspera for large file.
* [nucl_gb.accession2taxid](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz), Date: 17/10/2023
* [taxdump.tar.gz](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz), Date: 17/10/2023

### Step 2: Build the taxonomy database.
`kraken2-build --standard --threads 5 --db dirForDB`</br>

### Step 3: Extract nucleotide sequences of bacteria, archaea, and viruses.
`perl krakenDB.pl nodes.dmp names.dmp nucl_gb.accession2taxid nt.fna dirForDB/nt/library.fna`</br>

### Step 4: Edit the download_genomic_library.sh, to Avoid download again.

### Step 5: Build the genomic library.
`kraken2-build --download-library nt --db dirForDB --threads 10`</br>
`kraken2-build --build --db dirForDB --threads 10`</br>

## Citation
For usage of the tool, please cite the associated manuscript.
