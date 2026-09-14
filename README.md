# rvz — Retrievalizer
**rvz** is an abbreviation for **Retrievalizer**, a tool designed to rapidly retrieve FASTA files and their metadata from NCBI in a single command-line.

# What it does?
- Retrieve FASTA files and metadata from NCBI.
- Both Genome and Gene FASTA files can be downloaded.
- Both ```.fna``` and ```.fasta``` can be downloaded.
- Downloading multiple FASTA files concurrently using threads.
- Calculate and verify MD5 checksums to ensure that every downloaded file is complete and uncorrupted.
- Automatically retry downloads up to five times in the event of an MD5 checksum mismatch or server failure.

# How's it installed?
## Via Conda:
```
Conda (mamba) create -n rvz_env -c defaults -c conda-forge -c bioconda -c samrachhan11 rvz=1.0 -y
conda activate rvz_env
rvz -h
```
## Manual installation:
Please install Python >=3.18 along with biopython, openpyxl, and pandas in the path, and run with ```python3 rvz.py -h```
# How's it used?
```
usage: rvz -i list.txt -o dir -e user@gmail.com

Download NCBI FASTA and metadata

options:
  -h, --help       Show this help message and exit
  -i, --input      Input file with accessions (one per line in text file)
  -o, --outdir     Output directory
  -e, --email      Email for Entrez
  -t, --threads    Number of threads for parallel download (default=1)
  -m, --mode {both,only_genome,only_metadata} (default=both)
                   both — download FASTA and metadata; only_genome — download FASTA only; only_metadata — download metadata only
```
# What does it output?
It produces one directory, ```genomes``` containing the downloaded FASTA file, one log file called ```rvz.log``` and a metadata file named ```metadata.xlsx```.

## Output

`rvz` generates a metadata table containing information about each retrieved sequence and its corresponding FASTA file.

| Accession  | Organism                                                                                           | Assembly Name | Submitter | Submission Date | Geo Location | Collection Date | Host | Isolation Source      | Forward Primer Name | Forward Primer Sequence | Reverse Primer Name | Reverse Primer Sequence | Sequence Length | Record Date | Submitter (Authors)                                                                                                                                                    | Genome Size | Contig Count | GC (%) | Genome File                            |
| :--------- | :------------------------------------------------------------------------------------------------- | :------------ | :-------- | :-------------- | :----------- | :-------------- | :--- | :-------------------- | :------------------ | :---------------------- | :------------------ | :---------------------- | --------------: | :---------- | :--------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ----------: | -----------: | -----: | :------------------------------------- |
| PX550059.1 | Rattus rattus voucher BP0075_R2 cytochrome c oxidase subunit I (COX1) gene, partial cds; mitochondrial      | NA            | NA        | 2025/11/23      | Cambodia: Battambang        | 12-Jul-2021              | NA   | NA | Utyr                  | accyctgtcyttagatttacagtc                      | C1L705                  | acttcdgggtgnccraaraatca                      |             724 | 23-NOV-2025 |Han,S., Rahi,P., Heng,S., Gov,P., Heng,V., Khun,L., Hoem,T., Hul,V., Hak,S., Chiek,S., Banuls,A.-L., Ferdinand,S., Cheng,S. and Guillard,B.                                             |         724 |            1 |  41.58 | `./test_rvz1/genomes/PX550059.1.fasta` |
