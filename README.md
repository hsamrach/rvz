# rvz (Full name: Retrievalizer)
**rvz** is an abbreviation for **Retrievalizer**, a tool designed to rapidly retrieve FASTA files and their metadata from NCBI.

# What it does?
- Retrieve FASTA files and metadata from NCBI.
- Both Genome and Gene FASTA files can be downloaded.
- Downloading multiple FASTA files concurrently using threads.
- Calculate and verify MD5 checksums to ensure that every downloaded file is complete and uncorrupted.
- Automatically retry downloads up to five times in the event of an MD5 checksum mismatch or server failure.
