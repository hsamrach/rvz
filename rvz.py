#!/usr/bin/env python3
# All rights reserved. This script is provided "as is" without warranty of any kind.
# Please see the LICENSE file for details on usage and redistribution.
# Samrach HAN
# required: pip install pandas biopython openpyxl

import argparse
import os
import re
import sys
import logging
import gzip
import urllib.request
import urllib.error
import hashlib
import time
import random
import threading
from http.client import IncompleteRead

import pandas as pd
from Bio import Entrez, SeqIO
from concurrent.futures import ThreadPoolExecutor, as_completed

def parse_args():
    parser = argparse.ArgumentParser(usage="rvz -i list.txt -o dir -e user@gmail.com", description="Download genomes and metadata from NCBI")
    parser.add_argument("-i", "--input", metavar="txt", required=True, help="Input file with accessions (one per line in text file)")
    parser.add_argument("-o", "--outdir", metavar="dir", required=True, help="Output directory")
    parser.add_argument("-e", "--email", metavar="email", required=True, help="Email for Entrez")
    parser.add_argument("-t", "--threads", metavar="int", type=int, default=1,
        help="Number of threads for parallel download (default=1)")
    parser.add_argument("-m", "--mode",
        choices=["both", "only_genome", "only_metadata"], default="both",
        help="  both           — download FASTA and metadata\n"
             "  only_genome    — download FASTA only, skip metadata\n"
             "  only_metadata  — download metadata only")
    return parser.parse_args()

def setup_logging(outdir):
    log_file = os.path.join(outdir, "rvz.log")
    root = logging.getLogger()
    if root.handlers:
        root.handlers.clear()
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler(sys.stdout)
        ]
    )

def genome_valid(path):
    if not os.path.exists(path):
        return False
    try:
        opener = gzip.open if path.endswith(".gz") else open
        with opener(path, "rt") as handle:
            for _ in SeqIO.parse(handle, "fasta"):
                return True
        return False
    except Exception:
        return False

def compute_stats(path):
    total_len = 0
    contigs = 0
    gc = 0

    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            seq = str(record.seq).upper()
            total_len += len(seq)
            contigs += 1
            gc += seq.count("G") + seq.count("C")

    gc_percent = round((gc / total_len) * 100, 3) if total_len else 0
    return total_len, contigs, gc_percent

def compute_md5(filepath, chunk_size=8192):
    md5 = hashlib.md5()
    with open(filepath, "rb") as f:
        for chunk in iter(lambda: f.read(chunk_size), b""):
            md5.update(chunk)
    return md5.hexdigest()

_RETRIABLE_EXCEPTIONS = (
    IncompleteRead,
    ConnectionResetError,
    ConnectionAbortedError,
    urllib.error.URLError,
    OSError,
    TimeoutError,
)

_RETRIABLE_KEYWORDS = (
    "read failed",
    "eof",
    "incompleteread",
    "connection reset",
    "connection aborted",
    "timed out",
    "timeout",
    "broken pipe",
    "remote end closed",
    "429",
    "rate limit",
    "too many requests",
)


class NCBIRateLimiter:

    def __init__(self, min_interval: float = 0.34):
        self._min_interval = min_interval
        self._lock = threading.Lock()
        self._last = 0.0

    def wait(self):
        with self._lock:
            now = time.time()
            delta = self._last + self._min_interval - now
            if delta > 0:
                time.sleep(delta)
            self._last = time.time()

_ncbi_limiter = NCBIRateLimiter()


def _is_retriable_ncbi_error(exc: Exception) -> bool:
    msg = str(exc).lower()
    return any(kw in msg for kw in _RETRIABLE_KEYWORDS)


def _maybe_retry(attempt: int, retries: int, base_delay: float, exc: Exception) -> None:
    if attempt < retries - 1:
        wait = base_delay * (2 ** attempt) + random.uniform(0.5, 2.0)
        logging.warning(
            "Transient NCBI error (attempt %d/%d): %s — retrying in %.1fs",
            attempt + 1, retries, str(exc)[:200], wait,
        )
        time.sleep(wait)
    else:
        raise exc


def entrez_call(func, *args, retries: int = 5, base_delay: float = 3.0,
                read_func=None, **kwargs):

    if read_func is None:
        read_func = Entrez.read

    for attempt in range(retries):
        handle = None
        try:
            _ncbi_limiter.wait()
            handle = func(*args, **kwargs)
            result = read_func(handle)
            return result
        except _RETRIABLE_EXCEPTIONS as e:
            _maybe_retry(attempt, retries, base_delay, e)
        except Exception as e:
            if _is_retriable_ncbi_error(e):
                _maybe_retry(attempt, retries, base_delay, e)
            else:
                raise
        finally:

            if handle is not None:
                try:
                    handle.close()
                except Exception:
                    pass

    raise RuntimeError(f"entrez_call: exhausted {retries} retries for {func.__name__}")

def download_assembly(acc, genome_dir):

    result = entrez_call(Entrez.esearch, db="assembly", term=acc)
    time.sleep(0.2)

    if not result["IdList"]:
        raise ValueError("Assembly not found")

    uid = result["IdList"][0]

    summary = entrez_call(Entrez.esummary, db="assembly", id=uid, report="full")
    doc = summary["DocumentSummarySet"]["DocumentSummary"][0]

    ftp_path = doc.get("FtpPath_RefSeq") or doc.get("FtpPath_GenBank")
    if not ftp_path:
        raise ValueError("No FTP path available")

    fname = ftp_path.split("/")[-1] + "_genomic.fna.gz"
    url = f"{ftp_path}/{fname}"
    md5_url = f"{ftp_path}/md5checksums.txt"

    gz_path = os.path.join(genome_dir, f"{acc}.fna.gz")
    md5_file = os.path.join(genome_dir, f"{acc}_md5.txt")

    if genome_valid(gz_path):
        logging.info(f"{acc} already downloaded and valid")
        return gz_path, uid, doc

    logging.info(f"Fetching MD5 checksum for {acc}")
    try:
        urllib.request.urlretrieve(md5_url, md5_file)

        expected_md5 = None
        with open(md5_file) as f:
            for line in f:
                if fname in line:
                    expected_md5 = line.split()[0]
                    break
    finally:
        if os.path.exists(md5_file):
            os.remove(md5_file)

    if not expected_md5:
        raise ValueError("MD5 checksum not found in md5checksums.txt")

    for attempt in range(3):
        logging.info(f"Downloading assembly {acc} (attempt {attempt+1})")
        urllib.request.urlretrieve(url, gz_path)

        local_md5 = compute_md5(gz_path)

        if local_md5 == expected_md5:
            logging.info(f"MD5 verified {acc} successfully")
            break
        else:
            logging.warning(f"MD5 mismatch for {acc} — retrying download")
            os.remove(gz_path)
            time.sleep(2)
    else:
        raise ValueError(f"MD5 verification {acc} failed after 3 attempts")

    return gz_path, uid, doc

def download_nucleotide(acc, genome_dir):

    local_path = os.path.join(genome_dir, f"{acc}.fasta")

    if not genome_valid(local_path):
        logging.info(f"Downloading nucleotide {acc}")

        def _save_fasta(handle):
            content = handle.read()
            with open(local_path, "w") as out:
                out.write(content)

        entrez_call(
            Entrez.efetch, db="nuccore", id=acc,
            rettype="fasta", retmode="text",
            read_func=_save_fasta,
        )

    summary = entrez_call(Entrez.esummary, db="nuccore", id=acc)
    doc = summary[0]

    return local_path, acc, doc

def _parse_pcr_primers(raw: str) -> tuple:

    text = " ".join(raw.split())

    def _get(pattern):
        m = re.search(pattern, text)
        return m.group(1).strip() if m else "NA"

    fwd_name = _get(r'fwd_name:\s*([^,]+)')
    fwd_seq  = _get(r'fwd_seq:\s*([^,]+)')
    rev_name = _get(r'rev_name:\s*([^,]+)')
    rev_seq  = _get(r'rev_seq:\s*([^,]+)')

    return fwd_name, fwd_seq, rev_name, rev_seq


def extract_genbank_metadata(nuc_acc):

    country = collection_date = host = iso_source = "NA"
    seq_length = record_date = submitter_names = "NA"
    fwd_primer_name = fwd_primer_seq = rev_primer_name = rev_primer_seq = "NA"

    try:
        record = entrez_call(
            Entrez.efetch,
            db="nuccore", id=nuc_acc,
            rettype="gb", retmode="text",
            read_func=lambda h: SeqIO.read(h, "genbank"),
        )

        seq_length = len(record.seq)
        record_date = record.annotations.get("date", "NA")

        refs = record.annotations.get("references", [])
        if refs:
            submitter_names = ", ".join(refs[0].authors)

        for feature in record.features:
            if feature.type == "source":
                qualifiers = feature.qualifiers
                country = qualifiers.get(
                    "geo_loc_name",
                    qualifiers.get("country", ["NA"])
                )[0]
                collection_date = qualifiers.get(
                    "collection_date", ["NA"]
                )[0]
                host = qualifiers.get("host", ["NA"])[0]
                iso_source = qualifiers.get(
                    "isolation_source", ["NA"]
                )[0]

                pcr_list = qualifiers.get("PCR_primers", [])
                if pcr_list:
                    fwd_primer_name, fwd_primer_seq, \
                        rev_primer_name, rev_primer_seq = _parse_pcr_primers(pcr_list[0])

    except Exception as e:
        logging.warning(f"GenBank metadata failed for {nuc_acc}: {e}")

    return {
        "Geo Location": country,
        "Collection Date": collection_date,
        "Host": host,
        "Isolation Source": iso_source,
        "Forward Primer Name": fwd_primer_name,
        "Forward Primer Sequence": fwd_primer_seq,
        "Reverse Primer Name": rev_primer_name,
        "Reverse Primer Sequence": rev_primer_seq,
        "Sequence Length": seq_length,
        "Record Date": record_date,
        "Submitter (Authors)": submitter_names
    }

def download_genome(acc, args, genome_dir):
    start_time = time.time()
    try:
        logging.info(f"Downloading genome {acc}")

        if acc.startswith(("GCA_", "GCF_")):
            genome_path, uid, doc = download_assembly(acc, genome_dir)
        else:
            genome_path, uid, doc = download_nucleotide(acc, genome_dir)

        elapsed = time.time() - start_time
        logging.info(f"Genome download done for {acc}: {elapsed:.2f}s")
        return ("success", acc, genome_path, uid, doc)

    except Exception as e:
        logging.error(f"{acc} genome download failed: {e}")
        with open(os.path.join(args.outdir, "failed_download.txt"), "a") as f:
            f.write(f"{acc}\t{e}\n")
        return ("fail", acc, None, None, None)


def fetch_metadata(acc, args, genome_path, uid, doc,
                   retries: int = 4, base_delay: float = 5.0):

    for attempt in range(retries):
        start_time = time.time()
        try:
            if attempt:
                logging.info(f"Fetching metadata {acc} (retry {attempt}/{retries - 1})")
            else:
                logging.info(f"Fetching metadata {acc}")

            if acc.startswith(("GCA_", "GCF_")):

                if uid is None:
                    search_result = entrez_call(Entrez.esearch, db="assembly", term=acc)
                    time.sleep(0.2)
                    if not search_result["IdList"]:
                        raise ValueError(f"Assembly not found in NCBI: {acc}")
                    uid = search_result["IdList"][0]

                if doc is None:
                    summary = entrez_call(
                        Entrez.esummary, db="assembly", id=uid, report="full"
                    )
                    doc = summary["DocumentSummarySet"]["DocumentSummary"][0]
                    time.sleep(0.2)

                link_records = entrez_call(
                    Entrez.elink, dbfrom="assembly", db="nucleotide", id=uid
                )
                time.sleep(0.3)

                nuc_ids = []
                if link_records[0]["LinkSetDb"]:
                    nuc_ids = [
                        link["Id"]
                        for link in link_records[0]["LinkSetDb"][0]["Link"]
                    ]

                nucleotide_acc = "NA"
                if nuc_ids:
                    sum_result = entrez_call(
                        Entrez.esummary, db="nucleotide", id=nuc_ids[0]
                    )
                    time.sleep(0.3)
                    nucleotide_acc = sum_result[0]["AccessionVersion"]

                organism = doc.get("Organism", "NA")
                assembly_name = doc.get("AssemblyName", "NA")
                submitter = doc.get("SubmitterOrganization", "NA")
                date = doc.get("SubmissionDate", "NA")

            else:

                nucleotide_acc = uid if uid is not None else acc
                if doc is None:
                    doc_result = entrez_call(
                        Entrez.esummary, db="nuccore", id=nucleotide_acc
                    )
                    doc = doc_result[0]
                    time.sleep(0.2)
                organism = doc.get("Title", "NA")
                assembly_name = "NA"
                submitter = doc.get("Submitter", "NA")
                date = doc.get("CreateDate", "NA")

            gb_meta = extract_genbank_metadata(nucleotide_acc)

            if genome_path is not None:
                size, contigs, gc = compute_stats(genome_path)
            else:
                size, contigs, gc = "NA", "NA", "NA"

            new_row = {
                "Accession": acc,
                "Organism": organism,
                "Assembly Name": assembly_name,
                "Submitter": submitter,
                "Submission Date": date,
                **gb_meta,
                "Genome Size": size,
                "Contig Count": contigs,
                "GC (%)": gc,
                "Genome File": genome_path if genome_path is not None else "NA",
            }

            elapsed = time.time() - start_time
            logging.info(f"Metadata done for {acc}: {elapsed:.2f}s")
            time.sleep(0.5)
            return ("success", acc, new_row)

        except Exception as e:
            logging.error(
                f"{acc} metadata failed (attempt {attempt + 1}/{retries}): {e}"
            )
            if attempt < retries - 1:
                wait = base_delay * (2 ** attempt) + random.uniform(1.0, 3.0)
                logging.warning(f"Retrying {acc} metadata in {wait:.1f}s")
                time.sleep(wait)
            else:
                with open(os.path.join(args.outdir, "failed_metadata.txt"), "a") as f:
                    f.write(f"{acc}\t{e}\n")
                return ("fail", acc, None)

    return ("fail", acc, None)

def main():
    start_all = time.time()

    args = parse_args()
    Entrez.email = args.email
    Entrez.sleep_between_tries = 3
    Entrez.max_tries = 5

    os.makedirs(args.outdir, exist_ok=True)
    genome_dir = os.path.join(args.outdir, "genomes")
    if args.mode in ("both", "only_genome"):
        os.makedirs(genome_dir, exist_ok=True)

    setup_logging(args.outdir)

    excel_file = os.path.join(args.outdir, "metadata.xlsx")

    if os.path.exists(excel_file):
        df = pd.read_excel(excel_file)
        processed = set(df["Accession"])
        logging.info("Resume mode enabled")
    else:
        df = pd.DataFrame()
        processed = set()

    with open(args.input) as f:
        accessions = [x.strip() for x in f if x.strip()]

    logging.info(f"Mode: {args.mode} | Accessions: {len(accessions)}")

    downloaded = {}

    if args.mode in ("both", "only_genome"):
        logging.info("=== Downloading genomes ===")

        with ThreadPoolExecutor(max_workers=args.threads) as executor:
            futures = {
                executor.submit(download_genome, acc, args, genome_dir): acc
                for acc in accessions if acc not in processed
            }
            for future in as_completed(futures):
                status, acc, genome_path, uid, doc = future.result()
                if status == "success":
                    downloaded[acc] = (genome_path, uid, doc)

        logging.info(f"=== Finished: {len(downloaded)} genomes downloaded ===")

    elif args.mode == "only_metadata":
        for acc in accessions:
            if acc in processed:
                logging.info(f"Skipping {acc} (already in metadata.xlsx)")
                continue
            downloaded[acc] = (None, None, None)

        logging.info(f"=== {len(downloaded)} accessions queued for metadata fetch ===")

    if args.mode in ("both", "only_metadata"):
        logging.info("=== Fetching metadata ===")

        results = []
        lock = threading.Lock()

        def _fetch_and_collect(payload):
            acc, (genome_path, uid, doc) = payload
            status, acc, new_row = fetch_metadata(acc, args, genome_path, uid, doc)
            if status == "success" and new_row is not None:
                with lock:
                    results.append(new_row)

        with ThreadPoolExecutor(max_workers=args.threads) as executor:
            list(executor.map(_fetch_and_collect, downloaded.items()))

        if results:
            df = pd.concat([df, pd.DataFrame(results)], ignore_index=True)
            df.to_excel(excel_file, index=False)

    else:
        logging.info("Skipping (only_genome mode — no metadata fetched)")

    total_elapsed = time.time() - start_all
    logging.info(f"All done in {total_elapsed/60:.2f} minutes")


if __name__ == "__main__":
    main()
