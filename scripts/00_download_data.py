#!/usr/bin/env python3
"""
Fast Parallel Range Downloader for EBI FTP GWAS Summary Statistics.
Downloads large TSV files using 16 parallel HTTP byte-range connections.
"""

import os
import sys
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor

def download_chunk(url, start, end, part_file):
    expected_len = end - start + 1
    for attempt in range(5):
        try:
            req = urllib.request.Request(url)
            req.add_header('Range', f'bytes={start}-{end}')
            bytes_read = 0
            with urllib.request.urlopen(req, timeout=60) as resp, open(part_file, 'wb') as out:
                while True:
                    buf = resp.read(1024 * 1024)
                    if not buf:
                        break
                    out.write(buf)
                    bytes_read += len(buf)
            if bytes_read == expected_len:
                return True
            else:
                print(f"Short read ({bytes_read}/{expected_len}) for range {start}-{end}, retrying...")
        except Exception as e:
            print(f"Retry {attempt+1} for range {start}-{end}: {e}")
            time.sleep(2)
    return False

def parallel_download(url, output_path, num_threads=16):
    t0 = time.time()
    req = urllib.request.Request(url, method='HEAD')
    with urllib.request.urlopen(req) as resp:
        total_size = int(resp.headers['Content-Length'])
        accept_ranges = resp.headers.get('Accept-Ranges', 'none')
    
    print(f"Downloading {url} ({total_size / (1024*1024):.2f} MB), Accept-Ranges: {accept_ranges}")
    chunk_size = total_size // num_threads
    ranges = []
    part_files = []
    for i in range(num_threads):
        start = i * chunk_size
        end = total_size - 1 if i == num_threads - 1 else (i + 1) * chunk_size - 1
        part_file = f"{output_path}.part{i}"
        ranges.append((url, start, end, part_file))
        part_files.append(part_file)
        
    with ThreadPoolExecutor(max_workers=num_threads) as ex:
        results = list(ex.map(lambda r: download_chunk(*r), ranges))
        
    if not all(results):
        raise RuntimeError("One or more chunks failed to download.")
        
    print(f"Merging {num_threads} parts into {output_path}...")
    with open(output_path, 'wb') as outfile:
        for pf in part_files:
            with open(pf, 'rb') as infile:
                while True:
                    buf = infile.read(4 * 1024 * 1024)
                    if not buf:
                        break
                    outfile.write(buf)
            os.remove(pf)
            
    dt = time.time() - t0
    print(f"Completed {output_path} in {dt:.1f}s ({total_size / (1024*1024*dt):.2f} MB/s)")

if __name__ == "__main__":
    import hashlib
    os.makedirs("data", exist_ok=True)
    urls = [
        ("http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90502001-GCST90503000/GCST90502762/GCST90502762.tsv", "data/GCST90502762.tsv", "87185b7467ff79720da933a72fa2bf3e"),
        ("http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90502001-GCST90503000/GCST90502763/GCST90502763.tsv", "data/GCST90502763.tsv", "4e72255db1d27b3a4592a10c9430f11b")
    ]
    for url, out, expected_md5 in urls:
        if os.path.exists(out):
            h = hashlib.md5()
            with open(out, 'rb') as f:
                while chunk := f.read(8*1024*1024):
                    h.update(chunk)
            if h.hexdigest().lower() == expected_md5.lower():
                print(f"{out} already exists and MD5 matches ({expected_md5}). Skipping.")
                continue
        parallel_download(url, out, num_threads=16)
