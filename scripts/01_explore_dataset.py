#!/usr/bin/env python3
"""
Step 1: Dataset Exploration (Task 1)
Streams or reads the GWAS summary statistics datasets (GCST90502762 and GCST90502763),
computes structural metadata, missing value counts, chromosome coverage, numerical distributions,
and saves a local compressed copy and summary report.
"""

import os
import sys
import urllib.request
import pandas as pd
import numpy as np
import json

DATASETS = {
    "GCST90502762": {
        "name": "Prediabetes (Unadjusted for BMI)",
        "url": "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90502001-GCST90503000/GCST90502762/GCST90502762.tsv",
        "local_tsv": "data/GCST90502762.tsv"
    },
    "GCST90502763": {
        "name": "Prediabetes (BMI-adjusted)",
        "url": "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90502001-GCST90503000/GCST90502763/GCST90502763.tsv",
        "local_tsv": "data/GCST90502763.tsv"
    }
}

def explore_dataset(gcst_id, info):
    filepath = info["local_tsv"]
    if not os.path.exists(filepath):
        print(f"File {filepath} not found locally. Please download first.")
        return None

    print(f"Exploring {gcst_id} ({info['name']})...")
    
    # Read in chunks to keep memory usage low and fast
    total_rows = 0
    missing_counts = None
    chrom_counts = {}
    
    # Track min/max and quantiles via sampling or online stats
    beta_min, beta_max = float('inf'), float('-inf')
    se_min, se_max = float('inf'), float('-inf')
    eaf_min, eaf_max = float('inf'), float('-inf')
    p_min, p_max = float('inf'), float('-inf')
    
    invalid_p = 0
    invalid_se = 0
    invalid_eaf = 0
    missing_rsid = 0
    chr_pos_id_count = 0
    palindromic_count = 0
    indel_count = 0
    snv_count = 0
    maf_lt_001 = 0
    gw_sig_count = 0
    sugg_sig_count = 0

    dtypes = {
        'chromosome': str,
        'base_pair_location': 'Int64',
        'effect_allele': str,
        'other_allele': str,
        'beta': float,
        'standard_error': float,
        'effect_allele_frequency': float,
        'p_value': float,
        'rs_id': str
    }

    chunk_size = 1_000_000
    for chunk in pd.read_csv(filepath, sep='\t', dtype=dtypes, chunksize=chunk_size):
        total_rows += len(chunk)
        
        if missing_counts is None:
            missing_counts = chunk.isna().sum()
        else:
            missing_counts += chunk.isna().sum()
            
        # Chromosome counts
        vc = chunk['chromosome'].value_counts()
        for k, v in vc.items():
            chrom_counts[str(k)] = chrom_counts.get(str(k), 0) + int(v)
            
        # Ranges
        b_min, b_max = chunk['beta'].min(), chunk['beta'].max()
        s_min, s_max = chunk['standard_error'].min(), chunk['standard_error'].max()
        e_min, e_max = chunk['effect_allele_frequency'].min(), chunk['effect_allele_frequency'].max()
        pv_min, pv_max = chunk['p_value'].min(), chunk['p_value'].max()
        
        if pd.notna(b_min) and b_min < beta_min: beta_min = float(b_min)
        if pd.notna(b_max) and b_max > beta_max: beta_max = float(b_max)
        if pd.notna(s_min) and s_min < se_min: se_min = float(s_min)
        if pd.notna(s_max) and s_max > se_max: se_max = float(s_max)
        if pd.notna(e_min) and e_min < eaf_min: eaf_min = float(e_min)
        if pd.notna(e_max) and e_max > eaf_max: eaf_max = float(e_max)
        if pd.notna(pv_min) and pv_min < p_min: p_min = float(pv_min)
        if pd.notna(pv_max) and pv_max > p_max: p_max = float(pv_max)
        
        # Anomalies
        invalid_p += int(((chunk['p_value'] <= 0) | (chunk['p_value'] > 1)).sum())
        invalid_se += int(((chunk['standard_error'] <= 0) | (chunk['standard_error'] == np.inf)).sum())
        invalid_eaf += int(((chunk['effect_allele_frequency'] < 0) | (chunk['effect_allele_frequency'] > 1)).sum())
        
        # rs_id format check
        rs_na = chunk['rs_id'].isna()
        missing_rsid += int(rs_na.sum())
        chr_pos_id_count += int((~rs_na & ~chunk['rs_id'].str.startswith('rs', na=False)).sum())
        
        # Allele types (SNV vs Indel, Palindromic)
        ea = chunk['effect_allele'].fillna('')
        oa = chunk['other_allele'].fillna('')
        is_snv = (ea.str.len() == 1) & (oa.str.len() == 1)
        snv_count += int(is_snv.sum())
        indel_count += int((~is_snv).sum())
        
        # Palindromic (A/T, T/A, C/G, G/C)
        pairs = ea + oa
        palin = pairs.isin(['AT', 'TA', 'CG', 'GC'])
        palindromic_count += int(palin.sum())
        
        # MAF < 0.01
        eaf = chunk['effect_allele_frequency']
        maf = np.minimum(eaf, 1.0 - eaf)
        maf_lt_001 += int((maf < 0.01).sum())
        
        # Significance counts (raw)
        gw_sig_count += int((chunk['p_value'] < 5e-8).sum())
        sugg_sig_count += int((chunk['p_value'] < 1e-5).sum())

    summary = {
        "gcst_id": gcst_id,
        "trait_name": info["name"],
        "total_variants": total_rows,
        "columns": list(dtypes.keys()),
        "missing_values": {k: int(v) for k, v in missing_counts.items()},
        "chromosomes": chrom_counts,
        "ranges": {
            "beta": [beta_min, beta_max],
            "standard_error": [se_min, se_max],
            "effect_allele_frequency": [eaf_min, eaf_max],
            "p_value": [p_min, p_max]
        },
        "variant_types": {
            "biallelic_snv_length_1": snv_count,
            "indels_or_multinucleotide": indel_count,
            "palindromic_snps_AT_CG": palindromic_count,
            "non_rsid_identifiers_chr_pos": chr_pos_id_count,
            "maf_lt_0.01": maf_lt_001
        },
        "anomalies": {
            "invalid_p_value": invalid_p,
            "invalid_se": invalid_se,
            "invalid_eaf": invalid_eaf
        },
        "raw_significance_counts": {
            "genome_wide_p_lt_5e-8": gw_sig_count,
            "suggestive_p_lt_1e-5": sugg_sig_count
        }
    }
    return summary

if __name__ == "__main__":
    results = {}
    for gcst_id, info in DATASETS.items():
        res = explore_dataset(gcst_id, info)
        if res:
            results[gcst_id] = res
            
    os.makedirs("results", exist_ok=True)
    with open("results/task1_exploration_summary.json", "w") as f:
        json.dump(results, f, indent=2)
    print("Saved exploration summary to results/task1_exploration_summary.json")
