#!/usr/bin/env python3
"""
Step 2: Quality Assessment & Preprocessing (Task 2)
Evaluates data quality, applies rigorous and justified QC filters to both GWAS datasets
(GCST90502762 and GCST90502763), calculates Genomic Inflation Factors (lambda_GC) across
multiple MAF strata, verifies Wald vs SPA p-value consistency, and saves QC summary tables
plus filtered candidate subsets for fast downstream clumping and visualization.
"""

import os
import json
import numpy as np
import pandas as pd
from scipy.stats import chi2

DATASETS = {
    "GCST90502762": {
        "name": "Prediabetes (Unadjusted for BMI)",
        "local_tsv": "data/GCST90502762.tsv",
        "sig_out": "results/GCST90502762_qc_p1e-3.tsv"
    },
    "GCST90502763": {
        "name": "Prediabetes (BMI-adjusted)",
        "local_tsv": "data/GCST90502763.tsv",
        "sig_out": "results/GCST90502763_qc_p1e-3.tsv"
    }
}

def calc_lambda_gc(pvals):
    """Calculate genomic inflation factor lambda_GC from median p-value."""
    valid_p = pvals[(pvals > 0) & (pvals <= 1) & np.isfinite(pvals)]
    med_p = np.median(valid_p)
    obs_median_chi2 = chi2.isf(med_p, df=1)
    exp_median_chi2 = chi2.ppf(0.5, df=1)  # 0.454936423119572
    return float(obs_median_chi2 / exp_median_chi2), float(med_p)

def run_qc(gcst_id, info):
    print(f"\n--- Running Quality Assessment & QC on {gcst_id} ({info['name']}) ---")
    filepath = info["local_tsv"]
    
    dtypes = {
        'chromosome': 'int8',
        'base_pair_location': 'int32',
        'effect_allele': str,
        'other_allele': str,
        'beta': 'float64',
        'standard_error': 'float64',
        'effect_allele_frequency': 'float64',
        'p_value': 'float64',
        'rs_id': str
    }
    
    df = pd.read_csv(filepath, sep='\t', dtype=dtypes)
    n_raw = len(df)
    print(f"Loaded {n_raw:,} total variants.")
    
    # Compute derived columns
    df['maf'] = np.minimum(df['effect_allele_frequency'], 1.0 - df['effect_allele_frequency'])
    
    # Pre-QC lambda_GC
    lambda_raw, med_p_raw = calc_lambda_gc(df['p_value'].values)
    
    # Step-by-step QC flags
    # 1. Missing or invalid basic values (chr in 1..22, pos > 0, 0 < p <= 1, SE > 0, 0 < EAF < 1)
    valid_basic = (
        df['chromosome'].between(1, 22) &
        (df['base_pair_location'] > 0) &
        (df['p_value'] > 0) & (df['p_value'] <= 1) &
        (df['standard_error'] > 0) &
        (df['effect_allele_frequency'] > 0) & (df['effect_allele_frequency'] < 1) &
        df['beta'].notna() & df['rs_id'].notna()
    )
    n_fail_basic = int((~valid_basic).sum())
    
    # 2. Valid ACGT nucleotides only
    valid_alleles = (
        df['effect_allele'].str.fullmatch(r'[ACGT]+', na=False) &
        df['other_allele'].str.fullmatch(r'[ACGT]+', na=False) &
        (df['effect_allele'] != df['other_allele'])
    )
    n_fail_alleles = int((valid_basic & ~valid_alleles).sum())
    
    # 3. Extreme SAIGE numerical non-convergence outliers (|beta| >= 10 or SE >= 10)
    valid_convergence = (df['beta'].abs() < 10.0) & (df['standard_error'] < 10.0)
    n_fail_convergence = int((valid_basic & valid_alleles & ~valid_convergence).sum())
    
    # Inspect the extreme outliers before removing
    outliers_df = df[~valid_convergence]
    outlier_summary = {
        "count": len(outliers_df),
        "se_gt_1e6_count": int((outliers_df['standard_error'] > 1e6).sum()),
        "abs_beta_between_10_and_100": int(((outliers_df['beta'].abs() >= 10) & (outliers_df['standard_error'] <= 1e6)).sum()),
        "min_p_among_outliers": float(outliers_df['p_value'].min()) if len(outliers_df) > 0 else None
    }
    
    # 4. Duplicate variant coordinates+alleles or multi-allelic sites
    # Check exact duplicate (chr, pos, effect_allele, other_allele)
    is_dup_exact = df.duplicated(subset=['chromosome', 'base_pair_location', 'effect_allele', 'other_allele'], keep='first')
    n_dup_exact = int(is_dup_exact.sum())
    
    # Check multi-allelic same coordinate (chr, pos)
    is_multiallelic_pos = df.duplicated(subset=['chromosome', 'base_pair_location'], keep=False)
    n_multiallelic_rows = int(is_multiallelic_pos.sum())
    
    # Check duplicate rs_id (excluding chr:pos strings)
    is_rs = df['rs_id'].str.startswith('rs', na=False)
    is_dup_rsid = is_rs & df.duplicated(subset=['rs_id'], keep=False)
    n_dup_rsid = int(is_dup_rsid.sum())
    
    # 5. MAF check: how many variants have MAF < 0.01 or MAF == 0.01?
    n_maf_lt_001 = int((df['maf'] < 0.01).sum())
    n_maf_lt_005 = int((df['maf'] < 0.05).sum())
    
    # 6. Ambiguous palindromic SNPs (A/T or C/G with intermediate MAF > 0.40)
    pairs = df['effect_allele'] + df['other_allele']
    is_palindromic = pairs.isin(['AT', 'TA', 'CG', 'GC'])
    is_ambig_palindromic = is_palindromic & (df['maf'] > 0.40)
    n_palindromic = int(is_palindromic.sum())
    n_ambig_palindromic = int(is_ambig_palindromic.sum())
    
    # Primary QC filter (retaining clean, converged, non-duplicate variants with MAF >= 0.01)
    qc_pass = valid_basic & valid_alleles & valid_convergence & (~is_dup_exact) & (df['maf'] >= 0.01)
    df_qc = df[qc_pass].copy()
    n_qc_pass = len(df_qc)
    print(f"QC-passed variants: {n_qc_pass:,} (removed {n_raw - n_qc_pass:,})")
    
    # Strict Harmonisation-Ready QC count (also excluding ambiguous palindromic MAF > 0.40 and unmapped chr:pos IDs)
    qc_strict = qc_pass & (~is_ambig_palindromic) & is_rs
    n_qc_strict = int(qc_strict.sum())
    
    # Post-QC Genomic Inflation Factors (lambda_GC)
    lambda_qc_all, med_p_qc = calc_lambda_gc(df_qc['p_value'].values)
    lambda_common, _ = calc_lambda_gc(df_qc.loc[df_qc['maf'] >= 0.05, 'p_value'].values)
    lambda_lowfreq, _ = calc_lambda_gc(df_qc.loc[df_qc['maf'] < 0.05, 'p_value'].values)
    lambda_strict, _ = calc_lambda_gc(df.loc[qc_strict, 'p_value'].values)
    
    # Check consistency between Wald Z = beta / SE and reported SAIGE SPA p_value
    z_wald = df_qc['beta'] / df_qc['standard_error']
    p_wald = chi2.sf(z_wald**2, df=1)
    log10_p_rep = -np.log10(df_qc['p_value'])
    log10_p_wald = -np.log10(np.maximum(p_wald, 1e-300))
    corr_log10p = float(np.corrcoef(log10_p_rep, log10_p_wald)[0, 1])
    max_diff_gw = float(np.abs(log10_p_rep[df_qc['p_value'] < 1e-5] - log10_p_wald[df_qc['p_value'] < 1e-5]).max())
    
    # Save candidate subset (p < 1e-3) for fast clumping and fine analysis
    df_cand = df_qc[df_qc['p_value'] < 1e-3].copy()
    df_cand['odds_ratio'] = np.exp(df_cand['beta'])
    df_cand['ci_lower'] = np.exp(df_cand['beta'] - 1.96 * df_cand['standard_error'])
    df_cand['ci_upper'] = np.exp(df_cand['beta'] + 1.96 * df_cand['standard_error'])
    df_cand.sort_values('p_value', inplace=True)
    df_cand.to_csv(info["sig_out"], sep='\t', index=False)
    print(f"Saved {len(df_cand):,} variants with p < 1e-3 to {info['sig_out']}")
    
    # Also save sorted p-value quantiles for Q-Q plot generation in Task 4
    p_sorted = np.sort(df_qc['p_value'].values)
    n_total = len(p_sorted)
    # Keep top 50,000 smallest p-values + 50,000 evenly spaced across the rest for exact Q-Q plotting
    idx_top = np.arange(min(50000, n_total))
    idx_rest = np.linspace(50000, n_total - 1, 50000, dtype=int) if n_total > 50000 else np.array([], dtype=int)
    idx_keep = np.unique(np.concatenate([idx_top, idx_rest]))
    obs_log10p = -np.log10(p_sorted[idx_keep])
    exp_log10p = -np.log10((idx_keep + 0.5) / n_total)
    np.savez_compressed(
        f"results/{gcst_id}_qq_data.npz",
        obs_log10p=obs_log10p,
        exp_log10p=exp_log10p,
        n_total=n_total,
        lambda_gc=lambda_qc_all
    )
    
    summary = {
        "gcst_id": gcst_id,
        "trait_name": info["name"],
        "counts": {
            "raw_variants": n_raw,
            "failed_basic_range_or_missing": n_fail_basic,
            "failed_acgt_allele_check": n_fail_alleles,
            "failed_saige_convergence_beta_se_lt_10": n_fail_convergence,
            "exact_duplicate_variants": n_dup_exact,
            "maf_lt_0.01": n_maf_lt_001,
            "qc_passed_primary": n_qc_pass,
            "qc_removed_primary": n_raw - n_qc_pass,
            "multiallelic_coordinate_rows_flagged": n_multiallelic_rows,
            "duplicate_rsid_rows_flagged": n_dup_rsid,
            "palindromic_snps_total": n_palindromic,
            "ambiguous_palindromic_maf_gt_0.40_flagged": n_ambig_palindromic,
            "low_frequency_0.01_to_0.05": n_maf_lt_005,
            "common_variants_maf_ge_0.05": n_qc_pass - n_maf_lt_005,
            "qc_passed_strict_harmonisation_ready": n_qc_strict
        },
        "outlier_diagnostics": outlier_summary,
        "genomic_inflation": {
            "lambda_gc_raw": round(lambda_raw, 5),
            "lambda_gc_qc_all": round(lambda_qc_all, 5),
            "lambda_gc_common_maf_ge_0.05": round(lambda_common, 5),
            "lambda_gc_lowfreq_maf_lt_0.05": round(lambda_lowfreq, 5),
            "lambda_gc_strict_no_ambig_palindromic": round(lambda_strict, 5),
            "median_p_value_qc": round(med_p_qc, 6)
        },
        "statistical_consistency": {
            "pearson_r_reported_vs_wald_neglog10p": round(corr_log10p, 6),
            "max_abs_diff_neglog10p_at_p_lt_1e-5": round(max_diff_gw, 4)
        },
        "post_qc_significance_counts": {
            "genome_wide_p_lt_5e-8": int((df_qc['p_value'] < 5e-8).sum()),
            "suggestive_p_lt_1e-5": int((df_qc['p_value'] < 1e-5).sum()),
            "nominal_p_lt_1e-3": len(df_cand)
        }
    }
    return summary

if __name__ == "__main__":
    os.makedirs("results", exist_ok=True)
    qc_results = {}
    for gcst_id, info in DATASETS.items():
        qc_results[gcst_id] = run_qc(gcst_id, info)
        
    with open("results/task2_qc_summary.json", "w") as f:
        json.dump(qc_results, f, indent=2)
    print("\nSaved Task 2 QC summary to results/task2_qc_summary.json")
