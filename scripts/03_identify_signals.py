#!/usr/bin/env python3
"""
Step 3: Identification & Prioritization of Key Genetic Signals (Task 3)
Performs distance-based physical locus clumping (+/- 500 kb window on GRCh38)
and cross-dataset concordance analysis between GCST90502763 (BMI-adjusted) and
GCST90502762 (BMI-unadjusted).
"""

import os
import json
import numpy as np
import pandas as pd

def clump_loci(df_sig, window_bp=500_000):
    """
    Greedy physical clumping:
    1. Sort all candidate variants by p-value ascending.
    2. Take the most significant un-clumped variant as the Lead Variant of a new locus.
    3. Assign all remaining variants on the same chromosome within +/- window_bp to this locus.
    4. Repeat until all variants are assigned.
    """
    df_sorted = df_sig.sort_values('p_value', ascending=True).copy()
    assigned = set()
    loci = []
    
    for idx, row in df_sorted.iterrows():
        if idx in assigned:
            continue
        chrom = int(row['chromosome'])
        pos = int(row['base_pair_location'])
        
        # Find all variants within +/- window_bp on the same chromosome
        in_window = df_sorted[
            (df_sorted['chromosome'] == chrom) &
            (df_sorted['base_pair_location'] >= pos - window_bp) &
            (df_sorted['base_pair_location'] <= pos + window_bp) &
            (~df_sorted.index.isin(assigned))
        ]
        
        assigned.update(in_window.index.tolist())
        
        n_gw = int((in_window['p_value'] < 5e-8).sum())
        n_sugg = int((in_window['p_value'] < 1e-5).sum())
        
        # Risk-increasing orientation (if beta < 0, the other_allele is the risk-increasing allele)
        beta = float(row['beta'])
        se = float(row['standard_error'])
        eaf = float(row['effect_allele_frequency'])
        maf = float(row['maf'])
        
        if beta >= 0:
            risk_allele = row['effect_allele']
            non_risk_allele = row['other_allele']
            raf = eaf
            risk_beta = beta
        else:
            risk_allele = row['other_allele']
            non_risk_allele = row['effect_allele']
            raf = 1.0 - eaf
            risk_beta = -beta
            
        risk_or = float(np.exp(risk_beta))
        risk_or_ci_low = float(np.exp(risk_beta - 1.96 * se))
        risk_or_ci_high = float(np.exp(risk_beta + 1.96 * se))
        
        loci.append({
            "locus_rank": len(loci) + 1,
            "tier": "Tier 1 (Genome-Wide Significant, p < 5e-8)" if row['p_value'] < 5e-8 else "Tier 2 (Suggestive Significance, p < 1e-5)",
            "chromosome": chrom,
            "lead_bp_grch38": pos,
            "region_start_bp": int(in_window['base_pair_location'].min()),
            "region_end_bp": int(in_window['base_pair_location'].max()),
            "lead_rs_id": str(row['rs_id']),
            "effect_allele": str(row['effect_allele']),
            "other_allele": str(row['other_allele']),
            "effect_allele_frequency": round(eaf, 5),
            "maf": round(maf, 5),
            "beta": round(beta, 5),
            "standard_error": round(se, 5),
            "odds_ratio_effect_allele": round(float(row['odds_ratio']), 4),
            "ci_95_effect_allele": f"[{row['ci_lower']:.3f} - {row['ci_upper']:.3f}]",
            "risk_increasing_allele": f"{row['rs_id']}-{risk_allele}",
            "risk_allele_frequency": round(raf, 4),
            "risk_odds_ratio": round(risk_or, 4),
            "risk_or_95ci": f"[{risk_or_ci_low:.2f} - {risk_or_ci_high:.2f}]",
            "p_value": float(row['p_value']),
            "neg_log10_p": round(float(-np.log10(row['p_value'])), 3),
            "n_gw_sig_in_locus": n_gw,
            "n_suggestive_in_locus": n_sugg,
            "supporting_variants_rsids": in_window['rs_id'].tolist()[:15]
        })
        
    return loci

def main():
    df_unadj = pd.read_csv("results/GCST90502762_qc_p1e-3.tsv", sep='\t')
    df_bmi = pd.read_csv("results/GCST90502763_qc_p1e-3.tsv", sep='\t')
    
    # Filter to p < 1e-5 for suggestive + genome-wide clumping
    sugg_unadj = df_unadj[df_unadj['p_value'] < 1e-5].copy()
    sugg_bmi = df_bmi[df_bmi['p_value'] < 1e-5].copy()
    
    loci_bmi = clump_loci(sugg_bmi, window_bp=500_000)
    loci_unadj = clump_loci(sugg_unadj, window_bp=500_000)
    
    # Build lookup of all p < 1e-3 variants by (chromosome, base_pair_location) to cross-reference
    unadj_lookup = df_unadj.set_index(['chromosome', 'base_pair_location']).to_dict(orient='index')
    bmi_lookup = df_bmi.set_index(['chromosome', 'base_pair_location']).to_dict(orient='index')
    
    for loc in loci_bmi:
        key = (loc['chromosome'], loc['lead_bp_grch38'])
        if key in unadj_lookup:
            u = unadj_lookup[key]
            loc['unadjusted_p_value'] = float(u['p_value'])
            loc['unadjusted_beta'] = round(float(u['beta']), 5)
            loc['unadjusted_or'] = round(float(u['odds_ratio']), 4)
        else:
            loc['unadjusted_p_value'] = None
            loc['unadjusted_beta'] = None
            loc['unadjusted_or'] = None
            
    for loc in loci_unadj:
        key = (loc['chromosome'], loc['lead_bp_grch38'])
        if key in bmi_lookup:
            b = bmi_lookup[key]
            loc['bmi_adjusted_p_value'] = float(b['p_value'])
            loc['bmi_adjusted_beta'] = round(float(b['beta']), 5)
            loc['bmi_adjusted_or'] = round(float(b['odds_ratio']), 4)
        else:
            loc['bmi_adjusted_p_value'] = None
            loc['bmi_adjusted_beta'] = None
            loc['bmi_adjusted_or'] = None

    # Extract all genome-wide significant variants (p < 5e-8) in either study
    gw_bmi = df_bmi[df_bmi['p_value'] < 5e-8].copy()
    gw_unadj = df_unadj[df_unadj['p_value'] < 5e-8].copy()
    
    gw_merged = pd.merge(
        df_bmi[df_bmi['p_value'] < 5e-7],
        df_unadj[['chromosome', 'base_pair_location', 'beta', 'standard_error', 'p_value', 'odds_ratio', 'ci_lower', 'ci_upper']],
        on=['chromosome', 'base_pair_location'],
        how='outer',
        suffixes=('_bmi_adj', '_unadj')
    ).sort_values('p_value_bmi_adj')
    
    gw_merged.to_csv("results/task3_top_variants_comparison.tsv", sep='\t', index=False)
    pd.DataFrame(loci_bmi).to_csv("results/task3_clumped_loci_GCST90502763_bmi_adj.tsv", sep='\t', index=False)
    pd.DataFrame(loci_unadj).to_csv("results/task3_clumped_loci_GCST90502762_unadj.tsv", sep='\t', index=False)
    
    out = {
        "GCST90502763_bmi_adjusted_loci": loci_bmi,
        "GCST90502762_unadjusted_loci": loci_unadj
    }
    with open("results/task3_signals_summary.json", "w") as f:
        json.dump(out, f, indent=2)
        
    print(f"GCST90502763 (BMI-adjusted): {len(loci_bmi)} independent loci (p < 1e-5), of which {sum(1 for x in loci_bmi if x['p_value'] < 5e-8)} are genome-wide significant (p < 5e-8).")
    print(f"GCST90502762 (Unadjusted): {len(loci_unadj)} independent loci (p < 1e-5), of which {sum(1 for x in loci_unadj if x['p_value'] < 5e-8)} are genome-wide significant (p < 5e-8).")

if __name__ == "__main__":
    main()
