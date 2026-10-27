#!/usr/bin/env python3
"""
Step 6: Proposed Next Steps & Statistical Power Calculations (Task 6)
Computes exact sample sizes required for replication and multi-ancestry meta-analysis
at genome-wide significance (alpha = 5e-8) and Bonferroni-corrected replication
significance (alpha = 0.05 / 7 = 7.14e-3) for our top 7 prioritized prediabetes loci,
and outputs a structured roadmap to results/task6_proposed_next_steps.json.
"""

import json
import numpy as np
from scipy.stats import norm

PRIORITIZED_LOCI = [
    {
        "rank": 1,
        "gene": "GCK",
        "rs_id": "rs2908293",
        "chr_pos": "7:44169723",
        "raf": 0.2530,
        "maf": 0.2530,
        "beta_risk": 0.2235,
        "or_risk": 1.2504,
        "se": 0.0387,
        "p_bmi_adj": 7.68e-9,
        "tier": "Tier 1 (Genome-Wide Significant)",
        "axis": "Axis 1: Pancreatic Beta-Cell Glucose Sensing"
    },
    {
        "rank": 2,
        "gene": "HK1",
        "rs_id": "rs72805692",
        "chr_pos": "10:69339353",
        "raf": 0.9272,
        "maf": 0.0728,
        "beta_risk": 0.3769,
        "or_risk": 1.4578,
        "se": 0.0667,
        "p_bmi_adj": 1.59e-8,
        "tier": "Tier 1 (Genome-Wide Significant)",
        "axis": "Axis 2: Erythrocyte Glycolytic Flux & HbA1c Kinetics"
    },
    {
        "rank": 3,
        "gene": "PDE1C",
        "rs_id": "rs917778362",
        "chr_pos": "7:31805776",
        "raf": 0.0421,
        "maf": 0.0421,
        "beta_risk": 0.4787,
        "or_risk": 1.6140,
        "se": 0.0954,
        "p_bmi_adj": 5.53e-8,
        "tier": "Tier 2 (Borderline GW-Sig, Novel Qatari-Enriched)",
        "axis": "Axis 1: Beta-Cell cAMP Hydrolysis & Incretin Response"
    },
    {
        "rank": 4,
        "gene": "PLXNA4",
        "rs_id": "rs62468156",
        "chr_pos": "7:131987551",
        "raf": 0.2597,
        "maf": 0.2597,
        "beta_risk": 0.2039,
        "or_risk": 1.2262,
        "se": 0.0389,
        "p_bmi_adj": 1.62e-7,
        "tier": "Tier 2 (Suggestive)",
        "axis": "Candidate Neuro-Islet / Guidance Signaling"
    },
    {
        "rank": 5,
        "gene": "MTNR1B",
        "rs_id": "rs10830963",
        "chr_pos": "11:92975544",
        "raf": 0.2784,
        "maf": 0.2784,
        "beta_risk": 0.1965,
        "or_risk": 1.2171,
        "se": 0.0385,
        "p_bmi_adj": 3.43e-7,
        "tier": "Tier 2 (Canonical Glycemic Control)",
        "axis": "Axis 1: Beta-Cell Gi-Coupled cAMP Inhibition"
    },
    {
        "rank": 6,
        "gene": "G6PC2 / ABCB11",
        "rs_id": "rs853777",
        "chr_pos": "2:168955707",
        "raf": 0.7977,
        "maf": 0.2023,
        "beta_risk": 0.2188,
        "or_risk": 1.2446,
        "se": 0.0433,
        "p_bmi_adj": 4.57e-7,
        "tier": "Tier 2 (Canonical Glycemic Control)",
        "axis": "Axis 1: Beta-Cell Glucose-6-Phosphate Futile Cycle"
    },
    {
        "rank": 7,
        "gene": "FN3KRP / B3GNTL1 / METRNL",
        "rs_id": "rs12947062",
        "chr_pos": "17:82728494",
        "raf": 0.3074,
        "maf": 0.3074,
        "beta_risk": 0.1847,
        "or_risk": 1.2029,
        "se": 0.0374,
        "p_bmi_adj": 8.04e-7,
        "tier": "Tier 2 (Suggestive)",
        "axis": "Axis 2: Erythrocyte Fructosamine Deglycation (FN3KRP) & Insulin Sensitivity (METRNL)"
    }
]

def compute_power_and_sample_sizes():
    n_discovery = 10932
    phi = 2267 / 10932  # 0.207373

    z_gw = norm.ppf(1 - 5e-8 / 2)        # 5.4513
    z_rep = norm.ppf(1 - (0.05 / 7) / 2) # 2.6901
    z_80 = norm.ppf(0.80)                # 0.8416

    results = []
    for loc in PRIORITIZED_LOCI:
        v_emp = (loc["se"] ** 2) * n_discovery
        beta_disc = loc["beta_risk"]
        beta_wc = 0.85 * beta_disc

        n_gw80_disc = int(np.ceil(v_emp * ((z_gw + z_80) / beta_disc) ** 2))
        n_gw80_wc = int(np.ceil(v_emp * ((z_gw + z_80) / beta_wc) ** 2))
        n_rep80_wc = int(np.ceil(v_emp * ((z_rep + z_80) / beta_wc) ** 2))
        ncp_current = beta_disc / loc["se"]
        power_current_gw = float(1 - norm.cdf(z_gw - ncp_current))

        loc_res = dict(loc)
        loc_res.update({
            "current_wald_z": round(ncp_current, 2),
            "current_power_at_5e-8": round(power_current_gw, 3),
            "n_total_for_80pct_gw_sig_unshrunk": n_gw80_disc,
            "n_total_for_80pct_gw_sig_wc_shrunk": n_gw80_wc,
            "additional_n_beyond_qbb_for_gw_sig": max(0, n_gw80_disc - n_discovery),
            "independent_replication_n_for_80pct_bonferroni": n_rep80_wc,
            "independent_replication_cases_at_21pct_prevalence": int(np.ceil(n_rep80_wc * phi))
        })
        results.append(loc_res)

    out = {
        "discovery_cohort": {
            "name": "Qatar Biobank (QBB)",
            "n_total": n_discovery,
            "n_cases": 2267,
            "n_controls": 8665,
            "case_fraction": round(phi, 4),
            "effective_sample_size_neff": round(4 * phi * (1 - phi) * n_discovery, 1)
        },
        "power_and_replication_projections": results
    }

    with open("results/task6_proposed_next_steps.json", "w") as f:
        json.dump(out, f, indent=2)
    print("Saved results/task6_proposed_next_steps.json")

if __name__ == "__main__":
    compute_power_and_sample_sizes()
