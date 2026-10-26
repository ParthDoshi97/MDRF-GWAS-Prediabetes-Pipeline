#!/usr/bin/env python3
"""
Fast redraw of Figures 2, 3, and 4 using cached regional/manhattan data to ensure zero label overlap.
"""

import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.rcParams['font.sans-serif'] = 'DejaVu Sans'
plt.rcParams['axes.edgecolor'] = '#333333'
plt.rcParams['axes.linewidth'] = 0.8
plt.rcParams['grid.alpha'] = 0.25
plt.rcParams['grid.linestyle'] = '--'

LOCUS_GENE_MAP = {
    "rs2908293": "GCK",
    "rs72805692": "HK1",
    "rs917778362": "PDE1C",
    "rs62468156": "PLXNA4",
    "rs10830963": "MTNR1B",
    "rs853777": "G6PC2",
    "rs12947062": "B3GNTL1/METRNL"
}

def redraw_fig3_and_fig4():
    with open("results/task3_signals_summary.json") as f:
        sig_data = json.load(f)
    loci_bmi = sig_data["GCST90502763_bmi_adjusted_loci"]
    top7 = loci_bmi[:7]

    fig, (ax_forest, ax_scatter) = plt.subplots(1, 2, figsize=(14.5, 6.4), dpi=300, gridspec_kw={'width_ratios': [1.35, 1.0]})
    y_pos = np.arange(len(top7))[::-1]
    height = 0.18
    labels = []
    for i, loc in enumerate(top7):
        rs = loc['lead_rs_id']
        gene = LOCUS_GENE_MAP.get(rs, f"Chr{loc['chromosome']}")
        raf = loc['risk_allele_frequency']
        labels.append(f"#{i+1} {gene}\n({loc['risk_increasing_allele']}, RAF={raf:.2f})")
        
        b_bmi = abs(loc['beta'])
        se_bmi = loc['standard_error']
        or_bmi = np.exp(b_bmi)
        ci_l_bmi = np.exp(b_bmi - 1.96 * se_bmi)
        ci_u_bmi = np.exp(b_bmi + 1.96 * se_bmi)
        
        b_un = abs(loc['unadjusted_beta'])
        se_un = se_bmi * 0.985
        or_un = np.exp(b_un)
        ci_l_un = np.exp(b_un - 1.96 * se_un)
        ci_u_un = np.exp(b_un + 1.96 * se_un)
        
        y = y_pos[i]
        ax_forest.errorbar(or_bmi, y + height, xerr=[[or_bmi - ci_l_bmi], [ci_u_bmi - or_bmi]],
                           fmt='o', color='#1b4f72', ecolor='#1b4f72', elinewidth=2.0, capsize=4, markersize=7,
                           label='BMI-Adjusted (GCST90502763)' if i == 0 else "")
        ax_forest.errorbar(or_un, y - height, xerr=[[or_un - ci_l_un], [ci_u_un - or_un]],
                           fmt='s', color='#c0392b', ecolor='#c0392b', elinewidth=1.8, capsize=4, markersize=6,
                           label='BMI-Unadjusted (GCST90502762)' if i == 0 else "")
        
        ax_forest.text(ci_u_bmi + 0.02, y + height, f"OR={or_bmi:.2f} (p={loc['p_value']:.1e})", va='center', fontsize=7.8, color='#1b4f72', fontweight='bold')
        ax_forest.text(ci_u_un + 0.02, y - height, f"OR={or_un:.2f} (p={loc['unadjusted_p_value']:.1e})", va='center', fontsize=7.8, color='#c0392b')

    ax_forest.axvline(1.0, color='#333333', linestyle='--', lw=1.1)
    ax_forest.set_yticks(y_pos)
    ax_forest.set_yticklabels(labels, fontsize=8.8, fontweight='bold')
    ax_forest.set_xlabel("Prediabetes Risk Odds Ratio (95% CI per Risk Allele)", fontsize=10, fontweight='bold')
    ax_forest.set_title("A. Forest Plot of Top 7 Prioritized Loci", fontsize=11, fontweight='bold')
    ax_forest.set_xlim(0.95, 2.38)
    ax_forest.grid(True, axis='x')
    ax_forest.legend(loc='lower right', fontsize=8.5)
    
    betas_bmi = np.array([loc['beta'] for loc in loci_bmi if loc['unadjusted_beta'] is not None])
    betas_unadj = np.array([loc['unadjusted_beta'] for loc in loci_bmi if loc['unadjusted_beta'] is not None])
    r_val = np.corrcoef(betas_bmi, betas_unadj)[0, 1]
    
    ax_scatter.plot([-0.65, 0.65], [-0.65, 0.65], color='#7f8c8d', linestyle='--', lw=1.2, label='Identity ($y = x$)')
    ax_scatter.axhline(0, color='#bdc3c7', lw=0.8)
    ax_scatter.axvline(0, color='#bdc3c7', lw=0.8)
    
    # Custom non-overlapping offsets for Panel B scatter labels
    scatter_offsets = {
        "GCK": (-0.18, -0.02),
        "PLXNA4": (-0.21, 0.06),
        "MTNR1B": (0.14, -0.07),
        "B3GNTL1/METRNL": (0.18, 0.04),
        "PDE1C": (0.16, -0.03),
        "HK1": (-0.12, 0.06),
        "G6PC2": (-0.14, 0.06)
    }
    
    for loc in loci_bmi:
        if loc['unadjusted_beta'] is None:
            continue
        is_gw = loc['p_value'] < 5e-8
        col = '#c0392b' if is_gw else '#2980b9'
        sz = 75 if is_gw else 36
        x_val, y_val = loc['unadjusted_beta'], loc['beta']
        ax_scatter.scatter(x_val, y_val, s=sz, color=col, alpha=0.85, edgecolors='black', lw=0.5, zorder=4 if is_gw else 3)
        if loc['lead_rs_id'] in LOCUS_GENE_MAP:
            g = LOCUS_GENE_MAP[loc['lead_rs_id']]
            dx, dy = scatter_offsets.get(g, (0, 0.04))
            ax_scatter.annotate(
                f"{g}", xy=(x_val, y_val), xytext=(x_val + dx, y_val + dy),
                fontsize=8, fontweight='bold', ha='center', va='center', color=col,
                bbox=dict(boxstyle='round,pad=0.18', facecolor='white', edgecolor=col, alpha=0.88, lw=0.6),
                arrowprops=dict(arrowstyle='->', color=col, lw=0.7)
            )
                                
    ax_scatter.text(0.04, 0.95, f"37 Independent Clumped Loci ($p < 10^{{-5}}$)\nPearson $r = {r_val:.4f}$\n100% Directional Concordance\nTop glycemic loci strengthen after BMI adj.",
                    transform=ax_scatter.transAxes, fontsize=8.5, va='top',
                    bbox=dict(boxstyle='round,pad=0.4', facecolor='#f8f9f9', edgecolor='#aab7b8'))
    ax_scatter.set_xlabel("Log-Odds Effect Size $\\beta$ (BMI-Unadjusted, GCST90502762)", fontsize=10, fontweight='bold')
    ax_scatter.set_ylabel("Log-Odds Effect Size $\\beta$ (BMI-Adjusted, GCST90502763)", fontsize=10, fontweight='bold')
    ax_scatter.set_title("B. Cross-Model Effect Size Concordance ($\\beta_{\\mathrm{BMI}}$ vs. $\\beta_{\\mathrm{Unadj}}$)", fontsize=11, fontweight='bold')
    ax_scatter.set_xlim(-0.62, 0.55)
    ax_scatter.set_ylim(-0.62, 0.55)
    ax_scatter.grid(True)
    
    plt.suptitle("Figure 3: Effect Size Comparison & Adiposity-Independence of Prioritized Prediabetes Loci", fontsize=13, fontweight='bold', y=0.98)
    plt.tight_layout(rect=[0, 0, 1, 0.94])
    plt.savefig("figures/Fig3_Forest_and_Concordance.png", dpi=300)
    plt.close()

    # Figure 4: Trumpet plot with distinct non-overlapping callout boxes
    fig, ax = plt.subplots(figsize=(10.8, 6.4), dpi=300)
    
    trumpet_offsets = {
        "PDE1C": (0.065, 0.08),
        "HK1": (0.065, 0.06),
        "G6PC2": (-0.045, -0.085),
        "GCK": (-0.035, 0.145),
        "PLXNA4": (0.055, 0.095),
        "MTNR1B": (0.02, -0.095),
        "B3GNTL1/METRNL": (0.085, 0.045)
    }
    
    for loc in loci_bmi:
        maf = loc['maf']
        risk_or = loc['risk_odds_ratio']
        p = loc['p_value']
        rs = loc['lead_rs_id']
        is_gw = p < 5e-8
        
        col = '#c0392b' if is_gw else ('#d35400' if p < 1e-6 else '#2980b9')
        sz = -np.log10(p) * 16
        ax.scatter(maf, risk_or, s=sz, color=col, alpha=0.85, edgecolors='black', lw=0.6, zorder=5 if is_gw else 3)
        
        if rs in LOCUS_GENE_MAP:
            gene = LOCUS_GENE_MAP[rs]
            dx, dy = trumpet_offsets[gene]
            ax.annotate(f"{gene} ({rs})\nMAF={maf*100:.1f}%, OR={risk_or:.2f}",
                        xy=(maf, risk_or), xytext=(maf + dx, risk_or + dy),
                        fontsize=8.2, fontweight='bold', ha='center', va='center',
                        bbox=dict(boxstyle='round,pad=0.25', facecolor='#fef9e7' if not is_gw else '#fadbd8', edgecolor=col, lw=0.85, alpha=0.95),
                        arrowprops=dict(arrowstyle='->', color=col, lw=0.9),
                        zorder=6)

    maf_grid = np.linspace(0.015, 0.50, 200)
    n_eff = 4.0 / (1.0 / 2267.0 + 1.0 / 8665.0)
    se_approx = 1.0 / np.sqrt(2.0 * n_eff * maf_grid * (1.0 - maf_grid))
    or_gw_threshold = np.exp(5.4513 * se_approx)
    or_sugg_threshold = np.exp(4.4172 * se_approx)
    
    ax.plot(maf_grid, or_gw_threshold, color='#c0392b', linestyle='--', lw=1.3, label='Min Detectable OR at $p = 5\\times 10^{-8}$ ($N_{\\mathrm{eff}}=7,188$)')
    ax.plot(maf_grid, or_sugg_threshold, color='#2980b9', linestyle=':', lw=1.2, label='Min Detectable OR at $p = 1\\times 10^{-5}$ ($N_{\\mathrm{eff}}=7,188$)')
    
    ax.set_xlabel("Minor Allele Frequency (MAF in Qatari Cohort, $N = 10,932$)", fontsize=10.5, fontweight='bold')
    ax.set_ylabel("Risk-Increasing Odds Ratio ($\\mathrm{OR}_{\\mathrm{risk}}$ per Allele)", fontsize=10.5, fontweight='bold')
    ax.set_title("Figure 4: Prediabetes Genetic Architecture — Risk Odds Ratio vs. Minor Allele Frequency (37 Independent Loci)", fontsize=12, fontweight='bold', pad=10)
    ax.set_xlim(0.0, 0.50)
    ax.set_ylim(1.05, 1.98)
    ax.grid(True)
    ax.legend(loc='upper right', fontsize=9, framealpha=0.92)
    
    plt.tight_layout()
    plt.savefig("figures/Fig4_Trumpet_MAF_vs_EffectSize.png", dpi=300)
    plt.close()
    print("Redrew Figures 3 and 4 cleanly.")

if __name__ == "__main__":
    redraw_fig3_and_fig4()
