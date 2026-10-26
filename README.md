# GWAS Summary Statistics Analysis: Prediabetes in the Qatar Biobank Cohort (`GCST90502762` & `GCST90502763`)

This repository contains the end-to-end reproducible bioinformatics, statistical genetics, quality control, locus prioritization, visualization, and clinical/biological interpretation pipeline for **Use Case 1 – GWAS Summary Statistics Analysis**.

* **Original Study:** Aliyu et al., *"Genetic risk and polygenic risk score assessment of prediabetes and progression to type 2 diabetes"*, *Diabetes, Obesity and Metabolism* (2025); **PMID:** `40468860` | **DOI:** `10.1111/dom.16490`
* **GWAS Catalog Accessions:**
  * `GCST90502762`: Prediabetes (Unadjusted for BMI — adjusted for Age, Gender, `PC1–PC10`; $8,490,370$ autosomal variants)
  * `GCST90502763`: Prediabetes (BMI-Adjusted — adjusted for Age, Gender, **BMI**, `PC1–PC10`; $8,488,404$ autosomal variants)
* **Cohort:** Qatar Biobank (QBB) Discovery Cohort ($N = 10,932$: $2,267$ prediabetes cases and $8,665$ normoglycemic controls; Whole Genome Sequencing on **GRCh38**, analyzed with **SAIGE**).
* **Complete Report:** See [`Final_Scientific_Report.md`](Final_Scientific_Report.md) for the full scientific report covering Objective, Methodology, Results, Biological Interpretation, Proposed Next Steps, Limitations, and Conclusion.

---

## Repository Structure

```text
task1/
├── README.md                                           # Repository guide & execution instructions
├── Final_Scientific_Report.md                          # Complete scientific report (Tasks 1–7)
├── requirements.txt                                    # Python package dependencies
├── .gitignore                                          # Ignores 814 MB raw TSVs & temporary cache files
├── scripts/
│   ├── 00_download_data.py                             # Step 0: 16-thread HTTP byte-range downloader + MD5 check
│   ├── 01_explore_dataset.py                           # Task 1: Dataset structure, coverage, & anomaly scan
│   ├── 02_quality_control.py                           # Task 2: 6-step QC pipeline, SAIGE outlier removal, & lambda_GC
│   ├── 03_identify_signals.py                          # Task 3: +-500 kb clumping & multi-criteria prioritization
│   ├── 04_visualize_results.py                         # Task 4: Regional window extraction & initial plots
│   ├── 04b_polish_figures.py                           # Task 4: Publication-grade 300 DPI figures (Fig1–Fig5)
│   ├── 05_biological_annotation.py                     # Task 5: Ensembl VEP, GTEx v10, Open Targets, ClinVar, & Reactome
│   └── 06_proposed_next_steps.py                       # Task 6: Statistical power & replication sample size calculator
├── figures/
│   ├── Fig1_QQ_Plots_Comparison.png                    # MAF-stratified Q-Q plots (GCST90502762 vs GCST90502763)
│   ├── Fig2_Manhattan_Plots_Comparison.png             # Stacked genome-wide Manhattan comparison (GRCh38)
│   ├── Fig3_Forest_and_Concordance.png                 # Risk-standardized Forest plot & beta concordance scatter
│   ├── Fig4_Trumpet_MAF_vs_EffectSize.png              # MAF vs. Beta Trumpet plot with 50% and 80% power curves
│   └── Fig5_Regional_LocusZoom_GCK_HK1.png             # +-250 kb regional association plots for GCK, HK1, & PDE1C
├── results/
│   ├── task1_exploration_summary.json                  # Task 1 summary metrics
│   ├── task2_qc_summary.json                           # Task 2 QC metrics, lambda_GC, & Wald consistency
│   ├── GCST90502762_qc_p1e-3.tsv                       # QC-passed candidate variants (p < 1e-3, BMI-unadjusted)
│   ├── GCST90502763_qc_p1e-3.tsv                       # QC-passed candidate variants (p < 1e-3, BMI-adjusted)
│   ├── task3_signals_summary.json                      # Task 3 locus clumping & prioritization summary
│   ├── task3_clumped_loci_GCST90502763_bmi_adj.tsv     # 37 independent clumped loci (BMI-adjusted)
│   ├── task3_clumped_loci_GCST90502762_unadj.tsv       # 34 independent clumped loci (BMI-unadjusted)
│   ├── task3_top_variants_comparison.tsv               # Cross-model effect size comparison table
│   ├── task5_biological_interpretation.json            # Synthesized multi-database biological annotations
│   ├── reactome_enrichment.json                        # Reactome pathway over-representation results
│   └── task6_proposed_next_steps.json                  # Quantitative replication & meta-analysis sample size projections
└── data/
    ├── GCST90502762.tsv-meta.yaml                      # Official GWAS Catalog YAML metadata (BMI-unadjusted)
    ├── GCST90502763.tsv-meta.yaml                      # Official GWAS Catalog YAML metadata (BMI-adjusted)
    ├── GCST90502762.running.log                        # EBI harmonisation pipeline log (BMI-unadjusted)
    └── GCST90502763.running.log                        # EBI harmonisation pipeline log (BMI-adjusted)
```

---

## Reproducing the Entire Analysis Pipeline

### 1. Environment Setup
```bash
python -m venv .venv
# Windows PowerShell:
.venv\Scripts\Activate.ps1
# Linux / macOS:
# source .venv/bin/activate

pip install -r requirements.txt
```

### 2. Step-by-Step Execution
```bash
# Step 0: Download raw GWAS summary statistics (~1.63 GB total) and verify official EBI MD5 checksums
python scripts/00_download_data.py

# Task 1: Explore dataset schema, chromosome coverage, value ranges, and anomalies
python scripts/01_explore_dataset.py

# Task 2: Run 6-step Quality Control, filter 45 SAIGE non-convergence outliers, and compute lambda_GC
python scripts/02_quality_control.py

# Task 3: Perform +-500 kb locus clumping, cross-model comparison, and multi-criteria prioritization
python scripts/03_identify_signals.py

# Task 4: Extract regional windows and generate publication-grade 300 DPI visualizations (Fig1–Fig5)
python scripts/04_visualize_results.py
python scripts/04b_polish_figures.py

# Task 5: Synthesize Ensembl VEP, GTEx v10, Open Targets L2G, ClinVar, and Reactome annotations
python scripts/05_biological_annotation.py

# Task 6: Compute statistical power and required replication/meta-analysis sample sizes
python scripts/06_proposed_next_steps.py
```

---

## Pushing to Your GitLab Repository

This directory is initialized as a clean Git repository on branch `main` (with large raw data files `data/*.tsv` excluded via `.gitignore`). To push to your GitLab remote:

```bash
git remote add origin https://gitlab.com/<your-username>/<your-repo-name>.git
git branch -M main
git push -u origin main
```
