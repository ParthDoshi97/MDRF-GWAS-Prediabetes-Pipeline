# MDRF Use Case 1: GWAS Summary Statistics Analysis Pipeline

This repository contains the end-to-end Python bioinformatics pipeline for analyzing the **Qatar Biobank Prediabetes GWAS Summary Statistics** from the NHGRI-EBI GWAS Catalog:
* **`GCST90502762`**: Prediabetes (BMI-Unadjusted; $8,490,370$ variants)
* **`GCST90502763`**: Prediabetes (BMI-Adjusted; $8,488,404$ variants)
* **Original Study:** Aliyu et al., *Diabetes, Obesity and Metabolism* (2025) | **PMID:** `40468860` | **Genome Build:** `GRCh38`

---

## Folder Structure

```text
MDRF-GWAS-Prediabetes-Pipeline/
├── README.md                  # Quick-start guide and pipeline overview
├── requirements.txt           # Required Python libraries
├── .gitignore                 # Excludes raw data/ (~1.6 GB) and separate report files
├── scripts/                   # Step-by-step analysis scripts (Tasks 1–6)
│   ├── 00_download_data.py        # Step 0: Downloads GWAS TSV files & verifies MD5 checksums
│   ├── 01_explore_dataset.py      # Task 1: Explores schema, chromosomes, & value ranges
│   ├── 02_quality_control.py      # Task 2: Removes 45 SAIGE outliers & computes lambda_GC
│   ├── 03_identify_signals.py     # Task 3: +-500 kb LD clumping & locus prioritization
│   ├── 04_visualize_results.py    # Task 4: Extracts regional windows & builds base plots
│   ├── 04b_polish_figures.py      # Task 4: Generates 300 DPI publication figures (Fig1–Fig5)
│   ├── 05_biological_annotation.py# Task 5: Synthesizes VEP, GTEx, Open Targets, ClinVar, Reactome
│   └── 06_proposed_next_steps.py  # Task 6: Statistical power & replication sample size calculator
├── figures/                   # Generated 300 DPI visualizations (Task 4)
│   ├── Fig1_QQ_Plots_Comparison.png
│   ├── Fig2_Manhattan_Plots_Comparison.png
│   ├── Fig3_Forest_and_Concordance.png
│   ├── Fig4_Trumpet_MAF_vs_EffectSize.png
│   └── Fig5_Regional_LocusZoom_GCK_HK1.png
└── results/                   # Output summary tables (TSV and JSON) from Tasks 1–6
```

---

## How to Run the Pipeline

### 1. Install Dependencies
```bash
pip install -r requirements.txt
```

### 2. Run Each Step in Order
```bash
# Step 0: Download GCST90502762.tsv and GCST90502763.tsv into data/ and verify MD5 checksums
python scripts/00_download_data.py

# Task 1: Dataset Exploration (outputs results/task1_exploration_summary.json)
python scripts/01_explore_dataset.py

# Task 2: Quality Control & Genomic Inflation lambda_GC (outputs results/task2_qc_summary.json)
python scripts/02_quality_control.py

# Task 3: Identify & Prioritize Independent Loci (outputs results/task3_clumped_loci_*.tsv)
python scripts/03_identify_signals.py

# Task 4: Generate Visualizations (outputs figures/Fig1_*.png through Fig5_*.png)
python scripts/04_visualize_results.py
python scripts/04b_polish_figures.py

# Task 5: Clinical & Biological Annotation Synthesis (outputs results/task5_biological_interpretation.json)
python scripts/05_biological_annotation.py

# Task 6: Statistical Power & Replication Sample Size Projections (outputs results/task6_proposed_next_steps.json)
python scripts/06_proposed_next_steps.py
```

---

## Summary of Key Outputs

| Task | Script | Key Output File(s) | Main Result |
| :--- | :--- | :--- | :--- |
| **Task 1** | `01_explore_dataset.py` | `results/task1_exploration_summary.json` | ~8.49M autosomal variants (`Chr 1–22`, `GRCh38`), `0` missing values, $\text{MAF} \ge 0.01$. |
| **Task 2** | `02_quality_control.py` | `results/task2_qc_summary.json` | Removed $45$ SAIGE non-convergence outliers ($\text{SE} > 10^{13}$); $\lambda_{\text{GC}} = 1.0000$ in both models. |
| **Task 3** | `03_identify_signals.py` | `results/task3_clumped_loci_GCST90502763_bmi_adj.tsv` | $37$ independent loci ($p < 10^{-5}$); Top 7 prioritized: **`GCK`**, **`HK1`**, **`PDE1C`**, **`PLXNA4`**, **`MTNR1B`**, **`G6PC2`**, **`FN3KRP`/`METRNL`**. |
| **Task 4** | `04b_polish_figures.py` | `figures/Fig1`–`Fig5` (`.png`) | Q–Q, Stacked Manhattan, Forest & Concordance ($r = 0.9983$), Trumpet ($\text{MAF}$ vs. $\beta$), & Regional LocusZoom plots. |
| **Task 5** | `05_biological_annotation.py` | `results/task5_biological_interpretation.json` | Maps loci to **Axis 1** ($\beta$-cell glucose sensing & $\text{cAMP}$ secretion) and **Axis 2** (RBC glycolysis & $\text{HbA}_{1\text{c}}$ glycation). |
| **Task 6** | `06_proposed_next_steps.py` | `results/task6_proposed_next_steps.json` | $N \approx 7,500\text{--}7,750$ needed for $80\%$ Bonferroni replication power; $N_{\text{total}} \approx 17,800$ for $80\%$ GW discovery power. |
