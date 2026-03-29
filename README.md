# IRISE: Intersectional Regulation of Immune and Stress Endpoints

## Workflow Documentation

**Study:** The Substrate Is Not Neutral: Intersectional Psychoneuroendocrinology as a Framework for Neuroendocrine Regulation  
**Author:** Stephanie H. Cook, NYU School of Global Public Health  
**Data source:** All of Us Research Program, Controlled Tier CDR v8  
**Submitted:** Dirk Hellhammer Award 2026, International Society of Psychoneuroendocrinology  
**Last updated:** 2026-03-29 | Version: v2  

---

## Repository Contents

| File | Description |
|------|-------------|
| `IRISE_Complete_Analysis_v2.R` | Complete analysis script (v2, 2026-03-29). Run in RStudio within the All of Us Workbench environment. |
| `README.md` | This document. |

**Not included (data use agreement compliance):**  
No participant-level data, no `.RData` workspace files, no individual-level outputs. All analytic datasets remain inside the All of Us Researcher Workbench. Output CSVs containing only aggregated position-level statistics (N≥30 per cell) are safe for sharing and are included in the published supplement.

---

## Data Access

Analyses require access to the All of Us Research Program Controlled Tier.  
Apply at: https://www.researchallofus.org/register/  

CDR version used: **CDR v8** (`wb-silky-artichoke-2408.C2024Q3R8`)  
Analysis environment: **RStudio**, launched within the All of Us Researcher Workbench  

**Why RStudio and not the Workbench Python+R environment?** `lme4` could not be installed in Workbench 2.0. The pre-compiled binary failed and source compilation requires `cmake`, which is not available in that environment. RStudio was launched from within the Workbench and connects to the same BigQuery CDR and workspace bucket. All data access, credentialing, and DUA compliance requirements are identical. If you encounter `lme4` installation errors in the Python+R environment, switch to RStudio.

---

## How to Run

1. Log into the All of Us Researcher Workbench at https://workbench.researchallofus.org
2. Create a new workspace or open an existing one with Controlled Tier access
3. Upload `IRISE_Complete_Analysis_v2.R` to your workspace
4. Open a new Jupyter notebook or R terminal
5. Run the script in sequence from Section 0 through Section 23
6. All outputs save automatically to your workspace bucket

**Runtime:** Approximately 45-75 minutes. The bootstrap mediation (Section 17, 1,000 iterations) accounts for 20-30 minutes of that time.

**Critical configuration:** Update `WORKSPACE_BUCKET` in Section 0 to match your workspace bucket before running. The current script uses a hardcoded bucket specific to the original workspace.

---

## Analytic Design

### Sample
- N=47,059 adults from the All of Us Research Program
- Primary outcome: serum high-sensitivity C-reactive protein (hsCRP)
- Replication outcomes: glycated hemoglobin (HbA1c), type 2 diabetes prevalence
- Perceived stress: PSS-10 subsample N=23,059

### Design
Cross-sectional and descriptive throughout. The Basics survey was completed primarily in 2020–2021. EHR biomarkers are pulled retrospectively from the full clinical history dating to 1987. 53.6% of participants had hsCRP recorded before Basics completion. All findings describe how inflammatory burden is distributed across intersectional structural positions. No causal claims are made.

### One observation per participant
For all biomarkers (hsCRP, HbA1c, BMI), the script selects the most recent available measurement per participant using `slice(1)` after descending date sort. Each participant contributes exactly one measurement. No date-proximity restriction is applied between biomarkers or between biomarkers and survey date. The temporal gap distribution between hsCRP and Basics survey completion is reported in Table S5 of the supplement.

### Intersectional positions
27 positions defined by crossing:
- Race and ethnicity: Non-Hispanic White, Non-Hispanic Black, Hispanic or Latino, Other/Asian (5 categories, Asian and Other collapsed due to small cells)
- Gender identity: cisgender woman, cisgender man, gender minority (collapsed across race)
- Annual household income: low (<$35,000), middle ($35,000–$99,999), high (≥$100,000)

Positions with N<30 excluded from position-level reporting.

### Method
Multilevel Analysis of Individual Heterogeneity and Discriminatory Accuracy (MAIHDA). Individuals at level 1, intersectional positions at level 2. REML estimation throughout for VPC calculations. Three-level geographic extension adds three-digit zip code areas as level 3. Full model sequence in supplement Table S2 and Section S7.

---

## Script Structure

| Section | Content |
|---------|---------|
| 0 | Setup, package installation, BigQuery connection |
| 1 | hsCRP pull and processing |
| 2 | Demographics (race, gender, age) |
| 3 | Covariates (zip code, HbA1c, BMI, smoking, education) |
| 4 | Merge and missing data audit |
| 5 | Intersectional variable recoding |
| 6 | Position definitions |
| 7 | Analytic sample construction |
| 8 | Neighborhood deprivation index |
| 9 | Table 1 (sample characteristics) |
| 10 | MAIHDA models for hsCRP (Models 1–5) |
| 11 | Table 2 (position-level hsCRP rankings) |
| 12 | Table 3 (VPC decomposition) |
| 13 | Table 4 (fixed effects) |
| 14 | HbA1c and diabetes replication |
| 15 | Table 5 (replication summary) |
| 16 | PSS pathway analysis (cross-sectional descriptive) |
| 17 | Supplementary mediation (N=6,023 temporally restricted) |
| 18–22 | Figures 1–5 |
| 23 | Save all results |

---

## Output Files

The script generates the following output files. Aggregated CSVs (marked **public-safe**) contain only position-level statistics with N≥30 per cell and may be shared publicly. All other files remain in the Workbench bucket.

| File | Type | Public-safe |
|------|------|-------------|
| `IRISE_Table1_threecol.csv` | Sample characteristics | Yes |
| `IRISE_Table2_positions_hsCRP.csv` | Position-level hsCRP | Yes |
| `IRISE_Table3_VPC.csv` | VPC decomposition | Yes |
| `IRISE_Table4_fixed_effects.csv` | Fixed effects | Yes |
| `IRISE_Table5_replication.csv` | Replication summary | Yes |
| `IRISE_TableS3_PSS_pathway.csv` | PSS pathway results | Yes |
| `IRISE_TableS4_mediation.csv` | Mediation results | Yes |
| `IRISE_TableS5_temporal_ordering.csv` | Temporal ordering | Yes |
| `IRISE_PSS_by_position.csv` | PSS by position | Yes |
| `IRISE_position_diabetes.csv` | Diabetes by position | Yes |
| `IRISE_Figure1_caterpillar.png` | Caterpillar plot | Yes |
| `IRISE_Figure2_heatmap.png` | Heatmap | Yes |
| `IRISE_Figure3_VPC_waterfall.png` | VPC waterfall | Yes |
| `IRISE_Figure4_deprivation_scatter.png` | Deprivation scatter | Yes |
| `IRISE_Figure5_diabetes.png` | Diabetes dot plot | Yes |
| `IRISE_complete_workspace.RData` | All models and datasets | **No — Workbench only** |
| `IRISE_mediation_final.RData` | Mediation workspace | **No — Workbench only** |

---

## Key Analytic Decisions

All decisions are logged in the script using `log_decision()` and printed at the end of the run. Key decisions:

**hsCRP:** Winsorized at 99th percentile, log-transformed as log(hsCRP + 0.1). Standard in CRP literature. Winsorizing reduces influence of acute infection values.

**HbA1c:** Binary diabetes indicator (HbA1c ≥ 6.5%) is primary replication outcome due to bimodal distribution in mixed diabetic/non-diabetic sample.

**BMI:** Values outside 15–60 kg/m² excluded as implausible EHR errors. Winsorized at 99th percentile. Log-transformed.

**Smoking:** Binary ever/never (100+ lifetime cigarettes). Does not distinguish current from former smokers. Limitation noted in paper.

**Income:** Three categories. Low <$35,000, Middle $35,000–$99,999, High ≥$100,000. Based on All of Us survey response categories in CDR v8 full string format.

**Gender minority:** Collapsed across race due to small cells. Minimum position cell size 63 in adjusted sample.

**Missing data:** Complete case analysis. No multiple imputation. Missing ranges 3–11% across covariates.

**Mediation:** Supplementary only. Restricted to N=6,023 where PSS preceded hsCRP. Product of coefficients with 1,000 bootstrap iterations (Hayes, 2018). Framed as confirmatory sensitivity analysis with explicit selection bias caveat.

---

## Software

R version: 4.x (All of Us Workbench R environment)  
Key packages: `lme4`, `lmerTest`, `bigrquery`, `tidyverse`, `ggplot2`, `mediation`  
Deprivation index: geomarker-io 2023 ACS (https://github.com/geomarker-io/dep_index)

---

## Citation

Cook, S. H. (2026). The substrate is not neutral: Intersectional psychoneuroendocrinology as a framework for neuroendocrine regulation. *Psychoneuroendocrinology*.

---

## Contact

Stephanie H. Cook, PhD  
Associate Professor, NYU School of Global Public Health  
