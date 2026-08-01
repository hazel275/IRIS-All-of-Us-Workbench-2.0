# IRIS: Intersectional Regulation of Immune and Stress Endpoints

## Workflow Documentation

**Study:** The Intersectional Psychoneuroendocrinology Framework: Social and Geographic Patterning of Systemic Inflammation\
**Author:** Stephanie H. Cook, NYU School of Global Public Health\
**Data source:** All of Us Research Program, Controlled Tier\
**Submitted:** International Society of Psychoneuroendocrinology\
**Last updated:** 2026-08-01 \| Version: v9 (locked run 20260801_1815)

---

## Repository Contents

| File | Description |
|------|-------------|
| `code/pipeline_20260801.R` | Complete self-contained analysis pipeline (v9). Run in RStudio within the All of Us Workbench environment. |
| `README.md` | This document. |
| Aggregate result CSVs and diagnostic PNGs (`*_20260801_1815`) | Public-safe position-level and model-level outputs. |

**Not included (data use agreement compliance):**\
No participant-level data, no `.RData` or `.rds` workspace files, no individual-level outputs. All analytic datasets remain inside the All of Us Researcher Workbench. Output CSVs containing only aggregated position-level statistics (n>=27 per cell) are safe for sharing and are included in the published supplement.

---

## Data Access

Analyses require access to the All of Us Research Program Controlled Tier.\
Apply at: <https://www.researchallofus.org/register/>

CDR version used: **INSERT YOUR CDR RELEASE HERE**\
Analysis environment: **RStudio**, launched within the All of Us Researcher Workbench

**Why RStudio and not the Workbench Python+R environment?** `lme4` could not be installed in Workbench 2.0. The pre-compiled binary failed and source compilation requires `cmake`, which is not available in that environment. RStudio was launched from within the Workbench and connects to the same BigQuery CDR and workspace bucket. All data access, credentialing, and DUA compliance requirements are identical. If you encounter `lme4` installation errors in the Python+R environment, switch to RStudio.

---

## How to Run

1. Log into the All of Us Researcher Workbench at <https://workbench.researchallofus.org>
2. Create a new workspace or open an existing one with Controlled Tier access
3. Upload `pipeline_20260801.R` to your workspace
4. Open a new Jupyter notebook or R terminal
5. Set your environment values, then source the script
6. All outputs save automatically to your workspace bucket

**Runtime:** Approximately 20-40 minutes. The parametric bootstrap (Section 10, 1,000 iterations) accounts for several minutes of that time.

**Critical configuration:** Set your workspace values before running. The pipeline reads these from the environment and does not hardcode them:

- Workspace bucket: `INSERT YOUR BUCKET HERE`
- Google Cloud project: `INSERT YOUR PROJECT HERE`
- CDR release: `INSERT YOUR CDR RELEASE HERE`
- Deprivation source commit and checksum (see Software section)

---

## Analytic Design

### Sample

- N=14,601 adults in the primary sample (hs-CRP 0-10 mg/L) from the All of Us Research Program
- Primary outcome: serum high-sensitivity C-reactive protein (hs-CRP)
- Sensitivity samples: 0-20 mg/L (N=16,142), broad-range 0-200 mg/L (N=18,463), capped at 10 (N=18,463)
- Geographic sample: N=11,514 across 88 three-digit ZIP areas and 32 sites

### Design

Cross-sectional and associational throughout. The Basics survey was completed primarily in 2020-2021. hs-CRP is pulled retrospectively from the full clinical history. About 65% of participants had hs-CRP recorded before Basics completion, so the analysis is associational and cannot establish direction. All findings describe how systemic inflammatory burden is distributed across intersectional structural positions. No causal claims are made.

### One observation per participant

For hs-CRP, the pipeline selects the measurement closest to the Basics anchor within a five-year window, breaking ties by earliest date then lexicographic order. Each participant contributes exactly one measurement.

### Intersectional positions

27 positions defined by crossing race and ethnicity (White, Black, Hispanic, Other/Asian), gender (Man, Woman, Gender Minority collapsed across race due to small cells), and annual household income (Low, Middle, High). Positions with n<27 are not reported. Minimum position cell size is 27 in the primary sample.

### Method

Multilevel Analysis of Individual Heterogeneity and Discriminatory Accuracy (MAIHDA). Individuals at level 1, intersectional positions at level 2. REML estimation throughout for VPC calculations. A three-level geographic extension adds three-digit ZIP areas as level 3.

---

## Script Structure

| Section | Content |
|---------|---------|
| Connect and audit | Environment setup, BigQuery connection, audit-gate machinery |
| Anchor | Eligible population from four core surveys |
| hs-CRP | Closest-to-anchor selection with in-query site join |
| Covariates | Gender at source, race distinct aggregation, income, education, age, BMI, smoking |
| Positions | 27-position construction |
| Samples | Primary, 0-20, broad, capped |
| Models 1-4 | MAIHDA null through adjusted-additive |
| Axis-count | Structural-dimension contribution decomposition |
| Geographic | Three-level ZIP model, site adjustment, influential-area sensitivity |
| Robustness | Parametric bootstrap, cluster-robust inference, diagnostics |
| Export | Master results table and aggregate outputs |

---

## Output Files

The pipeline generates aggregated CSVs and figures. Files marked **public-safe** contain only position-level or model-level statistics with n>=27 per cell and may be shared publicly. All participant-level files remain in the Workbench bucket.

| File | Type | Public-safe |
|------|------|-------------|
| `MASTER_results_*.csv` | Master results table | Yes |
| `MASTER_geographic_results_*.csv` | Geographic model results | Yes |
| `M4_fixed_*.csv` | Fixed-effect coefficients per specification | Yes |
| `axis_count_*.csv` | Structural-dimension contributions | Yes |
| `race_mapping_*.csv`, `gender_mapping_*.csv` | Category mappings | Yes |
| `OLS_ZIP3_FE_sensitivity_*.csv`, `robust_fixed_CR2_*.csv` | Robustness and cluster-robust inference | Yes |
| `diag_resid_qq_*.png`, `diag_resid_fitted_*.png` | Residual diagnostics | Yes |
| `diag_ranef_position_qq_*.png`, `diag_ranef_zip3_qq_*.png` | Random-effect diagnostics | Yes |
| `samp_*.rds`, `geo.rds`, `closest.rds` | Individual-level model objects | **No, Workbench only** |

---

## Key Analytic Decisions

All decisions are logged in the pipeline and enforced by audit gates. Key decisions:

**hs-CRP:** Winsorized at the 99th percentile, log-transformed as log(hs-CRP + 0.1). Standard in CRP literature. Winsorizing reduces the influence of acute-infection values. Robustness evaluated across four outcome-range and capping specifications.

**Primary outcome identification:** hs-CRP is concept 3010156 (LOINC 30522-7), the true high-sensitivity assay, not standard CRP. Units harmonized from mg/dL to mg/L before range filtering.

**BMI:** Values outside plausible range excluded, then modeled on the log scale.

**Smoking:** Binary ever/never (100+ lifetime cigarettes). Does not distinguish current from former smokers. Limitation noted in the paper.

**Income:** Three categories, Low, Middle, High, based on the All of Us survey response categories.

**Gender:** Harmonized at source to three levels, Man, Woman, Gender Minority. Gender Minority collapsed across race due to small cells.

**Narrowed to inflammation:** HbA1c, diabetes, DHEA-S, and perceived stress were removed. The revised paper narrows the empirical demonstration to systemic inflammation indexed by hs-CRP, the outcome aligned most directly with the framework's inflammatory endpoint.

**Missing data:** Complete-case analysis. No multiple imputation.

**Geographic influence:** One low-inflammation ZIP3 area (782) drives a substantial share of the geographic variance. Both estimates, with and without the area, are reported.

---

## Software

R version 4.6.1 (All of Us Workbench, RStudio)\
Key packages: `lme4` 2.0.6, `lmerTest` 3.2.1, `clubSandwich`, `dplyr` 1.2.1\
Deprivation index: geomarker-io 2023 ACS deprivation index, <https://github.com/geomarker-io/dep_index>

Deprivation source is pinned for reproducibility:
- Commit SHA: bd897cedbe1da2bcd2649b0222d425e277dead38
- MD5 checksum: ba08b659aa5b25ecbbe80a7f208b3b00
- Aggregation: unweighted mean of ZCTA-level index within each 3-digit ZIP prefix

---

## Citation

Cook, S. H. (2026). The Intersectional Psychoneuroendocrinology Framework: Social and Geographic Patterning of Systemic Inflammation. *Psychoneuroendocrinology*.

---

## Contact

Stephanie H. Cook, PhD\
Associate Professor, NYU School of Global Public Health
