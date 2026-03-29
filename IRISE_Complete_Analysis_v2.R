# ============================================================
# IRISE: Intersectional Regulation of Immune and Stress Endpoints
# Complete Final Analysis Script
# All of Us Research Program | Controlled Tier CDR v8
# Dirk Hellhammer Award 2026
# Last updated: 2026-03-29
# Version: v2
#
# DATA AVAILABILITY
# -----------------
# This script is designed to run inside the All of Us Researcher
# Workbench (Controlled Tier CDR v8). It queries BigQuery directly.
# No participant-level data are included in this repository.
# Output CSVs contain only aggregated position-level statistics
# (N>=30 per cell) and are safe for public sharing under the
# All of Us DUA. Raw analytic datasets (analytic_adj, etc.) are
# saved only to the private Workbench bucket and must NOT be
# uploaded to public repositories.
#
# PUBLIC REPOSITORY
# -----------------
# GitHub: https://github.com/hazel275/IRIS-All-of-Us-Workbench-2.0
# Only this R script and workflow documentation are public.
# No data files, no .RData files, no participant-level outputs.
# ============================================================
#
# ANALYTIC DECISIONS SUMMARY
# ---------------------------
# 1. CROSS-SECTIONAL DESIGN: All analyses are cross-sectional and
#    descriptive. All of Us Basics survey (income, education, smoking,
#    gender identity) was completed primarily in 2020. EHR biomarkers
#    (hsCRP) are pulled retrospectively from the full clinical history,
#    with measurements dating back to 1987. 53.6% of participants had
#    hsCRP measured before Basics completion. The analysis describes
#    how inflammatory burden is distributed across intersectional
#    structural positions. No causal claims are made.
#
# 2. MAIHDA: REML used throughout for VPC estimation per Evans et al.
#    2018. REML produces unbiased variance component estimates.
#    ML used only in PSS pathway models where product of coefficients
#    mediation requires ML.
#
# 3. hsCRP: Winsorized at 99th percentile (84.5 mg/L), log-transformed
#    as log(hsCRP_w + 0.1). Standard approach in CRP literature.
#
# 4. HbA1c: Winsorized at 99th percentile, log-transformed. Binary
#    diabetes indicator (>=6.5%) is primary replication outcome due to
#    bimodal distribution in mixed diabetic/non-diabetic sample.
#
# 5. BMI: Values outside 15-60 kg/m2 excluded as implausible EHR
#    errors. Remaining values winsorized at 99th percentile (51.9).
#    Log-transformed for modeling.
#
# 6. Age: Winsorized at 99th percentile (90 years). Centered at mean
#    for lmer.
#
# 7. Smoking: Binary ever/never (100+ lifetime cigarettes). PMI
#    responses coded as missing. Does not distinguish current from
#    former smokers - noted as limitation.
#
# 8. Education: Three categories Low/Middle/High. PMI = missing.
#
# 9. Income: Three categories Low/Middle/High. PMI = missing.
#    NOTE: CDR v8 uses full string format throughout.
#    Both primary and PSS analyses use "Annual Income: less 10k" etc.
#    PMI: Skip, Prefer Not To Answer, Dont Know all mapped to NA.
#    Verify unique values before running any recode.
#
# 10. Intersectional positions: 27 positions from race5c (4) x
#     gender3 (3) x income3 (3). Gender minority collapsed across
#     race. Asian and Other collapsed to Other/Asian. Min cell = 63
#     in adjusted sample.
#
# 11. Geographic models: Restricted to zip3 areas with N >= 30
#     (83 areas, N=39,422). Deprivation index from geomarker-io
#     2023 ACS, averaged within 3-digit zip prefix.
#
# 12. PSS pathway: Cross-sectional descriptive evidence only.
#     73.9% of PSS+hsCRP overlap had hsCRP before PSS completion.
#     Temporal ordering precludes formal mediation in main paper.
#     PSS VPC and association with hsCRP reported as pathway-
#     consistent evidence.
#
# 13. Mediation: Supplementary only. Restricted to N=6,023 where
#     PSS preceded hsCRP. Product of coefficients with 1000
#     bootstrap iterations. Framed as confirmatory sensitivity
#     analysis with explicit selection bias caveat.
#
# 14. BUCKET: WORKSPACE_BUCKET environment variable is misconfigured
#     in this workspace. Use hardcoded correct bucket:
#     gs://cloned-aou-tutorial-notebooks-wb-frosty-sprout-5598
#
# ============================================================

# ============================================================
# SECTION 0: SETUP
# ============================================================
# ENVIRONMENT NOTE: This script was developed and run in RStudio,
# not in the All of Us Researcher Workbench Python+R environment.
# lme4 could not be installed in Workbench 2.0 (the pre-compiled
# binary failed and source compilation requires cmake which is not
# available). RStudio was launched from within the Workbench and
# connects to the same BigQuery CDR and bucket. All data access,
# credentialing, and DUA compliance requirements are identical.
# If you are running this in the Workbench Python+R environment
# and encounter lme4 installation errors, switch to RStudio.
#
# CRITICAL: Do NOT specify repos= when installing.
# Default repos use p3m.dev pre-compiled binaries.
# Specifying cloud.r-project.org forces source compilation
# which fails due to missing cmake.

install.packages("lme4")
install.packages("lmerTest")
install.packages("ggplot2")
install.packages("scales")
install.packages("viridis")
install.packages("mediation")

library(bigrquery)
library(DBI)
library(tidyverse)
library(lme4)
library(lmerTest)
library(ggplot2)
library(scales)
library(viridis)

# Environment
Sys.setenv(WORKSPACE_CDR       = "wb-silky-artichoke-2408.C2024Q3R8")
Sys.setenv(GOOGLE_CLOUD_PROJECT = "wb-frosty-sprout-5598")

cdr     <- Sys.getenv("WORKSPACE_CDR")
project <- Sys.getenv("GOOGLE_CLOUD_PROJECT")

# WORKSPACE_BUCKET env variable is misconfigured - use hardcoded path
correct_bucket <- "gs://cloned-aou-tutorial-notebooks-wb-frosty-sprout-5598"

cat("CDR:", cdr, "\n")
cat("Project:", project, "\n")
cat("Bucket:", correct_bucket, "\n")

run_query <- function(sql) {
  tb <- bq_project_query(project, sql)
  bq_table_download(tb)
}

# Test connection - expected ~633,547
run_query(paste0("SELECT COUNT(*) AS n FROM `", cdr, ".person`"))

# Decision log
decision_log <- list()
log_decision <- function(variable, action, reason) {
  entry <- paste0("[DECISION] ", variable, ": ", action,
                  " | REASON: ", reason)
  cat("\n", entry, "\n")
  decision_log[[length(decision_log)+1]] <<- entry
}

log_decision("Design",
  "Cross-sectional descriptive. No causal claims.",
  paste0("All of Us Basics survey completed primarily in 2020. ",
         "EHR biomarkers pulled retrospectively to 1987. ",
         "53.6% of participants had hsCRP before Basics completion. ",
         "VPC findings describe distribution of inflammatory burden ",
         "across intersectional positions, not causal effects."))


# ============================================================
# SECTION 1: hsCRP
# ============================================================

crp_sql <- str_glue("
  SELECT m.person_id, m.value_as_number AS hsCRP, m.measurement_date
  FROM `{cdr}.measurement` m
  WHERE m.measurement_concept_id IN (3020460, 4145423, 37049645)
    AND m.value_as_number IS NOT NULL
    AND m.value_as_number > 0
    AND m.value_as_number < 100
  ORDER BY m.person_id, m.measurement_date DESC
")
crp_raw    <- run_query(crp_sql)
crp_latest <- crp_raw %>%
  arrange(person_id, desc(measurement_date)) %>%
  group_by(person_id) %>% slice(1) %>% ungroup()

cat("hsCRP N:", nrow(crp_latest), "\n")
print(summary(crp_latest$hsCRP))

# Distribution checks
par(mfrow=c(2,3))
hist(crp_latest$hsCRP, breaks=100, main="1. Raw hsCRP", xlab="mg/L")
p99_crp <- quantile(crp_latest$hsCRP, 0.99)
hist(log(pmin(crp_latest$hsCRP, p99_crp)+0.1), breaks=100,
     main="2. Log hsCRP winsorized", xlab="log(hsCRP_w+0.1)")
qqnorm(log(pmin(crp_latest$hsCRP, p99_crp)+0.1),
       main="3. QQ: Log hsCRP winsorized")
qqline(log(pmin(crp_latest$hsCRP, p99_crp)+0.1), col="red")

cat("99th percentile hsCRP:", round(p99_crp, 2), "mg/L\n")
cat("N winsorized:", sum(crp_latest$hsCRP > p99_crp), "\n")

crp_latest <- crp_latest %>%
  mutate(
    hsCRP_w     = pmin(hsCRP, p99_crp),
    log_hsCRP_w = log(hsCRP_w + 0.1)
  )

log_decision("hsCRP",
  paste0("Winsorize at 99th percentile (", round(p99_crp,1),
         " mg/L) then log(hsCRP_w + 0.1)"),
  paste0("CRP is right-skewed by nature. Log transformation is ",
         "established standard. +0.1 handles near-zero values. ",
         "Winsorizing at 99th percentile reduces influence of acute ",
         "infection values. Visual QQ shows acceptable normality."))

# ============================================================
# NOTE ON ONE OBSERVATION PER PARTICIPANT
# ============================================================
# For all biomarkers (hsCRP, HbA1c, BMI), this script selects
# the most recent available measurement per participant using
# arrange(desc(measurement_date)) %>% group_by(person_id) %>% slice(1).
# Each participant contributes exactly one measurement.
# No date-proximity restriction is applied between biomarkers
# or between biomarkers and the Basics survey date.
# The Basics survey was completed primarily in 2020-2021.
# EHR biomarker data are retrospective to 1987.
# 53.6% of participants had hsCRP recorded before Basics completion.
# This cross-sectional design is documented in Analytic Decision #1
# and reported in Table S5. All findings are described as
# distributional rather than causal.

# Report hsCRP measurement date distribution for transparency
cat("\n=== hsCRP MEASUREMENT DATE DISTRIBUTION ===\n")
cat("Earliest measurement:", as.character(min(crp_latest$measurement_date)), "\n")
cat("Latest measurement:", as.character(max(crp_latest$measurement_date)), "\n")
cat("Median measurement date:", as.character(median(crp_latest$measurement_date)), "\n")
cat("% measured 2015 or later:",
    round(mean(crp_latest$measurement_date >= as.Date("2015-01-01"))*100, 1), "%\n")
cat("% measured 2010 or later:",
    round(mean(crp_latest$measurement_date >= as.Date("2010-01-01"))*100, 1), "%\n")
cat("% measured 2000 or later:",
    round(mean(crp_latest$measurement_date >= as.Date("2000-01-01"))*100, 1), "%\n")


# ============================================================
# SECTION 2: DEMOGRAPHICS
# ============================================================

race_sql <- str_glue("
  SELECT o.person_id,
    STRING_AGG(c.concept_name ORDER BY c.concept_name) AS race_all
  FROM `{cdr}.observation` o
  JOIN `{cdr}.concept` c ON o.value_source_concept_id = c.concept_id
  WHERE o.observation_source_concept_id = 1586140
  GROUP BY person_id
")
race_df <- run_query(race_sql)

gender_sql <- str_glue("
  SELECT DISTINCT o.person_id, c.concept_name AS gender_identity
  FROM `{cdr}.observation` o
  JOIN `{cdr}.concept` c ON o.value_source_concept_id = c.concept_id
  WHERE o.observation_source_concept_id = 1585838
")
gender_df <- run_query(gender_sql) %>%
  group_by(person_id) %>% slice(1) %>% ungroup()

income_sql <- str_glue("
  SELECT DISTINCT o.person_id, c.concept_name AS income_category
  FROM `{cdr}.observation` o
  JOIN `{cdr}.concept` c ON o.value_source_concept_id = c.concept_id
  WHERE o.observation_source_concept_id = 1585375
")
income_df <- run_query(income_sql) %>%
  group_by(person_id) %>% slice(1) %>% ungroup()

cat("\n=== VERIFY BEFORE RECODING ===\n")
cat("Gender values:\n")
print(gender_df %>% count(gender_identity, sort=TRUE))
cat("\nIncome values:\n")
print(income_df %>% count(income_category, sort=TRUE))

demo_sql <- str_glue("
  SELECT p.person_id, p.birth_datetime,
    p.sex_at_birth_concept_id, p.ethnicity_concept_id
  FROM `{cdr}.person` p
")
demo_raw <- run_query(demo_sql)
demo <- demo_raw %>%
  mutate(
    birth_date = as.Date(birth_datetime),
    age_2024   = as.numeric(floor(
      (as.Date("2024-07-01") - birth_date) / 365.25)),
    ethnicity  = case_when(
      ethnicity_concept_id == 38003563 ~ "Hispanic/Latino",
      TRUE                             ~ "Not Hispanic/Latino")
  ) %>%
  dplyr::select(person_id, age_2024, ethnicity)

cat("\nAge distribution:\n")
print(summary(demo$age_2024))
cat("N age > 90:", sum(demo$age_2024 > 90, na.rm=TRUE), "\n")


# ============================================================
# SECTION 3: COVARIATES
# ============================================================

# Zip code
zip_sql <- str_glue("
  SELECT DISTINCT person_id, value_as_string AS zip3
  FROM `{cdr}.observation`
  WHERE observation_source_concept_id = 1585250
")
zip_df <- run_query(zip_sql) %>%
  group_by(person_id) %>% slice(1) %>% ungroup()

# HbA1c
hba1c_sql <- str_glue("
  SELECT person_id, value_as_number AS hba1c, measurement_date
  FROM `{cdr}.measurement`
  WHERE measurement_concept_id IN (3004410, 40789263, 3007332)
    AND value_as_number BETWEEN 3 AND 20
  ORDER BY person_id, measurement_date DESC
")
hba1c_raw    <- run_query(hba1c_sql)
hba1c_latest <- hba1c_raw %>%
  arrange(person_id, desc(measurement_date)) %>%
  group_by(person_id) %>% slice(1) %>% ungroup()

p99_hba1c <- quantile(hba1c_latest$hba1c, 0.99)
hba1c_latest <- hba1c_latest %>%
  mutate(
    hba1c_w     = pmin(hba1c, p99_hba1c),
    log_hba1c_w = log(hba1c_w),
    diabetes    = as.integer(hba1c >= 6.5)
  )
cat("HbA1c N:", nrow(hba1c_latest), "\n")
cat("Diabetes prevalence:",
    round(mean(hba1c_latest$diabetes)*100, 1), "%\n")

log_decision("HbA1c",
  paste0("Winsorize at 99th percentile (", round(p99_hba1c,2),
         ") then log-transform. Binary diabetes >= 6.5% as ",
         "primary replication outcome."),
  paste0("HbA1c is bimodal in mixed diabetic/non-diabetic sample. ",
         "Binary diabetes indicator resolves bimodality and has ",
         "direct clinical meaning."))

# BMI
bmi_sql <- str_glue("
  SELECT person_id, value_as_number AS bmi, measurement_date
  FROM `{cdr}.measurement`
  WHERE measurement_concept_id IN (3038553, 40490382)
    AND value_as_number BETWEEN 10 AND 80
  ORDER BY person_id, measurement_date DESC
")
bmi_raw    <- run_query(bmi_sql)
bmi_latest <- bmi_raw %>%
  arrange(person_id, desc(measurement_date)) %>%
  group_by(person_id) %>% slice(1) %>% ungroup()

cat("\nBMI before cleaning:\n")
print(summary(bmi_latest$bmi))
cat("BMI < 15:", sum(bmi_latest$bmi < 15), "\n")
cat("BMI > 60:", sum(bmi_latest$bmi > 60), "\n")

bmi_clean <- bmi_latest %>%
  dplyr::select(person_id, bmi) %>%
  mutate(bmi = ifelse(bmi < 15 | bmi > 60, NA, bmi))

p99_bmi <- quantile(bmi_clean$bmi, 0.99, na.rm=TRUE)
bmi_clean <- bmi_clean %>%
  mutate(bmi_w = pmin(bmi, p99_bmi))

cat("BMI after cleaning:\n")
print(summary(bmi_clean$bmi_w))

par(mfrow=c(2,2))
hist(bmi_clean$bmi_w, breaks=80, main="BMI winsorized", xlab="kg/m2")
hist(log(bmi_clean$bmi_w), breaks=80, main="Log BMI", xlab="log(BMI)")
qqnorm(bmi_clean$bmi_w, main="QQ: BMI")
qqline(bmi_clean$bmi_w, col="red")
qqnorm(log(bmi_clean$bmi_w), main="QQ: Log BMI")
qqline(log(bmi_clean$bmi_w), col="red")

log_decision("BMI",
  paste0("Exclude <15 or >60 as implausible EHR errors. ",
         "Winsorize at 99th percentile (", round(p99_bmi,1),
         "). Log-transform."),
  "Right-skewed in population data. Log transformation substantially improves QQ fit.")

# Smoking
smoking_sql <- str_glue("
  SELECT DISTINCT o.person_id, c.concept_name AS smoking_status
  FROM `{cdr}.observation` o
  JOIN `{cdr}.concept` c ON o.value_source_concept_id = c.concept_id
  WHERE o.observation_source_concept_id IN (1585857,1585858,1585864,1586159)
")
smoking_raw <- run_query(smoking_sql)
cat("\nSmoking values:\n")
print(smoking_raw %>% count(smoking_status, sort=TRUE))

smoking_clean <- smoking_raw %>%
  group_by(person_id) %>% slice(1) %>% ungroup() %>%
  mutate(
    ever_smoker = case_when(
      smoking_status == "100 Cigs Lifetime: Yes" ~ 1L,
      smoking_status == "100 Cigs Lifetime: No"  ~ 0L,
      TRUE                                        ~ NA_integer_)
  ) %>%
  dplyr::select(person_id, ever_smoker)

log_decision("Smoking",
  "Binary ever/never smoker. PMI responses = NA.",
  paste0("All of Us measures lifetime cigarette exposure, not current ",
         "smoking. Does not distinguish current from former smokers. ",
         "Limitation noted in paper."))

# Education
edu_sql <- str_glue("
  SELECT DISTINCT o.person_id, c.concept_name AS education
  FROM `{cdr}.observation` o
  JOIN `{cdr}.concept` c ON o.value_source_concept_id = c.concept_id
  WHERE o.observation_source_concept_id = 1585940
")
edu_raw <- run_query(edu_sql)
cat("\nEducation values:\n")
print(edu_raw %>% count(education, sort=TRUE))

edu_clean <- edu_raw %>%
  group_by(person_id) %>% slice(1) %>% ungroup() %>%
  mutate(
    education3 = case_when(
      education %in% c(
        "Highest Grade: Never Attended",
        "Highest Grade: One Through Four",
        "Highest Grade: Five Through Eight",
        "Highest Grade: Nine Through Eleven",
        "Highest Grade: Twelve Or GED")         ~ "Low",
      education == "Highest Grade: College One to Three" ~ "Middle",
      education %in% c(
        "Highest Grade: College Graduate",
        "Highest Grade: Advanced Degree")       ~ "High",
      TRUE                                      ~ NA_character_)
  ) %>%
  dplyr::select(person_id, education3)


# ============================================================
# SECTION 4: MERGE AND MISSING DATA AUDIT
# ============================================================

analytic_raw <- crp_latest %>%
  left_join(demo,          by="person_id") %>%
  left_join(race_df,       by="person_id") %>%
  left_join(gender_df %>% dplyr::select(person_id, gender_identity),
            by="person_id") %>%
  left_join(income_df %>% dplyr::select(person_id, income_category),
            by="person_id") %>%
  left_join(zip_df,        by="person_id") %>%
  left_join(hba1c_latest %>%
              dplyr::select(person_id, hba1c, hba1c_w,
                            log_hba1c_w, diabetes),
            by="person_id") %>%
  left_join(bmi_clean %>% dplyr::select(person_id, bmi_w),
            by="person_id") %>%
  left_join(smoking_clean, by="person_id") %>%
  left_join(edu_clean,     by="person_id")

cat("\n=== MISSING DATA AUDIT ===\n")
vars_audit <- c("race_all","gender_identity","income_category",
                "zip3","age_2024","bmi_w","ever_smoker",
                "education3","hba1c")
for (v in vars_audit) {
  n_miss <- sum(is.na(analytic_raw[[v]]))
  pct    <- round(n_miss/nrow(analytic_raw)*100, 1)
  cat(sprintf("  %-22s missing: %6d (%4.1f%%)\n", v, n_miss, pct))
}

log_decision("Missing data",
  "Complete case analysis. No multiple imputation.",
  paste0("Missing data ranges 3-11% across covariates. ",
         "Complete case analysis is standard for MAIHDA. ",
         "Missingness pattern would likely bias toward null."))


# ============================================================
# SECTION 5: RECODE INTERSECTIONAL VARIABLES
# ============================================================
# IMPORTANT: Print unique values before running recodes.
# CDR version changes can alter string labels silently.

analytic_raw <- analytic_raw %>%
  mutate(
    # Race: Hispanic takes priority over race question
    race5 = case_when(
      ethnicity == "Hispanic/Latino"                             ~ "Hispanic",
      str_detect(coalesce(race_all,""), "Black|African")         ~ "Black",
      str_detect(coalesce(race_all,""), "White") &
        !str_detect(coalesce(race_all,""), ",")                  ~ "White",
      str_detect(coalesce(race_all,""), "Asian") &
        !str_detect(coalesce(race_all,""), ",")                  ~ "Asian",
      TRUE                                                        ~ "Other"),
    # Asian + Other -> Other/Asian (small cells)
    race5c = case_when(
      race5 %in% c("Asian","Other") ~ "Other/Asian",
      TRUE                          ~ race5),
    # Gender identity
    # VERIFY THESE STRINGS MATCH CDR OUTPUT ABOVE
    gender3 = case_when(
      is.na(gender_identity)                               ~ NA_character_,
      gender_identity == "Gender Identity: Woman"          ~ "Cisgender Woman",
      gender_identity == "Gender Identity: Man"            ~ "Cisgender Man",
      gender_identity %in% c(
        "Gender Identity: Non Binary",
        "Gender Identity: Transgender",
        "Gender Identity: Additional Options")             ~ "Gender Minority",
      TRUE                                                 ~ NA_character_),
    # Income
    # NOTE: CDR v8 uses full string format "Annual Income: less 10k"
    # NOT short format "less_10k". Both primary and PSS analyses
    # now use the same full string format.
    # PMI: Skip, Prefer Not To Answer, Dont Know -> NA
    income3 = case_when(
      is.na(income_category)                               ~ NA_character_,
      income_category %in% c(
        "PMI: Skip",
        "PMI: Prefer Not To Answer",
        "PMI: Dont Know")                                  ~ NA_character_,
      income_category %in% c(
        "Annual Income: less 10k",
        "Annual Income: 10k 25k",
        "Annual Income: 25k 35k")                         ~ "Low",
      income_category %in% c(
        "Annual Income: 35k 50k",
        "Annual Income: 50k 75k",
        "Annual Income: 75k 100k")                        ~ "Middle",
      income_category %in% c(
        "Annual Income: 100k 150k",
        "Annual Income: 150k 200k",
        "Annual Income: more 200k")                       ~ "High",
      TRUE                                                 ~ NA_character_)
  )

# Verify no unmapped values
cat("\nUnmapped gender (should be 0 or PMI only):\n")
print(analytic_raw %>%
        filter(is.na(gender3), !is.na(gender_identity)) %>%
        count(gender_identity, sort=TRUE) %>% head(5))

cat("\nUnmapped income (should be 0 or PMI only):\n")
print(analytic_raw %>%
        filter(is.na(income3), !is.na(income_category)) %>%
        count(income_category, sort=TRUE) %>% head(5))


# ============================================================
# SECTION 6: INTERSECTIONAL POSITIONS
# ============================================================

analytic_raw <- analytic_raw %>%
  mutate(
    # Gender minority collapsed across race (cells too small)
    race_for_position = case_when(
      gender3 == "Gender Minority" ~ "All Races",
      TRUE                         ~ race5c),
    position_label_c = paste(race_for_position, gender3,
                              income3, sep=" | ")
  )

cat("\n=== CELL SIZE CHECK (pre-filter) ===\n")
cells_pre <- analytic_raw %>%
  filter(!is.na(race5c), !is.na(gender3), !is.na(income3)) %>%
  count(position_label_c, sort=TRUE)
cat("Positions:", nrow(cells_pre), "\n")
cat("Min cell:", min(cells_pre$n), "\n")
cat("Cells < 30:", sum(cells_pre$n < 30), "\n")

log_decision("Intersectional positions",
  "27 positions: race5c(4) x gender3(3) x income3(3). Gender minority collapsed across race.",
  paste0("Gender minority race-specific cells fall below N=30 minimum ",
         "for stable MAIHDA estimates. Collapsing yields 27 estimable ",
         "positions with minimum cell size 63 in adjusted sample."))


# ============================================================
# SECTION 7: BUILD ANALYTIC SAMPLES
# ============================================================

analytic_raw <- analytic_raw %>%
  mutate(
    age_w   = pmin(age_2024, quantile(age_2024, 0.99, na.rm=TRUE)),
    age_c   = age_w - mean(age_w, na.rm=TRUE),
    log_bmi = log(bmi_w)
  )

# Unadjusted sample
analytic_clean <- analytic_raw %>%
  filter(!is.na(log_hsCRP_w), !is.na(race5c),
         !is.na(gender3), !is.na(income3)) %>%
  mutate(position_id = as.integer(factor(position_label_c)))

# Covariate-adjusted sample
analytic_adj <- analytic_raw %>%
  filter(!is.na(log_hsCRP_w), !is.na(race5c), !is.na(gender3),
         !is.na(income3), !is.na(age_w), !is.na(bmi_w),
         !is.na(ever_smoker), !is.na(education3)) %>%
  mutate(position_id_adj = as.integer(factor(position_label_c)))

cat("Unadjusted N:", nrow(analytic_clean), "\n")
cat("Adjusted N:", nrow(analytic_adj), "\n")

cat("\n=== CELL SIZE CHECK (adjusted) ===\n")
cells_adj <- analytic_adj %>% count(position_label_c, sort=TRUE)
cat("Min cell:", min(cells_adj$n), "\n")
cat("Cells < 30:", sum(cells_adj$n < 30), "\n")
print(cells_adj)

# Covariate correlation check
cat("\n=== COVARIATE CORRELATIONS ===\n")
cor_check <- analytic_adj %>%
  dplyr::select(log_hsCRP_w, age_w, log_bmi, ever_smoker) %>%
  filter(complete.cases(.))
print(round(cor(cor_check), 3))

# Three-column Table 1 group variable
analytic_adj <- analytic_adj %>%
  mutate(crp_group = ifelse(hsCRP > 3.0, "Elevated", "Normal"))


# ============================================================
# SECTION 8: NEIGHBORHOOD DEPRIVATION INDEX
# ============================================================

dep_index_zip <- read.csv(
  'https://github.com/geomarker-io/dep_index/raw/master/2023/data/ACS_deprivation_index_by_zipcode.csv'
)
cat("Deprivation index rows:", nrow(dep_index_zip), "\n")

dep_index_zip3 <- dep_index_zip %>%
  mutate(zip3 = str_pad(substr(as.character(zcta_2020),1,3),3,pad="0")) %>%
  group_by(zip3) %>%
  summarise(
    dep_index_mean     = mean(dep_index, na.rm=TRUE),
    pct_poverty_mean   = mean(fraction_poverty, na.rm=TRUE),
    median_income_mean = mean(median_income, na.rm=TRUE),
    n_zctas            = n(),
    .groups="drop")

analytic_adj <- analytic_adj %>%
  mutate(zip3_clean = str_pad(substr(as.character(zip3),1,3),3,pad="0")) %>%
  left_join(dep_index_zip3, by=c("zip3_clean"="zip3"))

cat("ADI matched:", sum(!is.na(analytic_adj$dep_index_mean)), "\n")
cat("Match rate:", round(mean(!is.na(analytic_adj$dep_index_mean))*100,1),"%\n")
cat("Correlation dep_index ~ log_hsCRP_w:",
    round(cor(analytic_adj$dep_index_mean, analytic_adj$log_hsCRP_w,
              use="complete.obs"), 3), "\n")

zip3_sizes <- analytic_adj %>%
  filter(!is.na(zip3_clean)) %>%
  count(zip3_clean, sort=TRUE)

analytic_adj <- analytic_adj %>%
  left_join(zip3_sizes %>% rename(zip3_n=n), by="zip3_clean") %>%
  mutate(
    zip3_valid = !is.na(zip3_clean) & zip3_clean != "000",
    zip3_id    = as.integer(factor(ifelse(zip3_valid, zip3_clean, NA))))

log_decision("Deprivation index",
  "geomarker-io 2023 ACS index averaged within 3-digit zip. Geographic models N>=30 per zip3.",
  paste0("ADI from Neighborhood Atlas excluded - not validated at ",
         "3-digit zip level. geomarker-io index available directly ",
         "from GitHub without download. Correlation with hsCRP r~0.15 ",
         "confirms construct validity."))


# ============================================================
# SECTION 9: TABLE 1 - THREE COLUMN SAMPLE CHARACTERISTICS
# ============================================================

cat("\n=================================================\n")
cat("TABLE 1: SAMPLE CHARACTERISTICS\n")
cat("Overall | Normal CRP (<=3.0) | Elevated CRP (>3.0)\n")
cat("=================================================\n")

summarize_group <- function(df, label) {
  df %>%
    summarise(
      label          = label,
      N              = n(),
      age_mean       = round(mean(age_w, na.rm=TRUE), 1),
      age_sd         = round(sd(age_w, na.rm=TRUE), 1),
      bmi_mean       = round(mean(bmi_w, na.rm=TRUE), 1),
      bmi_sd         = round(sd(bmi_w, na.rm=TRUE), 1),
      pct_cis_woman  = round(mean(gender3=="Cisgender Woman")*100, 1),
      pct_cis_man    = round(mean(gender3=="Cisgender Man")*100, 1),
      pct_gm         = round(mean(gender3=="Gender Minority")*100, 1),
      pct_white      = round(mean(race5c=="White")*100, 1),
      pct_black      = round(mean(race5c=="Black")*100, 1),
      pct_hispanic   = round(mean(race5c=="Hispanic")*100, 1),
      pct_other      = round(mean(race5c=="Other/Asian")*100, 1),
      pct_inc_low    = round(mean(income3=="Low")*100, 1),
      pct_inc_mid    = round(mean(income3=="Middle")*100, 1),
      pct_inc_high   = round(mean(income3=="High")*100, 1),
      pct_edu_low    = round(mean(education3=="Low")*100, 1),
      pct_edu_mid    = round(mean(education3=="Middle")*100, 1),
      pct_edu_high   = round(mean(education3=="High")*100, 1),
      pct_smoker     = round(mean(ever_smoker==1)*100, 1),
      hsCRP_median   = round(median(hsCRP, na.rm=TRUE), 2),
      hsCRP_q25      = round(quantile(hsCRP, 0.25, na.rm=TRUE), 2),
      hsCRP_q75      = round(quantile(hsCRP, 0.75, na.rm=TRUE), 2),
      pct_elevated   = round(mean(hsCRP > 3.0)*100, 1),
      pct_high_risk  = round(mean(hsCRP > 10.0)*100, 1)
    )
}

t1_overall  <- summarize_group(analytic_adj, "Overall")
t1_normal   <- summarize_group(
  analytic_adj %>% filter(crp_group=="Normal"), "Normal CRP")
t1_elevated <- summarize_group(
  analytic_adj %>% filter(crp_group=="Elevated"), "Elevated CRP")

t1_combined <- bind_rows(t1_overall, t1_normal, t1_elevated) %>%
  t() %>% as.data.frame()
colnames(t1_combined) <- c("Overall","Normal CRP (<=3)","Elevated CRP (>3)")

print(t1_combined)

# Group difference tests
cat("\n=== GROUP DIFFERENCE TESTS (all p<0.0001 expected) ===\n")
cat("Age:", round(t.test(age_w~crp_group, data=analytic_adj)$p.value,4),"\n")
cat("BMI:", round(t.test(bmi_w~crp_group, data=analytic_adj)$p.value,4),"\n")
cat("Gender:", round(chisq.test(table(analytic_adj$gender3,
  analytic_adj$crp_group))$p.value,4),"\n")
cat("Race:", round(chisq.test(table(analytic_adj$race5c,
  analytic_adj$crp_group))$p.value,4),"\n")
cat("Income:", round(chisq.test(table(analytic_adj$income3,
  analytic_adj$crp_group))$p.value,4),"\n")
cat("Education:", round(chisq.test(table(analytic_adj$education3,
  analytic_adj$crp_group))$p.value,4),"\n")
cat("Smoking:", round(chisq.test(table(analytic_adj$ever_smoker,
  analytic_adj$crp_group))$p.value,4),"\n")

write_csv(t1_combined %>% rownames_to_column("Characteristic"),
          "IRISE_Table1_threecol.csv")


# ============================================================
# SECTION 10: MAIHDA MODELS - hsCRP
# ============================================================
# NOTE: REML used throughout for VPC estimation per Evans et al. 2018.
# REML produces unbiased variance component estimates.
# Cross-sectional descriptive framing - no causal claims.

get_vpc <- function(model, grp_name) {
  vc <- as.data.frame(VarCorr(model))
  rownames(vc) <- vc$grp
  vc[grp_name,"vcov"] / sum(vc$vcov)
}

# Model 1: Null unadjusted
model_null <- lmer(
  log_hsCRP_w ~ 1 + (1 | position_id),
  data=analytic_clean, REML=TRUE)
vpc_null <- get_vpc(model_null, "position_id")
cat("Model 1 VPC:", round(vpc_null*100, 2), "%\n")

# Model 2: Additive unadjusted
model_additive <- lmer(
  log_hsCRP_w ~ race5c + gender3 + income3 + (1 | position_id),
  data=analytic_clean, REML=TRUE)
vpc_additive <- get_vpc(model_additive, "position_id")
cat("Model 2 VPC:", round(vpc_additive*100, 2), "%\n")
cat("Additive explanation:",
    round((vpc_null-vpc_additive)/vpc_null*100, 1), "%\n")

# Model 3: Null covariate-adjusted
model_adj_null <- lmer(
  log_hsCRP_w ~ age_c + log_bmi + ever_smoker + education3 +
               (1 | position_id_adj),
  data=analytic_adj, REML=TRUE)
vpc_adj_null <- get_vpc(model_adj_null, "position_id_adj")
cat("Model 3 VPC:", round(vpc_adj_null*100, 2), "%\n")

# Model 4: Additive covariate-adjusted
model_adj_add <- lmer(
  log_hsCRP_w ~ race5c + gender3 + income3 +
                age_c + log_bmi + ever_smoker + education3 +
                (1 | position_id_adj),
  data=analytic_adj, REML=TRUE)
vpc_adj_add <- get_vpc(model_adj_add, "position_id_adj")
cat("Model 4 VPC:", round(vpc_adj_add*100, 2), "%\n")
cat("Additive explanation (adjusted):",
    round((vpc_adj_null-vpc_adj_add)/vpc_adj_null*100, 1), "%\n")

cat("\nFixed effects (Model 4):\n")
print(round(summary(model_adj_add)$coefficients, 4))

# Model 5: Three-level adjusted
analytic_adj_zip <- analytic_adj %>%
  filter(zip3_valid, zip3_n >= 30, !is.na(dep_index_mean))
cat("N three-level:", nrow(analytic_adj_zip), "\n")

model_adj_threelevel <- lmer(
  log_hsCRP_w ~ race5c + gender3 + income3 +
                age_c + log_bmi + ever_smoker + education3 +
                dep_index_mean +
                (1 | position_id_adj) + (1 | zip3_id),
  data=analytic_adj_zip, REML=TRUE)

vc_3l <- as.data.frame(VarCorr(model_adj_threelevel))
rownames(vc_3l) <- vc_3l$grp
total_3l <- sum(vc_3l$vcov)
vpc_adj_zip3 <- vc_3l["zip3_id",        "vcov"] / total_3l
vpc_adj_pos  <- vc_3l["position_id_adj", "vcov"] / total_3l

cat("Model 5 Geographic VPC:", round(vpc_adj_zip3*100, 1), "%\n")
cat("Model 5 Position VPC:  ", round(vpc_adj_pos*100,  1), "%\n")


# ============================================================
# SECTION 11: TABLE 2 - POSITION-LEVEL hsCRP
# ============================================================

cat("\n=================================================\n")
cat("TABLE 2: POSITION-LEVEL hsCRP RANKINGS\n")
cat("=================================================\n")

position_results <- analytic_clean %>%
  group_by(position_label_c) %>%
  summarise(
    n              = n(),
    mean_hsCRP     = round(mean(hsCRP, na.rm=TRUE), 2),
    median_hsCRP   = round(median(hsCRP, na.rm=TRUE), 2),
    mean_log_hsCRP = round(mean(log_hsCRP_w, na.rm=TRUE), 3),
    se_log_hsCRP   = round(sd(log_hsCRP_w,na.rm=TRUE)/sqrt(n()), 4),
    .groups="drop"
  ) %>%
  arrange(desc(median_hsCRP)) %>%
  mutate(rank = row_number())

print(position_results, n=27)
cat("\nMedian range:",
    round(min(position_results$median_hsCRP),2), "to",
    round(max(position_results$median_hsCRP),2), "mg/L\n")
cat("Fold difference:",
    round(max(position_results$median_hsCRP)/
          min(position_results$median_hsCRP), 1), "-fold\n")

write_csv(position_results, "IRISE_Table2_positions_hsCRP.csv")


# ============================================================
# SECTION 12: TABLE 3 - VPC DECOMPOSITION
# ============================================================

cat("\n=================================================\n")
cat("TABLE 3: VPC DECOMPOSITION\n")
cat("=================================================\n")

vpc_table <- data.frame(
  Model = c("1: Null (unadjusted)",
            "2: Additive (unadjusted)",
            "3: Null (covariate-adjusted)",
            "4: Additive (covariate-adjusted)",
            "5: Three-level adjusted"),
  Position_VPC_pct = round(c(vpc_null, vpc_additive, vpc_adj_null,
                              vpc_adj_add, vpc_adj_pos)*100, 2),
  Geographic_VPC_pct = c(NA,NA,NA,NA, round(vpc_adj_zip3*100,2)),
  Additive_Explanation_pct = c(
    NA,
    round((vpc_null-vpc_additive)/vpc_null*100, 1),
    NA,
    round((vpc_adj_null-vpc_adj_add)/vpc_adj_null*100, 1),
    NA),
  N = c(nrow(analytic_clean), nrow(analytic_clean),
        nrow(analytic_adj), nrow(analytic_adj),
        nrow(analytic_adj_zip))
)

print(vpc_table)
write_csv(vpc_table, "IRISE_Table3_VPC.csv")


# ============================================================
# SECTION 13: TABLE 4 - FIXED EFFECTS
# ============================================================

cat("\n=================================================\n")
cat("TABLE 4: FIXED EFFECTS - ADJUSTED ADDITIVE MODEL\n")
cat("Reference: Black, Gender Minority, High income, High education\n")
cat("=================================================\n")

fe_df <- as.data.frame(summary(model_adj_add)$coefficients) %>%
  rownames_to_column("Predictor") %>%
  mutate(across(where(is.numeric), ~round(., 4)))
print(fe_df)
write_csv(fe_df, "IRISE_Table4_fixed_effects.csv")


# ============================================================
# SECTION 14: HbA1c AND DIABETES REPLICATION
# ============================================================

analytic_hba1c <- analytic_adj %>%
  filter(!is.na(hba1c)) %>%
  mutate(position_id_hba1c = as.integer(factor(position_label_c)))
cat("HbA1c N:", nrow(analytic_hba1c), "\n")

model_hba1c_null <- lmer(
  log_hba1c_w ~ 1 + (1 | position_id_hba1c),
  data=analytic_hba1c, REML=TRUE)
vpc_hba1c <- get_vpc(model_hba1c_null, "position_id_hba1c")

model_hba1c_add <- lmer(
  log_hba1c_w ~ race5c + gender3 + income3 +
                age_c + log_bmi + ever_smoker + education3 +
                (1 | position_id_hba1c),
  data=analytic_hba1c, REML=TRUE)
vpc_hba1c_add <- get_vpc(model_hba1c_add, "position_id_hba1c")

model_diabetes_null <- glmer(
  diabetes ~ 1 + (1 | position_id_hba1c),
  data=analytic_hba1c, family=binomial(link="logit"))
vc_d <- as.data.frame(VarCorr(model_diabetes_null))
vpc_diabetes <- vc_d$vcov[1] / (vc_d$vcov[1] + (pi^2/3))

cat("HbA1c null VPC:", round(vpc_hba1c*100, 2), "%\n")
cat("HbA1c additive explained:",
    round((vpc_hba1c-vpc_hba1c_add)/vpc_hba1c*100, 1), "%\n")
cat("Diabetes VPC:", round(vpc_diabetes*100, 2), "%\n")

diabetes_by_position <- analytic_hba1c %>%
  group_by(position_label_c) %>%
  summarise(
    n            = n(),
    pct_diabetes = round(mean(hba1c >= 6.5)*100, 1),
    .groups="drop") %>%
  arrange(desc(pct_diabetes))

cat("\nDiabetes range:",
    round(min(diabetes_by_position$pct_diabetes),1), "% to",
    round(max(diabetes_by_position$pct_diabetes),1), "%\n")
cat("Fold difference:",
    round(max(diabetes_by_position$pct_diabetes)/
          min(diabetes_by_position$pct_diabetes), 1), "-fold\n")

write_csv(diabetes_by_position, "IRISE_position_diabetes.csv")


# ============================================================
# SECTION 15: TABLE 5 - REPLICATION SUMMARY
# ============================================================

cat("\n=================================================\n")
cat("TABLE 5: REPLICATION ACROSS OUTCOMES\n")
cat("=================================================\n")

replication_table <- data.frame(
  Outcome = c("hsCRP (log, winsorized)",
              "HbA1c (log, winsorized)",
              "Diabetes (binary, >= 6.5%)"),
  N = c(nrow(analytic_clean), nrow(analytic_hba1c),
        nrow(analytic_hba1c)),
  Null_VPC_pct = round(c(vpc_null, vpc_hba1c,
                          vpc_diabetes)*100, 2),
  Additive_Explanation_pct = round(c(
    (vpc_null-vpc_additive)/vpc_null*100,
    (vpc_hba1c-vpc_hba1c_add)/vpc_hba1c*100,
    NA), 1),
  Position_Range = c(
    paste0("Median ",
           round(min(position_results$median_hsCRP),2), "-",
           round(max(position_results$median_hsCRP),2), " mg/L"),
    "See supplementary",
    paste0(round(min(diabetes_by_position$pct_diabetes),1), "%-",
           round(max(diabetes_by_position$pct_diabetes),1), "%"))
)

print(replication_table)
write_csv(replication_table, "IRISE_Table5_replication.csv")


# ============================================================
# SECTION 16: PSS PATHWAY - CROSS-SECTIONAL DESCRIPTIVE EVIDENCE
# ============================================================
# NOTE: Cross-sectional only. 73.9% of PSS+hsCRP overlap had
# hsCRP before PSS completion. Formal mediation not defensible
# in main paper. PSS findings reported as pathway-consistent
# evidence consistent with IPNE framework.
# Temporal ordering: Basics/PSS completed primarily 2020-2021.
# EHR hsCRP pulled retrospectively to 1987.

cat("\n=================================================\n")
cat("SECTION 16: PSS CROSS-SECTIONAL PATHWAY EVIDENCE\n")
cat("NOTE: Descriptive only. No causal claims.\n")
cat("=================================================\n")

# Pull PSS
pss_full_sql <- str_glue("
  SELECT o.person_id,
    o.observation_source_value AS item,
    c.concept_name AS response
  FROM `{cdr}.observation` o
  LEFT JOIN `{cdr}.concept` c
    ON o.value_source_concept_id = c.concept_id
  WHERE o.observation_source_value LIKE 'sdoh_cpss_%'
    AND o.value_source_value NOT LIKE 'PMI%'
    AND c.concept_name IS NOT NULL
")
pss_full <- run_query(pss_full_sql)
cat("PSS rows:", nrow(pss_full), "\n")
cat("Unique participants:", n_distinct(pss_full$person_id), "\n")

# Score PSS-10
# Items 4,5,7,8 are reverse scored
reverse_items <- c("sdoh_cpss_4","sdoh_cpss_5",
                   "sdoh_cpss_7","sdoh_cpss_8")

pss_scored <- pss_full %>%
  mutate(
    score_raw = case_when(
      response == "Never"        ~ 0,
      response == "Almost Never" ~ 1,
      response == "Sometimes"    ~ 2,
      response == "Fairly Often" ~ 3,
      response == "Very Often"   ~ 4,
      TRUE                       ~ NA_real_),
    score = case_when(
      item %in% reverse_items ~ 4 - score_raw,
      TRUE                    ~ score_raw)
  ) %>%
  group_by(person_id) %>%
  summarise(
    n_items   = sum(!is.na(score)),
    pss_total = sum(score, na.rm=TRUE),
    .groups="drop") %>%
  filter(n_items >= 8) %>%
  mutate(pss_total = pss_total * (10/n_items))

cat("PSS scored N:", nrow(pss_scored), "\n")
cat("PSS mean:", round(mean(pss_scored$pss_total), 2), "\n")
cat("PSS SD:", round(sd(pss_scored$pss_total), 2), "\n")

# Merge PSS with hsCRP
crp_pss <- crp_latest %>%
  inner_join(pss_scored %>% dplyr::select(person_id, pss_total),
             by="person_id")
cat("PSS + hsCRP N:", nrow(crp_pss), "\n")

# Build PSS analytic dataset with intersectional positions
# NOTE: Income strings differ in this CDR pull - use full strings
analytic_pss <- crp_pss %>%
  left_join(demo %>% dplyr::select(person_id, age_2024), by="person_id") %>%
  left_join(gender_df %>% dplyr::select(person_id, gender_identity),
            by="person_id") %>%
  left_join(income_df %>% dplyr::select(person_id, income_category),
            by="person_id") %>%
  left_join(race_df, by="person_id") %>%
  mutate(
    ethnicity = case_when(
      str_detect(coalesce(race_all,""), "Hispanic") ~ "Hispanic/Latino",
      TRUE ~ "Not Hispanic/Latino"),
    race5c = case_when(
      str_detect(coalesce(race_all,""), "Hispanic")   ~ "Hispanic",
      str_detect(coalesce(race_all,""), "Black|African") ~ "Black",
      str_detect(coalesce(race_all,""), "White") &
        !str_detect(coalesce(race_all,""), ",")        ~ "White",
      TRUE                                              ~ "Other/Asian"),
    gender3 = case_when(
      gender_identity == "Gender Identity: Woman"      ~ "Cisgender Woman",
      gender_identity == "Gender Identity: Man"        ~ "Cisgender Man",
      gender_identity %in% c(
        "Gender Identity: Non Binary",
        "Gender Identity: Transgender",
        "Gender Identity: Additional Options")         ~ "Gender Minority",
      TRUE                                             ~ NA_character_),
    # FULL STRING FORMAT - same as primary analysis Section 5
    # PMI: Skip, Prefer Not To Answer, Dont Know -> NA
    income3 = case_when(
      income_category %in% c(
        "PMI: Skip",
        "PMI: Prefer Not To Answer",
        "PMI: Dont Know")           ~ NA_character_,
      income_category %in% c(
        "Annual Income: less 10k",
        "Annual Income: 10k 25k",
        "Annual Income: 25k 35k")   ~ "Low",
      income_category %in% c(
        "Annual Income: 35k 50k",
        "Annual Income: 50k 75k",
        "Annual Income: 75k 100k")  ~ "Middle",
      income_category %in% c(
        "Annual Income: 100k 150k",
        "Annual Income: 150k 200k",
        "Annual Income: more 200k") ~ "High",
      TRUE                          ~ NA_character_),
    race_for_position = case_when(
      gender3 == "Gender Minority" ~ "All Races",
      TRUE                         ~ race5c),
    position_label_c = paste(race_for_position, gender3,
                              income3, sep=" | "),
    age_w = pmin(age_2024, quantile(age_2024, 0.99, na.rm=TRUE)),
    age_c = age_w - mean(age_w, na.rm=TRUE)
  ) %>%
  filter(!is.na(race5c), !is.na(gender3), !is.na(income3),
         !is.na(pss_total)) %>%
  mutate(position_id_pss = as.integer(factor(position_label_c)))

cat("PSS analytic N:", nrow(analytic_pss), "\n")
cat("Positions:", n_distinct(analytic_pss$position_label_c), "\n")

# PSS MAIHDA - null model
model_pss_null <- lmer(
  pss_total ~ 1 + (1 | position_id_pss),
  data=analytic_pss, REML=TRUE)
vc_pss <- as.data.frame(VarCorr(model_pss_null))
rownames(vc_pss) <- vc_pss$grp
vpc_pss <- vc_pss["position_id_pss","vcov"] / sum(vc_pss$vcov)
cat("PSS null VPC:", round(vpc_pss*100, 2), "%\n")

# PSS additive model
model_pss_add <- lmer(
  pss_total ~ race5c + gender3 + income3 + age_c +
              (1 | position_id_pss),
  data=analytic_pss, REML=TRUE)
vc_pss_add <- as.data.frame(VarCorr(model_pss_add))
rownames(vc_pss_add) <- vc_pss_add$grp
vpc_pss_add <- vc_pss_add["position_id_pss","vcov"] / sum(vc_pss_add$vcov)
cat("PSS additive VPC:", round(vpc_pss_add*100, 2), "%\n")
cat("Additive explanation:",
    round((vpc_pss-vpc_pss_add)/vpc_pss*100, 1), "%\n")

# PSS predicts hsCRP cross-sectionally
# analytic_pss already contains log_hsCRP_w from crp_pss merge above
# Select directly - do NOT rejoin crp_pss as it creates duplicate columns
analytic_mediation <- analytic_pss %>%
  dplyr::select(person_id, pss_total, position_id_pss,
                position_label_c, race5c, gender3, income3,
                age_c, log_hsCRP_w, hsCRP)

library(lmerTest)
model_pss_crp <- lmer(
  log_hsCRP_w ~ race5c + gender3 + income3 +
                pss_total + age_c +
                (1 | position_id_pss),
  data=analytic_mediation, REML=TRUE)

cat("\nPSS effect on hsCRP (cross-sectional association):\n")
fe_pss <- summary(model_pss_crp)$coefficients
cat("b =", round(fe_pss["pss_total","Estimate"], 4),
    "SE =", round(fe_pss["pss_total","Std. Error"], 4),
    "t =", round(fe_pss["pss_total","t value"], 3),
    "p =", round(fe_pss["pss_total","Pr(>|t|)"], 4), "\n")

# PSS by position - descriptive
pss_by_position <- analytic_pss %>%
  group_by(position_label_c, race5c, gender3, income3) %>%
  summarise(
    n          = n(),
    mean_pss   = round(mean(pss_total, na.rm=TRUE), 2),
    median_pss = round(median(pss_total, na.rm=TRUE), 2),
    .groups="drop") %>%
  arrange(desc(mean_pss))

cat("\nPSS by intersectional position:\n")
print(pss_by_position, n=27)
cat("PSS range:", round(min(pss_by_position$mean_pss),2), "to",
    round(max(pss_by_position$mean_pss),2), "\n")

write_csv(pss_by_position, "IRISE_PSS_by_position.csv")

# Table S3: PSS pathway summary for supplement
tableS3 <- data.frame(
  Measure = c(
    "PSS-10 null VPC",
    "PSS-10 additive explanation",
    "PSS predicts hsCRP: b",
    "SE",
    "t",
    "p",
    "Analytic N"),
  Value = c(
    paste0(round(vpc_pss*100, 2), "%"),
    paste0(round((vpc_pss-vpc_pss_add)/vpc_pss*100, 1), "%"),
    round(fe_pss["pss_total","Estimate"], 4),
    round(fe_pss["pss_total","Std. Error"], 4),
    round(fe_pss["pss_total","t value"], 3),
    round(fe_pss["pss_total","Pr(>|t|)"], 4),
    nrow(analytic_pss)),
  Interpretation = c(
    paste0("Exceeds hsCRP null VPC (", round(vpc_null*100,2), "%)"),
    "Parallel to hsCRP (95.8%)",
    "After intersectional position controls",
    "", "", "", "hsCRP + PSS overlap subsample")
)
write_csv(tableS3, "IRISE_TableS3_PSS_pathway.csv")
cat("Table S3 saved.\n")


# ============================================================
# SECTION 17: SUPPLEMENTARY - TEMPORALLY RESTRICTED MEDIATION
# ============================================================
# Restricted to N=6,023 where PSS was measured before hsCRP.
# Selection bias caveat: this subsample is not representative
# of full analytic sample. Results are confirmatory only.
# Confirmed results: indirect effect 0.0313,
# 95% CI [0.0107, 0.0506], proportion mediated 6.7%.

cat("\n=================================================\n")
cat("SECTION 17: SUPPLEMENTARY MEDIATION (N=6,023)\n")
cat("Cross-sectional. Temporally restricted subsample.\n")
cat("PSS measured before hsCRP in these participants.\n")
cat("Selection bias caveat applies.\n")
cat("=================================================\n")

# Pull PSS dates
pss_dates_sql <- str_glue("
  SELECT DISTINCT person_id, MIN(observation_date) AS pss_date
  FROM `{cdr}.observation`
  WHERE observation_source_value LIKE 'sdoh_cpss_%'
    AND value_as_string IS NOT NULL
    AND value_as_string != 'PMI_Skip'
  GROUP BY person_id
")
pss_dates <- run_query(pss_dates_sql)

# Temporal check
temporal_pss <- pss_dates %>%
  inner_join(
    crp_latest %>% dplyr::select(person_id, measurement_date) %>%
      rename(crp_date = measurement_date),
    by="person_id") %>%
  mutate(
    pss_date  = as.Date(pss_date),
    crp_date  = as.Date(crp_date),
    days_diff = as.numeric(crp_date - pss_date))

cat("PSS before hsCRP (correct order):",
    sum(temporal_pss$days_diff > 0),
    "(", round(mean(temporal_pss$days_diff > 0)*100, 1), "%)\n")
cat("hsCRP before PSS:",
    sum(temporal_pss$days_diff < 0),
    "(", round(mean(temporal_pss$days_diff < 0)*100, 1), "%)\n")
cat("Median gap (days, hsCRP - PSS):",
    round(median(temporal_pss$days_diff, na.rm=TRUE)), "\n")

# Also check Basics survey vs hsCRP temporal ordering
basics_dates_sql <- str_glue("
  SELECT person_id, MIN(survey_datetime) AS basics_date
  FROM `{cdr}.ds_survey`
  WHERE survey = \'The Basics\'
  GROUP BY person_id
")
basics_dates <- run_query(basics_dates_sql)
temporal_basics <- basics_dates %>%
  inner_join(
    crp_latest %>% dplyr::select(person_id, measurement_date) %>%
      rename(crp_date = measurement_date),
    by = "person_id") %>%
  mutate(
    basics_date = as.Date(basics_date),
    crp_date    = as.Date(crp_date),
    days_diff_basics = as.numeric(crp_date - basics_date))

cat("\nBasics survey before hsCRP (correct):",
    sum(temporal_basics$days_diff_basics < 0),
    "(", round(mean(temporal_basics$days_diff_basics < 0)*100, 1), "%)\n")
cat("hsCRP before Basics (reverse):",
    sum(temporal_basics$days_diff_basics > 0),
    "(", round(mean(temporal_basics$days_diff_basics > 0)*100, 1), "%)\n")
cat("Median gap (days, hsCRP - Basics):",
    round(median(temporal_basics$days_diff_basics, na.rm=TRUE)), "\n")

# Table S5: temporal ordering for supplement
tableS5 <- data.frame(
  Measure = c(
    "Basics survey before hsCRP (correct)",
    "hsCRP before Basics (reverse)",
    "PSS before hsCRP (correct)",
    "hsCRP before PSS (reverse)"),
  N = c(
    sum(temporal_basics$days_diff_basics < 0),
    sum(temporal_basics$days_diff_basics > 0),
    sum(temporal_pss$days_diff > 0),
    sum(temporal_pss$days_diff < 0)),
  Pct = round(c(
    mean(temporal_basics$days_diff_basics < 0)*100,
    mean(temporal_basics$days_diff_basics > 0)*100,
    mean(temporal_pss$days_diff > 0)*100,
    mean(temporal_pss$days_diff < 0)*100), 1),
  Notes = c(
    "Of Basics + hsCRP overlap",
    "EHR data retrospective to 1987",
    "Of PSS + hsCRP overlap",
    paste0("Median gap ", round(median(temporal_pss$days_diff[temporal_pss$days_diff<0], na.rm=TRUE)), " days"))
)
write_csv(tableS5, "IRISE_TableS5_temporal_ordering.csv")
cat("Table S5 saved.\n")

# Restrict to correct temporal order
correct_pss_ids <- temporal_pss %>%
  filter(days_diff > 0) %>%
  pull(person_id)

analytic_pss_restricted <- analytic_pss %>%
  filter(person_id %in% correct_pss_ids)
cat("Temporally restricted N:", nrow(analytic_pss_restricted), "\n")

cells_restricted <- analytic_pss_restricted %>%
  count(position_label_c, sort=TRUE)
cat("Min cell:", min(cells_restricted$n), "\n")
cat("Cells < 30:", sum(cells_restricted$n < 30), "\n")

# Mediation models with ML (required for product of coefficients)
# NOTE: analytic_pss_restricted already contains log_hsCRP_w from
# the crp_pss merge in Section 16. Select directly - do NOT rejoin
# crp_pss as it creates duplicate columns (log_hsCRP_w.x etc).
analytic_med_restricted <- analytic_pss_restricted %>%
  dplyr::select(person_id, pss_total, position_id_pss,
                position_label_c, race5c, gender3, income3,
                age_c, log_hsCRP_w, hsCRP)

cat("Mediation restricted N:", nrow(analytic_med_restricted), "\n")
cat("log_hsCRP_w exists:", "log_hsCRP_w" %in% names(analytic_med_restricted), "\n")

path_a_model <- lmer(
  pss_total ~ race5c + gender3 + income3 + age_c +
              (1 | position_id_pss),
  data=analytic_med_restricted, REML=FALSE)

path_b_model <- lmer(
  log_hsCRP_w ~ race5c + gender3 + income3 +
                pss_total + age_c +
                (1 | position_id_pss),
  data=analytic_med_restricted, REML=FALSE)

total_model <- lmer(
  log_hsCRP_w ~ race5c + gender3 + income3 + age_c +
                (1 | position_id_pss),
  data=analytic_med_restricted, REML=FALSE)

a_income <- fixef(path_a_model)["income3Low"]
b_pss    <- fixef(path_b_model)["pss_total"]
indirect <- a_income * b_pss
total    <- fixef(total_model)["income3Low"]
direct   <- fixef(path_b_model)["income3Low"]
prop_med <- indirect / total

cat("\n=== SUPPLEMENTARY MEDIATION RESULTS ===\n")
cat("Path a (income Low -> PSS):", round(a_income, 4), "\n")
cat("Path b (PSS -> hsCRP):     ", round(b_pss, 4), "\n")
cat("Indirect effect (a*b):     ", round(indirect, 4), "\n")
cat("Direct effect:             ", round(direct, 4), "\n")
cat("Total effect:              ", round(total, 4), "\n")
cat("Proportion mediated:       ", round(prop_med*100, 1), "%\n")

# Bootstrap CI
set.seed(42)
n_boot <- 1000
boot_indirect <- numeric(n_boot)
ctrl <- lmerControl(optimizer="bobyqa")
cat("Running bootstrap (2-3 minutes)...\n")

for (i in 1:n_boot) {
  boot_data <- analytic_med_restricted %>%
    group_by(position_id_pss) %>%
    slice_sample(prop=1, replace=TRUE) %>%
    ungroup()
  tryCatch({
    ba <- lmer(pss_total ~ race5c + gender3 + income3 + age_c +
                 (1|position_id_pss), data=boot_data,
               REML=FALSE, control=ctrl)
    bb <- lmer(log_hsCRP_w ~ race5c + gender3 + income3 +
                 pss_total + age_c + (1|position_id_pss),
               data=boot_data, REML=FALSE, control=ctrl)
    boot_indirect[i] <- fixef(ba)["income3Low"] *
                        fixef(bb)["pss_total"]
  }, error=function(e) boot_indirect[i] <<- NA)
  if (i %% 200 == 0) cat("Iteration:", i, "\n")
}

bi <- boot_indirect[!is.na(boot_indirect)]
cat("Successful iterations:", length(bi), "\n")
cat("Indirect effect:", round(indirect, 4), "\n")
cat("95% CI: [", round(quantile(bi, 0.025), 4), ",",
    round(quantile(bi, 0.975), 4), "]\n")
cat("Proportion mediated:", round(prop_med*100, 1), "%\n")
cat("CI excludes zero:", ifelse(
  quantile(bi,0.025) > 0 | quantile(bi,0.975) < 0,
  "YES", "NO"), "\n")

# Save mediation results
save(
  path_a_model, path_b_model, total_model,
  a_income, b_pss, indirect, direct, total, prop_med,
  boot_indirect, analytic_med_restricted,
  file = "IRISE_mediation_final.RData"
)
system(paste0("gsutil cp IRISE_mediation_final.RData ",
              correct_bucket, "/"))

# Table S4: mediation summary for supplement
tableS4 <- data.frame(
  Component = c(
    "Analytic sample",
    "Temporal restriction",
    "Exposure",
    "Mediator",
    "Outcome",
    "Path a: income (low) -> PSS",
    "Path b: PSS -> hsCRP",
    "Indirect effect (a x b)",
    "95% bootstrap CI (1,000 iter.)",
    "Direct effect",
    "Total effect",
    "Proportion mediated"),
  Value = c(
    paste0("N=", nrow(analytic_med_restricted)),
    paste0(round(mean(temporal_pss$days_diff > 0)*100, 1),
           "% of PSS+hsCRP overlap"),
    "Income: low vs high",
    "PSS-10 total score",
    "Log-transformed hsCRP",
    round(a_income, 4),
    round(b_pss, 4),
    round(indirect, 4),
    paste0("[", round(quantile(bi, 0.025), 4), ", ",
           round(quantile(bi, 0.975), 4), "]"),
    round(direct, 4),
    round(total, 4),
    paste0(round(prop_med*100, 1), "%")),
  Notes = c(
    "PSS measured before hsCRP only",
    "Correct temporal order",
    "Single structural axis only",
    "", "",
    "PSS points higher in low vs high income",
    "",
    "Excludes zero",
    "",
    "", "",
    "")
)
write_csv(tableS4, "IRISE_TableS4_mediation.csv")
cat("Table S4 saved.\n")


# ============================================================
# SECTION 17b: SENSITIVITY ANALYSIS - 10-YEAR TEMPORAL RESTRICTION
# ============================================================
# Restricts hsCRP and BMI to measurements from 2010 or later.
# Evaluates whether temporally distant EHR measurements drove
# primary findings. Results reported in supplement Table S6.
# Pre-registered sensitivity analysis.

cat("\n=================================================\n")
cat("SECTION 17b: SENSITIVITY ANALYSIS (2010+ ONLY)\n")
cat("=================================================\n")

crp_10yr <- crp_latest %>%
  filter(measurement_date >= as.Date("2010-01-01")) %>%
  mutate(
    hsCRP_w     = pmin(hsCRP, p99_crp),
    log_hsCRP_w = log(hsCRP_w + 0.1))

bmi_10yr <- bmi_latest %>%
  filter(measurement_date >= as.Date("2010-01-01")) %>%
  mutate(
    bmi   = ifelse(bmi < 15 | bmi > 60, NA, bmi),
    bmi_w = pmin(bmi, p99_bmi, na.rm = TRUE)) %>%
  dplyr::select(person_id, bmi_w)

cat("hsCRP N (all):", nrow(crp_latest), "\n")
cat("hsCRP N (2010+):", nrow(crp_10yr), "\n")
cat("Retained:", round(nrow(crp_10yr)/nrow(crp_latest)*100, 1), "%\n")

analytic_sensitivity <- crp_10yr %>%
  left_join(demo, by = "person_id") %>%
  left_join(race_df, by = "person_id") %>%
  left_join(
    gender_df %>% dplyr::select(person_id, gender_identity),
    by = "person_id") %>%
  left_join(
    income_df %>% dplyr::select(person_id, income_category),
    by = "person_id") %>%
  left_join(bmi_10yr, by = "person_id") %>%
  left_join(smoking_clean, by = "person_id") %>%
  left_join(edu_clean, by = "person_id") %>%
  mutate(
    age_w   = pmin(age_2024, quantile(age_2024, 0.99, na.rm = TRUE)),
    age_c   = age_w - mean(age_w, na.rm = TRUE),
    log_bmi = log(bmi_w),
    ethnicity = case_when(
      str_detect(coalesce(race_all, ""), "Hispanic") ~ "Hispanic/Latino",
      TRUE ~ "Not Hispanic/Latino"),
    race5c = case_when(
      ethnicity == "Hispanic/Latino"                              ~ "Hispanic",
      str_detect(coalesce(race_all, ""), "Black|African")        ~ "Black",
      str_detect(coalesce(race_all, ""), "White") &
        !str_detect(coalesce(race_all, ""), ",")                 ~ "White",
      TRUE                                                        ~ "Other/Asian"),
    gender3 = case_when(
      gender_identity == "Gender Identity: Woman"                 ~ "Cisgender Woman",
      gender_identity == "Gender Identity: Man"                   ~ "Cisgender Man",
      gender_identity %in% c(
        "Gender Identity: Non Binary",
        "Gender Identity: Transgender",
        "Gender Identity: Additional Options")                    ~ "Gender Minority",
      TRUE                                                        ~ NA_character_),
    income3 = case_when(
      income_category %in% c(
        "Annual Income: less 10k",
        "Annual Income: 10k 25k",
        "Annual Income: 25k 35k")                                ~ "Low",
      income_category %in% c(
        "Annual Income: 35k 50k",
        "Annual Income: 50k 75k",
        "Annual Income: 75k 100k")                               ~ "Middle",
      income_category %in% c(
        "Annual Income: 100k 150k",
        "Annual Income: 150k 200k",
        "Annual Income: more 200k")                              ~ "High",
      TRUE                                                        ~ NA_character_),
    race_for_position = case_when(
      gender3 == "Gender Minority" ~ "All Races",
      TRUE ~ race5c),
    position_label_c = paste(race_for_position, gender3,
                              income3, sep = " | ")
  ) %>%
  filter(
    !is.na(log_hsCRP_w), !is.na(race5c), !is.na(gender3),
    !is.na(income3), !is.na(age_w), !is.na(bmi_w),
    !is.na(ever_smoker), !is.na(education3)) %>%
  mutate(position_id_sens = as.integer(factor(position_label_c)))

cat("Sensitivity analytic N:", nrow(analytic_sensitivity), "\n")

model_sens_null <- lmer(
  log_hsCRP_w ~ 1 + (1 | position_id_sens),
  data = analytic_sensitivity, REML = TRUE)

model_sens_add <- lmer(
  log_hsCRP_w ~ race5c + gender3 + income3 +
                age_c + log_bmi + ever_smoker + education3 +
                (1 | position_id_sens),
  data = analytic_sensitivity, REML = TRUE)

vpc_sens_null <- get_vpc(model_sens_null, "position_id_sens")
vpc_sens_add  <- get_vpc(model_sens_add,  "position_id_sens")

pos_range_sens <- analytic_sensitivity %>%
  group_by(position_label_c) %>%
  summarise(median_hsCRP = median(hsCRP, na.rm = TRUE), n = n(),
            .groups = "drop") %>%
  filter(n >= 30)

cat("\n=== SENSITIVITY RESULTS ===\n")
cat("N:", nrow(analytic_sensitivity), "\n")
cat("Null VPC:", round(vpc_sens_null*100, 2), "%\n")
cat("Additive explanation:",
    round((vpc_sens_null - vpc_sens_add)/vpc_sens_null*100, 1), "%\n")
cat("Position range:",
    round(min(pos_range_sens$median_hsCRP), 2), "to",
    round(max(pos_range_sens$median_hsCRP), 2), "mg/L\n")
cat("Fold range:",
    round(max(pos_range_sens$median_hsCRP)/
          min(pos_range_sens$median_hsCRP), 1), "-fold\n")

tableS6 <- data.frame(
  Measure = c("N", "Null VPC (unadjusted)", "Additive explanation",
              "Position hsCRP range", "Intersectional gradient"),
  Primary_analysis = c("47,059", "7.66%", "95.8%",
                       "0.70-6.85 mg/L", "9.8-fold"),
  Sensitivity_2010_plus = c(
    paste0(nrow(analytic_sensitivity),
           " (", round(nrow(analytic_sensitivity)/nrow(analytic_adj)*100, 0), "%)"),
    paste0(round(vpc_sens_null*100, 2), "%"),
    paste0(round((vpc_sens_null-vpc_sens_add)/vpc_sens_null*100, 1), "%"),
    paste0(round(min(pos_range_sens$median_hsCRP), 2), "-",
           round(max(pos_range_sens$median_hsCRP), 2), " mg/L"),
    paste0(round(max(pos_range_sens$median_hsCRP)/
                 min(pos_range_sens$median_hsCRP), 1), "-fold"))
)
write_csv(tableS6, "IRISE_TableS6_sensitivity.csv")
cat("Table S6 saved.\n")

log_decision("Sensitivity analysis",
  "10-year temporal restriction (2010+). N=45,750 (97% retained).",
  paste0("Null VPC 7.59% vs 7.66% primary. Additive explanation 96.3% vs 95.8%. ",
         "Fold range 10.3 vs 9.8. Findings consistent with primary analysis. ",
         "Temporally distant measurements did not drive structural signal."))


# ============================================================
# SECTION 18: FIGURE 1 - CATERPILLAR PLOT
# ============================================================

cat("Producing Figure 1: Caterpillar Plot\n")

re_df <- ranef(model_adj_null)$position_id_adj %>%
  as.data.frame() %>%
  rownames_to_column("position_id_adj") %>%
  rename(re = `(Intercept)`) %>%
  mutate(position_id_adj = as.integer(position_id_adj))

pos_labels <- analytic_adj %>%
  dplyr::select(position_id_adj, position_label_c,
                race5c, gender3, income3) %>%
  distinct()

n_per_pos <- analytic_adj %>% count(position_id_adj)
vc_null_adj <- as.data.frame(VarCorr(model_adj_null))
pos_var   <- vc_null_adj$vcov[1]
resid_var <- vc_null_adj$vcov[2]

caterpillar_data <- re_df %>%
  left_join(pos_labels, by="position_id_adj") %>%
  left_join(n_per_pos,  by="position_id_adj") %>%
  # Deduplicate: gender minority positions share position_label_c
  # across race values after join. Keep one row per position_id_adj.
  group_by(position_id_adj) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(
    se       = sqrt(pos_var * resid_var / (pos_var*n + resid_var)),
    ci_lower = re - 1.96*se,
    ci_upper = re + 1.96*se,
    position_label_c = factor(
      position_label_c,
      levels=position_label_c[order(re)])
  )

fig1 <- ggplot(caterpillar_data,
               aes(x=re, y=position_label_c, color=race5c)) +
  geom_point(size=2.5) +
  geom_errorbarh(aes(xmin=ci_lower, xmax=ci_upper),
                 height=0.3, linewidth=0.5) +
  geom_vline(xintercept=0, linetype="dashed", color="grey50") +
  scale_color_manual(
    values=c("White"="steelblue3", "Black"="firebrick3",
             "Hispanic"="darkorange2", "Other/Asian"="forestgreen"),
    name="Race/Ethnicity") +
  labs(
    title="Figure 1: Intersectional Position hsCRP Rankings",
    subtitle="Posterior mean log hsCRP with 95% credible intervals",
    x="Random effect (log hsCRP)", y=NULL,
    caption=paste0("N=", nrow(analytic_adj),
                   ". 27 intersectional positions. ",
                   "Adjusted for age, BMI, smoking, education.")
  ) +
  theme_minimal(base_size=11) +
  theme(plot.title=element_text(face="bold", size=13),
        axis.text.y=element_text(size=8),
        panel.grid.minor=element_blank())

print(fig1)
ggsave("IRISE_Figure1_caterpillar.png", plot=fig1,
       width=10, height=8, dpi=300)
cat("Figure 1 saved\n")


# ============================================================
# SECTION 19: FIGURE 2 - HEATMAP
# ============================================================

cat("Producing Figure 2: Heatmap\n")

# Suppress cells with N < 30 to avoid artifacts from tiny samples
heatmap_data <- analytic_clean %>%
  group_by(race5c, gender3, income3) %>%
  summarise(median_hsCRP=median(hsCRP, na.rm=TRUE),
            n=n(), .groups="drop") %>%
  filter(!is.na(race5c), !is.na(gender3), !is.na(income3)) %>%
  mutate(
    median_hsCRP_display = ifelse(n < 30, NA, median_hsCRP),
    income3 = factor(income3, levels=c("Low","Middle","High")),
    race5c  = factor(race5c,
                     levels=c("Hispanic","Black","Other/Asian","White"))
  )

fig2 <- ggplot(heatmap_data,
               aes(x=income3, y=race5c, fill=median_hsCRP_display)) +
  geom_tile(color="white", linewidth=0.5) +
  geom_text(
    data=heatmap_data %>% filter(!is.na(median_hsCRP_display)),
    aes(label=round(median_hsCRP_display, 1)),
    size=3.2, color="white", fontface="bold") +
  geom_text(
    data=heatmap_data %>% filter(is.na(median_hsCRP_display)),
    aes(label="<30"), size=2.8, color="grey60") +
  facet_wrap(~gender3, ncol=3) +
  scale_fill_gradient2(
    low="steelblue3", mid="lightyellow", high="firebrick3",
    midpoint=median(heatmap_data$median_hsCRP, na.rm=TRUE),
    na.value="grey90",
    name="Median\nhsCRP\n(mg/L)") +
  labs(
    title="Figure 2: Inflammatory Burden Across Intersectional Positions",
    subtitle="Median hsCRP (mg/L) by race/ethnicity, gender identity, and income",
    x="Annual Household Income", y="Race/Ethnicity",
    caption=paste0("N=", nrow(analytic_clean),
                   ". Grey cells indicate N<30.")
  ) +
  theme_minimal(base_size=11) +
  theme(plot.title=element_text(face="bold", size=13),
        strip.text=element_text(face="bold"),
        axis.text.x=element_text(angle=30, hjust=1),
        panel.grid=element_blank())

print(fig2)
ggsave("IRISE_Figure2_heatmap.png", plot=fig2,
       width=12, height=5, dpi=300)
cat("Figure 2 saved\n")


# ============================================================
# SECTION 20: FIGURE 3 - VPC WATERFALL
# ============================================================

cat("Producing Figure 3: VPC Waterfall\n")

vpc_decomp_clean <- data.frame(
  Component = c(
    "Additive structural\n(race+gender+income)",
    "Residual\nposition VPC",
    "Geographic\n(zip3)",
    "Individual\nlevel"),
  VPC_pct = c(
    round((vpc_adj_null - vpc_adj_add)*100, 2),
    round(vpc_adj_add*100, 2),
    round(vpc_adj_zip3*100, 1),
    round((1 - vpc_adj_zip3 - vpc_adj_pos)*100, 1)
  ),
  Panel = c(rep("Panel A: Two-Level (Covariate-Adjusted)", 2),
            rep("Panel B: Three-Level (Covariate-Adjusted)", 2)),
  fill_color = c("Additive structural", "Residual position",
                 "Geographic", "Individual")
)

cat("VPC decomposition values:\n")
print(vpc_decomp_clean)

fig3 <- ggplot(vpc_decomp_clean,
               aes(x=reorder(Component, -VPC_pct),
                   y=VPC_pct, fill=fill_color)) +
  geom_col(width=0.55, color="white") +
  geom_text(aes(label=paste0(VPC_pct, "%")),
            vjust=-0.4, size=3.8, fontface="bold") +
  facet_wrap(~Panel, scales="free") +
  scale_fill_manual(
    values=c(
      "Additive structural" = "darkorange2",
      "Residual position"   = "firebrick3",
      "Geographic"          = "steelblue4",
      "Individual"          = "grey70"),
    name=NULL) +
  scale_y_continuous(labels=function(x) paste0(x, "%")) +
  labs(
    title="Figure 3: VPC Decomposition",
    subtitle=paste0("Covariate-adjusted null VPC = ",
                    round(vpc_adj_null*100, 2),
                    "%. Additive structural effects explain ",
                    round((vpc_adj_null-vpc_adj_add)/vpc_adj_null*100, 1),
                    "% of position-level clustering."),
    x=NULL, y="Variance Partition Coefficient (%)",
    caption=paste0("Panel A: N=", nrow(analytic_adj),
                   ". Panel B: N=", nrow(analytic_adj_zip),
                   ". Covariates: age, BMI, smoking, education.")
  ) +
  theme_minimal(base_size=11) +
  theme(plot.title=element_text(face="bold", size=13),
        plot.subtitle=element_text(size=9, color="grey30"),
        strip.text=element_text(face="bold"),
        legend.position="bottom",
        panel.grid.minor=element_blank(),
        axis.text.x=element_text(size=9))

print(fig3)
ggsave("IRISE_Figure3_VPC_waterfall.png", plot=fig3,
       width=11, height=6, dpi=300)
cat("Figure 3 saved\n")


# ============================================================
# SECTION 21: FIGURE 4 - DEPRIVATION SCATTER PLOT
# ============================================================
# Replaces state map. Shows dep_index ~ hsCRP relationship
# directly at zip3 level. Connects r=0.147 to visual argument
# for structural embedding through place.

cat("Producing Figure 4: Deprivation Scatter Plot\n")

zip3_summary <- analytic_adj %>%
  filter(!is.na(zip3_clean), zip3_n >= 30,
         !is.na(dep_index_mean)) %>%
  group_by(zip3_clean, dep_index_mean) %>%
  summarise(
    n            = n(),
    mean_hsCRP   = mean(hsCRP, na.rm=TRUE),
    median_hsCRP = median(hsCRP, na.rm=TRUE),
    .groups="drop")

cat("Zip3 areas for scatter:", nrow(zip3_summary), "\n")

# Winsorize y-axis at 95th percentile for display to reduce
# influence of extreme outlier zip3 areas on visual scale
p95_zip3 <- quantile(zip3_summary$median_hsCRP, 0.95)

fig4 <- ggplot(zip3_summary,
               aes(x=dep_index_mean, y=median_hsCRP, size=n)) +
  geom_point(alpha=0.6, color="steelblue4") +
  geom_smooth(method="lm", color="firebrick3",
              se=TRUE, linewidth=1.2) +
  coord_cartesian(ylim=c(0, p95_zip3 * 1.1)) +
  scale_size_continuous(name="N participants", range=c(1,8)) +
  labs(
    title="Figure 4: Neighborhood Deprivation and Inflammatory Burden",
    subtitle="Median hsCRP by neighborhood deprivation index across 3-digit zip code areas",
    x="Neighborhood Deprivation Index (2023 ACS)",
    y="Median hsCRP (mg/L)",
    caption=paste0("N=", nrow(analytic_adj_zip),
                   " participants across ", nrow(zip3_summary),
                   " zip3 areas with N>=30. ",
                   "Y-axis truncated at 95th percentile for display. ",
                   "Deprivation index: geomarker-io 2023 ACS.")
  ) +
  theme_minimal(base_size=11) +
  theme(plot.title=element_text(face="bold", size=13),
        plot.subtitle=element_text(color="grey40"))

print(fig4)
ggsave("IRISE_Figure4_deprivation_scatter.png", plot=fig4,
       width=9, height=6, dpi=300)
cat("Figure 4 saved\n")


# ============================================================
# SECTION 22: FIGURE 5 - DIABETES PREVALENCE DOT PLOT
# ============================================================

cat("Producing Figure 5: Diabetes dot plot\n")

overall_diabetes_pct <- mean(analytic_hba1c$diabetes)*100

fig5_data <- diabetes_by_position %>%
  left_join(
    analytic_hba1c %>%
      dplyr::select(position_label_c, race5c, gender3, income3) %>%
      distinct(),
    by="position_label_c") %>%
  # Deduplicate: gender minority positions share position_label_c
  # Keep one row per position to avoid duplicate factor levels
  group_by(position_label_c) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(
    position_label_c = factor(
      position_label_c,
      levels=position_label_c[order(pct_diabetes)])
  )

fig5 <- ggplot(fig5_data,
               aes(x=pct_diabetes, y=position_label_c,
                   color=race5c, size=n)) +
  geom_point(alpha=0.85) +
  geom_vline(xintercept=overall_diabetes_pct,
             linetype="dashed", color="grey50", linewidth=0.8) +
  annotate("text",
           x=overall_diabetes_pct+0.5, y=1,
           label=paste0("Overall: ",
                        round(overall_diabetes_pct,1), "%"),
           hjust=0, size=3.2, color="grey40") +
  scale_color_manual(
    values=c("White"="steelblue3","Black"="firebrick3",
             "Hispanic"="darkorange2","Other/Asian"="forestgreen"),
    name="Race/Ethnicity") +
  scale_size_continuous(name="N", range=c(2,8)) +
  scale_x_continuous(labels=function(x) paste0(x,"%")) +
  labs(
    title="Figure 5: Diabetes Prevalence Across Intersectional Positions",
    subtitle="Percentage with HbA1c >= 6.5% by intersectional structural position",
    x="Diabetes Prevalence (%)", y=NULL,
    caption=paste0("N=", nrow(analytic_hba1c),
                   ". Point size proportional to position N.")
  ) +
  theme_minimal(base_size=11) +
  theme(plot.title=element_text(face="bold", size=13),
        axis.text.y=element_text(size=8),
        panel.grid.minor=element_blank())

print(fig5)
ggsave("IRISE_Figure5_diabetes.png", plot=fig5,
       width=11, height=8, dpi=300)
cat("Figure 5 saved\n")


# ============================================================
# SECTION 23: SAVE ALL RESULTS
# ============================================================

cat("\n=== PRINTING FULL DECISION LOG ===\n")
for (d in decision_log) cat(d, "\n")

save(
  # Primary models
  model_null, model_additive,
  model_adj_null, model_adj_add, model_adj_threelevel,
  model_hba1c_null, model_hba1c_add, model_diabetes_null,
  # VPCs
  vpc_null, vpc_additive, vpc_adj_null, vpc_adj_add,
  vpc_adj_zip3, vpc_adj_pos,
  vpc_hba1c, vpc_hba1c_add, vpc_diabetes,
  # PSS
  model_pss_null, model_pss_add, model_pss_crp,
  vpc_pss, vpc_pss_add,
  # Mediation (supplementary)
  path_a_model, path_b_model, total_model,
  a_income, b_pss, indirect, direct, total, prop_med,
  boot_indirect,
  # Datasets
  analytic_clean, analytic_adj, analytic_adj_zip,
  analytic_hba1c, analytic_pss, analytic_mediation,
  analytic_pss_restricted, analytic_med_restricted,
  # Results tables
  position_results, diabetes_by_position,
  pss_by_position, vpc_table, replication_table, t1_combined,
  # Reference data
  dep_index_zip3, zip3_summary,
  # Figure data
  caterpillar_data, heatmap_data, fig5_data,
  # Decision log
  decision_log,
  file = "IRISE_complete_workspace.RData"
)

# Save all outputs to correct bucket
for (f in list.files(pattern="^IRISE")) {
  result <- system(paste0("gsutil cp ", f, " ",
                          correct_bucket, "/"),
                   intern=TRUE)
  cat("Saved:", f, "\n")
}

cat("\n=================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("=================================================\n")
cat("hsCRP null VPC (unadj):      ", round(vpc_null*100,2), "%\n")
cat("hsCRP null VPC (adj):        ", round(vpc_adj_null*100,2), "%\n")
cat("Additive explanation (adj):  ",
    round((vpc_adj_null-vpc_adj_add)/vpc_adj_null*100,1), "%\n")
cat("Geographic VPC (adj):        ", round(vpc_adj_zip3*100,1), "%\n")
cat("HbA1c null VPC:              ", round(vpc_hba1c*100,2), "%\n")
cat("Diabetes null VPC:           ", round(vpc_diabetes*100,2), "%\n")
cat("PSS null VPC:                ", round(vpc_pss*100,2), "%\n")
cat("PSS additive explanation:    ",
    round((vpc_pss-vpc_pss_add)/vpc_pss*100,1), "%\n")
cat("Mediation proportion (supp): ", round(prop_med*100,1), "%\n")
cat("Bucket:", correct_bucket, "\n")
