# #############################################################################
# IRISE v9 — DEFINITIVE ANALYTIC PIPELINE (fully self-contained, runnable)
#
# Manuscript: PNEC-D-26-00480 "The Substrate Is Not Neutral"
# Author: S. H. Cook.  Platform: All of Us Controlled Tier v9 (C2025Q4R6),
#         Verily Workbench 2.0.  Statistical approach: MAIHDA (lme4).
#
# WHAT THIS FILE IS
#   The consolidated pipeline that reproduces the locked v9 results from the raw
#   CDR with NO dependency on prior in-memory objects. Gender is harmonized once
#   at source; the geographic inputs (ZIP, site, deprivation) are rebuilt from
#   the CDR (site via measurement_ext.src_id joined on measurement_id); audit
#   targets are ENFORCED gates that halt on mismatch; caches carry provenance.
#
# COMPLETION-PASS CHANGES (vs prior FINAL):
#   [1] Stage 5 geographic inputs rebuilt FROM CDR (no `samp` dependency; no
#       silent skip). Geographic numbers RE-LOCKED on this self-contained build.
#   [2] Section 10 residual robustness is LIVE code (parametric bootstrap + robust
#       SE + diagnostics), not a template. Its CI is locked on first run.
#   [3] AUDIT TARGETS are enforced via check_audit() -> stops on out-of-tolerance.
#   [4] Cache files carry a PROVENANCE stamp (CDR + key params); load refuses a
#       cache built under different specs.
#   [5] Every arbitrary slice(1) replaced with a DETERMINISTIC, documented rule.
#
# CONVENTIONS
#   - Aggregate in BigQuery; pull only per-person rows into R.
#   - AUDIT gates: absolute tolerance TOL_VPC on VPC%, TOL_FOLD on gradients,
#     TOL_N (relative) on sample sizes.
#   - Style: outcome = chronic low-grade inflammation (hs-CRP 0-10 primary).
# #############################################################################

suppressMessages({
  library(dplyr); library(tidyr); library(stringr)
  library(lme4); library(lmerTest); library(digest)
})

# ============================================================================
# SECTION 0 — CONNECT, CONSTANTS, AUDIT MACHINERY
# ============================================================================
stopifnot(exists("run_query"), exists("cdr"))

Q_RACE   <- 1586140L; Q_GENDER <- 1585838L
Q_INCOME <- 1585375L; Q_EDU    <- 1585940L
HS_CID   <- 3010156L  # hs-CRP LOINC 30522-7 (true high-sensitivity)
HBA1C_CID<- 3004410L; BMI_CID <- 3038553L
ZIP_CID  <- 1585250L  # 3-digit ZIP (Basics)

WIN_YRS   <- 5
HS_PRIMARY_HI <- 10; HS_S20_HI <- 20; HS_BROAD_HI <- 200
CELL_MIN_REPORT <- 30
GEO_MIN_ZIP     <- 30   # min persons per ZIP for geographic sample

# expected CDR (provenance). Set to NA to disable the hard CDR check.
EXPECTED_CDR <- "wb-silky-artichoke-2408.C2025Q4R6"

# tolerances for enforced audit gates
TOL_VPC  <- 0.10   # absolute %, on VPC percentages
TOL_FOLD <- 0.20   # absolute, on fold gradients
TOL_N    <- 0.01   # relative, on sample sizes (1%)

ROOT <- "~/IRISE_v9"
EXP  <- file.path(ROOT, "exports", "definitive")
dir.create(EXP, showWarnings = FALSE, recursive = TRUE)
stamp <- format(Sys.time(), "%Y%m%d_%H%M")
LOG <- file(file.path(EXP, paste0("PIPELINE_", stamp, ".txt")), open = "wt")
say <- function(...) { l <- paste0(...); cat(l, "\n"); writeLines(l, LOG) }
prt <- function(df, n = 40) for (i in 1:min(n, nrow(df))) say("   ", paste(format(df[i, ]), collapse = "  "))

# ---- ENFORCED audit gate ---------------------------------------------------
AUDIT_FAILURES <- character(0)
check_audit <- function(label, observed, target, tol, kind = c("abs","rel")) {
  kind <- match.arg(kind)
  diff <- if (kind == "abs") abs(observed - target) else abs(observed - target)/target
  ok <- is.finite(diff) && diff <= tol
  say(sprintf("   AUDIT %-38s obs=%.4f target=%.4f %s(tol=%.3g) -> %s",
              label, observed, target, kind, tol, ifelse(ok, "PASS", "**FAIL**")))
  if (!ok) AUDIT_FAILURES <<- c(AUDIT_FAILURES, sprintf("%s (obs=%.4f target=%.4f)", label, observed, target))
  invisible(ok)
}
check_identity <- function(label, observed, target) {
  ok <- identical(as.character(observed), as.character(target))
  say(sprintf("   AUDIT %-38s obs=%s target=%s -> %s",
              label, as.character(observed), as.character(target),
              ifelse(ok, "PASS", "**FAIL**")))
  if (!ok) AUDIT_FAILURES <<- c(
    AUDIT_FAILURES,
    sprintf("%s (obs=%s target=%s)", label, observed, target)
  )
  invisible(ok)
}
finalize_audits <- function() {
  if (length(AUDIT_FAILURES) == 0) { say("\nALL AUDIT GATES PASSED."); return(invisible(TRUE)) }
  say("\n**AUDIT GATES FAILED**:"); for (f in AUDIT_FAILURES) say("   - ", f)
  stop("Audit gate(s) failed; pipeline did not reproduce locked results within tolerance.")
}

# ---- provenance-stamped cache ---------------------------------------------
save_cache <- function(obj, path, spec) {
  attr(obj, "provenance") <- list(cdr = cdr, spec = spec,
    spec_hash = digest(spec), built = as.character(Sys.time()))
  saveRDS(obj, path); obj
}
load_cache <- function(path, spec) {
  if (!file.exists(path)) return(NULL)
  obj <- readRDS(path); pv <- attr(obj, "provenance")
  if (is.null(pv) || !identical(pv$cdr, cdr) || !identical(pv$spec_hash, digest(spec))) {
    say("   cache provenance mismatch at ", path, " -> rebuilding.")
    return(NULL)
  }
  say("   cache OK (", path, ", built ", pv$built, ")"); obj
}

say("IRISE v9 DEFINITIVE PIPELINE  ", stamp, "  cdr=", cdr)
if (!is.na(EXPECTED_CDR) && !identical(cdr, EXPECTED_CDR))
  stop(sprintf("CDR mismatch: cdr='%s' but EXPECTED_CDR='%s'. Set EXPECTED_CDR=NA to override.", cdr, EXPECTED_CDR))

# ============================================================================
# DETERMINISTIC single-response helper
# For one-per-person survey items: take the EARLIEST observation_date response;
# break ties by lexicographic concept_name. Replaces arbitrary slice(1).
# ============================================================================
one_response_sql <- function(qid, valcol) {
  sprintf("
  WITH r AS (
    SELECT o.person_id, %s, o.observation_date AS d, c.concept_name AS cname
    FROM `%s.observation` o JOIN `%s.concept` c ON o.value_source_concept_id=c.concept_id
    WHERE o.observation_source_concept_id=%d),
  pick AS (SELECT person_id, val,
      ROW_NUMBER() OVER (PARTITION BY person_id ORDER BY d ASC, cname ASC) rn FROM r)
  SELECT person_id, val FROM pick WHERE rn=1", valcol, cdr, cdr, qid)
}

# ============================================================================
# SECTION 1 — ANCHOR
# ============================================================================
anchor <- run_query(sprintf("
  SELECT person_id, MIN(observation_date) AS anchor
  FROM `%s.observation` WHERE observation_source_concept_id IN (%d,%d,%d,%d)
  GROUP BY person_id HAVING COUNT(DISTINCT observation_source_concept_id)=4
", cdr, Q_RACE, Q_GENDER, Q_INCOME, Q_EDU))
say("SECTION 1 anchor people: ", nrow(anchor))
stopifnot(nrow(anchor) > 500000)

# ============================================================================
# SECTION 2 — hs-CRP CLOSEST-TO-ANCHOR (cached w/ provenance)
# ============================================================================
CLOSEST_SPEC <- list(concept=HS_CID, win=WIN_YRS, units=c(8751,4122414,8840,4121396),
                     harmonize="mgdL_x10_to_mgL", rule="closest_abs_then_earliest",
                     site="src_id_of_selected_measurement")
closest_path <- file.path(ROOT, "closest.rds")
closest <- load_cache(closest_path, CLOSEST_SPEC)
if (is.null(closest)) {
  closest <- run_query(sprintf("
  WITH a AS (SELECT person_id, MIN(observation_date) AS anchor FROM `%s.observation`
     WHERE observation_source_concept_id IN (%d,%d,%d,%d) GROUP BY person_id
     HAVING COUNT(DISTINCT observation_source_concept_id)=4),
  v0 AS (
    SELECT person_id, measurement_id, measurement_date AS mdate,
      CASE WHEN unit_concept_id IN (8840,4121396) THEN value_as_number*10
           ELSE value_as_number END AS mgL
    FROM `%s.measurement`
    WHERE measurement_concept_id=%d AND value_as_number IS NOT NULL
      AND unit_concept_id IN (8751,4122414,8840,4121396)
  ),
  v AS (
    SELECT v0.person_id, v0.measurement_id, v0.mdate, v0.mgL,
      ARRAY_AGG(e.src_id IGNORE NULLS ORDER BY e.src_id LIMIT 1)[SAFE_OFFSET(0)] AS site
    FROM v0
    LEFT JOIN `%s.measurement_ext` e USING(measurement_id)
    GROUP BY v0.person_id, v0.measurement_id, v0.mdate, v0.mgL
  ),
  w AS (SELECT v.person_id, v.measurement_id, v.mdate, v.site, v.mgL,
      ABS(DATE_DIFF(v.mdate,a.anchor,DAY)) AS abs_days,
      DATE_DIFF(v.mdate,a.anchor,DAY) AS signed_days
      FROM v JOIN a USING(person_id)
      WHERE v.mgL>0 AND ABS(DATE_DIFF(v.mdate,a.anchor,DAY))<=%d*365),
  pick AS (SELECT person_id, measurement_id, site, mgL, signed_days,
      ROW_NUMBER() OVER (PARTITION BY person_id ORDER BY abs_days ASC, mdate ASC, measurement_id ASC) rn FROM w)
  SELECT person_id, CAST(measurement_id AS STRING) AS measurement_id, site, mgL, signed_days FROM pick WHERE rn=1
  ", cdr, Q_RACE,Q_GENDER,Q_INCOME,Q_EDU, cdr, HS_CID, cdr, WIN_YRS))
  closest <- save_cache(closest, closest_path, CLOSEST_SPEC)
}
say("SECTION 2 closest: ", nrow(closest))
# person_id MUST be unique (one selected measurement per person).
stopifnot(!anyDuplicated(closest$person_id))
# measurement_id should be unique; if not, report before deciding (do not crash
# blindly). A duplicate here would indicate a measurement_ext fan-out or an id
# representation issue; site is already joined in-query so this is diagnostic.
dup_mid <- sum(duplicated(closest$measurement_id))
if (dup_mid > 0) {
  say("   WARNING: ", dup_mid, " duplicate measurement_id values in closest.")
  say("   (site already joined in-query; person_id is unique, so analysis is safe.)")
}
check_audit("closest_N", nrow(closest), 24729, TOL_N, "rel")
check_audit("closest_pct_gt10", 100*mean(closest$mgL>10), 23.8, 0.5, "abs")

# ============================================================================
# SECTION 3 — COVARIATES FROM CDR (deterministic rules; gender at source)
# ============================================================================
# CRITICAL: STRING_AGG DISTINCT — repeated identical race responses (e.g. two
# "White" rows in one sitting) must NOT create "White,White" and trip the
# no-comma single-race test. DISTINCT collapses exact repeats; a comma then
# genuinely indicates multiple DIFFERENT race categories (multiracial).
race_df <- run_query(sprintf("
  SELECT o.person_id, STRING_AGG(DISTINCT c.concept_name ORDER BY c.concept_name) AS race_all,
         COUNT(DISTINCT c.concept_name) AS n_race_cats
  FROM `%s.observation` o JOIN `%s.concept` c ON o.value_source_concept_id=c.concept_id
  WHERE o.observation_source_concept_id=%d GROUP BY o.person_id", cdr,cdr,Q_RACE))
eth_df <- run_query(sprintf("
  SELECT person_id, CASE WHEN ethnicity_concept_id=38003563 THEN 'Hispanic/Latino'
         ELSE 'Not' END AS ethnicity FROM `%s.person`", cdr))
income_df <- run_query(one_response_sql(Q_INCOME, "c.concept_name AS val")) %>% rename(income_category=val)
edu_df    <- run_query(one_response_sql(Q_EDU,    "c.concept_name AS val")) %>% rename(education=val)

# GENDER HARMONIZED AT SOURCE
map_gender <- function(x){
  x0 <- str_squish(as.character(x))
  dplyr::case_when(
    is.na(x0) | x0=="" ~ NA_character_,
    str_detect(x0, regex("prefer not|unknown|skip|none of these|no matching|don't know|pmi", ignore_case=TRUE)) ~ NA_character_,
    str_detect(x0, regex("woman", ignore_case=TRUE)) ~ "Woman",
    str_detect(x0, regex("\\bman\\b|cisgender man", ignore_case=TRUE)) ~ "Man",
    str_detect(x0, regex("non.?binary|transgender|additional|genderqueer|two.?spirit|other", ignore_case=TRUE)) ~ "Gender Minority",
    TRUE ~ NA_character_)
}
gender_df <- run_query(one_response_sql(Q_GENDER, "c.concept_name AS val")) %>%
  rename(raw_gender=val) %>% mutate(gender3=map_gender(raw_gender))
gmap <- gender_df %>% count(raw_gender, gender3, name="n") %>% arrange(desc(n))
say("SECTION 3b gender mapping:"); prt(gmap, 12)
write.csv(gmap, file.path(EXP, paste0("gender_mapping_", stamp, ".csv")), row.names=FALSE)
stopifnot(setequal(na.omit(unique(gmap$gender3)), c("Man","Woman","Gender Minority")))

smoke_df <- run_query(sprintf("
  WITH r AS (SELECT o.person_id, c.concept_name AS s, o.observation_date AS d
    FROM `%s.observation` o JOIN `%s.concept` c ON o.value_source_concept_id=c.concept_id
    WHERE o.observation_source_concept_id IN (1585857,1585858,1585864,1586159)),
  pick AS (SELECT person_id, s, ROW_NUMBER() OVER (PARTITION BY person_id ORDER BY d ASC, s ASC) rn FROM r)
  SELECT person_id, s FROM pick WHERE rn=1", cdr,cdr)) %>%
  mutate(ever_smoker=case_when(s=="100 Cigs Lifetime: Yes"~1L, s=="100 Cigs Lifetime: No"~0L, TRUE~NA_integer_)) %>%
  select(person_id, ever_smoker)

bmi_df <- run_query(sprintf("
  WITH a AS (SELECT person_id, MIN(observation_date) AS anchor FROM `%s.observation`
    WHERE observation_source_concept_id IN (%d,%d,%d,%d) GROUP BY person_id
    HAVING COUNT(DISTINCT observation_source_concept_id)=4),
  v AS (SELECT person_id, measurement_date AS mdate, value_as_number AS bmi
    FROM `%s.measurement` WHERE measurement_concept_id=%d AND value_as_number BETWEEN 10 AND 80),
  p AS (SELECT v.person_id, v.bmi, ROW_NUMBER() OVER
    (PARTITION BY v.person_id ORDER BY ABS(DATE_DIFF(v.mdate,a.anchor,DAY)) ASC, v.mdate ASC) rn
    FROM v JOIN a USING(person_id))
  SELECT person_id, bmi FROM p WHERE rn=1
  ", cdr,Q_RACE,Q_GENDER,Q_INCOME,Q_EDU, cdr, BMI_CID)) %>%
  mutate(bmi=ifelse(bmi<15|bmi>60,NA,bmi))
bmi_df$log_bmi <- log(pmin(bmi_df$bmi, quantile(bmi_df$bmi,0.99,na.rm=TRUE)))

age_df <- run_query(sprintf("
  WITH a AS (SELECT person_id, MIN(observation_date) AS anchor FROM `%s.observation`
    WHERE observation_source_concept_id IN (%d,%d,%d,%d) GROUP BY person_id
    HAVING COUNT(DISTINCT observation_source_concept_id)=4)
  SELECT p.person_id, DATE_DIFF(a.anchor, DATE(p.birth_datetime), DAY)/365.25 AS age_anchor
  FROM `%s.person` p JOIN a USING(person_id)
  ", cdr,Q_RACE,Q_GENDER,Q_INCOME,Q_EDU, cdr)) %>%
  mutate(age_anchor=ifelse(age_anchor<18|age_anchor>120,NA,age_anchor))
age_df$age_w <- pmin(age_df$age_anchor, quantile(age_df$age_anchor,0.99,na.rm=TRUE))
age_df$age_c <- age_df$age_w - mean(age_df$age_w, na.rm=TRUE)
age_df <- age_df %>% select(person_id, age_c)

ha_df <- run_query(sprintf("
  WITH a AS (SELECT person_id, MIN(observation_date) AS anchor FROM `%s.observation`
    WHERE observation_source_concept_id IN (%d,%d,%d,%d) GROUP BY person_id
    HAVING COUNT(DISTINCT observation_source_concept_id)=4),
  v AS (SELECT person_id, measurement_date AS mdate, value_as_number AS pct
    FROM `%s.measurement` WHERE measurement_concept_id=%d AND unit_concept_id=8554
      AND value_as_number BETWEEN 3 AND 20),
  p AS (SELECT v.person_id, v.pct, ROW_NUMBER() OVER
    (PARTITION BY v.person_id ORDER BY ABS(DATE_DIFF(v.mdate,a.anchor,DAY)) ASC, v.mdate ASC) rn
    FROM v JOIN a USING(person_id) WHERE ABS(DATE_DIFF(v.mdate,a.anchor,DAY))<=%d*365)
  SELECT person_id, pct AS hba1c FROM p WHERE rn=1
  ", cdr,Q_RACE,Q_GENDER,Q_INCOME,Q_EDU, cdr, HBA1C_CID, WIN_YRS)) %>%
  mutate(hba1c_w=pmin(hba1c,quantile(hba1c,0.99,na.rm=TRUE)),
         log_hba1c_w=log(hba1c_w), diabetes=as.integer(hba1c>=6.5))

# ZIP (deterministic), SITE (carried on the selected hs-CRP record),
# DEPRIVATION (external, cached)
zip_df <- run_query(sprintf("
  WITH r AS (SELECT person_id, value_as_string AS zip3, observation_date AS d
    FROM `%s.observation` WHERE observation_source_concept_id=%d AND value_as_string IS NOT NULL),
  pick AS (SELECT person_id, zip3, ROW_NUMBER() OVER
    (PARTITION BY person_id ORDER BY d ASC, zip3 ASC) rn FROM r)
  SELECT person_id, zip3 FROM pick WHERE rn=1", cdr, ZIP_CID)) %>%
  mutate(zip3 = substr(str_trim(as.character(zip3)), 1, 3))

say("SECTION 3d zip people=", nrow(zip_df),
    " selected-measurement site nonmissing=", sum(!is.na(closest$site)))

# DEPRIVATION: pin to a specific commit (not a moving branch), record checksum
# and retrieval date, and DOCUMENT the ZIP3 aggregation rule.
# NOTE ON AGGREGATION: the current rule is an UNWEIGHTED mean of the ADI-style
# index across 5-digit ZCTAs within each 3-digit prefix. This gives each ZCTA
# equal weight regardless of population. If a population-weighted ZIP3 index is
# intended, replace with a weighted mean using ZCTA population. Flagged for a
# deliberate decision, not left implicit.
DEP_COMMIT <- Sys.getenv("IRISE_DEP_COMMIT")
EXPECTED_DEP_MD5 <- Sys.getenv("IRISE_DEP_MD5")
# CONFIRMATION-RUN ALLOWANCE: if IRISE_DEP_ALLOW_MASTER=1, permit the 'master'
# branch (working source) and skip the checksum requirement. Use ONLY for a
# confirmation run; for the LOCKED submission analysis, set a real SHA + md5.
DEP_ALLOW_MASTER <- Sys.getenv("IRISE_DEP_ALLOW_MASTER") == "1"
if (!nzchar(DEP_COMMIT)) DEP_COMMIT <- "master"
if (!DEP_ALLOW_MASTER && DEP_COMMIT %in% c("master","main"))
  stop("Set IRISE_DEP_COMMIT to a commit SHA, or IRISE_DEP_ALLOW_MASTER=1 for a confirmation run.")
if (!DEP_ALLOW_MASTER && !nzchar(EXPECTED_DEP_MD5))
  stop("Set IRISE_DEP_MD5 to the MD5 checksum of the locked deprivation CSV.")
if (DEP_COMMIT %in% c("cached","")) DEP_COMMIT <- "master"  # never build a fake-commit URL
DEP_URL <- sprintf("https://raw.githubusercontent.com/geomarker-io/dep_index/%s/2023/data/ACS_deprivation_index_by_zipcode.csv", DEP_COMMIT)
DEP_SPEC <- list(source="geomarker_dep_index_2023", url=DEP_URL, commit=DEP_COMMIT,
                 agg="unweighted_mean_zcta_within_zip3")
dep_path <- file.path(ROOT, "dep_by_zip3.rds")
dep_df <- load_cache(dep_path, DEP_SPEC)
if (is.null(dep_df)) {
  raw_path <- file.path(ROOT, "dep_raw_source.csv")
  dep_df <- tryCatch({
    download.file(DEP_URL, raw_path, quiet=TRUE)
    chk <- tools::md5sum(raw_path)[[1]]
    if (nzchar(EXPECTED_DEP_MD5) && !identical(unname(chk), EXPECTED_DEP_MD5))
      stop("Deprivation source checksum does not match IRISE_DEP_MD5.")
    say("   DEP source md5=", chk, " retrieved=", as.character(Sys.Date()), " commit=", DEP_COMMIT)
    di <- read.csv(raw_path)
    out <- di %>% mutate(zip3=substr(sprintf("%05d", as.integer(zcta_2020)),1,3)) %>%
      group_by(zip3) %>% summarise(deprivation_index=mean(dep_index,na.rm=TRUE),
                                   n_zcta=n(), .groups="drop")
    attr(out, "dep_source") <- list(url=DEP_URL, md5=chk, retrieved=as.character(Sys.Date()), commit=DEP_COMMIT)
    out
  }, error=function(e){ say("   DEP fetch failed: ", conditionMessage(e)); NULL })
  if (!is.null(dep_df)) dep_df <- save_cache(dep_df, dep_path, DEP_SPEC)
}
if (is.null(dep_df)) stop("Deprivation index unavailable (network). Cannot build geographic inputs self-contained.")

# ============================================================================
# SECTION 4 — POSITION VARIABLES + 27-POSITION RECONSTRUCTION
# ============================================================================
covars <- anchor %>% select(person_id) %>%
  left_join(race_df,"person_id") %>% left_join(eth_df,"person_id") %>%
  left_join(gender_df %>% select(person_id, gender3),"person_id") %>%
  left_join(income_df,"person_id") %>% left_join(edu_df,"person_id") %>%
  left_join(smoke_df,"person_id") %>% left_join(bmi_df %>% select(person_id,log_bmi),"person_id") %>%
  left_join(age_df,"person_id") %>% left_join(ha_df,"person_id") %>%
  left_join(zip_df,"person_id") %>%
  left_join(dep_df,"zip3") %>%
  mutate(
    # single-race = exactly one DISTINCT substantive category (n_race_cats==1);
    # multiple distinct categories -> multiracial -> Other (unless Black/Hispanic
    # per priority hierarchy). Uses n_race_cats, not comma detection.
    race5 = case_when(
      ethnicity=="Hispanic/Latino" ~ "Hispanic",
      str_detect(coalesce(race_all,""),"Black|African") ~ "Black",
      str_detect(coalesce(race_all,""),"White") & n_race_cats==1 ~ "White",
      str_detect(coalesce(race_all,""),"Asian") & n_race_cats==1 ~ "Asian",
      TRUE ~ "Other"),
    race5c = ifelse(race5 %in% c("Asian","Other"),"Other/Asian",race5),
    income3 = case_when(is.na(income_category)~NA_character_,
      str_detect(income_category,"less 10k|10k 25k|25k 35k")~"Low",
      str_detect(income_category,"35k 50k|50k 75k|75k 100k")~"Middle",
      str_detect(income_category,"100k|150k|200k|more")~"High", TRUE~NA_character_),
    education3 = case_when(
      education %in% c("Highest Grade: Never Attended","Highest Grade: One Through Four",
        "Highest Grade: Five Through Eight","Highest Grade: Nine Through Eleven",
        "Highest Grade: Twelve Or GED")~"Low",
      education=="Highest Grade: College One to Three"~"Middle",
      education %in% c("Highest Grade: College Graduate","Highest Grade: Advanced Degree")~"High",
      TRUE~NA_character_),
    gender3 = factor(gender3, levels=c("Man","Woman","Gender Minority")),
    race_for_pos = ifelse(as.character(gender3)=="Gender Minority","All Races", race5c),
    position = ifelse(!is.na(race_for_pos)&!is.na(gender3)&!is.na(income3),
                      paste(race_for_pos, as.character(gender3), income3, sep=" | "), NA))
stopifnot(!anyDuplicated(covars$person_id))
stopifnot(nlevels(covars$gender3)==3)
# RACE validation gate + exported raw->analytic race mapping
stopifnot(setequal(unique(na.omit(covars$race5c)), c("White","Black","Hispanic","Other/Asian")))
race_map <- covars %>% count(race_all, n_race_cats, ethnicity, race5c, name="n") %>% arrange(desc(n))
write.csv(race_map, file.path(EXP, paste0("race_mapping_", stamp, ".csv")), row.names=FALSE)
say("SECTION 4 race categories: ", paste(sort(unique(na.omit(covars$race5c))), collapse=", "))
say("   race5c frequencies:")
prt(covars %>% count(race5c) %>% arrange(desc(n)), 6)
np_all <- dplyr::n_distinct(covars$position[!is.na(covars$position)])
say("SECTION 4 gender levels=", nlevels(covars$gender3), " positions=", np_all)
stopifnot(np_all == 27)   # EXACT: a deficient reconstruction must fail

# ============================================================================
# SECTION 5 — FOUR ANALYTIC SAMPLES + GATES + AUDITS
# ============================================================================
CC_VARS <- c("log_hsCRP_w","race5c","gender3","income3","age_c","log_bmi","ever_smoker","education3")
mk_sample <- function(hi, tag) {
  d <- closest %>% filter(mgL>0, mgL<=hi) %>%
    mutate(hsCRP=mgL, hsCRP_w=pmin(mgL, quantile(mgL,0.99)), log_hsCRP_w=log(hsCRP_w+0.1)) %>%
    select(person_id, site, hsCRP, hsCRP_w, log_hsCRP_w, hs_signed=signed_days) %>%
    inner_join(covars, "person_id")
  cc <- d %>% filter(if_all(all_of(CC_VARS), ~ !is.na(.)))
  cc$position_id <- as.integer(factor(cc$position))
  psz <- cc %>% count(position)
  say(sprintf("SECTION 5 [%s] N=%d positions=%d min_cell=%d", tag, nrow(cc), n_distinct(cc$position), min(psz$n)))
  stopifnot(n_distinct(cc$position)==27)
  cc
}
samp_primary <- mk_sample(HS_PRIMARY_HI, "PRIMARY 0-10")
samp_s20     <- mk_sample(HS_S20_HI,     "SENS 0-20")
samp_broad   <- mk_sample(HS_BROAD_HI,   "BROAD 0-200")
samp_capped  <- samp_broad %>% mutate(hsCRP_cap=pmin(hsCRP,10), log_hsCRP_cap=log(hsCRP_cap+0.1))
# EXACT sample-size gates (locked CDR + deterministic pipeline => exact match)
check_audit("primary_N", nrow(samp_primary), 14601, 5, "abs")
check_audit("s20_N",     nrow(samp_s20),     16142, 5, "abs")
check_audit("broad_N",   nrow(samp_broad),   18463, 5, "abs")
for (nm in c("samp_primary","samp_s20","samp_broad"))
  save_cache(get(nm), file.path(EXP, paste0(nm,"_",stamp,".rds")), list(sample=nm, closest=CLOSEST_SPEC))

# ============================================================================
# SECTION 6 — CORE MAIHDA MODELS 1-4 + AUDIT GATES
# ============================================================================
ctrl <- lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=3e5))
vc <- function(f){v<-as.data.frame(VarCorr(f));setNames(v$vcov,v$grp)}
maihda <- function(d, yvar, tag) {
  d$position_id <- as.integer(factor(d$position))
  m1<-lmer(reformulate("(1|position_id)", yvar), data=d, REML=TRUE, control=ctrl)
  m2<-lmer(reformulate(c("race5c","gender3","income3","(1|position_id)"), yvar), data=d, REML=TRUE, control=ctrl)
  m3<-lmer(reformulate(c("age_c","log_bmi","ever_smoker","education3","(1|position_id)"), yvar), data=d, REML=TRUE, control=ctrl)
  m4<-lmer(reformulate(c("race5c","gender3","income3","age_c","log_bmi","ever_smoker","education3","(1|position_id)"), yvar), data=d, REML=TRUE, control=ctrl)
  V<-sapply(list(m1,m2,m3,m4), function(m) vc(m)["position_id"])
  tot<-sapply(list(m1,m2,m3,m4), function(m) sum(vc(m)))
  vpc1<-100*V[1]/tot[1]; vpc4<-100*V[4]/tot[4]; addexpl<-100*(V[3]-V[4])/V[3]
  say(""); say("== MAIHDA ", tag, " N=", nrow(d), " ==")
  say(sprintf("   M1 null VPC=%.3f%%  M4 adj VPC=%.3f%%  adj additive expl=%.1f%%  M4 singular=%s",
      vpc1, vpc4, addexpl, isSingular(m4,1e-4)))
  fe<-as.data.frame(summary(m4)$coefficients); fe$term<-rownames(fe)
  ci<-tryCatch(confint(m4,method="Wald"),error=function(e)NULL)
  write.csv(cbind(fe, ci_lo=if(!is.null(ci))ci[fe$term,1] else NA, ci_hi=if(!is.null(ci))ci[fe$term,2] else NA),
            file.path(EXP,paste0("M4_fixed_",gsub("[^0-9]","",tag),"_",stamp,".csv")))
  rawvar <- if(yvar=="log_hsCRP_cap") "hsCRP_cap" else "hsCRP"
  gr <- d %>% group_by(position) %>% summarise(m=median(.data[[rawvar]]), n=n(), .groups="drop") %>%
    filter(n>=CELL_MIN_REPORT) %>% arrange(m)
  grad<-max(gr$m)/min(gr$m)
  say(sprintf("   DESCRIPTIVE gradient %.1f-fold  lo=%s hi=%s", grad, gr$position[1], gr$position[nrow(gr)]))
  list(vpc1=vpc1, vpc4=vpc4, addexpl=addexpl, grad=grad, m4=m4, N=nrow(d))
}
r_prim  <- maihda(samp_primary, "log_hsCRP_w", "PRIMARY 0-10")
r_s20   <- maihda(samp_s20,     "log_hsCRP_w", "SENS 0-20")
r_broad <- maihda(samp_broad,   "log_hsCRP_w", "BROAD 0-200")
check_audit("primary_null_VPC",  r_prim$vpc1,   3.08, 0.15, "abs")
check_audit("primary_adj_VPC",   r_prim$vpc4,   0.30, 0.15, "abs")
check_audit("primary_add_expl",  r_prim$addexpl,83.6, 2.0,     "abs")
check_audit("primary_gradient",  r_prim$grad,   2.5,  TOL_FOLD,"abs")
# sensitivity-specification gates (0-20 and broad-range)
check_audit("s20_null_VPC",   r_s20$vpc1,    3.97, 0.15, "abs")
check_audit("s20_adj_VPC",    r_s20$vpc4,    0.26, 0.15, "abs")
check_audit("s20_add_expl",   r_s20$addexpl, 86.4, 2.0,     "abs")
check_audit("s20_gradient",   r_s20$grad,    3.3,  TOL_FOLD,"abs")
check_audit("broad_null_VPC", r_broad$vpc1,  4.77, 0.15, "abs")
check_audit("broad_adj_VPC",  r_broad$vpc4,  0.14, 0.15, "abs")
check_audit("broad_add_expl", r_broad$addexpl,94.1,2.0,     "abs")
check_audit("broad_gradient", r_broad$grad,  4.0,  TOL_FOLD,"abs")

# ============================================================================
# SECTION 7 — CAPPED-AT-10 + LOSS LADDER + >10-by-position + AUDIT
# ============================================================================
say(""); say("### SECTION 7 SAMPLE-LOSS LADDER")
say(sprintf("   closest=%d  >10=%d (%.1f%%)  <=10=%d  broad_cc=%d  primary_cc=%d",
    nrow(closest), sum(closest$mgL>10), 100*mean(closest$mgL>10), sum(closest$mgL<=10),
    nrow(samp_broad), nrow(samp_primary)))
byp <- samp_broad %>% mutate(over10=hsCRP>10) %>% group_by(position) %>%
  summarise(n=n(), n_over10=sum(over10), pct_over10=round(100*mean(over10),1), .groups="drop") %>% arrange(desc(pct_over10))
write.csv(byp, file.path(EXP,paste0("over10_by_position_",stamp,".csv")), row.names=FALSE)
r_capped <- maihda(samp_capped, "log_hsCRP_cap", "CAPPED at 10")
check_audit("capped_null_VPC", r_capped$vpc1, 4.56, 0.15, "abs")
check_audit("capped_adj_VPC",  r_capped$vpc4, 0.14, 0.15, "abs")
check_audit("capped_add_expl", r_capped$addexpl, 92.9, 2.0, "abs")

# ============================================================================
# SECTION 8 — SECONDARY axis-count
# ============================================================================
say(""); say("### SECTION 8 AXIS-COUNT (primary)")
d <- samp_primary; d$position_id <- as.integer(factor(d$position))
nullv <- {c<-vc(lmer(log_hsCRP_w~1+(1|position_id),data=d,REML=TRUE,control=ctrl)); c["position_id"]/sum(c)}
sets <- list(race="race5c",gender="gender3",income="income3",
  `race+gender`=c("race5c","gender3"),`race+income`=c("race5c","income3"),
  `gender+income`=c("gender3","income3"),`all three`=c("race5c","gender3","income3"))
ax <- do.call(rbind, lapply(names(sets), function(nm){
  f<-lmer(reformulate(c(sets[[nm]],"(1|position_id)"),"log_hsCRP_w"),data=d,REML=TRUE,control=ctrl)
  v<-vc(f); data.frame(set=nm, resid_VPC=100*v["position_id"]/sum(v),
                       explains=100*(nullv - v["position_id"]/sum(v))/nullv)}))
prt(ax,7); write.csv(ax, file.path(EXP,paste0("axis_count_",stamp,".csv")), row.names=FALSE)

# ============================================================================
# SECTION 9 — STAGE 5 GEOGRAPHIC (SELF-CONTAINED; geo inputs from CDR)
# ============================================================================
g0 <- samp_primary %>%
  filter(!is.na(zip3), zip3!="000", zip3!="", !is.na(deprivation_index), !is.na(site))
zsz <- g0 %>% count(zip3, name="zn")
geo <- g0 %>% left_join(zsz,"zip3") %>% filter(zn>=GEO_MIN_ZIP) %>% mutate(deprivation=deprivation_index)
geo$position_id <- as.integer(factor(geo$position))
save_cache(geo, file.path(ROOT,"geo.rds"), list(build="cdr_selfcontained", min_zip=GEO_MIN_ZIP))
say(""); say(sprintf("### SECTION 9 GEO N=%d zips=%d sites=%d", nrow(geo), n_distinct(geo$zip3), n_distinct(geo$site)))
FEg <- "race5c+gender3+income3+age_c+log_bmi+ever_smoker+education3+deprivation"
m5  <- lmer(as.formula(paste("log_hsCRP_w ~",FEg,"+ (1|position_id)+(1|zip3)")), data=geo, REML=TRUE, control=ctrl)
m5s <- lmer(as.formula(paste("log_hsCRP_w ~",FEg,"+ site + (1|position_id)+(1|zip3)")), data=geo, REML=TRUE, control=ctrl)
c5<-vc(m5); cs<-vc(m5s); zip_vpc <- 100*c5["zip3"]/sum(c5)
say(sprintf("   M5 zip VPC=%.3f%%  zip var=%.5f  position VPC=%.3f%%  singular=%s",
    zip_vpc, c5["zip3"], 100*c5["position_id"]/sum(c5), isSingular(m5,1e-4)))
say(sprintf("   +site(fixed): zip var %.5f -> %.5f (area variation REMAINS after site means)", c5["zip3"], cs["zip3"]))
re<-ranef(m5)$zip3[,1]; names(re)<-rownames(ranef(m5)$zip3); ext<-names(which.max(abs(re)))
m5x<-lmer(as.formula(paste("log_hsCRP_w ~",FEg,"+ (1|position_id)+(1|zip3)")), data=geo %>% filter(zip3!=ext), REML=TRUE, control=ctrl)
cx<-vc(m5x); infl_pct<-100*(cx["zip3"]-c5["zip3"])/c5["zip3"]
say(sprintf("   extreme zip %s BLUP=%.3f; refit w/o: zip var %.5f -> %.5f (%+.0f%%)", ext, re[ext], c5["zip3"], cx["zip3"], infl_pct))
prof<-tryCatch({pp<-profile(m5,which="theta_",signames=FALSE);confint(pp)},error=function(e)NULL)
zip_var_ci <- c(NA,NA)
if(!is.null(prof)){zr<-grep("zip3",rownames(prof)); zip_var_ci<-c(prof[zr,1]^2, prof[zr,2]^2)
  say(sprintf("   zip VAR CI [%.5f,%.5f] (profile)", zip_var_ci[1], zip_var_ci[2]))}
# GEOGRAPHIC audit gates (all locked; a mismatch triggers INVESTIGATION, not relock).
# Targets from the overlap-build; if the CDR-rebuilt geo diverges, STOP and reconcile.
check_audit("geo_N",            nrow(geo),              11514, 20,     "abs")
check_audit("geo_n_zips",       n_distinct(geo$zip3),   88,    3,      "abs")
check_audit("geo_n_sites",      n_distinct(geo$site),   32,    3,      "abs")
check_audit("geo_zip_VPC",      zip_vpc,                3.38,  0.20,   "abs")
check_audit("geo_zip_var",      c5["zip3"],             0.02662, 0.003,"abs")
check_audit("geo_position_var", c5["position_id"],      0.00401, 0.002,"abs")
check_audit("geo_site_adj_zip_var", cs["zip3"],         0.04930, 0.006,"abs")
check_audit("geo_infl_pct",     infl_pct,               -57.0, 6.0,    "abs")
# This is an audit target, not a forced analytic result. If correcting site
# assignment changes the influential ZIP3, exports are still written and the
# pipeline stops at final audit review rather than overwriting the new result.
check_identity("geo_influential_zip3", ext, "782")

# ============================================================================
# SECTION 10 — RESIDUAL ROBUSTNESS (LIVE)
# ============================================================================
say(""); say("### SECTION 10 RESIDUAL ROBUSTNESS (live)")
# WHAT EACH PIECE ESTABLISHES (stated precisely):
#  - Parametric bootstrap: uncertainty in the ZIP variance UNDER the fitted
#    Gaussian mixed model. It is a parametric uncertainty analysis; it does NOT
#    remedy or test robustness to the observed non-normal residual pattern.
#  - CR2 cluster-robust (sandwich) SEs: fixed-effect inference that does not
#    rely on the Gaussian residual assumption. This is the residual-nonnormality
#    safeguard for the FIXED effects (not the variance component).
if (!requireNamespace("clubSandwich", quietly=TRUE))
  stop("Package 'clubSandwich' is required for Section 10. Install: install.packages('clubSandwich').")
set.seed(1)
bf <- function(fm){ v<-as.data.frame(VarCorr(fm)); v$sdcor[v$grp=="zip3"] }
NSIM <- 1000
bb <- bootMer(m5, bf, nsim=NSIM, type="parametric", use.u=FALSE)
v <- bb$t[,1]; ok <- sum(is.finite(v)); q_sd <- quantile(v[is.finite(v)], c(.025,.975))
zip_sd_boot_ci <- as.numeric(q_sd)
say(sprintf("   zip SD parametric-bootstrap CI [%.4f, %.4f] (var [%.5f,%.5f]) success=%d/%d",
    q_sd[1], q_sd[2], q_sd[1]^2, q_sd[2]^2, ok, NSIM))
saveRDS(bb, file.path(EXP, paste0("zip_bootstrap_", stamp, ".rds")))
check_audit("zip_bootstrap_success", ok, NSIM, 0, "abs")
# AUDIT: lock on first run. Uncomment with observed values after the first run:
# check_audit("zip_sd_boot_lo", q_sd[1], <LOCK_LO>, 0.01, "abs")
# check_audit("zip_sd_boot_hi", q_sd[2], <LOCK_HI>, 0.01, "abs")
# CR2 on the mixed model fails for CROSS-CLASSIFIED effects (position x zip3).
# Fall back to an OLS marginal model with zip3 cluster-robust (CR2) SEs, which
# answers whether the fixed-effect pattern persists under geographic clustering.
ct <- tryCatch({
  vcr <- clubSandwich::vcovCR(m5, type="CR2", cluster=geo$zip3)
  clubSandwich::coef_test(m5, vcov=vcr, test="Satterthwaite")
}, error=function(e){
  say("   clubSandwich on cross-classified m5 not supported (", conditionMessage(e), ")")
  say("   -> OLS marginal model with zip3 CR2 cluster-robust SEs (sensitivity).")
  lm_fe <- lm(as.formula(paste("log_hsCRP_w ~", FEg)), data=geo)
  vcr <- clubSandwich::vcovCR(lm_fe, type="CR2", cluster=geo$zip3)
  clubSandwich::coef_test(lm_fe, vcov=vcr, test="Satterthwaite")
})
write.csv(ct, file.path(EXP, paste0("robust_fixed_CR2_", stamp, ".csv")))
say("   CR2 zip3-clustered fixed-effect tests written (OLS marginal if m5 unsupported).")

# A separate OLS sensitivity removes both random intercepts and therefore does
# NOT estimate the M5 mixed-model fixed effects with alternate standard errors.
# It answers whether the fixed-effect pattern persists in a marginal model.
ols_sens <- lm(as.formula(paste("log_hsCRP_w ~", FEg, "+ factor(zip3)")), data=geo)
ols_tab <- as.data.frame(summary(ols_sens)$coefficients)
ols_tab$term <- rownames(ols_tab)
write.csv(ols_tab, file.path(EXP, paste0("OLS_ZIP3_FE_sensitivity_", stamp, ".csv")),
          row.names=FALSE)
say("   OLS ZIP3 fixed-effects sensitivity written (different estimand from M5).")

# Publication-sized diagnostics. The random-effect plots are exported
# separately because the ZIP3 and position effects have very different scales.
png(file.path(EXP, paste0("diag_resid_qq_M5_", stamp, ".png")),
    width=1400, height=1400, res=180)
qqnorm(resid(m5), main="Conditional residual Q-Q (M5)", pch=16, cex=.45)
qqline(resid(m5), col="red", lwd=2)
dev.off()

png(file.path(EXP, paste0("diag_resid_fitted_M5_", stamp, ".png")),
    width=1400, height=1400, res=180)
plot(fitted(m5), resid(m5), pch=16, cex=.35,
     xlab="Fitted values", ylab="Conditional residuals",
     main="Conditional residuals versus fitted values (M5)")
abline(h=0, col="red", lwd=2)
dev.off()

zip_re <- ranef(m5)$zip3[, 1]
pos_re <- ranef(m5)$position_id[, 1]
png(file.path(EXP, paste0("diag_ranef_zip3_qq_M5_", stamp, ".png")),
    width=1400, height=1400, res=180)
qqnorm(zip_re, main="Random-effect Q-Q: three-digit ZIP area (M5)",
       pch=16, col="darkgreen")
qqline(zip_re, col="red", lwd=2)
dev.off()

png(file.path(EXP, paste0("diag_ranef_position_qq_M5_", stamp, ".png")),
    width=1400, height=1400, res=180)
qqnorm(pos_re, main="Random-effect Q-Q: intersectional position (M5)",
       pch=16, col="navy")
qqline(pos_re, col="red", lwd=2)
dev.off()
say("   M5 residual and random-effect diagnostics saved.")

# ============================================================================
# SECTION 11 — MASTER RESULTS TABLE (authoritative source of truth) + MANIFEST
# One table with sample sizes, M1/M4 position + residual variances, VPCs,
# additive explanation, gradients, endpoints, and geographic results.
# ============================================================================
grab <- function(r, d){
  dd<-d; dd$position_id<-as.integer(factor(dd$position))
  rawvar <- if("hsCRP_cap" %in% names(dd)) "hsCRP_cap" else "hsCRP"
  m1<-lmer(reformulate("(1|position_id)", if("log_hsCRP_cap" %in% names(dd)) "log_hsCRP_cap" else "log_hsCRP_w"),
           data=dd, REML=TRUE, control=ctrl)
  gr <- dd %>% group_by(position) %>% summarise(m=median(.data[[rawvar]]),n=n(),.groups="drop") %>%
    filter(n>=CELL_MIN_REPORT) %>% arrange(m)
  c1<-vc(m1); c4<-vc(r$m4)
  data.frame(
    spec=NA, N=r$N,
    M1_pos_var=unname(c1["position_id"]), M1_resid_var=unname(c1["Residual"]),
    null_VPC=r$vpc1,
    M4_pos_var=unname(c4["position_id"]), M4_resid_var=unname(c4["Residual"]),
    adj_VPC=r$vpc4, adj_additive_expl=r$addexpl,
    gradient_fold=r$grad, lowest_pos=gr$position[1], highest_pos=gr$position[nrow(gr)],
    lowest_median=gr$m[1], highest_median=gr$m[nrow(gr)])
}
master <- bind_rows(
  cbind(spec="primary_0_10",  grab(r_prim,  samp_primary)[-1]),
  cbind(spec="sens_0_20",     grab(r_s20,   samp_s20)[-1]),
  cbind(spec="broad_0_200",   grab(r_broad, samp_broad)[-1]),
  cbind(spec="capped_at_10",  grab(r_capped,samp_capped)[-1]))
master_path <- file.path(EXP, paste0("MASTER_results_", stamp, ".csv"))
write.csv(master, master_path, row.names=FALSE)
geo_master <- data.frame(
  spec="geographic_M5", N=nrow(geo), n_zip3=n_distinct(geo$zip3),
  n_sites=n_distinct(geo$site), zip_VPC=zip_vpc,
  zip_variance=unname(c5["zip3"]), position_variance=unname(c5["position_id"]),
  residual_variance=unname(c5["Residual"]), site_adjusted_zip_variance=unname(cs["zip3"]),
  influential_zip3=ext, exclusion_change_percent=infl_pct,
  profile_zip_variance_lo=zip_var_ci[1], profile_zip_variance_hi=zip_var_ci[2],
  parametric_boot_zip_variance_lo=zip_sd_boot_ci[1]^2,
  parametric_boot_zip_variance_hi=zip_sd_boot_ci[2]^2)
geo_master_path <- file.path(EXP, paste0("MASTER_geographic_results_", stamp, ".csv"))
write.csv(geo_master, geo_master_path, row.names=FALSE)
say(""); say("### MASTER RESULTS TABLE -> ", master_path)
for(i in 1:nrow(master)) say(sprintf("   %-14s N=%s null=%s adj=%s addexpl=%s grad=%s",
  master$spec[i], master$N[i],
  ifelse(is.na(master$null_VPC[i]),"-",sprintf("%.3f",master$null_VPC[i])),
  sprintf("%.3f",master$adj_VPC[i]),
  ifelse(is.na(master$adj_additive_expl[i]),"-",sprintf("%.1f",master$adj_additive_expl[i])),
  ifelse(is.na(master$gradient_fold[i]),"-",sprintf("%.1f",master$gradient_fold[i]))))

# ============================================================================
# EXPORT MANIFEST
# ============================================================================
say(""); say("### EXPORT MANIFEST (", EXP, ")")
for (f in list.files(EXP, pattern=stamp)) say("   ", f)
finalize_audits()   # STOPS here if any enforced gate failed
close(LOG)
LOG <- NULL
cat("\n\nDEFINITIVE PIPELINE COMPLETE (self-contained). Exports in", EXP, "\n")
cat("Geo rebuilt from CDR; Section 10 executed; audit gates enforced.\n")
