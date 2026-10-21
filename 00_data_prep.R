# ==============================================================================
# 00_data_prep.R
#
# Builds the analysis-ready panel used throughout this replication package,
# starting ONLY from the three raw ACT Study extracts provided to the
# research team (data/raw/ACT_296_{person,visit,event}_table_2025_12_01.csv).
#
# No pre-built "event_all_*.dta" file is read anywhere in this package: this
# script IS the data-construction step. (During development, the output of
# this script was checked against the frozen event_all_120125.dta that was
# used to produce the submitted manuscript -- see 00b_validate_data_prep.R --
# but that file is not a dependency of the pipeline itself.)
#
# Output: data/analysis_panel.rds (one row per subject x ACT visit)
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(lubridate)
  library(stringr)
})

raw_dir <- file.path("data", "raw")

# ------------------------------------------------------------------------
# 1) Load raw extracts
# ------------------------------------------------------------------------
visit <- read_csv(file.path(raw_dir, "ACT_296_visit_table_2025_12_01.csv"),
                   guess_max = 50000, show_col_types = FALSE)
person <- read_csv(file.path(raw_dir, "ACT_296_person_table_2025_12_01.csv"),
                    guess_max = 10000, show_col_types = FALSE)
event <- read_csv(file.path(raw_dir, "ACT_296_event_table_2025_12_01.csv"),
                   guess_max = 50000, show_col_types = FALSE)

names(visit)  <- tolower(names(visit))
names(person) <- tolower(names(person))
names(event)  <- tolower(names(event))

# ------------------------------------------------------------------------
# 2) Visit-level core panel
#
# The raw visit table arrives from the ACT Data Coordinating Center with
# CASI item/subscale/total scores, CDR domain scores, DSM-IV dementia
# diagnosis, and APOE genotype already derived -- these are NOT recomputed
# here, only renamed/retained.
# ------------------------------------------------------------------------
visit_core <- visit %>%
  transmute(
    subject          = subject,
    visit            = visit,
    visit_date       = as.Date(visit_date, format = "%m/%d/%Y"),
    visit_year       = visit_year,
    visit_age        = visit_age,
    days_since_baseline = days_since_baseline,
    cohort           = cohort,

    # demographics (repeated at every visit; collapsed to one value/subject below)
    sex_c            = sex_c,
    race_c           = race_c,
    hispanic_c       = hispanic_c,
    raceeth_c        = raceeth_c,
    education        = education,
    college          = college,
    educat4          = educat4,

    # LASSO dementia-risk-score features (see 00c_lasso_risk_score.R) --
    # raw columns, collapsed to baseline (first/mean/last pre-vax visit) there
    cesd20_score     = cesd20_score,
    iadl_sum         = iadl_sum,
    income           = income,
    basicneeds_met   = basicneeds_met,
    exercise_reg     = exercise_reg,
    diabetes         = diabetes,
    chf              = chf,
    mi               = mi,
    hypertension     = hypertension,
    stroke           = stroke,
    sleepapnea       = sleepapnea,
    rheumarth        = rheumarth,
    lupus            = lupus,
    othautoimm       = othautoimm,
    organbmt         = organbmt,
    chemo            = chemo,
    othimmsup        = othimmsup,
    smoke_bl         = smoke_bl,
    smoke_pack_bl    = smoke_pack_bl,
    alcohol_past_bl  = alcohol_past_bl,
    avesys           = avesys,
    avedia           = avedia,
    bmi              = bmi,
    gait_time        = gait_time,
    grip_strength    = grip_strength,

    # CASI total + 9 domain subscales
    casi_valid       = casi_valid,
    casi_total       = casi_total,
    casi_att         = casi_att,
    casi_mmc         = casi_mmc,
    casi_ori         = casi_ori,
    casi_ltm         = casi_ltm,
    casi_stm         = casi_stm,
    casi_lang        = casi_lang,
    casi_draw        = casi_draw,
    casi_flu         = casi_flu,
    casi_ataj        = casi_ataj,
    remote_casi_status = remote_casi_status,

    # CDR global + 6 domains
    cdr_total        = cdr_total,
    cdr_memory       = cdr_memory,
    cdr_orient       = cdr_orient,
    cdr_judge        = cdr_judge,
    cdr_comm         = cdr_comm,
    cdr_home         = cdr_home,
    cdr_pers         = cdr_pers,

    dsmiv_dx         = dsmiv_dx,

    apoe_genotype    = apoe_genotype,
    apoe_any4        = apoe_any4,

    # frailty / ADL (used for Figure 4 panel C and Table S "frailty" split)
    adl_walk = adl_walk, adl_bed = adl_bed, adl_feed = adl_feed,
    adl_dress = adl_dress, adl_bath = adl_bath, adl_toilet = adl_toilet,
    adl_sum  = adl_sum,
    adl_flag = adl_flag,

    # family history (Figure 4 panels A/B)
    firstdeg_memprob = firstdeg_memprob,
    firstdeg_dem     = firstdeg_dem,
    firstdeg_ad      = firstdeg_ad,
    firstdeg_vasc    = firstdeg_vasc,

    # self-reported shingles history at/before baseline (if present on form)
    frm51 = frm51,

    # CESD depression score, exercise covariate (used in several HTE/
    # robustness splits and the LASSO risk-score feature set)
    cesd20_score = cesd20_score,
    exercise_reg = exercise_reg
  )

# ------------------------------------------------------------------------
# 3) Person-level lookback comorbidities (lb_*) -- baseline covariates used
#    for Charlson-style indices / "adjust for time-varying health" check
# ------------------------------------------------------------------------
person_core <- person  # already one row per subject; keep as-is

# ------------------------------------------------------------------------
# 4) Event table -> subject-level dated variables
#
# `event` is a category label, `days_to_event` is days from the subject's
# baseline visit. Multiple rows per (subject, event) occur for repeat doses
# / repeat diagnoses; some rows are exact duplicates. We rank distinct event
# days ascending to recover 1st/2nd occurrence dates.
# ------------------------------------------------------------------------
baseline_date <- visit_core %>%
  filter(days_since_baseline == 0) %>%
  distinct(subject, .keep_all = TRUE) %>%
  select(subject, baseline_date = visit_date)

nth_event_date <- function(ev_name, n_occ = 1) {
  event %>%
    filter(event == ev_name) %>%
    distinct(subject, days_to_event) %>%
    arrange(subject, days_to_event) %>%
    group_by(subject) %>%
    slice(seq_len(min(n(), n_occ))) %>%
    mutate(occ = row_number()) %>%
    ungroup() %>%
    filter(occ == n_occ) %>%
    inner_join(baseline_date, by = "subject") %>%
    transmute(subject, date = baseline_date + days_to_event)
}

make_dose_dates <- function(ev_name, prefix, n_doses = 2) {
  out <- NULL
  for (k in 1:n_doses) {
    d <- nth_event_date(ev_name, k) %>% rename(!!paste0(prefix, k, "_date") := date)
    out <- if (is.null(out)) d else full_join(out, d, by = "subject")
  }
  out
}

zost_dates  <- make_dose_dates("vax_zost", "zost", 2)
shin_dates  <- make_dose_dates("vax_shin", "shin", 2)
pneu_dates  <- make_dose_dates("vax_pneu", "pneu", 1)
td_dates    <- make_dose_dates("vax_td",   "td",   1)
flu_dates   <- make_dose_dates("vax_flu",  "flu",  1)
dtap_dates  <- make_dose_dates("vax_dtap", "dtap", 1)

# herpes zoster (shingles), herpes simplex, PHN ("pa" = pain), "ph" = post-herpetic
hz_dates <- make_dose_dates("ep_hz", "hz", 2)   # hz1_date, hz2_date (first/second shingles episode)
hs_dates <- event %>% filter(event == "ep_hs") %>% distinct(subject, days_to_event) %>%
  inner_join(baseline_date, by = "subject") %>%
  group_by(subject) %>% summarise(hs1_date = min(baseline_date + days_to_event), .groups = "drop")

# antivirals (acyclovir/famciclovir/valacyclovir) -- any use, first date
antiviral_dates <- event %>%
  filter(event %in% c("av_acycl", "av_famci", "av_valac")) %>%
  distinct(subject, days_to_event) %>%
  inner_join(baseline_date, by = "subject") %>%
  group_by(subject) %>%
  summarise(anti1_date = min(baseline_date + days_to_event), .groups = "drop")

# CCW chronic-condition first-onset dates (absorbing indicators downstream)
ccw_events <- event %>% filter(str_starts(event, "ccw_")) %>% distinct(event) %>% pull(event)
ccw_dates <- event %>%
  filter(event %in% ccw_events) %>%
  distinct(subject, event, days_to_event) %>%
  inner_join(baseline_date, by = "subject") %>%
  group_by(subject, event) %>%
  summarise(onset_date = min(baseline_date + days_to_event), .groups = "drop") %>%
  pivot_wider(names_from = event, values_from = onset_date, names_glue = "{event}_onset")

# IP hospitalization-category dates (non-absorbing, year-flow indicators downstream)
ip_events <- event %>% filter(str_starts(event, "ip_")) %>% distinct(event) %>% pull(event)
ip_dates <- event %>%
  filter(event %in% ip_events) %>%
  distinct(subject, event, days_to_event) %>%
  inner_join(baseline_date, by = "subject") %>%
  mutate(adm_date = baseline_date + days_to_event, adm_year = year(adm_date)) %>%
  select(subject, event, adm_year)

# ------------------------------------------------------------------------
# 5) Collapse demographics to one value per subject (first non-missing,
#    chronological visit order) -- mirrors Stata's (firstnm) collapse
#    convention seen in the labels of event_all_120125.dta.
# ------------------------------------------------------------------------
firstnm <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA)
  x[1]
}

demog <- visit_core %>%
  arrange(subject, visit_date) %>%
  group_by(subject) %>%
  summarise(
    cohort    = firstnm(cohort),
    sex_c     = firstnm(sex_c),
    race_c    = firstnm(race_c),
    hispanic_c = firstnm(hispanic_c),
    raceeth_c = firstnm(raceeth_c),
    education = firstnm(education),
    college   = firstnm(college),
    educat4   = firstnm(educat4),
    firstdeg_memprob = firstnm(firstdeg_memprob),
    firstdeg_dem     = firstnm(firstdeg_dem),
    firstdeg_ad      = firstnm(firstdeg_ad),
    firstdeg_vasc    = firstnm(firstdeg_vasc),
    .groups = "drop"
  ) %>%
  mutate(
    gender = sex_c,                             # 1 = male, 2 = female (ACT coding)
    # Matches HTE.do exactly: `gen byte white = (raceeth_c==7)`,
    # `gen byte college = (educat4==3)`
    white   = as.integer(raceeth_c == 7),
    college = as.integer(educat4 == 3),
  )

# ------------------------------------------------------------------------
# 6) Assemble the final subject-visit panel
# ------------------------------------------------------------------------
panel <- visit_core %>%
  select(-cohort, -sex_c, -race_c, -hispanic_c, -raceeth_c, -education, -college, -educat4,
         -firstdeg_memprob, -firstdeg_dem, -firstdeg_ad, -firstdeg_vasc) %>%
  left_join(demog, by = "subject") %>%
  left_join(zost_dates, by = "subject") %>%
  left_join(shin_dates, by = "subject") %>%
  left_join(pneu_dates, by = "subject") %>%
  left_join(td_dates, by = "subject") %>%
  left_join(flu_dates, by = "subject") %>%
  left_join(dtap_dates, by = "subject") %>%
  left_join(hz_dates, by = "subject") %>%
  left_join(hs_dates, by = "subject") %>%
  left_join(antiviral_dates, by = "subject") %>%
  left_join(person_core, by = "subject")

# ------------------------------------------------------------------------
# 6b) Drop visits recorded after a subject's study-exit event (freeze,
#     withdrawal, death, or enrollment end). A small number of raw visit
#     rows occur after a subject's documented exit date; these are
#     excluded from the analysis dataset, consistent with the dataset used
#     to produce the submitted manuscript.
# ------------------------------------------------------------------------
panel <- panel %>%
  mutate(end_days = pmin(days_to_freeze, days_to_withdraw,
                          days_to_death, days_to_enrend, na.rm = TRUE)) %>%
  filter(is.na(end_days) | days_since_baseline <= end_days) %>%
  select(-end_days)

# ------------------------------------------------------------------------
# 7) IP hospitalization flow indicators: 1 in the calendar year of any
#    admission in that category, else 0 (non-absorbing)
# ------------------------------------------------------------------------
for (ev in ip_events) {
  yrs <- ip_dates %>% filter(event == ev) %>% distinct(subject, adm_year)
  flagged <- panel %>%
    select(subject, visit_year) %>%
    inner_join(yrs, by = c("subject", "visit_year" = "adm_year")) %>%
    mutate(flag = 1L) %>%
    distinct(subject, visit_year, flag)
  colname <- ev  # e.g. "ip_myoca"
  panel <- panel %>%
    left_join(flagged, by = c("subject", "visit_year")) %>%
    mutate(!!colname := if_else(!is.na(flag), 1L, 0L)) %>%
    select(-flag)
}

# ------------------------------------------------------------------------
# 8) CCW chronic-condition absorbing indicators: 1 from first onset date
#    onward, else 0
# ------------------------------------------------------------------------
onset_cols <- grep("_onset$", names(ccw_dates), value = TRUE)
panel <- panel %>% left_join(ccw_dates, by = "subject")
for (col in onset_cols) {
  base <- sub("_onset$", "", col)
  panel[[base]] <- as.integer(!is.na(panel[[col]]) & panel$visit_date >= panel[[col]])
}
panel <- panel %>% select(-all_of(onset_cols))

# ------------------------------------------------------------------------
# 8b) Charlson Comorbidity Index (CCI), time-varying, excluding dementia --
#     negative control for Figure S8 (per Lilac Zhao's negativecontroloutcomes.R,
#     shared 2026-10-15: cci_excludedementia = max(cci - dementia, 0), then
#     min-max normalized). That script reads "cci"/"dementia" as precomputed
#     fields with no documented weighting formula; the cci_* raw event codes
#     map onto the standard 17-item Charlson/Deyo weights unambiguously by
#     abbreviation (confirmed against the standard clinical definition):
#       weight 1: cci_myoc (MI), cci_chv (CHF), cci_pvd (PVD), cci_cvd (CVD),
#                 cci_dem (dementia), cci_cpd (chronic pulmonary disease),
#                 cci_rheu (rheumatologic/connective tissue), cci_pep (peptic
#                 ulcer), cci_livm (mild liver disease), cci_dian (diabetes,
#                 no complications)
#       weight 2: cci_diaw (diabetes with complications), cci_hemi
#                 (hemiplegia/paraplegia), cci_ren (renal disease), cci_mal
#                 (any malignancy)
#       weight 3: cci_livs (moderate/severe liver disease)
#       weight 6: cci_met (metastatic solid tumor), cci_aids (AIDS)
#     Each component is absorbing (1 from first onset date onward, as for
#     CCW), and cci = sum of present components' weights at each visit.
# ------------------------------------------------------------------------
cci_weights <- c(cci_myoc = 1, cci_chv = 1, cci_pvd = 1, cci_cvd = 1, cci_dem = 1,
                  cci_cpd = 1, cci_rheu = 1, cci_pep = 1, cci_livm = 1, cci_dian = 1,
                  cci_diaw = 2, cci_hemi = 2, cci_ren = 2, cci_mal = 2,
                  cci_livs = 3, cci_met = 6, cci_aids = 6)
cci_events <- event %>% filter(event %in% names(cci_weights)) %>% distinct(event) %>% pull(event)
cci_dates <- event %>%
  filter(event %in% cci_events) %>%
  distinct(subject, event, days_to_event) %>%
  inner_join(baseline_date, by = "subject") %>%
  group_by(subject, event) %>%
  summarise(onset_date = min(baseline_date + days_to_event), .groups = "drop") %>%
  pivot_wider(names_from = event, values_from = onset_date, names_glue = "{event}_onset")

panel <- panel %>% left_join(cci_dates, by = "subject")
cci_onset_cols <- grep("_onset$", names(cci_dates), value = TRUE)
cci_flag_cols <- character(0)
for (col in cci_onset_cols) {
  comp <- sub("_onset$", "", col)
  flagcol <- paste0(comp, "_flag")
  panel[[flagcol]] <- as.integer(!is.na(panel[[col]]) & panel$visit_date >= panel[[col]])
  cci_flag_cols <- c(cci_flag_cols, flagcol)
}
panel <- panel %>% select(-all_of(cci_onset_cols))

present_flags <- intersect(paste0(names(cci_weights), "_flag"), cci_flag_cols)
present_weights <- cci_weights[sub("_flag$", "", present_flags)]
panel$cci <- as.numeric(as.matrix(panel[present_flags]) %*% present_weights)
panel$dementia_flag <- panel$cci_dem_flag  # component already computed above
panel$cci_excludedementia_raw <- pmax(panel$cci - cci_weights["cci_dem"] * panel$dementia_flag, 0)
cci_min <- min(panel$cci_excludedementia_raw, na.rm = TRUE)
cci_max <- max(panel$cci_excludedementia_raw, na.rm = TRUE)
panel$cci_excludedementia <- (panel$cci_excludedementia_raw - cci_min) / (cci_max - cci_min)
panel <- panel %>% select(-all_of(cci_flag_cols), -dementia_flag, -cci_excludedementia_raw)

# ------------------------------------------------------------------------
# 9) Baseline (pre-vaccination) carried-forward variables: b_casi_*,
#    b_visit_age, b_visit_year -- taken at each subject's FIRST ACT visit
#    (i.e., study baseline, not vaccination-relative baseline).
# ------------------------------------------------------------------------
baseline_vals <- panel %>%
  arrange(subject, visit_date) %>%
  group_by(subject) %>%
  slice(1) %>%
  ungroup() %>%
  select(subject,
         b_visit_year = visit_year, b_visit_age = visit_age,
         b_casi_total = casi_total, b_casi_att = casi_att, b_casi_mmc = casi_mmc,
         b_casi_ori = casi_ori, b_casi_ltm = casi_ltm, b_casi_stm = casi_stm,
         b_casi_lang = casi_lang, b_casi_draw = casi_draw, b_casi_flu = casi_flu,
         b_casi_ataj = casi_ataj)

panel <- panel %>% left_join(baseline_vals, by = "subject")

# ------------------------------------------------------------------------
# 10) Baseline rate of cognition decline (Figure 4 Panel D / HTE variable):
#     simple difference between a subject's first two PRE-VACCINATION ACT
#     visits (validated construction -- see replication package README and
#     Figure4D_Table_Quintile_FINAL.docx for why this replaces the original,
#     circular "prior cognitive decline" construction). Positive = CASI
#     dropped (declined) between visit 0 and visit 1.
# ------------------------------------------------------------------------
decline_dat <- panel %>%
  filter(!is.na(zost1_date), !is.na(casi_total)) %>%
  mutate(event_time = pmax(pmin(floor(as.numeric(visit_date - zost1_date) / (2*365)) * 2, 12), -14)) %>%
  filter(event_time < 0) %>%
  arrange(subject, visit_date) %>%
  group_by(subject) %>%
  slice_head(n = 2) %>%
  filter(n() == 2) %>%
  summarise(decline01 = first(casi_total) - last(casi_total), .groups = "drop")

panel <- panel %>%
  left_join(decline_dat, by = "subject") %>%
  mutate(fastdec = as.integer(decline01 >= 0))

# ------------------------------------------------------------------------
# 11) Save
# ------------------------------------------------------------------------
dir.create("data", showWarnings = FALSE)
saveRDS(panel, file.path("data", "analysis_panel.rds"))

cat(sprintf("Built analysis_panel.rds: %d rows, %d subjects, %d columns\n",
            nrow(panel), n_distinct(panel$subject), ncol(panel)))
