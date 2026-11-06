# ==============================================================================
# 28_table_s3_generalizability.R
#
# Table S3: ACT ZVL-vaccinated analysis cohort vs. US adults 65+. Ported from
# Xi/v0926/revision_xi.R PART 5 (reviewer generalizability comment), on the
# frozen data/raw/event_all_120125.dta -- "White" here uses raceeth_c (not
# race_c, which Table 1/S2 use), matching that script exactly. US-population
# figures are external published reference statistics (Census/NHANES/CDC/
# NCHS/etc.), not derived from this package's own data.
# ==============================================================================
suppressPackageStartupMessages({ library(dplyr); library(haven) })

dat <- read_dta("data/raw/event_all_120125.dta") %>%
  mutate(zost1_date = as.Date(zost1_date, origin = "1960-01-01")) %>%
  filter(!is.na(zost1_date)) %>%
  mutate(visit_age_num = as.numeric(visit_age))

bl <- dat %>% arrange(subject, visit_year) %>% group_by(subject) %>% slice(1) %>% ungroup()
Nb  <- nrow(bl)
pct <- function(x) round(100 * mean(x, na.rm = TRUE), 1)
gv  <- function(v) if (v %in% names(bl)) bl[[v]] else rep(NA, Nb)

gen_tab <- tibble::tribble(
  ~characteristic,                              ~ACT_cohort,                                   ~US_65plus,        ~ref,
  "N (ZVL-vaccinated analysis cohort)",         as.character(Nb),                             "55.8 million",   "[1]",
  "Age at enrolment, mean (years)",             sprintf("%.1f", mean(gv("visit_age_num"), na.rm=TRUE)), "-",     "-",
  "Age 65-74, %",                               sprintf("%.1f", pct(gv("visit_age_num") < 75)),        "59.3",   "[1]",
  "Age 75-84, %",                               sprintf("%.1f", pct(gv("visit_age_num") >= 75 & gv("visit_age_num") < 85)), "29.3", "[1]",
  "Age 85+, %",                                 sprintf("%.1f", pct(gv("visit_age_num") >= 85)),       "11.4",   "[1]",
  "Female, %",                                  sprintf("%.1f", pct(gv("gender") == 2)),               "55.2",   "[1]",
  "White, %",                                   sprintf("%.1f", pct(gv("raceeth_c") == 7)),            "76.6 (White alone); ~72-75 non-Hispanic White", "[1]",
  "Black, %",                                   sprintf("%.1f", pct(gv("raceeth_c") == 2)),            "9.2",    "[1]",
  "Asian, %",                                   sprintf("%.1f", pct(gv("raceeth_c") == 3)),            "4.5",    "[1]",
  "Hispanic or Latino, %",                      "<2 (not separately coded in ACT)",                   "8.8",    "[1]",
  "Bachelor's degree or higher, %",             sprintf("%.1f", pct(gv("educat4") == 3)),              "~32",    "[2]",
  "Less than high school, %",                   sprintf("%.1f", pct(gv("educat4") == 0)),              "~13",    "[2]",
  "Education, mean (years)",                    sprintf("%.1f", mean(gv("education"), na.rm=TRUE)),    "-",      "-",
  "Diagnosed diabetes, %",                      sprintf("%.1f", pct(gv("diabetes") == 1)),             "24.4 (diagnosed); 29.2 total", "[3]",
  "Hypertension, %",                            sprintf("%.1f", pct(gv("hypertension") == 1)),         "74.5 (measured/treated, 60+)", "[4]",
  "Stroke history, %",                          sprintf("%.1f", pct(gv("stroke") == 1)),               "~6-8",   "[9]",
  "Any cancer history, %",                      sprintf("%.1f", pct(gv("cancer") == 1)),               "~15-19", "[10]",
  "COPD, %",                                    sprintf("%.1f", pct(gv("copd") == 1)),                 "~11-13", "[11]",
  "Obesity (BMI >= 30), %",                     sprintf("%.1f", pct(gv("bmi") >= 30)),                 "42.8 (60+)", "[5]",
  "Body-mass index, mean",                      sprintf("%.1f", mean(gv("bmi"), na.rm=TRUE)),          "~28-29", "[5]",
  "Any ADL limitation, %",                      sprintf("%.1f", pct(gv("adl_flag") == 1)),             "~20-25 (any disability ~33)", "[12]",
  "APOE-e4 carrier (>=1 e4 allele), %",         sprintf("%.1f", pct(gv("apoe_any4") > 0)),             "~17.5-26.9 (US White, by European descent)", "[6]",
  "Baseline CASI Total Score, mean",            sprintf("%.1f", mean(gv("casi_total"), na.rm=TRUE)),   "n/a (dementia-free at entry)", "-",
  "Alzheimer's dementia (age 65+)",             "0 at entry (eligibility)",                            "10.9",   "[7]",
  "Life expectancy at 65 (years)",              "-",                                                  "18.9",   "[8]",
  "In a nursing / skilled-nursing facility, %", "0 at entry (eligibility)",                            "~2.5",   "[1]",
  "Health-plan coverage, %",                    "100 (continuous KPWA enrolment >=2 y required)",      "~99 (Medicare)", "-"
)

dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)
readr::write_csv(gen_tab, "output/tables/TableS3_generalizability.csv")
print(gen_tab, n = Inf)
