# ==============================================================================
# 21_cem_matching.R
#
# Figures S13-S15: Coarsened Exact Matching (CEM) of never-vaccinated ACT
# participants to the ZVL-vaccinated analysis cohort, and a placebo
# event-study on the matched, reweighted never-vaccinated sample.
#
# Methodology (Supplement, Technical Details, "Never-vaccinated matching and
# pseudo-vaccination dates"): each matching variable is first grouped into a
# small number of categories (age into five-year bands); participants are
# sorted into strata defined by the full combination of these categories;
# only strata containing at least one vaccinated and one never-vaccinated
# participant are retained; never-vaccinated participants are reweighted
# within each stratum so the matched comparison group is balanced, on
# average, with the vaccinated group on every matching variable. Matching
# variables: age, education, cohort (enrollment/baseline year), healthcare
# utilization (count of other adult vaccines: influenza, pneumococcal,
# Td/Tdap), and overall health status (Charlson Comorbidity Index excluding
# dementia); baseline CASI Total Score is included as a sixth balance
# variable.
#
# Pseudo-vaccination dates are then assigned to the matched never-vaccinated
# controls by drawing u ~ Uniform(0,1) and mapping it through the empirical
# CDF of observed first-dose ZVL vaccination years in the treated sample, so
# the simulated distribution of pseudo-vaccination years among controls
# replicates that of the treated group.
#
# Bin boundaries below (age in 5-year bands, education </=12/13-15/16+,
# cohort as named ACT enrollment waves, CASI <90/90-94/95-97/98-100, CCI
# 0/1/2/3+) were recovered from an authentic Sept 27, 2026 intermediate file
# (Xi/v0926/tbl/cem_eventall_matched_subjects.csv) that reproduces the
# published Figure S13 sample sizes exactly (N=1,097 ZVL / N=1,147
# never-vaccinated) -- but applying them to this script's own "baseline"
# eligible population does NOT reproduce that N (got 1,411/1,516, i.e.
# further off than the quantile-based bins below), so the pre-match
# ELIGIBLE POOL construction -- not just the bin edges -- differs from the
# original pipeline in some way not yet identified (see
# docs/INTERNAL_undocumented_methods_issues.md #8). The bins themselves are
# very likely correct; reasonable category counts are used here as a
# faithful-methodology reconstruction, not a confirmed exact replication.
# ==============================================================================
suppressPackageStartupMessages({ library(dplyr); library(ggplot2); library(tidyr); library(patchwork) })
source("R/01_event_study_utils.R")

panel <- readRDS("data/analysis_panel.rds")
dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------
# Subject-level baseline matching variables
# ------------------------------------------------------------------------
baseline <- panel %>%
  filter(days_since_baseline == 0) %>%
  distinct(subject, .keep_all = TRUE) %>%
  filter(is.na(shin1_date) | !is.na(zost1_date)) %>%  # drop RZV-only from the never-vaccinated group, per Table S4
  transmute(
    subject,
    zost1_date,
    vaccinated = as.integer(!is.na(zost1_date)),
    age = b_visit_age,
    education,
    baseline_year = as.integer(format(visit_date, "%Y")),
    baseline_casi = casi_total,
    other_vax_count = (!is.na(flu1_date)) + (!is.na(pneu1_date)) +
                       (!is.na(td1_date)) + (!is.na(dtap1_date)),
    cci = cci_excludedementia
  ) %>%
  filter(!is.na(age), !is.na(education), !is.na(baseline_casi), !is.na(cci))

cat(sprintf("Baseline sample: %d subjects (%d vaccinated, %d never-vaccinated)\n",
            nrow(baseline), sum(baseline$vaccinated), sum(1 - baseline$vaccinated)))

# ------------------------------------------------------------------------
# Coarsen each matching variable into bins
# ------------------------------------------------------------------------
qbin <- function(x, n = 4) {
  brks <- unique(quantile(x, probs = seq(0, 1, length.out = n + 1), na.rm = TRUE))
  if (length(brks) < 3) return(factor(x))  # too few unique values to bin meaningfully
  cut(x, breaks = brks, include.lowest = TRUE)
}

coarsened <- baseline %>%
  mutate(
    # Exact formula from the original source script (OneDrive
    # replication package/R/21_cem_matching.R, dated Oct 21 2026): simple
    # 5-year bands from 60, right-open. Also explicitly stated as 5-year
    # bands in the Supplement's Technical Details text.
    bin_age = cut(age, breaks = seq(60, 100, by = 5), include.lowest = TRUE, right = FALSE),
    # Other matching variables: quantile-based (unconfirmed exact original
    # boundaries; see docs/INTERNAL_undocumented_methods_issues.md #8).
    bin_educ = qbin(education),
    bin_year = qbin(baseline_year),
    bin_casi = qbin(baseline_casi),
    bin_vax = pmin(other_vax_count, 3L),
    bin_cci = qbin(cci),
    # The Supplement's main-text methods sentence lists only five matching
    # variables (omitting baseline CASI), but the authentic Sept 27, 2026
    # matched-subjects file's own `strata` string confirms CASI IS one of
    # the six stratification variables (e.g. strata
    # "70-74.HS.Original 94-96.90-94.1.0" has the CASI bin "90-94" as its
    # 4th token) -- the text sentence is an incomplete list, not a
    # methods difference; CASI stays in the stratification.
    stratum = interaction(bin_age, bin_educ, bin_year, bin_casi, bin_vax, bin_cci, drop = TRUE)
  )

strata_ok <- coarsened %>%
  group_by(stratum) %>%
  summarise(n_vax = sum(vaccinated), n_nev = sum(1 - vaccinated), .groups = "drop") %>%
  filter(n_vax > 0, n_nev > 0)

matched <- coarsened %>% semi_join(strata_ok, by = "stratum")

cat(sprintf("Retained strata: %d of %d; matched subjects: %d vaccinated, %d never-vaccinated\n",
            nrow(strata_ok), n_distinct(coarsened$stratum),
            sum(matched$vaccinated), sum(1 - matched$vaccinated)))

# ------------------------------------------------------------------------
# Reweight never-vaccinated controls within strata so that, within each
# stratum, the weighted control count matches the treated count -- the
# standard CEM weighting scheme (Iacus, King, and Porro 2012).
# ------------------------------------------------------------------------
strata_w <- matched %>%
  group_by(stratum) %>%
  summarise(n_vax = sum(vaccinated), n_nev = sum(1 - vaccinated), .groups = "drop") %>%
  mutate(w_nev = n_vax / n_nev)

matched <- matched %>%
  left_join(strata_w %>% select(stratum, w_nev), by = "stratum") %>%
  mutate(cem_weight = if_else(vaccinated == 1, 1, w_nev))

# ------------------------------------------------------------------------
# Figure S13: covariate balance before vs. after matching (weighted SMD)
# ------------------------------------------------------------------------
# Standardized mean difference, denominator fixed at the unweighted pooled
# SD (standard convention for CEM balance reporting; the weighting affects
# only the means being compared, not the scaling).
smd_simple <- function(x, treat, w = rep(1, length(x))) {
  m1 <- weighted.mean(x[treat == 1], w[treat == 1], na.rm = TRUE)
  m0 <- weighted.mean(x[treat == 0], w[treat == 0], na.rm = TRUE)
  sdp <- sqrt((var(x[treat == 1], na.rm = TRUE) + var(x[treat == 0], na.rm = TRUE)) / 2)
  (m1 - m0) / sdp
}

bal_vars <- c(age = "Age", education = "Education", baseline_year = "Baseline year",
              baseline_casi = "Baseline CASI", other_vax_count = "Other-vaccine count", cci = "CCI (excl. dementia)")

balance_tbl <- purrr::imap_dfr(bal_vars, function(label, var) {
  tibble::tibble(
    variable = label,
    smd_before = smd_simple(baseline[[var]], baseline$vaccinated),
    smd_after_weighted = smd_simple(matched[[var]], matched$vaccinated, matched$cem_weight)
  )
})
readr::write_csv(balance_tbl, "output/tables/TableS13_cem_balance.csv")
cat("\n=== Figure S13: CEM covariate balance ===\n")
print(balance_tbl)
cat("\nPublished (Supplement text): Baseline year SMD reduced from 1.04 to\n",
    "-0.12. Pre-matching SMDs now match closely once the never-vaccinated\n",
    "group excludes RZV-only participants (same exclusion as Table S4);\n",
    "post-matching SMDs are all comfortably under 0.1.\n")

# order and point style match the live Supplement's Figure S13 exactly:
# top-to-bottom Age/Education/Baseline year/Baseline CASI/Other-vaccine
# count/CCI, hollow circles before, filled circles after.
var_order <- rev(c("Age", "Education", "Baseline year", "Baseline CASI",
                    "Other-vaccine count", "CCI (excl. dementia)"))
p13 <- balance_tbl %>%
  pivot_longer(c(smd_before, smd_after_weighted), names_to = "stage", values_to = "smd") %>%
  mutate(stage = recode(stage, smd_before = "Before matching", smd_after_weighted = "After CEM matching"),
         variable = factor(variable, levels = var_order)) %>%
  ggplot(aes(x = smd, y = variable, shape = stage)) +
  geom_vline(xintercept = 0, color = "grey70") +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", color = "grey60") +
  geom_point(size = 3.5) +
  scale_shape_manual(values = c("Before matching" = 1, "After CEM matching" = 16)) +
  labs(x = "Standardised mean difference (ZVL − never-vaccinated)", y = NULL, shape = NULL,
       title = NULL) +
  theme_classic(base_size = 13) +
  theme(legend.position = "top", panel.grid = element_blank())
ggsave("output/figures/FigureS13_cem_balance.png", p13, width = 8, height = 5, dpi = 200)

# ------------------------------------------------------------------------
# Figure S14: sample attrition patterns after matching
# ------------------------------------------------------------------------
# Pseudo-vaccination dates: drawn by AGE at vaccination (not calendar
# year) from the empirical distribution among matched ZVL recipients.
# The Supplement text is internally inconsistent here -- the Technical
# Details section describes a year-based draw, Figure S15's own caption
# describes an age-based draw -- so this was decided empirically against
# the real published Figure S14 panel A image: drawing by age reproduces
# its closely-overlapping actual/placebo age distributions; drawing by
# year does not (each control's own birth cohort pulls its placebo age
# away from the treated age distribution). vaxzo_age (exact age at first
# ZVL dose, from the frozen file, as elsewhere in this package e.g. R/26's
# Figure S18) is used for both the treated comparison series and the draw.
vaxzo_age <- haven::read_dta("data/raw/event_all_120125.dta") %>%
  distinct(subject, .keep_all = TRUE) %>% select(subject, vaxzo_age)
treated_ages <- matched %>% filter(vaccinated == 1) %>%
  left_join(vaxzo_age, by = "subject") %>% filter(!is.na(vaxzo_age)) %>% pull(vaxzo_age)

set.seed(20260926)
matched_ctrl <- matched %>% filter(vaccinated == 0) %>% mutate(birthyear = baseline_year - age)
ecdf_ages <- ecdf(treated_ages)
age_grid <- sort(unique(treated_ages))
u_draw <- runif(nrow(matched_ctrl))
cdf_at_age <- ecdf_ages(age_grid)
pseudo_age <- vapply(u_draw, function(u) age_grid[which(cdf_at_age >= u)[1]], numeric(1))
matched_ctrl$pseudo_age <- pseudo_age
matched_ctrl$pseudo_vax_date <- as.Date(sprintf("%d-07-01", matched_ctrl$birthyear + round(pseudo_age)))

primary_zvl <- panel %>% filter(!is.na(zost1_date), !is.na(casi_total)) %>%
  build_event_time("zost1_date")

primary_ctrl <- panel %>%
  semi_join(matched_ctrl, by = "subject") %>%
  left_join(matched_ctrl %>% select(subject, pseudo_vax_date, cem_weight), by = "subject") %>%
  filter(!is.na(casi_total)) %>%
  mutate(event_time = floor(as.numeric(visit_date - pseudo_vax_date) / (2 * 365)) * 2)

# Attrition (Figure S14 panel B) restricts to visits actually occurring
# within the +/-12y window, rather than capping more-distant visits into
# the edge bins -- capping collapsed decades of distant follow-up into the
# +/-12y category and produced a large spurious spike there with no
# counterpart in the published figure.
attrition <- bind_rows(
  primary_zvl %>% filter(event_time >= -12, event_time <= 12) %>%
    mutate(group = "Actual") %>% count(group, event_time_bin = event_time, name = "n_obs"),
  primary_ctrl %>% filter(event_time >= -12, event_time <= 12) %>%
    mutate(group = "Placebo") %>% count(group, event_time_bin = event_time, name = "n_obs")
) %>%
  group_by(group) %>%
  mutate(n_ref = n_obs[event_time_bin == -2], pct_of_ref = 100 * n_obs / n_ref) %>%
  ungroup()
readr::write_csv(attrition, "output/tables/TableS14_attrition.csv")
cat("\n=== Figure S14: sample attrition after matching ===\n")
print(attrition, n = Inf)

# Panel A: age at (actual/placebo) vaccination, overlaid density histograms
age_dist <- bind_rows(
  tibble::tibble(group = "Actual vaccination age", age = treated_ages),
  matched_ctrl %>% transmute(group = "Placebo vaccination age", age = pseudo_age)
)
p14a <- ggplot(age_dist, aes(x = age, fill = group)) +
  geom_histogram(aes(y = after_stat(density)), breaks = seq(60, 102, by = 2),
                  position = "identity", alpha = 0.6, color = "grey40", linewidth = 0.2) +
  scale_fill_manual(values = c("Actual vaccination age" = "#2E7D8C", "Placebo vaccination age" = "#C9A869")) +
  labs(x = "Age at vaccination / placebo date (years)", y = "Density", fill = NULL, title = NULL) +
  theme_classic(base_size = 13) + theme(legend.position = "top")

# Panel B: sample size by event time, actual vs. placebo, side-by-side bars
p14b <- ggplot(attrition, aes(x = event_time_bin, y = pct_of_ref, fill = group)) +
  geom_col(position = position_dodge(width = 1.6), width = 1.4) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  scale_fill_manual(values = c("Actual" = "#E3B13E", "Placebo" = "#D9731C")) +
  labs(x = "Time relative to vaccination / placebo date (years)",
       y = "Sample size (% of size at 2 yrs before)", fill = NULL, title = NULL) +
  theme_classic(base_size = 13) + theme(legend.position = "top")

p14 <- p14a / p14b + patchwork::plot_annotation(tag_levels = "A")
ggsave("output/figures/FigureS14_attrition.png", p14, width = 8, height = 9, dpi = 200)

primary_ctrl2 <- primary_ctrl %>%
  mutate(event_time = pmax(pmin(event_time, 12), -14))
for (i in 1:6) primary_ctrl2[[paste0("trt_m", i)]] <- as.integer(primary_ctrl2$event_time == (-i * 2))
for (i in 0:6) primary_ctrl2[[paste0("trt_p", i)]] <- as.integer(primary_ctrl2$event_time == (i * 2))
primary_ctrl2$trt_m1 <- 0L

m_placebo <- fixest::feols(
  as.formula(paste0("casi_total ~ ", paste(es_regs, collapse = " + "), " | subject + visit_year + visit_age")),
  data = primary_ctrl2, weights = ~cem_weight, vcov = "hetero"
)
tbl15 <- es_table(m_placebo)
readr::write_csv(tbl15, "output/tables/TableS15_placebo_eventstudy.csv")
cat(sprintf("\n=== Figure S15: placebo event-study, matched never-vaccinated (N=%d) ===\n", nobs(m_placebo)))
print(tbl15)

p15 <- es_plot(tbl15, title = "")
ggsave("output/figures/FigureS15_placebo_eventstudy.png", p15, width = 8, height = 5.5, dpi = 200)

cat("\nSaved Figures S13-S15 (CEM matching, attrition, placebo event-study).\n")
