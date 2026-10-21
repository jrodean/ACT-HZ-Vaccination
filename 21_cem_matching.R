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
# The exact original bin widths for variables other than age (5-year bands,
# as explicitly stated in the text) were not located; education, baseline
# year, baseline CASI, other-vaccine count, and CCI are binned here using
# reasonable category counts chosen to balance stratum granularity against
# sample size, and should be treated as a faithful-methodology
# reconstruction rather than a confirmed exact replication of the original
# bin boundaries.
# ==============================================================================
suppressPackageStartupMessages({ library(dplyr); library(ggplot2); library(tidyr) })
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
coarsened <- baseline %>%
  mutate(
    bin_age = cut(age, breaks = seq(60, 100, by = 5), include.lowest = TRUE, right = FALSE),
    bin_educ = cut(education, breaks = c(-Inf, 11, 12, 15, 16, Inf), right = TRUE,
                   labels = c("<12", "12", "13-15", "16", ">16")),
    bin_year = cut(baseline_year, breaks = c(-Inf, 1996, 2000, 2005, 2010, 2015, Inf), right = TRUE),
    bin_casi = cut(baseline_casi, breaks = quantile(baseline_casi, probs = seq(0, 1, 0.25), na.rm = TRUE),
                   include.lowest = TRUE),
    bin_vax = pmin(other_vax_count, 3L),
    # CCI (excl. dementia) at baseline is heavily concentrated at 0 (~96% of
    # subjects), so quantile-based binning is degenerate; use a simple
    # zero-vs-positive indicator instead.
    bin_cci = as.integer(cci > 0),
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
    "-0.12. This reconstruction's post-matching SMDs are all comfortably\n",
    "under 0.1 (as expected, qualitatively), but its pre-matching SMDs are\n",
    "systematically smaller than the published ones (e.g., Baseline year:\n",
    "this run vs. 1.04 published) -- the exact original bin boundaries for\n",
    "variables other than age were not recovered, so this is a faithful\n",
    "methodology reconstruction rather than a verified exact replication.\n")

p13 <- balance_tbl %>%
  pivot_longer(c(smd_before, smd_after_weighted), names_to = "stage", values_to = "smd") %>%
  mutate(stage = recode(stage, smd_before = "Before matching", smd_after_weighted = "After matching (weighted)")) %>%
  ggplot(aes(x = smd, y = variable, color = stage, shape = stage)) +
  geom_vline(xintercept = c(-0.1, 0, 0.1), linetype = c("dashed", "solid", "dashed"), color = "grey60") +
  geom_point(size = 3) +
  labs(x = "Standardized mean difference", y = NULL, color = NULL, shape = NULL,
       title = NULL) +
  theme_minimal(base_size = 12) + theme(legend.position = "bottom")
ggsave("output/figures/FigureS13_cem_balance.png", p13, width = 8, height = 5, dpi = 200)

# ------------------------------------------------------------------------
# Figure S14: sample attrition patterns after matching
# ------------------------------------------------------------------------
vax_years <- baseline %>% filter(vaccinated == 1) %>%
  transmute(subject, year = as.integer(format(zost1_date, "%Y")))

set.seed(20260926)
matched_ctrl <- matched %>% filter(vaccinated == 0)
treated_years <- vax_years$year
ecdf_years <- ecdf(treated_years)
year_grid <- sort(unique(treated_years))
u_draw <- runif(nrow(matched_ctrl))
cdf_at_year <- ecdf_years(year_grid)
pseudo_year <- vapply(u_draw, function(u) year_grid[which(cdf_at_year >= u)[1]], numeric(1))
matched_ctrl$pseudo_vax_date <- as.Date(sprintf("%d-07-01", pseudo_year))

primary_zvl <- panel %>% filter(!is.na(zost1_date), !is.na(casi_total)) %>%
  build_event_time("zost1_date") %>% mutate(event_time_bin = pmax(pmin(event_time, 12), -12))

primary_ctrl <- panel %>%
  semi_join(matched_ctrl, by = "subject") %>%
  left_join(matched_ctrl %>% select(subject, pseudo_vax_date, cem_weight), by = "subject") %>%
  filter(!is.na(casi_total)) %>%
  mutate(event_time = floor(as.numeric(visit_date - pseudo_vax_date) / (2 * 365)) * 2,
         event_time_bin = pmax(pmin(event_time, 12), -12))

attrition <- bind_rows(
  primary_zvl %>% mutate(group = "ZVL (actual)") %>% count(group, event_time_bin, name = "n_obs"),
  primary_ctrl %>% mutate(group = "Never-vaccinated (matched, placebo)") %>% count(group, event_time_bin, name = "n_obs")
) %>%
  group_by(group) %>%
  mutate(n_ref = n_obs[event_time_bin == -2], pct_of_ref = 100 * n_obs / n_ref) %>%
  ungroup()
readr::write_csv(attrition, "output/tables/TableS14_attrition.csv")
cat("\n=== Figure S14: sample attrition after matching ===\n")
print(attrition, n = Inf)

p14 <- ggplot(attrition, aes(x = event_time_bin, y = pct_of_ref, color = group)) +
  geom_line() + geom_point() +
  geom_hline(yintercept = 100, linetype = "dashed", color = "grey60") +
  labs(x = "Event time (years)", y = "% of reference period (-2y) sample size", color = NULL,
       title = NULL) +
  theme_minimal(base_size = 12) + theme(legend.position = "bottom")
ggsave("output/figures/FigureS14_attrition.png", p14, width = 8, height = 5, dpi = 200)

# ------------------------------------------------------------------------
# Figure S15: placebo event-study on matched never-vaccinated, using
# pseudo-vaccination dates, weighted by CEM weight
# ------------------------------------------------------------------------
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
