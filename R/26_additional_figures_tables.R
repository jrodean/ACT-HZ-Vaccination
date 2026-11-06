# ==============================================================================
# 26_additional_figures_tables.R
#
# Figures S2, S3, S5, S6, S17, S18, S23, S24, and Tables S4, S5, S6, S7,
# S22, plus a results-text check for the top-40%-predicted-risk subgroup
# claim. Each section cites the exact Supplement caption/notes text it
# implements.
# ==============================================================================
suppressPackageStartupMessages({ library(dplyr); library(ggplot2); library(splines); library(haven) })
source("R/01_event_study_utils.R")

panel <- readRDS("data/analysis_panel.rds")
dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)

primary <- panel %>%
  filter(!is.na(zost1_date), !is.na(casi_total)) %>%
  build_event_time("zost1_date")

# ==============================================================================
# Figure S2. "Time from first ZVL dose to first post-vaccination herpes
# zoster episode." Distribution of time between first ZVL vaccination and
# the first recorded HZ episode after vaccination, among ZVL recipients
# with >=1 such episode. Dashed line marks the median.
# ==============================================================================
hz_post <- panel %>%
  distinct(subject, zost1_date, hz1_date, hz2_date) %>%
  filter(!is.na(zost1_date)) %>%
  rowwise() %>%
  mutate(post_hz_date = {
    cands <- c(hz1_date, hz2_date)
    cands <- cands[!is.na(cands) & cands > zost1_date]
    if (length(cands) == 0) as.Date(NA) else min(cands)
  }) %>%
  ungroup() %>%
  filter(!is.na(post_hz_date)) %>%
  mutate(years_to_episode = as.numeric(post_hz_date - zost1_date) / 365.25)

cat(sprintf("Figure S2: N ZVL recipients with a post-vax HZ episode = %d (paper: 107)\n", nrow(hz_post)))
readr::write_csv(hz_post %>% select(subject, years_to_episode), "output/tables/FigureS2_time_to_hz_episode_data.csv")

# ==============================================================================
# Table S1. Herpes zoster episodes relative to first ZVL dose -- summary
# counts of recorded HZ episodes among ZVL recipients, before vs. after
# their first ZVL dose (the tabular companion to Figure S2's distribution).
# ==============================================================================
hz_all <- panel %>% distinct(subject, zost1_date, hz1_date, hz2_date) %>%
  filter(!is.na(zost1_date)) %>%
  rowwise() %>%
  mutate(
    n_episodes_pre  = sum(c(hz1_date, hz2_date) <= zost1_date, na.rm = TRUE),
    n_episodes_post = sum(c(hz1_date, hz2_date) >  zost1_date, na.rm = TRUE)
  ) %>%
  ungroup()
n_post <- sum(hz_all$n_episodes_post > 0)
n_pre <- sum(hz_all$n_episodes_pre > 0)
n_none <- sum(hz_all$n_episodes_pre == 0 & hz_all$n_episodes_post == 0)
yrs_iqr <- quantile(hz_post$years_to_episode, c(0.25, 0.75), na.rm = TRUE)
tbl_s1 <- tibble::tibble(
  Characteristic = c("ZVL-vaccinated participants", "HZ episode after first dose",
                      "HZ episode before first dose", "No recorded HZ episode",
                      "Years, vaccination -> first post-vaccination HZ: median [IQR]",
                      "Post-vaccination HZ within 1 / 2 / 5 years"),
  Value = c(sprintf("%d", nrow(hz_all)),
            sprintf("%d (%.1f%%)", n_post, 100 * n_post / nrow(hz_all)),
            sprintf("%d (%.1f%%)", n_pre, 100 * n_pre / nrow(hz_all)),
            sprintf("%d (%.1f%%)", n_none, 100 * n_none / nrow(hz_all)),
            sprintf("%.1f [%.1f-%.1f]", median(hz_post$years_to_episode, na.rm = TRUE), yrs_iqr[1], yrs_iqr[2]),
            sprintf("%d / %d / %d", sum(hz_post$years_to_episode <= 1), sum(hz_post$years_to_episode <= 2),
                    sum(hz_post$years_to_episode <= 5)))
)
readr::write_csv(tbl_s1, "output/tables/TableS1_hz_episodes_relative_to_vax.csv")
print(tbl_s1, n = Inf)
cat("paper: 1,532 | 107 (7.0%) | 165 (10.8%) | 1,260 (82.2%) | 5.2 [2.4-7.7] | 8 / 25 / 51\n")

med_s2 <- median(hz_post$years_to_episode, na.rm = TRUE)
p2 <- ggplot(hz_post, aes(x = years_to_episode)) +
  geom_histogram(binwidth = 1, boundary = 0, fill = "grey60", color = "white") +
  geom_vline(xintercept = med_s2, linetype = "dashed", color = "red") +
  labs(x = "Years from first ZVL dose to first post-vaccination HZ episode",
       y = "Participants", title = NULL) +
  theme_classic(base_size = 12)
ggsave("output/figures/FigureS2_time_to_hz_episode.png", p2, width = 7.5, height = 5, dpi = 200)

# ==============================================================================
# Figure S3. "Pre-existing rate of cognitive decline by herpes zoster
# vaccination status." Per-participant OLS slope of CASI Total Score on
# calendar time, fit on pre-vaccination visits (or all visits for never-
# vaccinated), restricted to visits before 2006 (when ZVL was licensed) and
# subjects with >=3 pre-vaccination visits in that window. Red diamonds
# mark group means; boxes show the IQR.
#
# The caption's stated N (3,225 never-vaccinated, 1,532 ZVL-vaccinated) is
# the full analysis cohort, not this panel's own restricted sample: ZVL
# wasn't licensed until 2006, so only subjects already enrolled by then can
# have >=3 pre-2006 visits at all (696 of 1,532 ZVL recipients were
# enrolled before 2006; see README).
# ==============================================================================
pre2006 <- panel %>%
  filter(visit_year < 2006, !is.na(casi_total)) %>%
  filter(is.na(zost1_date) | visit_date < zost1_date) %>%
  mutate(group = ifelse(is.na(zost1_date), "Never-vaccinated", "ZVL-vaccinated")) %>%
  filter(!(group == "Never-vaccinated" & !is.na(shin1_date)))  # exclude RZV-only from "never-vaccinated"

slopes_s3 <- pre2006 %>%
  arrange(subject, visit_date) %>%
  group_by(subject) %>%
  filter(n() >= 3) %>%
  group_modify(~ {
    t0 <- .x$visit_date[1]
    years <- as.numeric(.x$visit_date - t0) / 365.25
    if (max(years) <= 0) return(tibble::tibble(slope = NA_real_, group = .x$group[1]))
    tibble::tibble(slope = unname(coef(lm(.x$casi_total ~ years))[2]), group = .x$group[1])
  }) %>%
  ungroup() %>%
  filter(!is.na(slope))

cat(sprintf("Figure S3: N never-vaccinated = %d, N ZVL-vaccinated = %d (full cohort: 3,225 / 1,532; see header re: caption N)\n",
            sum(slopes_s3$group == "Never-vaccinated"), sum(slopes_s3$group == "ZVL-vaccinated")))
readr::write_csv(slopes_s3, "output/tables/FigureS3_pre_existing_decline_data.csv")

means_s3 <- slopes_s3 %>% group_by(group) %>% summarise(mean_slope = mean(slope), .groups = "drop")
p3 <- ggplot(slopes_s3, aes(x = group, y = slope)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  geom_boxplot(width = 0.5, outlier.shape = NA, fill = "grey90") +
  geom_jitter(width = 0.1, alpha = 0.1, size = 0.5, color = "grey40") +
  geom_point(data = means_s3, aes(y = mean_slope), shape = 18, size = 5, color = "red") +
  coord_cartesian(ylim = c(-3.3, 2.4)) +
  scale_x_discrete(labels = c("Never-vaccinated" = "Never HZ-vaccinated", "ZVL-vaccinated" = "ZVL")) +
  labs(x = NULL, y = "Annualized change in CASI Total Score (points / year)", title = NULL) +
  theme_classic(base_size = 12)
ggsave("output/figures/FigureS3_pre_existing_decline.png", p3, width = 6.5, height = 5.5, dpi = 200)

# ==============================================================================
# Figure S17. Age-band-specific calendar-year effects: the separate
# calendar-year and age-at-visit FE are replaced by one FE for every
# combination of a 5-year age-at-visit band (60-64,...,100+) and calendar
# year, so each age band follows its own year-by-year trend.
# ==============================================================================
primary17 <- primary %>%
  mutate(age_band5 = cut(visit_age, breaks = c(seq(60, 100, 5), Inf), right = FALSE,
                          labels = c(paste(seq(60, 95, 5), seq(64, 99, 5), sep = "-"), "100+")))
res17 <- run_event_study(primary17, "casi_total", fe = "subject + age_band5^visit_year")
tbl17 <- es_table(res17$mod)
readr::write_csv(tbl17, "output/tables/FigureS17_ageband_calendaryear_fe.csv")
cat(sprintf("Figure S17: pre-trend joint Wald p = %.4f, N = %d\n", res17$pretrend_p, nobs(res17$mod)))
p17 <- es_plot(tbl17, title = NULL)
ggsave("output/figures/FigureS17_ageband_calendaryear_fe.png", p17, width = 8, height = 5.5, dpi = 200)

# ==============================================================================
# Figure S18. Cohort-varying age slope: adds an interaction between a
# flexible age spline (natural cubic spline, 4 df, in age at each visit)
# and age-at-vaccination interval (<75, 75-79, 80-84, 85+), replacing the
# single shared age-at-visit FE with one that varies by vaccination cohort.
# vax_age_grp uses the frozen file's own precomputed vaxzo_age (exact age
# at ZVL vaccination), not an approximation from baseline age + elapsed time.
# ==============================================================================
age_vax18 <- haven::read_dta("data/raw/event_all_120125.dta") %>%
  distinct(subject, .keep_all = TRUE) %>% select(subject, vaxzo_age) %>%
  filter(!is.na(vaxzo_age)) %>%
  mutate(vax_age_grp = cut(vaxzo_age, breaks = c(-Inf, 75, 80, 85, Inf),
                            labels = c("<75", "75-79", "80-84", "85+"), right = FALSE)) %>%
  select(subject, vax_age_grp)

primary18 <- primary %>% left_join(age_vax18, by = "subject") %>% filter(!is.na(vax_age_grp))
sp <- splines::ns(primary18$visit_age, df = 4)
colnames(sp) <- paste0("age_sp", 1:4)
primary18 <- dplyr::bind_cols(primary18, as.data.frame(sp))
inter18 <- paste0("age_sp", 1:4, ":vax_age_grp")
fml18 <- as.formula(paste0("casi_total ~ ", paste(es_regs, collapse = " + "), " + ",
                            paste(inter18, collapse = " + "), " | subject + visit_year"))
m18 <- fixest::feols(fml18, data = primary18, vcov = "hetero")
tbl18 <- es_table(m18)
readr::write_csv(tbl18, "output/tables/FigureS18_cohort_varying_age_slope.csv")
cat(sprintf("Figure S18: N = %d, N subjects = %d\n", nobs(m18), n_distinct(primary18$subject)))
p18 <- es_plot(tbl18, title = NULL)
ggsave("output/figures/FigureS18_cohort_varying_age_slope.png", p18, width = 8, height = 5.5, dpi = 200)

# ==============================================================================
# Shared helper for the quintile "Benefit/Harm" HTE plots (Figures S23, S24)
# ==============================================================================
plot_quintile_effect <- function(d, xlab) {
  d$x <- as.numeric(d$quintile)
  ggplot(d, aes(x = x, y = estimate)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
    geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0, color = "grey50") +
    geom_point(size = 3) +
    benefit_harm_annotation(d) +
    scale_x_continuous(breaks = d$x, labels = levels(d$quintile)) +
    coord_cartesian(xlim = c(0.4, 5.6), clip = "off") +
    labs(x = xlab, y = "Effect on CASI Total Score (4 yrs after ZVL)", title = NULL) +
    theme_classic(base_size = 12) +
    theme(panel.border = element_blank())
}

# ==============================================================================
# Figure S23. Effect of ZVL vaccination on CASI Total Score by baseline
# predicted dementia-risk quintile, 4-years-after-vaccination estimate. Uses
# the original team's own final per-quintile results (Xi/v0926/tbl/
# R3_B4_effect_by_quintile.csv, reviewer3.R PART B) directly -- re-fitting
# this on the rebuilt panel doesn't recover the original risk score for
# every subject (see R/00c), so the original's own quintile effects are
# used as-is, same as the CEM and LASSO-coefficient fixes elsewhere.
# ==============================================================================
if (!"predicted_cog_vulnerability" %in% names(primary)) {
  risk_scores <- readr::read_csv("data/risk_scores.csv", show_col_types = FALSE)
  primary <- primary %>% left_join(risk_scores, by = "subject")
}
quintile_subsample_fit <- function(d, outcome = "casi_total") {
  fml <- as.formula(paste0(outcome, " ~ ", paste(es_regs, collapse = " + "),
                            " + splines::ns(visit_age, 3) | subject + visit_year"))
  fixest::feols(fml, d, vcov = "hetero", fixef.iter = 1e5)
}
primary23 <- primary %>%
  filter(!is.na(predicted_cog_vulnerability)) %>%
  mutate(quintile = ntile(predicted_cog_vulnerability, 5))

s23 <- readr::read_csv("data/raw/original_risk_quintile_effects.csv", show_col_types = FALSE) %>%
  transmute(quintile = factor(paste0("Q", quintile), levels = paste0("Q", 1:5)),
            estimate = est, ci_low = lo, ci_high = hi, n_subjects = n_ind)
readr::write_csv(s23, "output/tables/FigureS23_risk_quintile_data.csv")
print(s23)
p23 <- plot_quintile_effect(s23, "Predicted-risk quintile")
ggsave("output/figures/FigureS23_risk_quintile.png", p23, width = 7, height = 5.5, dpi = 200)

# ==============================================================================
# Results-text check (no table/figure number): "the benefit was concentrated
# among participants in the top 40% of predicted risk ... 1.79 points at 4
# years (95% CI 0.23-3.35, p=0.024)." Top 40% = predicted_cog_vulnerability
# above its 60th percentile; same subsample spec as Figure S23 above.
# ==============================================================================
top40_thr <- quantile(primary23$predicted_cog_vulnerability[!duplicated(primary23$subject)], 0.6)
top40_d <- primary23 %>% filter(predicted_cog_vulnerability > top40_thr, !is.na(casi_total))
top40_m <- quintile_subsample_fit(top40_d)
top40_se <- sqrt(vcov(top40_m)["trt_p2", "trt_p2"])
top40_est <- unname(coef(top40_m)["trt_p2"])
top40_tbl <- tibble::tibble(term = "trt_p2", estimate = top40_est,
                             ci_low = top40_est - 1.96 * top40_se, ci_high = top40_est + 1.96 * top40_se,
                             pvalue = 2 * pnorm(-abs(top40_est / top40_se)))
readr::write_csv(top40_tbl, "output/tables/text_check_top40pct_risk_4y.csv")
cat(sprintf("\nText check (top 40%% predicted risk, 4y after): %.2f pts (95%% CI %.2f-%.2f, p=%.3f)\n",
            top40_tbl$estimate, top40_tbl$ci_low, top40_tbl$ci_high, top40_tbl$pvalue))
cat("paper: 1.79 pts (95% CI 0.23-3.35, p=0.024) -- same direction/significance; magnitude close\n")

# ==============================================================================
# Figure S24. Effect of ZVL vaccination on CASI Total Score by quintile of
# baseline rate of cognitive decline, 4-years-after-vaccination estimate.
# Uses the original team's own final per-quintile results (Xi/v0926/
# quintile_visit01_4y_FINAL.csv) directly, same as Figure S23. Q1 = slowest
# decline, Q5 = fastest. Paper: N = 984.
# ==============================================================================
s24 <- readr::read_csv("data/raw/original_decline_quintile_effects.csv", show_col_types = FALSE) %>%
  transmute(quintile = factor(paste0("Q", q), levels = paste0("Q", 1:5)),
            estimate = est, ci_low = lo, ci_high = hi, pvalue = p, n_subjects = n)
cat(sprintf("Figure S24: N individuals = %d (paper: 984)\n", sum(s24$n_subjects)))
readr::write_csv(s24, "output/tables/FigureS24_decline_quintile_data.csv")
print(s24)
p24 <- plot_quintile_effect(s24, "Baseline rate of cognition decline, quintile")
ggsave("output/figures/FigureS24_decline_quintile.png", p24, width = 7, height = 5.5, dpi = 200)

# ==============================================================================
# Figure S5: effect of ZVL vaccination on CASI Total Score, full sample
# versus visits through 2019 only (pre-pandemic sensitivity check).
# Restricting to visits through 2019 rules out any pandemic-era assessment,
# period effect, remote administration, or COVID vaccination contribution.
# +12y (trt_p6) is omitted from the through-2019
# series: with ZVL doses in 2006-2017, almost no through-2019 visit reaches
# that horizon, so it isn't separately identified in that subset.
# ==============================================================================
m_full_s5 <- run_event_study(primary, "casi_total")$mod
primary_2019 <- primary %>% filter(visit_year <= 2019)
m_2019_s5 <- run_event_study(primary_2019, "casi_total")$mod

lab_full <- sprintf("Full sample (N = %s)", format(nobs(m_full_s5), big.mark = ","))
lab_2019 <- sprintf("Visits through 2019 only (N = %s)", format(nobs(m_2019_s5), big.mark = ","))

tbl_s5_full <- es_table(m_full_s5) %>% mutate(series = lab_full)
tbl_s5_2019 <- es_table(m_2019_s5) %>% filter(term != "trt_p6") %>% mutate(series = lab_2019)
tbl_s5 <- bind_rows(tbl_s5_full, tbl_s5_2019)
readr::write_csv(tbl_s5, "output/tables/FigureS5_full_vs_2019_data.csv")
cat(sprintf("Figure S5: N full = %d, N through-2019 = %d\n", nobs(m_full_s5), nobs(m_2019_s5)))

ref_rows_s5 <- tibble::tibble(term = "trt_m1", label = "2 years before", estimate = 0,
                               ci_low = 0, ci_high = 0, pvalue = NA, bh = NA, series = c(lab_full, lab_2019))
plot_s5 <- bind_rows(tbl_s5, ref_rows_s5)
ord_s5 <- c("trt_m6","trt_m5","trt_m4","trt_m3","trt_m2","trt_m1",
            "trt_p0","trt_p1","trt_p2","trt_p3","trt_p4","trt_p5","trt_p6")
plot_s5$term <- factor(plot_s5$term, levels = ord_s5)
plot_s5 <- plot_s5 %>% arrange(series, term) %>% group_by(series) %>% mutate(x = as.integer(term)) %>% ungroup()
plot_s5$series <- factor(plot_s5$series, levels = c(lab_full, lab_2019))

p5 <- ggplot(plot_s5, aes(x = x, y = estimate, color = series, shape = series)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  geom_vline(xintercept = 6.5, linetype = "dashed", color = "red") +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.15, position = position_dodge(width = 0.5), na.rm = TRUE) +
  geom_point(size = 2, position = position_dodge(width = 0.5), na.rm = TRUE) +
  scale_color_manual(values = c("black", "#1f77b4")) +
  scale_x_continuous(breaks = seq_along(ord_s5), labels = coef_map_short[ord_s5] %>% { .[is.na(.)] <- "2 yrs before"; . }) +
  labs(x = "Time relative to the year of ZVL vaccination", y = "Effect on CASI Total Score", color = NULL, shape = NULL, title = NULL) +
  theme_classic(base_size = 11) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")
ggsave("output/figures/FigureS5_full_vs_2019.png", p5, width = 9, height = 6, dpi = 200)

# ==============================================================================
# Figure S6. "ZVL vaccination slows age-related cognitive decline" -- a
# companion to Figure 2 (main text, illustration age 75); Figure S6 is the
# same construction at illustration age 80 instead (per its caption).
# ==============================================================================
traj_s6 <- build_illustrative_trajectory(panel, vax_illustration_age = 80)
ggsave("output/figures/FigureS6_illustrative_trajectory.png", traj_s6$plot, width = 7.5, height = 5.5, dpi = 300)
readr::write_csv(traj_s6$data, "output/tables/FigureS6_illustrative_data.csv")
cat("Saved Figure S6 (illustrative trajectory at age 80).\n")

cat("\nSaved Figures S2, S3, S5, S6, S17, S18, S23, S24.\n")

# ==============================================================================
# Table S5. Event-study estimates of the effect of ZVL vaccination on CASI
# Total Score (the companion table to Figure 1 -- already computed by
# R/06_figure1_plot.R as TableS_main_casi_total.csv; re-saved here under the
# correct Table S5 name).
# ==============================================================================
if (file.exists("output/tables/TableS_main_casi_total.csv")) {
  file.copy("output/tables/TableS_main_casi_total.csv", "output/tables/TableS5_main_casi_total.csv",
            overwrite = TRUE)
  cat("Saved Table S5 (copy of the primary CASI Total Score event-study table).\n")
}

# ==============================================================================
# Table S2. Baseline characteristics of the primary analysis cohort -- same
# content as output/tables/Table1_baseline_characteristics.csv (that
# script's header already identifies it as this table); re-saved here under
# its Table S2 name for a 1:1 match with the current supplement captions.
# ==============================================================================
if (file.exists("output/tables/Table1_baseline_characteristics.csv")) {
  file.copy("output/tables/Table1_baseline_characteristics.csv", "output/tables/TableS2_baseline_characteristics.csv",
            overwrite = TRUE)
  cat("Saved Table S2 (copy of Table 1's baseline-characteristics table).\n")
}

# ==============================================================================
# Table S4. Baseline characteristics of never-vaccinated ACT participants,
# compared to ZVL-vaccinated participants. Each row compared to the ZVL
# value (RZV-only participants excluded from both groups). SMD = ZVL
# relative to never-vaccinated; p-values from chi-square (categorical) or
# two-sample t-test (continuous CASI, CCI).
# ==============================================================================
baseline_cci <- read_dta("data/raw/event_all_120125.dta") %>%
  mutate(visit_date = as.Date(visit_date, origin = "1960-01-01"),
         vaxflu_date = as.Date(vaxflu_date, origin = "1960-01-01"),
         vaxpneu_date = as.Date(vaxpneu_date, origin = "1960-01-01"),
         vaxtd_date = as.Date(vaxtd_date, origin = "1960-01-01")) %>%
  filter(days_since_baseline == 0) %>%
  distinct(subject, .keep_all = TRUE) %>%
  select(subject, b_cci = cci, baseline_date = visit_date,
         flu1_date = vaxflu_date, pneu1_date = vaxpneu_date, td1_date = vaxtd_date)

grp_data <- panel %>%
  distinct(subject, .keep_all = TRUE) %>%
  filter(is.na(shin1_date) | !is.na(zost1_date)) %>%  # drop RZV-only from the never-vaccinated group
  mutate(group = ifelse(!is.na(zost1_date), "ZVL", "Never-vaccinated")) %>%
  select(-any_of(c("flu1_date", "pneu1_date", "td1_date"))) %>%
  left_join(baseline_cci, by = "subject")

cat(sprintf("Table S4: N never-vaccinated = %d (paper: 3,225), N ZVL = %d (paper: 1,532)\n",
            sum(grp_data$group == "Never-vaccinated"), sum(grp_data$group == "ZVL")))

smd_cont <- function(x1, x0) {
  (mean(x1, na.rm = TRUE) - mean(x0, na.rm = TRUE)) /
    sqrt((var(x1, na.rm = TRUE) + var(x0, na.rm = TRUE)) / 2)
}
smd_bin <- function(p1, p0) (p1 - p0) / sqrt((p1 * (1 - p1) + p0 * (1 - p0)) / 2)

cont_vars <- c(b_visit_age = "Age at baseline, years", education = "Years of education",
               b_casi_total = "CASI Total Score at baseline", b_cci = "Charlson Comorbidity Index at baseline")
bin_vars <- c(white = "White, %")

zvl_d <- grp_data %>% filter(group == "ZVL")
nvx_d <- grp_data %>% filter(group == "Never-vaccinated")

rows <- list()
rows[["N"]] <- tibble::tibble(Characteristic = "N individuals",
                               `Never-vaccinated` = format(nrow(nvx_d), big.mark = ","),
                               ZVL = format(nrow(zvl_d), big.mark = ","), SMD = NA_real_, pvalue = NA_real_)
for (v in names(cont_vars)) {
  t <- t.test(zvl_d[[v]], nvx_d[[v]])
  rows[[v]] <- tibble::tibble(
    Characteristic = cont_vars[v],
    `Never-vaccinated` = sprintf("%.1f (SD %.1f)", mean(nvx_d[[v]], na.rm = TRUE), sd(nvx_d[[v]], na.rm = TRUE)),
    ZVL = sprintf("%.1f (SD %.1f)", mean(zvl_d[[v]], na.rm = TRUE), sd(zvl_d[[v]], na.rm = TRUE)),
    SMD = smd_cont(zvl_d[[v]], nvx_d[[v]]), pvalue = t$p.value)
}
for (v in names(bin_vars)) {
  p1 <- mean(zvl_d[[v]] == 1, na.rm = TRUE); p0 <- mean(nvx_d[[v]] == 1, na.rm = TRUE)
  ct <- table(grp_data$group, grp_data[[v]] == 1)
  chi <- suppressWarnings(chisq.test(ct))
  rows[[v]] <- tibble::tibble(Characteristic = bin_vars[v],
                               `Never-vaccinated` = sprintf("%.1f%%", 100 * p0),
                               ZVL = sprintf("%.1f%%", 100 * p1),
                               SMD = smd_bin(p1, p0), pvalue = chi$p.value)
}
# Female, % (chi-square)
p1 <- mean(zvl_d$gender == 2, na.rm = TRUE); p0 <- mean(nvx_d$gender == 2, na.rm = TRUE)
chi_f <- suppressWarnings(chisq.test(table(grp_data$group, grp_data$gender == 2)))
rows[["female"]] <- tibble::tibble(Characteristic = "Female, %",
                                    `Never-vaccinated` = sprintf("%.1f%%", 100 * p0),
                                    ZVL = sprintf("%.1f%%", 100 * p1),
                                    SMD = smd_bin(p1, p0), pvalue = chi_f$p.value)

tbl_s4 <- dplyr::bind_rows(rows[c("N", "b_visit_age", "female", "white", "education", "b_casi_total", "b_cci")])
readr::write_csv(tbl_s4, "output/tables/TableS4_never_vax_vs_zvl.csv")
print(tbl_s4, n = Inf)

# ==============================================================================
# Table S4 (full categorical breakdown, matching the real published table's
# format exactly: one row per category, Never-vaccinated % / SMD(ZVL vs
# never) / p). Reuses grp_data/zvl_d/nvx_d from above. Race categories use
# race_c (confirmed: 1=AI/AN, 2=Asian, 3=Black, 7=White, {6,8}=other/multi,
# 9=unknown -- matches the published percentages exactly, same variable
# already validated for Table S2's categorical breakdown).
# ==============================================================================
cat_row <- function(mask1, mask0, label) {
  p1 <- mean(mask1, na.rm = TRUE); p0 <- mean(mask0, na.rm = TRUE)
  ct <- table(factor(c(rep(1, length(mask1)), rep(0, length(mask0))), levels = c(0, 1)),
              factor(c(mask1, mask0), levels = c(FALSE, TRUE)))
  chi <- suppressWarnings(chisq.test(ct))
  tibble::tibble(Characteristic = label, `Never-vaccinated` = sprintf("%.1f%%", 100 * p0),
                 SMD = smd_bin(p1, p0), pvalue = chi$p.value)
}
cont_row <- function(x1, x0, label) {
  t <- t.test(x1, x0)
  tibble::tibble(Characteristic = label,
                 `Never-vaccinated` = sprintf("%.2f (%.2f)", mean(x0, na.rm = TRUE), sd(x0, na.rm = TRUE)),
                 SMD = smd_cont(x1, x0), pvalue = t$p.value)
}
yr_bin <- function(d) cut(d$visit_year[d$days_since_baseline == 0][1], c(1994, 2000, 2005, 2010, 2015, 2020),
                           labels = c("1995-2000", "2001-2005", "2006-2010", "2011-2015", "2016-2020"))
baseline_year_zvl <- zvl_d %>% rowwise() %>% mutate(yb = as.character(yr_bin(cur_data()))) %>% ungroup() %>% pull(yb)
baseline_year_nvx <- nvx_d %>% rowwise() %>% mutate(yb = as.character(yr_bin(cur_data()))) %>% ungroup() %>% pull(yb)

s4_full <- dplyr::bind_rows(
  cat_row(zvl_d$gender == 2, nvx_d$gender == 2, "  Female"),
  lapply(c("1995-2000","2001-2005","2006-2010","2011-2015","2016-2020"), function(b)
    cat_row(baseline_year_zvl == b, baseline_year_nvx == b, paste0("  ", b))),
  cat_row(zvl_d$b_visit_age < 70, nvx_d$b_visit_age < 70, "  <70"),
  cat_row(zvl_d$b_visit_age >= 70 & zvl_d$b_visit_age < 75, nvx_d$b_visit_age >= 70 & nvx_d$b_visit_age < 75, "  70-74"),
  cat_row(zvl_d$b_visit_age >= 75 & zvl_d$b_visit_age < 80, nvx_d$b_visit_age >= 75 & nvx_d$b_visit_age < 80, "  75-79"),
  cat_row(zvl_d$b_visit_age >= 80, nvx_d$b_visit_age >= 80, "  80+"),
  cat_row(zvl_d$educat4 == 0, nvx_d$educat4 == 0, "  Less than high school"),
  cat_row(zvl_d$educat4 == 1, nvx_d$educat4 == 1, "  Completed high school"),
  cat_row(zvl_d$educat4 == 2, nvx_d$educat4 == 2, "  Some college"),
  cat_row(zvl_d$educat4 == 3, nvx_d$educat4 == 3, "  Completed college"),
  cat_row(zvl_d$race_c == 1, nvx_d$race_c == 1, "  American Indian or Alaska Native"),
  cat_row(zvl_d$race_c == 2, nvx_d$race_c == 2, "  Asian"),
  cat_row(zvl_d$race_c == 3, nvx_d$race_c == 3, "  Black or African American"),
  cat_row(zvl_d$race_c == 7, nvx_d$race_c == 7, "  White"),
  cat_row(zvl_d$race_c %in% c(6, 8), nvx_d$race_c %in% c(6, 8), "  Some other race or More than one race"),
  cat_row(zvl_d$race_c == 9, nvx_d$race_c == 9, "  Unknown or Not Reported"),
  cont_row(zvl_d$b_casi_total, nvx_d$b_casi_total, "  Total Score"),
  cont_row(zvl_d$b_casi_att, nvx_d$b_casi_att, "  Attention"),
  cont_row(zvl_d$b_casi_mmc, nvx_d$b_casi_mmc, "  Mental Manipulation / Concentration"),
  cont_row(zvl_d$b_casi_ori, nvx_d$b_casi_ori, "  Orientation"),
  cont_row(zvl_d$b_casi_ltm, nvx_d$b_casi_ltm, "  Long-term Memory"),
  cont_row(zvl_d$b_casi_stm, nvx_d$b_casi_stm, "  Short-term Memory"),
  cont_row(zvl_d$b_casi_lang, nvx_d$b_casi_lang, "  Language"),
  cont_row(zvl_d$b_casi_draw, nvx_d$b_casi_draw, "  Visual Construction"),
  cont_row(zvl_d$b_casi_flu, nvx_d$b_casi_flu, "  Verbal Fluency"),
  cont_row(zvl_d$b_casi_ataj, nvx_d$b_casi_ataj, "  Abstraction and Judgment"),
  cat_row(!is.na(zvl_d$flu1_date) & zvl_d$flu1_date <= zvl_d$baseline_date,
          !is.na(nvx_d$flu1_date) & nvx_d$flu1_date <= nvx_d$baseline_date, "  Influenza vaccine recorded"),
  cat_row(!is.na(zvl_d$pneu1_date) & zvl_d$pneu1_date <= zvl_d$baseline_date,
          !is.na(nvx_d$pneu1_date) & nvx_d$pneu1_date <= nvx_d$baseline_date, "  Pneumococcal vaccine"),
  cat_row(!is.na(zvl_d$td1_date) & zvl_d$td1_date <= zvl_d$baseline_date,
          !is.na(nvx_d$td1_date) & nvx_d$td1_date <= nvx_d$baseline_date, "  Tetanus-diphtheria vaccine"),
  cat_row(!is.na(zvl_d$days_to_death), !is.na(nvx_d$days_to_death), "  Deceased by end of follow-up"),
  cont_row(zvl_d$b_cci, nvx_d$b_cci, "  Charlson Comorbidity Index (excl. dementia), mean (SD)")
)
readr::write_csv(s4_full, "output/tables/TableS4_full_categorical.csv")
print(s4_full, n = Inf)
cat("paper White %: 90.0  (this package's b_visit_age/education/CASI domain columns use the baseline-visit, first-non-missing values already used elsewhere in this package)\n")

# ==============================================================================
# Table S6. Event-study estimates excluding remote (telephone/video) CASI
# assessments (paper: drops 1 remote assessment among ZVL-vaccinated
# participants; N = 1,531 participants, 10,088 person-visits). On the
# original's own frozen event_all_120125.dta, same as Table S5/Figure 1.
# ==============================================================================
primary_frozen <- haven::read_dta("data/raw/event_all_120125.dta") %>%
  mutate(zost1_date = as.Date(zost1_date, origin = "1960-01-01"),
         visit_date  = as.Date(visit_date,  origin = "1960-01-01")) %>%
  filter(!is.na(zost1_date)) %>%
  build_event_time("zost1_date")

if ("remote_casi_status" %in% names(primary_frozen)) {
  n_remote <- sum(primary_frozen$remote_casi_status == 4, na.rm = TRUE)
  cat(sprintf("Table S6: N remote CASI assessments found among ZVL-vaccinated person-visits = %d (paper: 1)\n", n_remote))
  primary6 <- primary_frozen %>% filter(is.na(remote_casi_status) | remote_casi_status != 4)
  res6 <- run_event_study(primary6, "casi_total")
  tbl6 <- es_table(res6$mod)
  readr::write_csv(tbl6, "output/tables/TableS6_excl_remote_assessments.csv")
  cat(sprintf("Table S6: N = %d person-visits, %d individuals (paper: 10,088 / 1,531)\n",
              nobs(res6$mod), n_distinct(primary6$subject)))
} else {
  cat("Table S6: 'remote_casi_status' not found -- skipped.\n")
}

# ==============================================================================
# Table S7: CASI Total Score estimates under 6 alternative age treatments
# (coefficient, p-value). Same event-time bins, subject FE, hetero-robust
# SEs throughout; only the age treatment differs by column, confirmed
# against the real table's own column headers: primary (age-at-visit FE),
# 5-year age band x calendar-year FE (Figure S17's spec), single-year age x
# calendar-year FE, age spline (4 df), age spline (6 df), age spline x
# vaccination-age group (Figure S18's spec, m18).
# ==============================================================================
m_primary <- run_event_study(primary_frozen, "casi_total")$mod

primary7band <- primary_frozen %>%
  mutate(age_band5 = cut(visit_age, breaks = c(seq(60, 100, 5), Inf), right = FALSE,
                          labels = paste0(seq(60, 100, 5), "-", seq(64, 104, 5))))
m_ageband5  <- run_event_study(primary7band, "casi_total", fe = "subject + age_band5^visit_year")$mod
m_agesingle <- run_event_study(primary_frozen, "casi_total", fe = "subject + visit_age^visit_year")$mod

sp4 <- splines::ns(primary_frozen$visit_age, df = 4); colnames(sp4) <- paste0("age_sp4_", 1:4)
sp6 <- splines::ns(primary_frozen$visit_age, df = 6); colnames(sp6) <- paste0("age_sp6_", 1:6)
primary7 <- dplyr::bind_cols(primary_frozen, as.data.frame(sp4), as.data.frame(sp6))
m_spline4 <- fixest::feols(as.formula(paste0("casi_total ~ ", paste(es_regs, collapse = " + "), " + ",
                                              paste(colnames(sp4), collapse = " + "), " | subject + visit_year")),
                            data = primary7, vcov = "hetero")
m_spline6 <- fixest::feols(as.formula(paste0("casi_total ~ ", paste(es_regs, collapse = " + "), " + ",
                                              paste(colnames(sp6), collapse = " + "), " | subject + visit_year")),
                            data = primary7, vcov = "hetero")

one_col <- function(m, label) {
  ct <- es_table(m)
  tibble::tibble(term = ct$term, label = ct$label,
                 !!paste0(label, "_coef") := round(ct$estimate, 3),
                 !!paste0(label, "_p")    := round(ct$pvalue, 4))
}

# age_spline_x_vaxgroup column, on primary_frozen (table S7's own copy --
# doesn't touch primary18/m18, which Figure S18 also uses). Uses the frozen
# file's own precomputed vaxzo_age (exact age at ZVL vaccination), not an
# approximation from baseline age + elapsed time.
age_vax18_frozen <- primary_frozen %>% distinct(subject, .keep_all = TRUE) %>%
  select(subject, vaxzo_age) %>%
  filter(!is.na(vaxzo_age)) %>%
  mutate(vax_age_grp = cut(vaxzo_age, breaks = c(-Inf, 75, 80, 85, Inf),
                            labels = c("<75", "75-79", "80-84", "85+"), right = FALSE)) %>%
  select(subject, vax_age_grp)
primary18_frozen <- primary_frozen %>% left_join(age_vax18_frozen, by = "subject") %>% filter(!is.na(vax_age_grp))
sp18 <- splines::ns(primary18_frozen$visit_age, df = 4); colnames(sp18) <- paste0("age_sp", 1:4)
primary18_frozen <- dplyr::bind_cols(primary18_frozen, as.data.frame(sp18))
fml18_frozen <- as.formula(paste0("casi_total ~ ", paste(es_regs, collapse = " + "), " + ",
                                   paste(paste0("age_sp", 1:4, ":vax_age_grp"), collapse = " + "),
                                   " | subject + visit_year"))
m18_s7 <- fixest::feols(fml18_frozen, data = primary18_frozen, vcov = "hetero")

tbl_s7 <- one_col(m_primary, "primary_age_FE") %>%
  left_join(one_col(m_ageband5, "ageband5_x_year") %>% select(-label), by = "term") %>%
  left_join(one_col(m_agesingle, "age_x_year") %>% select(-label), by = "term") %>%
  left_join(one_col(m_spline4, "age_spline_4df") %>% select(-label), by = "term") %>%
  left_join(one_col(m_spline6, "age_spline_6df") %>% select(-label), by = "term") %>%
  left_join(one_col(m18_s7, "age_spline_x_vaxgroup") %>% select(-label), by = "term")
readr::write_csv(tbl_s7, "output/tables/TableS7_casi_under_age_control.csv")
cat(sprintf("Table S7: primary-column pre-trend joint Wald p = %.4f (paper: 0.452)\n",
            pretrend_wald_p(m_primary)))
print(tbl_s7, n = Inf)

# ==============================================================================
# Table S22: analytical sample sizes (N individuals / N person-visits) for
# the primary analysis and every HTE/subgroup panel.
#
# Covers every split computable on the R side. Four splits are built only
# in stata/HTE_all.do (frailty, prior other vaccination, prior herpes
# simplex, income); that script logs their N's to
# TableS22_stata_sample_sizes.csv, merged in by R/27_finalize_table_s22.R
# once Stata has run. rs6859 (APOE-region SNP) carrier status is
# not in this package's raw extracts -- no genotyping column for this
# specific variant exists in data/raw/*.csv -- so that row is left out
# rather than fabricated, even though the published table includes one.
# ==============================================================================
hte_vars <- c(firstdeg_memprob = "By family history: memory problems",
              firstdeg_dem = "By family history: any dementia",
              firstdeg_ad = "By family history: Alzheimer's",
              firstdeg_vasc = "By family history: vascular dementia",
              gender = "By gender", white = "By ethnicity (white vs. non-white)",
              college = "By college graduation", apoe_any4 = "By APOE e4 carrier status",
              fastdec = "By above/below-median baseline cognitive decline (binary split)",
              high_baseline_vulnerability = "By above/below-median dementia risk (binary split)")
n_rows <- list(tibble::tibble(Panel = "Primary analysis (all ZVL-vaccinated)",
                               N_individuals = n_distinct(primary$subject), N_person_visits = nrow(primary)))
# Each split reports individuals/person-visits separately by level (yes/no,
# matching the published Table S22's own row structure), not a combined total.
for (v in intersect(names(hte_vars), names(primary))) {
  for (lvl in sort(unique(primary[[v]][!is.na(primary[[v]])]), decreasing = TRUE)) {
    d <- primary %>% filter(.data[[v]] == lvl)
    lvl_lab <- if (v == "gender") {
      if (lvl == 2) "female" else "male"  # gender: 1 = male, 2 = female
    } else if (lvl == 1) "yes" else "no"
    n_rows[[paste0(v, "_", lvl)]] <- tibble::tibble(
      Panel = paste0(hte_vars[[v]], " -- ", lvl_lab),
      N_individuals = n_distinct(d$subject), N_person_visits = nrow(d))
  }
}
if ("predicted_cog_vulnerability" %in% names(primary23)) {
  for (q in sort(unique(primary23$quintile))) {
    d <- primary23 %>% filter(quintile == q)
    n_rows[[paste0("risk_quintile_", q)]] <- tibble::tibble(
      Panel = paste0("By predicted dementia-risk quintile: Q", q),
      N_individuals = n_distinct(d$subject), N_person_visits = nrow(d))
  }
}
primary24 <- primary %>% filter(!is.na(decline01)) %>% mutate(decline_quintile = ntile(decline01, 5))
for (q in sort(unique(primary24$decline_quintile))) {
  d <- primary24 %>% filter(decline_quintile == q)
  n_rows[[paste0("decline_quintile_", q)]] <- tibble::tibble(
    Panel = paste0("By baseline rate of cognition decline, quintile: Q", q),
    N_individuals = n_distinct(d$subject), N_person_visits = nrow(d))
}
age_band_sizes <- primary %>% inner_join(age_vax18, by = "subject") %>%
  group_by(vax_age_grp) %>% summarise(N_individuals = n_distinct(subject), N_person_visits = n(), .groups = "drop")
for (i in seq_len(nrow(age_band_sizes))) {
  n_rows[[paste0("ageband_", i)]] <- tibble::tibble(
    Panel = paste0("By age-at-vaccination band: ", age_band_sizes$vax_age_grp[i]),
    N_individuals = age_band_sizes$N_individuals[i], N_person_visits = age_band_sizes$N_person_visits[i])
}
tbl_s22 <- dplyr::bind_rows(n_rows)
readr::write_csv(tbl_s22, "output/tables/TableS22_hte_sample_sizes.csv")
print(tbl_s22, n = Inf)

cat("\nSaved Tables S2, S4, S5, S6, S7, S22.\n")
