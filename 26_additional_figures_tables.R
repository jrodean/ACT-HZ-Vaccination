# ==============================================================================
# 26_additional_figures_tables.R
#
# Figures S2, S3, S17, S18, S23, S24, and Tables S4, S5, S6, S7, S22 --
# identified as missing from the package by directly re-extracting the
# current Figure S#/Table S# captions and methodology notes from the
# Supplement's own XML (word/document.xml), since the supplement's
# figure/table numbering changed again after this package's last audit.
# Each section below cites the exact caption/notes text it implements.
# ==============================================================================
suppressPackageStartupMessages({ library(dplyr); library(ggplot2); library(splines) })
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
tbl_s1 <- tibble::tibble(
  Characteristic = c("N ZVL-vaccinated participants", "N with >=1 HZ episode before first ZVL dose",
                      "N with >=1 HZ episode after first ZVL dose",
                      "N with an HZ episode both before and after"),
  Value = c(nrow(hz_all),
            sum(hz_all$n_episodes_pre > 0),
            sum(hz_all$n_episodes_post > 0),
            sum(hz_all$n_episodes_pre > 0 & hz_all$n_episodes_post > 0))
)
readr::write_csv(tbl_s1, "output/tables/TableS1_hz_episodes_relative_to_vax.csv")
print(tbl_s1, n = Inf)

med_s2 <- median(hz_post$years_to_episode, na.rm = TRUE)
p2 <- ggplot(hz_post, aes(x = years_to_episode)) +
  geom_histogram(binwidth = 1, boundary = 0, fill = "grey60", color = "white") +
  geom_vline(xintercept = med_s2, linetype = "dashed", color = "red") +
  labs(x = "Years from first ZVL dose to first post-vaccination HZ episode",
       y = "N participants", title = NULL) +
  theme_minimal()
ggsave("output/figures/FigureS2_time_to_hz_episode.png", p2, width = 7.5, height = 5, dpi = 200)

# ==============================================================================
# Figure S3. "Pre-existing rate of cognitive decline by herpes zoster
# vaccination status." Per-participant OLS slope of CASI Total Score on
# calendar time, fit on pre-vaccination visits (or all visits for never-
# vaccinated), restricted to visits before 2006 (when ZVL was licensed) and
# subjects with >=3 pre-vaccination visits in that window. Red diamonds
# mark group means; boxes show the IQR. Paper: N = 3,225 never-vaccinated,
# 1,532 ZVL-vaccinated.
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

cat(sprintf("Figure S3: N never-vaccinated = %d (paper: 3,225), N ZVL-vaccinated = %d (paper: 1,532)\n",
            sum(slopes_s3$group == "Never-vaccinated"), sum(slopes_s3$group == "ZVL-vaccinated")))
readr::write_csv(slopes_s3, "output/tables/FigureS3_pre_existing_decline_data.csv")

means_s3 <- slopes_s3 %>% group_by(group) %>% summarise(mean_slope = mean(slope), .groups = "drop")
p3 <- ggplot(slopes_s3, aes(x = group, y = slope)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, fill = "grey90") +
  geom_jitter(width = 0.1, alpha = 0.05, size = 0.5) +
  geom_point(data = means_s3, aes(y = mean_slope), shape = 18, size = 5, color = "red") +
  labs(x = NULL, y = "Pre-vaccination annualized rate of change in CASI Total Score", title = NULL) +
  theme_minimal()
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
# ==============================================================================
baseline_date18 <- panel %>% filter(days_since_baseline == 0) %>%
  distinct(subject, .keep_all = TRUE) %>% select(subject, baseline_date = visit_date, b_visit_age)
age_vax18 <- panel %>% distinct(subject, .keep_all = TRUE) %>% select(subject, zost1_date) %>%
  left_join(baseline_date18, by = "subject") %>%
  mutate(age_at_zost1 = b_visit_age + as.numeric(zost1_date - baseline_date) / 365.25) %>%
  filter(!is.na(age_at_zost1)) %>%
  mutate(vax_age_grp = cut(age_at_zost1, breaks = c(-Inf, 75, 80, 85, Inf),
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
  ggplot(d, aes(x = quintile, y = estimate)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
    geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0, color = "grey50") +
    geom_point(size = 3) +
    annotate("segment", x = 0.55, xend = 0.55, y = 0, yend = max(d$ci_high) * 0.3,
             arrow = arrow(length = unit(0.15, "cm")), color = "grey40") +
    annotate("segment", x = 0.55, xend = 0.55, y = 0, yend = min(d$ci_low) * 0.3,
             arrow = arrow(length = unit(0.15, "cm")), color = "grey40") +
    annotate("text", x = 0.55, y = max(d$ci_high) * 0.45, label = "Benefit", color = "grey40", size = 3.2) +
    annotate("text", x = 0.55, y = min(d$ci_low) * 0.45, label = "Harm", color = "grey40", size = 3.2) +
    coord_cartesian(xlim = c(0.4, 5.6), clip = "off") +
    labs(x = xlab, y = "Effect on CASI Total Score (4 yrs after ZVL)", title = NULL) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank())
}

# ==============================================================================
# Figure S23. Effect of ZVL vaccination on CASI Total Score by baseline
# predicted dementia-risk quintile (LASSO score from 00c), focusing on the
# 4-years-after-vaccination estimate (trt_p2).
# ==============================================================================
if (!"predicted_cog_vulnerability" %in% names(primary)) {
  risk_scores <- readr::read_csv("data/risk_scores.csv", show_col_types = FALSE)
  primary <- primary %>% left_join(risk_scores, by = "subject")
}
primary23 <- primary %>%
  filter(!is.na(predicted_cog_vulnerability)) %>%
  mutate(quintile = ntile(predicted_cog_vulnerability, 5))

s23_rows <- lapply(1:5, function(q) {
  d <- primary23 %>% filter(quintile == q)
  res <- run_event_study(d, "casi_total")
  tbl <- es_table(res$mod)
  r <- tbl[tbl$term == "trt_p2", ]
  tibble::tibble(quintile = q, estimate = r$estimate, ci_low = r$ci_low, ci_high = r$ci_high,
                 pvalue = r$pvalue, n_subjects = n_distinct(d$subject))
})
s23 <- dplyr::bind_rows(s23_rows) %>% mutate(quintile = factor(quintile, labels = paste0("Q", 1:5)))
readr::write_csv(s23, "output/tables/FigureS23_risk_quintile_data.csv")
print(s23)
p23 <- plot_quintile_effect(s23, "Predicted-risk quintile")
ggsave("output/figures/FigureS23_risk_quintile.png", p23, width = 7, height = 5.5, dpi = 200)

# ==============================================================================
# Figure S24. Effect of ZVL vaccination on CASI Total Score by quintile of
# baseline rate of cognitive decline -- computed as the slope of CASI Total
# Score on calendar time using each participant's FIRST TWO pre-vaccination
# ACT visits (Q1 = slowest decline, Q5 = fastest). Paper: N = 984.
# ==============================================================================
decline_slope <- panel %>%
  filter(!is.na(zost1_date), visit_date < zost1_date, !is.na(casi_total)) %>%
  arrange(subject, visit_date) %>%
  group_by(subject) %>%
  filter(n() >= 2) %>%
  slice_head(n = 2) %>%
  filter(n() == 2) %>%
  summarise(
    slope = unname(diff(casi_total)) / (unname(diff(as.numeric(visit_date))) / 365.25),
    .groups = "drop"
  ) %>%
  filter(is.finite(slope))

cat(sprintf("Figure S24: N with a baseline decline slope = %d (paper: 984)\n", nrow(decline_slope)))

primary24 <- primary %>% inner_join(decline_slope, by = "subject") %>%
  mutate(quintile = ntile(-slope, 5))  # Q1 = slowest decline (largest/least-negative slope), Q5 = fastest

s24_rows <- lapply(1:5, function(q) {
  d <- primary24 %>% filter(quintile == q)
  res <- run_event_study(d, "casi_total")
  tbl <- es_table(res$mod)
  r <- tbl[tbl$term == "trt_p2", ]
  tibble::tibble(quintile = q, estimate = r$estimate, ci_low = r$ci_low, ci_high = r$ci_high,
                 pvalue = r$pvalue, n_subjects = n_distinct(d$subject))
})
s24 <- dplyr::bind_rows(s24_rows) %>% mutate(quintile = factor(quintile, labels = paste0("Q", 1:5)))
readr::write_csv(s24, "output/tables/FigureS24_decline_quintile_data.csv")
print(s24)
p24 <- plot_quintile_effect(s24, "Baseline rate of cognition decline, quintile")
ggsave("output/figures/FigureS24_decline_quintile.png", p24, width = 7, height = 5.5, dpi = 200)

# ==============================================================================
# Figure S5. "Effect of ZVL vaccination on the CASI Total Score, full
# sample versus visits through 2019 only." Original source: Xi/v0926/
# reviewer3.R, PART D ("pre-pandemic sensitivity: censor all analyses at
# end-2019") -- re-estimates the primary ZVL event study restricted to
# visits through 2019, so no pandemic-era assessment, period effect, remote
# administration, or COVID vaccination can contribute. +12y (trt_p6) is
# omitted from the through-2019 series: with ZVL doses in 2006-2017, almost
# no through-2019 visit reaches that horizon, so it is not separately
# identified in that subset.
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
# Figure S6. "ZVL vaccination slows age-related cognitive decline" -- the
# identical caption to Figure 2 (main text); the Supplement shows the same
# illustrative trajectory figure again as its own numbered item. Reuses
# Figure 2's output directly (R/05_descriptive_figures.R) rather than
# recomputing it.
# ==============================================================================
if (file.exists("output/figures/Figure2_illustrative_trajectory.png")) {
  file.copy("output/figures/Figure2_illustrative_trajectory.png", "output/figures/FigureS6_illustrative_trajectory.png",
            overwrite = TRUE)
  cat("Saved Figure S6 (copy of Figure 2's illustrative trajectory).\n")
}

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
baseline_cci <- panel %>% filter(days_since_baseline == 0) %>%
  distinct(subject, .keep_all = TRUE) %>% select(subject, b_cci = cci)

grp_data <- panel %>%
  distinct(subject, .keep_all = TRUE) %>%
  filter(is.na(shin1_date) | !is.na(zost1_date)) %>%  # drop RZV-only from the never-vaccinated group
  mutate(group = ifelse(!is.na(zost1_date), "ZVL", "Never-vaccinated")) %>%
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
# Table S6. Event-study estimates excluding remote (telephone/video) CASI
# assessments (paper: drops 1 remote assessment among ZVL-vaccinated
# participants; N = 1,531 participants, 10,088 person-visits).
# ==============================================================================
if ("remote_casi_status" %in% names(panel)) {
  # Matches the original definition in Xi/v0926/reviewer3.R ([C3] remote vs
  # in-person CASI): remote = !is.na(remote_casi_status), i.e. any non-
  # missing status code (not a specific value) marks a remote assessment.
  n_remote <- sum(!is.na(primary$remote_casi_status), na.rm = TRUE)
  cat(sprintf("Table S6: N remote CASI assessments found among ZVL-vaccinated person-visits = %d (paper: 1)\n", n_remote))
  primary6 <- primary %>% filter(is.na(remote_casi_status))
  res6 <- run_event_study(primary6, "casi_total")
  tbl6 <- es_table(res6$mod)
  readr::write_csv(tbl6, "output/tables/TableS6_excl_remote_assessments.csv")
  cat(sprintf("Table S6: N = %d person-visits, %d individuals (paper: 10,088 / 1,531)\n",
              nobs(res6$mod), n_distinct(primary6$subject)))
} else {
  cat("Table S6: 'remote_casi_status' not found in the analysis panel -- skipped.\n")
}

# ==============================================================================
# Table S7. CASI Total Score estimates under control of age (coefficient;
# p-value). Same event-time bins, individual FE, heteroskedasticity-robust
# SEs throughout; only the treatment of age differs across columns. The
# exact column set is not individually labelled in the extracted
# supplement text beyond "as labelled in each column header" -- this
# reconstructs the natural progression of age controls used elsewhere in
# this revision (no age control -> linear age -> primary age FE -> age
# spline -> age spline x vaccination-age group, i.e. Figure S18's spec).
# ==============================================================================
m_noage   <- fixest::feols(as.formula(paste0("casi_total ~ ", paste(es_regs, collapse = " + "), " | subject + visit_year")),
                            data = primary, vcov = "hetero")
m_linage  <- fixest::feols(as.formula(paste0("casi_total ~ ", paste(es_regs, collapse = " + "), " + visit_age | subject + visit_year")),
                            data = primary, vcov = "hetero")
m_primary <- run_event_study(primary, "casi_total")$mod
sp7 <- splines::ns(primary$visit_age, df = 4); colnames(sp7) <- paste0("age_sp", 1:4)
primary7 <- dplyr::bind_cols(primary, as.data.frame(sp7))
m_spline  <- fixest::feols(as.formula(paste0("casi_total ~ ", paste(es_regs, collapse = " + "), " + ",
                                              paste(colnames(sp7), collapse = " + "), " | subject + visit_year")),
                            data = primary7, vcov = "hetero")

one_col <- function(m, label) {
  ct <- es_table(m)
  tibble::tibble(term = ct$term, label = ct$label,
                 !!paste0(label, "_coef") := round(ct$estimate, 3),
                 !!paste0(label, "_p")    := round(ct$pvalue, 4))
}
tbl_s7 <- one_col(m_noage, "no_age_control") %>%
  left_join(one_col(m_linage, "linear_age") %>% select(-label), by = "term") %>%
  left_join(one_col(m_primary, "age_FE_primary") %>% select(-label), by = "term") %>%
  left_join(one_col(m_spline, "age_spline_4df") %>% select(-label), by = "term") %>%
  left_join(one_col(m18, "age_spline_x_vaxgroup") %>% select(-label), by = "term")
readr::write_csv(tbl_s7, "output/tables/TableS7_casi_under_age_control.csv")
cat(sprintf("Table S7: primary-column pre-trend joint Wald p = %.4f (paper: 0.452)\n",
            pretrend_wald_p(m_primary)))
print(tbl_s7, n = Inf)

# ==============================================================================
# Table S22. Analytical sample sizes for heterogeneity analyses -- N
# individuals / N person-visits contributing to the primary analysis and to
# every HTE/subgroup panel in the manuscript and Supplement.
# ==============================================================================
# Variables available on the R side of the panel (Stata-only HTE splits
# built inside stata/HTE_all.do from raw columns not carried over to the R
# panel under the same name -- e.g. high_income, any_vax_pre, prior_hsv,
# frail_pre, shingles_followup -- are not duplicated here; their sample
# sizes are reported directly by that script's own log/console output).
hte_vars <- c(firstdeg_memprob = "By family history: memory problems",
              firstdeg_dem = "By family history: any dementia",
              gender = "By gender", white = "By ethnicity (white vs. non-white)",
              college = "By college graduation", apoe_any4 = "By APOE e4 carrier status",
              fastdec = "By above/below-median baseline cognitive decline (binary split)",
              high_baseline_vulnerability = "By above/below-median dementia risk (binary split)")
n_rows <- list(tibble::tibble(Panel = "Primary analysis (all ZVL-vaccinated)",
                               N_individuals = n_distinct(primary$subject), N_person_visits = nrow(primary)))
for (v in intersect(names(hte_vars), names(primary))) {
  d <- primary %>% filter(!is.na(.data[[v]]))
  n_rows[[v]] <- tibble::tibble(Panel = hte_vars[[v]],
                                 N_individuals = n_distinct(d$subject), N_person_visits = nrow(d))
}
if ("predicted_cog_vulnerability" %in% names(primary23)) {
  n_rows[["risk_quintile"]] <- tibble::tibble(Panel = "By predicted dementia-risk quintile",
                                               N_individuals = n_distinct(primary23$subject),
                                               N_person_visits = nrow(primary23))
}
n_rows[["decline_quintile"]] <- tibble::tibble(Panel = "By baseline rate of cognition decline, quintile",
                                                N_individuals = n_distinct(primary24$subject),
                                                N_person_visits = nrow(primary24))
tbl_s22 <- dplyr::bind_rows(n_rows)
readr::write_csv(tbl_s22, "output/tables/TableS22_hte_sample_sizes.csv")
print(tbl_s22, n = Inf)

cat("\nSaved Tables S2, S4, S5, S6, S7, S22.\n")
