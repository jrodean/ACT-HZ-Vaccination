# ==============================================================================
# 17_figS38_age_at_vaccination.R
#
# Figure S38, "Effect of ZVL vaccination on the CASI Total Score by age at
# vaccination." Continuous interaction model per the supplement's
# "Interaction effect between age at vaccination and vaccine effect"
# section:
#   Y_it = subject FE + year FE + age FE
#          + sum_tau delta_tau 1{event_time=tau}
#          + sum_tau beta_tau  1{event_time=tau} * c_i
#          + sum_tau gamma_tau 1{event_time=tau} * c_i^2 + e_it
#   where c_i = (age at ZVL vaccination - 80) / 5 (5-year units, centered at 80)
# Age at vaccination is time-invariant so its main effect/square are absorbed
# by subject FE; only interactions with event time are identified (matches
# the supplement text exactly). Effect tau years after vaccination for a
# person vaccinated at age a: delta_tau + beta_tau*c(a) + gamma_tau*c(a)^2,
# with the 95% CI via the delta method on the 3-coefficient linear combination.
# Plotted for tau = 2, 4, 6, 8 years after, over age 66-92 (5th-95th pctile
# of age at ZVL vaccination in this cohort), matching the published figure.
# ==============================================================================
suppressPackageStartupMessages({ library(dplyr); library(ggplot2) })
source("R/01_event_study_utils.R")

panel <- readRDS("data/analysis_panel.rds")
baseline_date <- panel %>% filter(days_since_baseline == 0) %>%
  distinct(subject, .keep_all = TRUE) %>% select(subject, baseline_date = visit_date, b_visit_age)

age_vax <- panel %>% distinct(subject, .keep_all = TRUE) %>% select(subject, zost1_date) %>%
  left_join(baseline_date, by = "subject") %>%
  mutate(age_at_zost1 = b_visit_age + as.numeric(zost1_date - baseline_date) / 365.25) %>%
  filter(!is.na(age_at_zost1)) %>% select(subject, age_at_zost1)

primary <- panel %>%
  filter(!is.na(zost1_date), !is.na(casi_total)) %>%
  build_event_time("zost1_date") %>%
  left_join(age_vax, by = "subject") %>%
  mutate(c_age = (age_at_zost1 - 80) / 5, c_age2 = c_age^2)

inter1 <- paste0(es_regs, ":c_age")
inter2 <- paste0(es_regs, ":c_age2")
fml <- as.formula(paste0("casi_total ~ ", paste(es_regs, collapse = " + "), " + ",
                          paste(inter1, collapse = " + "), " + ", paste(inter2, collapse = " + "),
                          " | subject + visit_year + visit_age"))
m <- fixest::feols(fml, data = primary, vcov = "hetero")
cat(sprintf("N rows=%d, N subjects=%d\n", nobs(m), n_distinct(primary$subject)))

V <- vcov(m); B <- coef(m)
age_grid <- seq(quantile(age_vax$age_at_zost1, .05), quantile(age_vax$age_at_zost1, .95), length.out = 60)
taus <- c(trt_p1 = 2, trt_p2 = 4, trt_p3 = 6, trt_p4 = 8)

pred_rows <- list()
for (term in names(taus)) {
  d_name <- term; b_name <- paste0(term, ":c_age"); g_name <- paste0(term, ":c_age2")
  if (!all(c(d_name, b_name, g_name) %in% names(B))) { cat("missing terms for", term, "\n"); next }
  idx <- c(d_name, b_name, g_name)
  Vb <- V[idx, idx]
  for (a in age_grid) {
    c_a <- (a - 80) / 5
    cvec <- c(1, c_a, c_a^2)
    est <- sum(cvec * B[idx])
    se <- sqrt(as.numeric(t(cvec) %*% Vb %*% cvec))
    pred_rows[[length(pred_rows) + 1]] <- tibble::tibble(
      years_after = taus[term], age = a, estimate = est,
      ci_low = est - 1.96 * se, ci_high = est + 1.96 * se
    )
  }
}
pred <- bind_rows(pred_rows) %>%
  mutate(panel_label = factor(paste(years_after, "yrs after ZVL"),
                               levels = paste(c(2, 4, 6, 8), "yrs after ZVL")))

dir.create("output/tables", showWarnings = FALSE)
dir.create("output/figures", showWarnings = FALSE)
readr::write_csv(pred, "output/tables/FigureS38_age_at_vaccination_data.csv")

p <- ggplot(pred, aes(x = age, y = estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  geom_line(color = "black", linewidth = 0.7) +
  geom_line(aes(y = ci_low), linetype = "dashed", color = "grey40") +
  geom_line(aes(y = ci_high), linetype = "dashed", color = "grey40") +
  facet_wrap(~panel_label, nrow = 1) +
  labs(x = "Age at ZVL vaccination", y = "Effect on CASI Total Score") +
  theme_bw(base_size = 18) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(color = "black", size = 16),
    panel.border = element_rect(color = "black", fill = NA),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14, color = "black")
  )
ggsave("output/figures/FigureS38_age_at_vaccination.png", p, width = 14, height = 4.5, dpi = 600)
cat("saved output/figures/FigureS38_age_at_vaccination.png\n")
