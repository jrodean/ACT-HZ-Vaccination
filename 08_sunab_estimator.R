# ==============================================================================
# 08_sunab_estimator.R
#
# Figure S20. Sun & Abraham (2021) interaction-weighted event-study
# estimator, as a robustness check against heterogeneous treatment-effect
# dynamics under staggered ZVL adoption timing. The published figure was
# implemented via Stata's `eventstudyinteract`; this is an R equivalent via
# fixest::sunab(), using the design stated in the supplement caption:
#   - cohort = calendar year of first ZVL dose
#   - the last-treated cohort (vaccinated in 2017) serves as the control
#     cohort, so observations from 2017 onward are excluded
#   - event-time window: -12 to +8 years (asymmetric, per the caption)
# ==============================================================================
suppressPackageStartupMessages({ library(dplyr); library(fixest); library(ggplot2) })
source("R/01_event_study_utils.R")

panel <- readRDS("data/analysis_panel.rds")
dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)

primary <- panel %>%
  filter(!is.na(zost1_date), !is.na(casi_total)) %>%
  mutate(
    vax_year = lubridate::year(zost1_date),
    event_time = floor(as.numeric(visit_date - zost1_date) / (2 * 365)) * 2,
    rel = pmin(pmax(event_time, -12), 8)
  ) %>%
  filter(vax_year <= 2017)

m_sa <- fixest::feols(casi_total ~ sunab(vax_year, rel, ref.p = -2) | subject + visit_year + visit_age,
                       data = primary, vcov = "hetero")

ct <- broom::tidy(m_sa, conf.int = TRUE)
ct$rel <- suppressWarnings(as.numeric(sub(".*rel::(-?[0-9]+).*", "\\1", ct$term)))
ct <- ct[!is.na(ct$rel), ]
sat <- ct %>% transmute(rel, estimate, ci_low = conf.low, ci_high = conf.high, pvalue = p.value) %>%
  arrange(rel)

readr::write_csv(sat, "output/tables/TableS20_sunab_estimator.csv")
cat("=== Sun-Abraham (interaction-weighted) event study ===\n")
print(sat, n = Inf)

agg <- tryCatch(fixest::aggregate(m_sa, agg = "att"), error = function(e) NULL)
if (!is.null(agg)) cat(sprintf("\nOverall ATT (post-period average): %.3f (SE %.3f)\n", agg[1,1], agg[1,2]))

p <- ggplot(sat, aes(x = rel, y = estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high), alpha = 0.15) +
  geom_line() + geom_point(size = 2) +
  labs(x = "Years relative to ZVL vaccination", y = "Effect on CASI Total Score (Sun-Abraham IW)",
       title = NULL) +
  theme_minimal()
ggsave("output/figures/FigureS20_sunab_estimator.png", p, width = 7.5, height = 5, dpi = 300)

cat("\nSaved Sun-Abraham table and figure.\n")
