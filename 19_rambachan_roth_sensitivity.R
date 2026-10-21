# ==============================================================================
# 19_rambachan_roth_sensitivity.R
#
# Figure S19. Sensitivity of the ZVL effect on the CASI Total Score to
# violations of parallel pre-trends, using the relative-magnitudes approach
# of Rambachan and Roth (2023), via the HonestDiD package.
#
# For each post-period of interest (4 years after, 8 years after), this
# reports how large a violation of parallel pre-trends (relative to the
# largest observed pre-period deviation, Mbar) would be required before the
# confidence set for the treatment effect includes zero. Mbar = 0
# corresponds to the standard parallel-trends assumption; Mbar = 1 allows
# a post-period violation up to the size of the largest pre-period
# deviation; larger Mbar values allow proportionally larger violations.
# ==============================================================================
suppressPackageStartupMessages({ library(dplyr); library(ggplot2); library(HonestDiD); library(patchwork) })
source("R/01_event_study_utils.R")

panel <- readRDS("data/analysis_panel.rds")
dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)

primary <- panel %>%
  filter(!is.na(zost1_date), !is.na(casi_total)) %>%
  build_event_time("zost1_date")

res <- run_event_study(primary, "casi_total")
b <- coef(res$mod)[es_regs]
S <- vcov(res$mod)[es_regs, es_regs]

# es_regs = 5 pre-periods (trt_m6..trt_m2) + 7 post-periods (trt_p0..trt_p6),
# matching numPrePeriods/numPostPeriods below. Post-period index within
# es_regs[6:12]: trt_p0=1, trt_p1=2, trt_p2=3 (4 years after), trt_p3=4,
# trt_p4=5 (8 years after), trt_p5=6, trt_p6=7.
Mbarvec <- c(0, 0.5, 1, 1.5, 2)
targets <- list(
  list(label = "4y", idx = 3L, name = "4 years after"),
  list(label = "8y", idx = 5L, name = "8 years after")
)

plot_sensitivity <- function(res_df, point_est, y_lab) {
  d <- data.frame(x = as.numeric(res_df$Mbar), lb = as.numeric(res_df$lb), ub = as.numeric(res_df$ub))
  d <- d[order(d$x), ]
  xb <- unique(d$x)
  xlabs <- ifelse(xb == 0, "Parallel\ntrends", formatC(xb, format = "f", digits = 2))
  ggplot(d, aes(x = x)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.8, colour = "grey60") +
    geom_hline(yintercept = point_est, linetype = "dotted", linewidth = 0.5, colour = "grey45") +
    geom_linerange(aes(ymin = lb, ymax = ub), linewidth = 0.6, linetype = "dashed", colour = "grey60", na.rm = TRUE) +
    geom_point(aes(y = lb), shape = 16, size = 2.4, colour = "black", na.rm = TRUE) +
    geom_point(aes(y = ub), shape = 16, size = 2.4, colour = "black", na.rm = TRUE) +
    scale_x_continuous(breaks = xb, labels = xlabs, expand = expansion(mult = c(0.15, 0.1))) +
    labs(x = expression(paste("Relative magnitude of the pre-period violation, ", bar(M))), y = y_lab) +
    theme_classic(base_size = 13) +
    theme(panel.grid = element_blank())
}

panels <- list()
for (tg in targets) {
  l_vec <- HonestDiD::basisVector(index = tg$idx, size = 7L)
  rm_res <- HonestDiD::createSensitivityResults_relativeMagnitudes(
    betahat = b, sigma = S, numPrePeriods = 5L, numPostPeriods = 7L,
    l_vec = l_vec, Mbarvec = Mbarvec
  )
  pt_est <- as.numeric(b[5 + tg$idx])
  readr::write_csv(rm_res, sprintf("output/tables/TableS19_honestdid_%s.csv", tg$label))
  cat(sprintf("\n=== Rambachan-Roth sensitivity, %s (point est = %.3f) ===\n", tg$name, pt_est))
  print(rm_res)

  # breakdown Mbar: smallest value at which the CI starts to include zero
  excl <- (rm_res$lb > 0) | (rm_res$ub < 0)
  bd <- if (any(excl & rm_res$Mbar > 0)) max(rm_res$Mbar[excl & rm_res$Mbar > 0]) else NA_real_
  cat(sprintf("Breakdown Mbar (largest Mbar for which CI still excludes zero): %s\n",
              ifelse(is.na(bd), "none within range tested", sprintf("%.2f", bd))))

  panels[[tg$label]] <- plot_sensitivity(rm_res, pt_est, y_lab = paste0("Effect on CASI Total Score, ", tg$name))
}

# Figure S19 is a single two-panel figure: 4 years after (left), 8 years
# after (right) -- matches the published caption exactly.
p19 <- panels[["4y"]] + panels[["8y"]]
ggsave("output/figures/FigureS19_honestdid.png", p19, width = 13, height = 5.5, dpi = 200)

cat("\nSaved Figure S19 (Rambachan-Roth sensitivity) tables and figure.\n")
