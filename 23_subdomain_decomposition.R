# ==============================================================================
# 23_subdomain_decomposition.R
#
# Figure S21. Decomposition of CASI Total Score effects by subdomain.
#
# The CASI Total Score is the exact sum of its 9 domain subscores (verified
# directly against the analysis panel: correlation = 1, row sums match
# exactly). This script re-estimates the primary event-study specification
# separately for each domain subscore and shows, at each post-vaccination
# horizon, how much of the total effect each domain contributes, as a
# stacked bar chart.
#
# Depends on R/03_domains_and_dementia_outcomes.R having been run first
# (reads its Table S8-S16 output).
# ==============================================================================
suppressPackageStartupMessages({ library(dplyr); library(ggplot2); library(readr) })

dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)

domain_files <- c(
  "Short-Term Memory" = "output/tables/TableS8_domain_casi_stm.csv",
  "Orientation" = "output/tables/TableS9_domain_casi_ori.csv",
  "Long-Term Memory" = "output/tables/TableS10_domain_casi_ltm.csv",
  "Mental Manipulation" = "output/tables/TableS11_domain_casi_mmc.csv",
  "Language" = "output/tables/TableS12_domain_casi_lang.csv",
  "Verbal Fluency" = "output/tables/TableS13_domain_casi_flu.csv",
  "Attention" = "output/tables/TableS14_domain_casi_att.csv",
  "Visual Construction" = "output/tables/TableS15_domain_casi_draw.csv",
  "Abstraction and Judgment" = "output/tables/TableS16_domain_casi_ataj.csv"
)

missing <- domain_files[!file.exists(domain_files)]
if (length(missing)) {
  stop("Missing domain table(s) -- run R/03_domains_and_dementia_outcomes.R first: ",
       paste(names(missing), collapse = ", "))
}

domain_tbl <- purrr::imap_dfr(domain_files, function(path, label) {
  read_csv(path, show_col_types = FALSE) %>% mutate(domain = label)
})

post_terms <- c(trt_p1 = "2 years after", trt_p2 = "4 years after", trt_p3 = "6 years after",
                 trt_p4 = "8 years after", trt_p5 = "10 years after", trt_p6 = "12 years after")

plot_tbl <- domain_tbl %>%
  filter(term %in% names(post_terms)) %>%
  mutate(horizon = factor(post_terms[term], levels = unname(post_terms)))

readr::write_csv(plot_tbl, "output/tables/TableS21_subdomain_decomposition.csv")

# Sanity check: at each horizon, do the 9 domain estimates sum to
# approximately the primary CASI Total Score estimate?
total_check <- plot_tbl %>% group_by(horizon) %>% summarise(sum_domains = sum(estimate), .groups = "drop")
cat("=== Sum of 9 domain estimates by horizon (sanity check vs. Figure 1's total) ===\n")
print(total_check)

p21 <- ggplot(plot_tbl, aes(x = horizon, y = estimate, fill = domain)) +
  geom_col(position = "stack") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  labs(x = "Time relative to ZVL vaccination", y = "Effect on CASI Total Score", fill = "CASI subdomain",
       title = NULL) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave("output/figures/FigureS21_subdomain_decomposition.png", p21, width = 9, height = 6, dpi = 200)

cat("\nSaved Figure S21 (subdomain decomposition).\n")
