# Export the ROC figures (Figure 1, S2, S3, S5) as PNGs ------------------------
#
#   Rscript export_roc_figures.R
#
# Redraws the curves from the held-out predictions saved by
# 02_ensemble_models.Rmd (audit/ensemble_predictions.csv), so the figures come
# from the same run as the AUCs in audit/ensemble_aucs.csv without refitting.
# Each panel is labelled with its AUC (95% CI) from that file.
# Figure S4 (Bangladesh, day-of-admission symptoms) is not reproducible here.

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(cowplot)
})
source(here("R", "ensemble.R"))

preds <- read.csv(here("audit", "ensemble_predictions.csv"))
aucs <- read.csv(here("audit", "ensemble_aucs.csv"))
out_dir <- here("figures")
dir.create(out_dir, showWarnings = FALSE)

# ROC points and AUC labels for one analysis
panel <- function(pop, test, window, colour, title = NULL) {
  p <- filter(preds, population == pop, reference_test == test, time_window == window)
  a <- filter(aucs, population == pop, reference_test == test, time_window == window)
  stopifnot(nrow(a) == 2, length(unique(p$predictor_set)) == 2)
  roc_df <- p %>%
    group_by(predictor_set) %>%
    group_modify(~ roc_points(.x$infect_bin, .x$pred)) %>%
    ungroup() %>%
    as.data.frame()
  lab <- with(a, sprintf("%s: AUC %.2f (%.2f, %.2f)", predictor_set, auc, lower, upper))
  plot_roc(roc_df, colour, title) +
    annotate("text", x = 1, y = c(0.12, 0.05), label = lab, hjust = 1, size = 3) +
    coord_equal()
}

save_png <- function(plot, name, width) {
  file <- file.path(out_dir, paste0(name, ".png"))
  ggsave(file, plot, width = width, height = 5, dpi = 300, bg = "white")
  message("Wrote ", file)
}

save_png(plot_grid(panel("South Sudan", "IgM", "<=30 days", "red", "A"),
                   panel("Bangladesh", "IgM", "<=30 days", "orange", "B")),
         "figure_1_roc_30days", 10)
save_png(plot_grid(panel("South Sudan", "IgM", "<=14 days", "red", "A"),
                   panel("Bangladesh", "IgM", "<=14 days", "orange", "B")),
         "figure_S2_roc_14days", 10)
save_png(panel("South Sudan", "PCR", "<=30 days", "red"), "figure_S3_roc_pcr", 5.5)
save_png(panel("South Sudan children", "IgM", "<=30 days", "red"), "figure_S5_roc_children", 5.5)
