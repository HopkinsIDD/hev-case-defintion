# Ensemble (SuperLearner) models and ROC curves ------------------------------
#
# The ensemble is a flexible benchmark for how well signs/symptoms can separate
# HEV-positive from HEV-negative acute jaundice (within this learner library). Discrimination is estimated
# with leave-one-out cross-validation (LOOCV): refit on n-1 patients, predict
# the held-out patient, repeat for every patient, then compute the AUC of the
# held-out predictions.

library(SuperLearner)
library(pROC)
library(ggplot2)

SL_LIBRARY <- list(
  "SL.ranger",   # random forest
  "SL.gam",      # generalized additive model
  "SL.xgboost",  # extreme gradient boosting
  "SL.glmnet"    # lasso
)

#' LOOCV predicted probabilities from the SuperLearner ensemble.
#'
#' @param dat data with `infect_bin`; every other column is a predictor.
#' @param cores > 1 runs held-out fits in parallel (macOS/Linux). Results are
#'   then reproducible for a given `cores` value but differ slightly from the
#'   sequential run, because random numbers are drawn in a different order.
loocv_predictions <- function(dat, cores = 1) {
  X <- dat[setdiff(names(dat), "infect_bin")]
  y <- dat$infect_bin
  fit_one <- function(i) {
    fit <- SuperLearner(Y = y[-i], X = X[-i, , drop = FALSE],
                        family = binomial(), SL.library = SL_LIBRARY)
    as.numeric(predict(fit, newdata = X[i, , drop = FALSE])$pred)
  }
  if (cores > 1) {
    unlist(parallel::mclapply(seq_len(nrow(X)), fit_one, mc.cores = cores, mc.set.seed = TRUE))
  } else {
    vapply(seq_len(nrow(X)), function(i) {
      if (i %% 50 == 0) message("LOOCV fit ", i, " of ", nrow(X))
      fit_one(i)
    }, numeric(1))
  }
}

#' AUC with a 2,000-sample stratified percentile bootstrap CI, plus ROC curve
#' points, all from one pROC object.
#'
#' The direction is fixed: higher predicted probability = HEV positive
#' (controls = 0, cases = 1). pROC's default would otherwise pick whichever
#' direction gives AUC >= 0.5, hiding a model that ranks patients backwards.
auc_summary <- function(y, pred, predictor_set, boot.n = 2000) {
  roc_obj <- roc(y, pred, levels = c(0, 1), direction = "<", quiet = TRUE)
  ci <- ci.auc(roc_obj, method = "bootstrap", boot.n = boot.n, boot.stratified = TRUE,
               progress = "none")
  list(
    auc = data.frame(predictor_set = predictor_set, auc = as.numeric(auc(roc_obj)),
                     lower = ci[1], upper = ci[3], row.names = NULL),
    roc = data.frame(predictor_set = predictor_set,
                     fpr = 1 - roc_obj$specificities, tpr = roc_obj$sensitivities)
  )
}

#' Fit symptom-only and symptom+demographic models on one dataset.
#'
#' @param population,reference_test,time_window identify the analysis; together
#'   with `predictor_set` they form a unique key in the saved outputs.
#' @return list of `auc` (one row per predictor set), `roc` (curve points) and
#'   `pred` (held-out LOOCV predictions, one row per patient and predictor set).
run_ensembles <- function(dat, demographics, population, reference_test, time_window, cores = 1) {
  models <- list(
    "Signs/symptoms" = dat[setdiff(names(dat), demographics)],
    "Signs/symptoms plus demographics" = dat
  )
  id <- data.frame(population = population, reference_test = reference_test,
                   time_window = time_window)
  res <- lapply(names(models), function(m) {
    message("Fitting: ", population, " | ", reference_test, " | ", time_window, " | ", m)
    d <- models[[m]]
    pred <- loocv_predictions(d, cores)
    out <- auc_summary(d$infect_bin, pred, m)
    out$pred <- data.frame(predictor_set = m, row = seq_along(pred), infect_bin = d$infect_bin,
                           pred = pred)
    out
  })
  lapply(c(auc = "auc", roc = "roc", pred = "pred"), function(part)
    cbind(id, do.call(rbind, lapply(res, `[[`, part)), row.names = NULL))
}

#' ROC plot in the style of Figure 1.
plot_roc <- function(roc_df, colour = "red", title = NULL) {
  roc_df <- roc_df[order(roc_df$predictor_set, roc_df$fpr, roc_df$tpr), ]
  ggplot(roc_df, aes(fpr, tpr, alpha = predictor_set)) +
    geom_path(colour = colour, linewidth = 1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    scale_alpha_manual(values = c("Signs/symptoms" = 0.5, "Signs/symptoms plus demographics" = 1)) +
    labs(x = "1 - Specificity", y = "Sensitivity", alpha = NULL, title = title) +
    theme_minimal() +
    theme(legend.position = "bottom")
}
