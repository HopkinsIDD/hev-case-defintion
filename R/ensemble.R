# Ensemble (SuperLearner) models and ROC curves ------------------------------
#
# The ensemble is an "upper bound" on how well signs/symptoms can separate
# HEV-positive from HEV-negative acute jaundice. Discrimination is estimated
# with leave-one-out cross-validation (LOOCV): refit on n-1 patients, predict
# the held-out patient, repeat for every patient, then compute the AUC of the
# held-out predictions.

library(SuperLearner)
library(pROC)
library(ROCR)
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

#' AUC with a 2,000-sample percentile bootstrap CI, plus ROC curve points.
auc_summary <- function(y, pred, model) {
  roc_obj <- roc(y, pred, quiet = TRUE)
  ci <- ci.auc(roc_obj, method = "bootstrap", boot.n = 2000, progress = "none")
  perf <- performance(prediction(pred, y), measure = "tpr", x.measure = "fpr")
  list(
    auc = data.frame(model = model, auc = as.numeric(auc(roc_obj)), lower = ci[1], upper = ci[3],
                     row.names = NULL),
    roc = data.frame(model = model, fpr = perf@x.values[[1]], tpr = perf@y.values[[1]])
  )
}

#' Fit symptom-only and symptom+demographic models on one dataset.
run_ensembles <- function(dat, demographics, label = "", cores = 1) {
  models <- list(
    "Signs/symptoms" = dat[setdiff(names(dat), demographics)],
    "Signs/symptoms plus demographics" = dat
  )
  res <- lapply(names(models), function(m) {
    message("Fitting: ", label, " ", m)
    d <- models[[m]]
    auc_summary(d$infect_bin, loocv_predictions(d, cores), m)
  })
  list(
    auc = cbind(population = label, do.call(rbind, lapply(res, `[[`, "auc"))),
    roc = cbind(population = label, do.call(rbind, lapply(res, `[[`, "roc")))
  )
}

#' ROC plot in the style of Figure 1.
plot_roc <- function(roc_df, colour = "red", title = NULL) {
  ggplot(roc_df, aes(fpr, tpr, alpha = model)) +
    geom_line(colour = colour, linewidth = 1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    scale_alpha_manual(values = c("Signs/symptoms" = 0.5, "Signs/symptoms plus demographics" = 1)) +
    labs(x = "1 - Specificity", y = "Sensitivity", alpha = NULL, title = title) +
    theme_minimal() +
    theme(legend.position = "bottom")
}
