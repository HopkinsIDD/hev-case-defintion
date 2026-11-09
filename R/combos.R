# Exhaustive search over alternative case definitions ------------------------
#
# A case definition is a set of symptoms plus a rule "acute jaundice and at
# least `min_k` of these symptoms". We evaluate every non-empty subset of the
# available symptoms (2^p - 1 definitions: 32,767 for SSD, 131,071 for BD).
#
# Instead of building one data-frame column per definition (slow, and needed
# batch files on disk), we use matrix algebra: if X is the patient x symptom
# 0/1 matrix and D is the symptom x definition 0/1 matrix, then X %*% D counts,
# for every patient and definition, how many of the definition's symptoms the
# patient has.

#' All non-empty symptom subsets, in the same order as the original scripts
#' (all singles, then all pairs, ...), so definition ids are comparable.
enumerate_case_defs <- function(symptoms) {
  unlist(lapply(seq_along(symptoms), function(k) combn(symptoms, k, simplify = FALSE)),
         recursive = FALSE)
}

# Symptom x definition indicator matrix
def_matrix <- function(defs, symptoms) {
  D <- matrix(0L, length(symptoms), length(defs), dimnames = list(symptoms, NULL))
  for (j in seq_along(defs)) D[defs[[j]], j] <- 1L
  D
}

#' Evaluate every case definition.
#'
#' @param dat data with `infect_bin` and 0/1 symptom columns.
#' @param symptoms symptoms to combine.
#' @param min_k number of listed symptoms required (1 = "at least one").
#' @param weights optional row weights (e.g. post-stratification weights).
#' @param chunk number of definitions processed at once (limits memory use).
#' @return one row per definition with 2x2 counts (weighted if `weights` given)
#'   and sensitivity, specificity, PPV, NPV.
evaluate_case_defs <- function(dat, symptoms, min_k = 1, weights = NULL, chunk = 10000) {
  defs <- enumerate_case_defs(symptoms)
  X <- as.matrix(dat[symptoms])
  y <- dat$infect_bin
  w <- if (is.null(weights)) rep(1, nrow(dat)) else weights

  counts <- lapply(split(seq_along(defs), ceiling(seq_along(defs) / chunk)), function(ix) {
    meets <- (X %*% def_matrix(defs[ix], symptoms)) >= min_k   # patients x definitions
    cbind(tp = colSums(meets * (w * y)), fp = colSums(meets * (w * (1 - y))))
  })
  counts <- do.call(rbind, counts)

  pos <- sum(w * y); neg <- sum(w * (1 - y))
  data.frame(
    id = seq_along(defs),
    case_def = vapply(defs, paste, "", collapse = " + "),
    n_symptoms = lengths(defs),
    tp = counts[, "tp"], fn = pos - counts[, "tp"],
    fp = counts[, "fp"], tn = neg - counts[, "fp"]
  ) %>%
    mutate(
      sensitivity = tp / (tp + fn),
      specificity = tn / (tn + fp),
      ppv = tp / (tp + fp),
      npv = tn / (tn + fn),
      informative = sensitivity + specificity > 1
    )
}

#' Pick the alternative case definition(s) to report.
#'
#' Among informative definitions with sensitivity AND specificity above
#' `threshold`, take the highest specificity; break ties by highest
#' sensitivity, then fewest symptoms.
#'
#' Without the threshold, the "best" definition is a single rare symptom with
#' near-100% specificity and near-0% sensitivity, so the threshold matters and
#' must be reported in the methods.
select_case_defs <- function(results, threshold, n = 1) {
  results %>%
    filter(informative, sensitivity > threshold, specificity > threshold) %>%
    arrange(desc(specificity), desc(sensitivity), n_symptoms) %>%
    head(n)
}

#' Add Wilson CIs (unweighted analyses) to selected definitions.
add_wilson_ci <- function(selected) {
  ci <- function(x, n) t(mapply(prop_ci, x, n))
  for (m in list(c("sensitivity", "tp", "fn"), c("specificity", "tn", "fp"),
                 c("ppv", "tp", "fp"), c("npv", "tn", "fn"))) {
    x <- selected[[m[2]]]; cis <- ci(x, x + selected[[m[3]]])
    selected[[paste0(m[1], "_lower")]] <- cis[, "lower"]
    selected[[paste0(m[1], "_upper")]] <- cis[, "upper"]
  }
  selected
}

#' Percentile bootstrap CIs for weighted (joint) analyses.
#'
#' Resamples patients *within each study population* so every resample keeps
#' the original population sizes, matching the fixed 50/50 weights. The same
#' resamples are used for every definition.
bootstrap_case_defs <- function(dat, symptoms, min_k, ids, n_boot = 1000, strata = dat$dataset) {
  defs <- enumerate_case_defs(symptoms)[ids]
  meets <- (as.matrix(dat[symptoms]) %*% def_matrix(defs, symptoms)) >= min_k
  y <- dat$infect_bin
  by_stratum <- split(seq_len(nrow(dat)), strata)

  # n x n_boot matrix: how many times each patient appears in each resample
  times <- vapply(seq_len(n_boot), function(b) {
    idx <- unlist(lapply(by_stratum, function(i) i[sample.int(length(i), replace = TRUE)]))
    tabulate(idx, nbins = nrow(dat))
  }, numeric(nrow(dat)))
  W <- dat$weight * times

  tp <- crossprod(meets, W * y)                       # definitions x resamples
  fn <- matrix(colSums(W * y), nrow(tp), n_boot, byrow = TRUE) - tp
  fp <- crossprod(meets, W * (1 - y))
  tn <- matrix(colSums(W * (1 - y)), nrow(fp), n_boot, byrow = TRUE) - fp

  q <- function(m) t(apply(m, 1, quantile, c(0.025, 0.975), na.rm = TRUE))
  ci <- list(sensitivity = q(tp / (tp + fn)), specificity = q(tn / (tn + fp)),
             ppv = q(tp / (tp + fp)), npv = q(tn / (tn + fn)))
  out <- data.frame(id = ids)
  for (m in names(ci)) {
    out[[paste0(m, "_lower")]] <- ci[[m]][, 1]
    out[[paste0(m, "_upper")]] <- ci[[m]][, 2]
  }
  out
}

#' Format selected definitions like manuscript Table 5.
format_selected <- function(selected) {
  f <- function(m) sprintf("%.0f%% (%.0f, %.0f)", 100 * selected[[m]],
                           100 * selected[[paste0(m, "_lower")]], 100 * selected[[paste0(m, "_upper")]])
  data.frame(
    case_def = selected$case_def,
    tp = round(selected$tp, 1), fn = round(selected$fn, 1),
    fp = round(selected$fp, 1), tn = round(selected$tn, 1),
    sensitivity = f("sensitivity"), specificity = f("specificity"),
    ppv = f("ppv"), npv = f("npv")
  )
}
