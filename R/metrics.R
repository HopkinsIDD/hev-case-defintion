# Diagnostic accuracy of binary case definitions ---------------------------

#' Proportion with a 95% Wilson score interval (prop.test without continuity
#' correction). Returns NAs when the denominator is 0 instead of erroring.
prop_ci <- function(x, n) {
  if (n == 0) return(c(est = NA_real_, lower = NA_real_, upper = NA_real_))
  ci <- prop.test(x, n, conf.level = 0.95, correct = FALSE)$conf.int
  c(est = x / n, lower = ci[1], upper = ci[2])
}

#' 2x2 counts plus sensitivity, specificity, PPV and NPV (with Wilson CIs) for
#' one case definition.
#'
#' @param case_def 0/1 (or logical) vector: does the patient meet the definition?
#' @param infect_bin 0/1 vector: laboratory-confirmed hepatitis E.
case_def_performance <- function(case_def, infect_bin) {
  case_def <- as.logical(case_def)
  infected <- infect_bin == 1
  tp <- sum(case_def & infected, na.rm = TRUE)
  fn <- sum(!case_def & infected, na.rm = TRUE)
  fp <- sum(case_def & !infected, na.rm = TRUE)
  tn <- sum(!case_def & !infected, na.rm = TRUE)

  ci <- rbind(
    sensitivity = prop_ci(tp, tp + fn),
    specificity = prop_ci(tn, tn + fp),
    ppv         = prop_ci(tp, tp + fp),
    npv         = prop_ci(tn, tn + fn)
  )
  data.frame(
    tp = tp, fn = fn, fp = fp, tn = tn,
    metric = rownames(ci), est = ci[, "est"], lower = ci[, "lower"], upper = ci[, "upper"],
    row.names = NULL
  )
}

#' One-row summary formatted like the manuscript tables: "71% (61, 80)".
format_performance <- function(perf) {
  fmt <- sprintf("%.0f%% (%.0f, %.0f)", 100 * perf$est, 100 * perf$lower, 100 * perf$upper)
  out <- data.frame(tp = perf$tp[1], fn = perf$fn[1], fp = perf$fp[1], tn = perf$tn[1])
  out[perf$metric] <- as.list(fmt)
  out
}

# Existing international case definitions (WHO hepatitis E outbreak toolbox).
# Every patient has acute jaundice, so each reduces to "any of these symptoms".
# Malaise is recorded as fatigue (SSD) / drowsiness (BD). Epigastric pain and
# joint pain were not collected in Bangladesh; right upper quadrant tenderness
# was not collected in either population.
INTERNATIONAL_DEFS <- list(
  ssd = list(
    MSF = c("fatigue", "epigastric_pain", "nausea", "loss_appetite"),
    WHO = c("fatigue", "fever"),
    IMC = c("fatigue", "loss_appetite", "fever", "abdomen_pain", "joint_pains")
  ),
  bd = list(
    MSF = c("drowsy_ever", "nausea_ever", "appetite_ever"),
    WHO = c("drowsy_ever", "fever_ever"),
    IMC = c("drowsy_ever", "appetite_ever", "fever_ever", "abdomen_ever")
  )
)

#' Sensitivity/specificity table for the international case definitions.
#' `dat` must have 0/1 symptom columns (load with `binary = TRUE`).
international_table <- function(dat, defs) {
  rows <- lapply(names(defs), function(nm) {
    meets <- rowSums(dat[defs[[nm]]]) >= 1
    cbind(case_definition = nm, format_performance(case_def_performance(meets, dat$infect_bin)))
  })
  do.call(rbind, rows)
}
