# Manuscript audit -------------------------------------------------------------
#
# Recomputes every number in the manuscript and supplement that can be derived
# from the analytic datasets, and compares it with the published value in
# audit/published_values.csv (one row per number or claim).
#
#   Rscript audit/audit.R
#
# Takes about two minutes. The ensemble AUCs are read from
# audit/ensemble_aucs.csv, which 02_ensemble_models.Rmd writes (the models take
# about an hour to refit).
#
# Each published value gets a status:
#   PASS         recomputed value matches after rounding to the published precision
#                (within `tolerance` for bootstrap CIs and AUCs)
#   FAIL         recomputed value differs, or the claim is false
#   NOT IN DATA  the number needs data that is not in this repository
#   NOT RUN      the AUC cache is missing
# Rounding: a value counts as matching if it matches under either round-half-up
# or round-half-to-even, so 22.5% may be printed as 22% or 23%.
#
# A second set of checks uses only the published numbers: estimates inside
# their CIs, n (%) consistent with the denominator, 2x2 margins equal to the
# group sizes, and point estimates consistent with the published counts.
#
# Results: audit/audit_results.csv and audit/consistency_results.csv.
#
# To update after the manuscript changes, edit `published` in
# audit/published_values.csv. Every id there must be computed below or be
# marked as not in the data.

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
})
source(here("R", "data.R"))
source(here("R", "metrics.R"))
source(here("R", "combos.R"))

published <- read.csv(here("audit", "published_values.csv"), colClasses = "character",
                      na.strings = character())
published$tolerance <- as.numeric(published$tolerance)
pub <- setNames(published$published, published$id)

computed <- list()        # id -> numeric vector (percentages on 0-100 scale), logical or string
not_run <- character()    # ids that need a cache that does not exist
put <- function(id, value) computed[[id]] <<- value

ssd <- load_ssd(binary = TRUE)
ssd_children <- load_ssd(children = TRUE, binary = TRUE)
bd <- load_bd(binary = TRUE)
joint <- load_joint()

# Table 3, Table S3: patient characteristics ----------------------------------

SSD_ROWS <- c(fever = "fever", loss_of_appetite = "loss_appetite", nausea = "nausea",
              diarrhoea = "diarrhea", vomiting = "vomiting", headache = "head",
              abdominal_pain = "abdomen_pain", altered_mental_state = "mental",
              fatigue_drowsy = "fatigue", convulsion = "convulsion",
              epigastric_pain = "epigastric_pain", cough = "cough", joint_pains = "joint_pains",
              itch = "itch", bleeding = "bleeding")
BD_ROWS <- c(fever = "fever_ever", loss_of_appetite = "appetite_ever", nausea = "nausea_ever",
             diarrhoea = "diarrhea_ever", vomiting = "vomit_ever", headache = "head_ever",
             abdominal_pain = "abdomen_ever", altered_mental_state = "mental_ever",
             fatigue_drowsy = "drowsy_ever", convulsion = "convulsion_ever",
             melaena = "melaena_ever", constipation = "constipate_ever",
             unconsciousness = "unconsc_ever", distended_abdomen = "exam_dist_abdomen",
             oedema = "exam_oedema", dehydration = "exam_dehydrate", liver_palpable = "exam_liver_palp")

delay_band <- function(days) cut(days, c(-Inf, 7, 14, Inf),
                                 labels = c("delay_lt1w", "delay_1_2w", "delay_2w_1m"))
n_pct <- function(x) c(sum(x), 100 * mean(x))
med_iqr <- function(x) unname(quantile(x, c(0.5, 0.25, 0.75)))

describe <- function(dat, prefix, female, days, rows, age_breaks, delays = FALSE) {
  groups <- list(all = dat, neg = dat[dat$infect_bin == 0, ], pos = dat[dat$infect_bin == 1, ])
  for (g in names(groups)) {
    d <- groups[[g]]
    id <- function(key) paste(prefix, key, g, sep = ".")
    put(id("N"), nrow(d))
    put(id("female"), n_pct(female(d)))
    put(id("age_median"), med_iqr(d$age))
    band <- cut(d$age, c(age_breaks, Inf), labels = names(age_breaks), right = FALSE)
    for (b in names(age_breaks)) put(id(b), n_pct(band == b))
    put(id("days_median"), med_iqr(d[[days]]))
    if (delays) {   # delay bands are only tabulated in Table S3
      delay <- delay_band(d[[days]])
      for (b in levels(delay)) put(id(b), n_pct(delay == b))
    }
    for (key in names(rows)) put(id(key), n_pct(d[[rows[[key]]]] == 1))
  }
}
ADULT_AGES <- c(age_14_19 = 14, age_20_39 = 20, age_40_59 = 40, age_60plus = 60)
describe(ssd, "T3.ssd", function(d) d$id_sex == "Female", "days_since_ajs_baseline", SSD_ROWS, ADULT_AGES)
describe(bd, "T3.bd", function(d) d$sex == "2", "days_since_ajs", BD_ROWS, ADULT_AGES)
describe(ssd_children, "S3.children", function(d) d$id_sex == "Female", "days_since_ajs_baseline",
         SSD_ROWS, c(age_0_4 = 0, age_5_10 = 5, age_11_13 = 11), delays = TRUE)

# Table S2: symptoms by delay band ---------------------------------------------

S2_ROWS <- c(fever = "fever_ever", loss_of_appetite = "appetite_ever", nausea = "nausea_ever",
             diarrhea = "diarrhea_ever", vomiting = "vomit_ever", headache = "head_ever",
             abdominal_pain = "abdomen_ever", altered_mental_state = "mental_ever",
             fatigue_drowsiness = "drowsy_ever", convulsion = "convulsion_ever")
ssd_name <- c(SSD_TO_JOINT, convulsion_ever = "convulsion")    # BD name -> SSD name
for (key in names(S2_ROWS)) {
  for (site in c("ssd", "bd")) {
    d <- if (site == "ssd") ssd else bd
    var <- if (site == "ssd") ssd_name[[S2_ROWS[[key]]]] else S2_ROWS[[key]]
    days <- if (site == "ssd") d$days_since_ajs_baseline else d$days_since_ajs
    pct <- 100 * tapply(d[[var]], delay_band(days), mean)
    for (b in names(pct)) put(paste("S2", site, key, b, sep = "."), pct[[b]])
  }
}

# Table 4, Table S4: international case definitions ----------------------------

intl <- list()   # site -> data frame of sens/spec per definition (for claims in the text)
international <- function(dat, defs, prefix, aj = TRUE) {
  pos <- sum(dat$infect_bin == 1); neg <- sum(dat$infect_bin == 0)
  if (aj) {
    counts <- c(tp = pos, fp = neg, fn = 0, tn = 0)
    for (f in names(counts)) put(paste(prefix, "AJ", f, sep = "."), counts[[f]])
    put(paste(prefix, "AJ.sens", sep = "."), 100)
    put(paste(prefix, "AJ.spec", sep = "."), 0)
  }
  do.call(rbind, lapply(names(defs), function(nm) {
    perf <- case_def_performance(rowSums(dat[defs[[nm]]]) >= 1, dat$infect_bin)
    for (f in c("tp", "fp", "fn", "tn")) put(paste(prefix, nm, f, sep = "."), perf[[f]][1])
    ci <- function(m) 100 * unlist(perf[perf$metric == m, c("est", "lower", "upper")])
    put(paste(prefix, nm, "sens", sep = "."), unname(ci("sensitivity")))
    put(paste(prefix, nm, "spec", sep = "."), unname(ci("specificity")))
    data.frame(def = nm, sens = ci("sensitivity")[1] / 100, spec = ci("specificity")[1] / 100)
  }))
}
intl$ssd <- international(ssd, INTERNATIONAL_DEFS$ssd, "T4.ssd")
intl$bd <- international(bd, INTERNATIONAL_DEFS$bd, "T4.bd")
intl$children <- international(ssd_children, INTERNATIONAL_DEFS$ssd, "S4.children", aj = FALSE)

# Table 5, Table S5, Table S6: selected alternative definitions ----------------

sorted_def <- function(x) vapply(strsplit(x, " + ", fixed = TRUE), function(s) paste(sort(s), collapse = " + "), "")
put_selected <- function(sel, prefix, ci) {
  put(paste0(prefix, ".def"), sorted_def(sel$case_def))
  for (f in c("tp", "fp", "fn", "tn")) put(paste(prefix, f, sep = "."), sel[[f]])
  metrics <- c(sens = "sensitivity", spec = "specificity", ppv = "ppv", npv = "npv")
  for (m in names(metrics)) {
    full <- metrics[[m]]
    put(paste(prefix, m, sep = "."),
        100 * c(sel[[full]], ci[[paste0(full, "_lower")]], ci[[paste0(full, "_upper")]]))
  }
}

POPULATIONS <- list(
  `South Sudan` = list(dat = ssd, symptoms = SSD_SYMPTOMS, prefix = "T5.ssd"),
  `South Sudan children` = list(dat = ssd_children, symptoms = SSD_SYMPTOMS, prefix = "S5.children"),
  Bangladesh = list(dat = bd, symptoms = BD_SYMPTOMS, prefix = "T5.bd")
)
selected <- list()
for (i in seq_len(nrow(SELECTION_THRESHOLDS))) {
  th <- SELECTION_THRESHOLDS[i, ]
  if (th$population == "Combined") next
  p <- POPULATIONS[[th$population]]
  sel <- select_case_defs(evaluate_case_defs(p$dat, p$symptoms, th$min_k), th$threshold) %>% add_wilson_ci()
  prefix <- paste0(p$prefix, ".k", th$min_k)
  put_selected(sel, prefix, sel)
  selected[[prefix]] <- sel
}

joint_results <- lapply(c(k1 = 1, k2 = 2), function(k)
  evaluate_case_defs(joint, JOINT_SYMPTOMS, min_k = k, weights = joint$weight))
set.seed(123)
for (k in 1:2) {
  th <- SELECTION_THRESHOLDS %>% filter(population == "Combined", min_k == k)
  sel <- select_case_defs(joint_results[[k]], th$threshold)
  ci <- bootstrap_case_defs(joint, JOINT_SYMPTOMS, k, sel$id)
  put_selected(sel, paste0("S6.joint.k", k), ci)
}

# Ensemble AUCs (cached by 02_ensemble_models.Rmd) -----------------------------

auc_ids <- c(
  text.abstract.auc_ssd = "South Sudan, <=30 days|Signs/symptoms",
  text.abstract.auc_bd = "Bangladesh, <=30 days|Signs/symptoms",
  text.results.auc_ssd_symptoms = "South Sudan, <=30 days|Signs/symptoms",
  text.results.auc_bd_symptoms = "Bangladesh, <=30 days|Signs/symptoms",
  text.results.auc_ssd_demog = "South Sudan, <=30 days|Signs/symptoms plus demographics",
  text.results.auc_bd_demog = "Bangladesh, <=30 days|Signs/symptoms plus demographics",
  text.results.auc_children_symptoms = "South Sudan children|Signs/symptoms",
  text.results.auc_children_demog = "South Sudan children|Signs/symptoms plus demographics"
)
auc_file <- here("audit", "ensemble_aucs.csv")
if (file.exists(auc_file)) {
  aucs <- read.csv(auc_file)
  key <- paste(aucs$population, aucs$model, sep = "|")
  auc_of <- function(k) unlist(aucs[key == k, c("auc", "lower", "upper")])
  for (id in names(auc_ids)) put(id, unname(auc_of(auc_ids[[id]])))
  diff_14 <- sapply(c("South Sudan", "Bangladesh"), function(s)
    sapply(c("Signs/symptoms", "Signs/symptoms plus demographics"), function(m)
      auc_of(paste0(s, ", <=14 days|", m))[1] - auc_of(paste0(s, ", <=30 days|", m))[1]))
  put("text.results.auc_14day_similar", all(abs(diff_14) <= 0.02))
} else {
  not_run <- c(names(auc_ids), "text.results.auc_14day_similar")
}

# Numbers and claims in the text -----------------------------------------------

prev <- function(d) c(100 * mean(d$infect_bin), sum(d$infect_bin), nrow(d))
put("text.abstract.prev_ssd", prev(ssd)); put("text.results.prev_ssd", prev(ssd))
put("text.abstract.prev_bd", prev(bd));   put("text.results.prev_bd", prev(bd))

adult_intl <- rbind(intl$ssd, intl$bd)
for (where in c("abstract", "discussion")) {
  put(paste0("text.", where, ".intl_sens_range"), 100 * range(adult_intl$sens))
  put(paste0("text.", where, ".intl_spec_range"), 100 * range(adult_intl$spec))
}

# Abstract: does any definition built from the shared symptoms reach >60% in both populations?
both_above_60 <- any(sapply(1:2, function(k) {
  by_site <- lapply(split(joint, joint$dataset), evaluate_case_defs, symptoms = JOINT_SYMPTOMS, min_k = k)
  ok <- lapply(by_site, function(r) r$sensitivity > 0.6 & r$specificity > 0.6)
  any(ok[[1]] & ok[[2]])
}))
put("text.abstract.none_above_60_both", !both_above_60)

pct_pos <- function(d, is_male) 100 * c(mean(d$infect_bin[is_male]), mean(d$infect_bin[!is_male]))
put("text.results.male_female_ssd", pct_pos(ssd, ssd$id_sex == "Male"))
put("text.results.male_female_bd", pct_pos(bd, bd$sex == "1"))
put("text.results.median_delay", c(median(ssd$days_since_ajs_baseline), median(bd$days_since_ajs)))
put("text.limitations.median_delay", c(median(bd$days_since_ajs), median(ssd$days_since_ajs_baseline)))

# IgG is not in the data, so check the text against Table 3 instead
pct_in <- function(s) as.numeric(sub(".*\\((.*)\\).*", "\\1", s))
put("text.results.igg_neg", c(pct_in(pub[["T3.ssd.igg.neg"]]), pct_in(pub[["T3.bd.igg.neg"]])))

put("text.results.msf_spec_ssd", computed[["T4.ssd.MSF.spec"]])
put("text.results.msf_spec_bd", computed[["T4.bd.MSF.spec"]])
put("text.results.imc_extremes", all(sapply(intl[c("ssd", "bd")], function(x)
  x$def[which.max(x$sens)] == "IMC" && x$def[which.min(x$spec)] == "IMC")))
put("text.results.all_misinformative_bd", all(intl$bd$sens + intl$bd$spec < 1))
put("text.results.msf_only_spec40_ssd",
    identical(intl$ssd$def[intl$ssd$sens + intl$ssd$spec > 1 & intl$ssd$spec > 0.4], "MSF"))
who <- function(x) x[x$def == "WHO", ]
put("text.results.children_who",
    who(intl$children)$sens > who(intl$ssd)$sens && who(intl$children)$spec < who(intl$ssd)$spec)

sel_metric <- function(prefixes, m) unlist(lapply(selected[prefixes], `[[`, m))
max_intl_spec <- function(site) max(intl[[site]]$spec)
put("text.results.t5_more_specific",
    min(sel_metric(c("T5.ssd.k1", "T5.ssd.k2"), "specificity")) > max_intl_spec("ssd") &&
      min(sel_metric(c("T5.bd.k1", "T5.bd.k2"), "specificity")) > max_intl_spec("bd"))
below <- function(prefixes, limit)
  all(c(sel_metric(prefixes, "sensitivity"), sel_metric(prefixes, "specificity")) < limit)
put("text.results.t5_below_70_ssd", below(c("T5.ssd.k1", "T5.ssd.k2"), 0.7))
put("text.results.t5_below_60_bd", below(c("T5.bd.k1", "T5.bd.k2"), 0.6))
put("text.results.s5_below_65", below(c("S5.children.k1", "S5.children.k2"), 0.65))
put("text.results.s5_ppv_range", 100 * range(sel_metric(c("S5.children.k1", "S5.children.k2"), "ppv")))
put("text.results.t5_ssd_ppv", 100 * selected[["T5.ssd.k1"]]$ppv)
put("text.results.joint_none_above_60",
    !any(sapply(joint_results, function(r) any(r$sensitivity > 0.6 & r$specificity > 0.6))))
put("text.discussion.npv_example", computed[["S6.joint.k1.npv"]])

# Compare -----------------------------------------------------------------------

NUM <- "[0-9]+(\\.[0-9]+)?"
numbers <- function(s) as.numeric(regmatches(s, gregexpr(NUM, s))[[1]])
decimals <- function(s) {
  m <- regmatches(s, gregexpr(NUM, s))[[1]]
  ifelse(grepl(".", m, fixed = TRUE), nchar(sub(".*\\.", "", m)), 0)
}
round_half_up <- function(x, d) floor(x * 10^d + 0.5 + 1e-9) / 10^d

compare <- function(id, published_value, tol) {
  if (id %in% not_run) return(list(status = "NOT RUN", computed = ""))
  value <- computed[[id]]
  if (is.null(value)) return(list(status = "NOT IN DATA", computed = ""))
  if (is.logical(value) || is.character(value)) {
    shown <- as.character(value)
    return(list(status = if (identical(shown, published_value)) "PASS" else "FAIL", computed = shown))
  }
  p <- numbers(published_value)
  if (length(p) != length(value)) {
    return(list(status = "FAIL", computed = paste(round(value, 1), collapse = ", ")))
  }
  d <- decimals(published_value)
  half_up <- round_half_up(value, d)
  half_even <- round(value, d)
  ok <- all(pmin(abs(half_up - p), abs(half_even - p)) <= tol + 1e-9)
  # Show the recomputed numbers in the published format
  shown <- published_value
  regmatches(shown, gregexpr(NUM, shown)) <- list(sprintf("%.*f", as.integer(d), half_up))
  list(status = if (ok) "PASS" else "FAIL", computed = shown)
}

unused <- setdiff(names(computed), published$id)
if (length(unused)) warning("Computed but not in published_values.csv: ", paste(unused, collapse = ", "))

res <- lapply(seq_len(nrow(published)), function(i)
  compare(published$id[i], published$published[i], published$tolerance[i]))
results <- published %>%
  mutate(computed = vapply(res, `[[`, "", "computed"),
         status = vapply(res, `[[`, "", "status")) %>%
  select(id, source, item, published, computed, status, note)
# Blank cells in a published table are failures whatever the data say
results$status[results$published == "" & results$status == "PASS"] <- "FAIL"

# Internal consistency of the published numbers ---------------------------------

checks <- list()
check <- function(id, what, ok) checks[[length(checks) + 1]] <<- data.frame(id = id, check = what, ok = ok)
matches_pct <- function(pct, x, n) {
  exact <- 100 * x / n; d <- decimals(as.character(pct))
  min(abs(round_half_up(exact, d) - pct), abs(round(exact, d) - pct)) < 1e-9
}

for (i in seq_len(nrow(published))) {
  id <- published$id[i]; s <- published$published[i]
  # estimate (lower, upper): the estimate must lie inside its CI
  for (m in regmatches(s, gregexpr(paste0(NUM, "%? \\(", NUM, ", ", NUM, "\\)"), s))[[1]]) {
    x <- numbers(m)
    check(id, sprintf("%s: estimate inside CI", m), x[2] <= x[1] && x[1] <= x[3])
  }
  # x% (a of b): percentage matches the counts
  for (m in regmatches(s, gregexpr(paste0(NUM, "% \\(", NUM, " of ", NUM, "\\)"), s))[[1]]) {
    x <- numbers(m)
    check(id, sprintf("%s: percentage matches counts", m), matches_pct(x[1], x[2], x[3]))
  }
}

# n (%) cells in Tables 3 and S3 against the column denominators
descr <- published %>% filter(source %in% c("T3", "S3"), grepl("^[0-9]+ ?\\([0-9.]+\\)$", published))
for (i in seq_len(nrow(descr))) {
  parts <- strsplit(descr$id[i], ".", fixed = TRUE)[[1]]
  if (parts[3] %in% c("age_median", "days_median")) next
  n <- as.numeric(pub[[paste(parts[1], parts[2], "N", parts[4], sep = ".")]])
  x <- numbers(descr$published[i])
  check(descr$id[i], sprintf("%s: %% consistent with N = %d", descr$published[i], n), matches_pct(x[2], x[1], n))
}

# 2x2 tables: margins equal group sizes; estimates match counts
group_n <- function(site, g) {
  tab <- if (site == "children") "S3" else "T3"
  as.numeric(pub[[paste(tab, site, "N", g, sep = ".")]])
}
two_by_two <- unique(sub("\\.(tp|fp|fn|tn)$", "", grep("^(T4|S4|T5|S5)\\..*\\.tp$", published$id, value = TRUE)))
for (prefix in two_by_two) {
  site <- strsplit(prefix, ".", fixed = TRUE)[[1]][2]
  n <- sapply(c("tp", "fp", "fn", "tn"), function(f) as.numeric(pub[[paste(prefix, f, sep = ".")]]))
  check(prefix, sprintf("tp + fn = %d IgM positive", group_n(site, "pos")), n[["tp"]] + n[["fn"]] == group_n(site, "pos"))
  check(prefix, sprintf("fp + tn = %d IgM negative", group_n(site, "neg")), n[["fp"]] + n[["tn"]] == group_n(site, "neg"))
  ratios <- list(sens = c("tp", "fn"), spec = c("tn", "fp"), ppv = c("tp", "fp"), npv = c("tn", "fn"))
  for (m in names(ratios)) {
    est_id <- paste(prefix, m, sep = ".")
    if (!est_id %in% names(pub)) next
    a <- n[[ratios[[m]][1]]]; b <- n[[ratios[[m]][2]]]
    check(est_id, sprintf("%s = %s/(%s + %s) from published counts", m, ratios[[m]][1], ratios[[m]][1], ratios[[m]][2]),
          matches_pct(numbers(pub[[est_id]])[1], a, a + b))
  }
}

# Patient flow in the Results
flow <- numbers(pub[["text.results.flow"]])
check("text.results.flow", "South Sudan: 524 - 68 - 35 = 421", flow[1] - flow[3] - flow[5] == as.numeric(pub[["T3.ssd.N.all"]]))
check("text.results.flow", "Bangladesh: 1925 - 280 - 20 = 1625", flow[2] - flow[4] - flow[6] == as.numeric(pub[["T3.bd.N.all"]]))

consistency <- do.call(rbind, checks)

# Report --------------------------------------------------------------------------

write.csv(results, here("audit", "audit_results.csv"), row.names = FALSE)
write.csv(consistency, here("audit", "consistency_results.csv"), row.names = FALSE)

cat("\n== Published numbers vs recomputed values\n\n")
print(as.data.frame.matrix(table(factor(results$source, unique(results$source)),
                                 factor(results$status, c("PASS", "FAIL", "NOT IN DATA", "NOT RUN")))))
cat("\n== FAIL\n\n")
fails <- results %>% filter(status == "FAIL")
for (i in seq_len(nrow(fails))) {
  cat(sprintf("%-32s %s\n%32s published: %s\n%32s computed:  %s\n", fails$id[i], fails$item[i], "",
              if (fails$published[i] == "") "(blank)" else fails$published[i], "", fails$computed[i]))
}
cat("\n== Internal consistency of published numbers:", sum(consistency$ok), "of", nrow(consistency), "pass\n\n")
bad <- consistency %>% filter(!ok)
for (i in seq_len(nrow(bad))) cat(sprintf("%-32s %s\n", bad$id[i], bad$check[i]))
cat("\nNot in the repository's data:", sum(results$status == "NOT IN DATA"), "numbers (see note column).\n")
cat("Full results: audit/audit_results.csv, audit/consistency_results.csv\n")
