# Data loading and symptom definitions -------------------------------------
#
# All analyses use the three analytic datasets in the repo root. Each row is
# one patient with acute jaundice (onset <= 30 days before the visit).
#
# Outcome contract: every loader returns `infect_bin` (0/1) as the only
# laboratory column, so all other columns can safely be used as predictors.
#   - Adults, South Sudan: the RDS holds `baseline_igm_result` and
#     `baseline_pcr_result`; `load_ssd(outcome = )` picks one (IgM by default,
#     PCR for the sensitivity analysis) and drops both raw columns.
#   - Children, South Sudan, and Bangladesh: the RDS holds `infect_bin` (IgM).
#
# Symptoms are stored differently in each file:
#   South Sudan (SSD): factors with levels "No"/"Yes"
#   Bangladesh  (BD):  factors with levels "0"/"1"
# Never compare these factors with `== 1` directly: for the SSD data that is
# always FALSE. Use `binary = TRUE` to get 0/1 integer symptom columns instead.

library(dplyr)
library(here)

# Signs/symptoms available in each study population
SSD_SYMPTOMS <- c(
  "vomiting", "mental", "fever", "abdomen_pain", "cough", "head",
  "epigastric_pain", "fatigue", "loss_appetite", "joint_pains",
  "diarrhea", "nausea", "itch", "convulsion", "bleeding"
)

BD_SYMPTOMS <- c(
  "fever_ever", "nausea_ever", "vomit_ever", "appetite_ever",
  "abdomen_ever", "diarrhea_ever", "constipate_ever", "melaena_ever",
  "head_ever", "convulsion_ever", "unconsc_ever", "drowsy_ever",
  "mental_ever", "exam_liver_palp", "exam_dehydrate", "exam_oedema",
  "exam_dist_abdomen"
)

# Low-prevalence SSD symptoms excluded from the ensemble models (Figure 1)
SSD_RARE_SYMPTOMS <- c("bleeding", "itch", "convulsion")

# Symptoms collected in both populations, using Bangladesh names. Fatigue (SSD)
# and drowsiness (BD) are both used for malaise.
SSD_TO_JOINT <- c(
  vomit_ever = "vomiting", mental_ever = "mental", fever_ever = "fever",
  abdomen_ever = "abdomen_pain", head_ever = "head",
  appetite_ever = "loss_appetite", diarrhea_ever = "diarrhea",
  nausea_ever = "nausea", drowsy_ever = "fatigue"
)
JOINT_SYMPTOMS <- names(SSD_TO_JOINT)

# Demographic columns (dropped for symptom-only models)
SSD_DEMOGRAPHICS <- c("id_sex", "age", "days_since_ajs_baseline")
BD_DEMOGRAPHICS  <- c("sex", "age", "days_since_ajs")

#' Recode factor symptom columns to 0/1 integers.
to_binary <- function(dat, symptoms) {
  dat %>% mutate(across(all_of(symptoms), ~ as.integer(as.character(.) %in% c("Yes", "1"))))
}

# Column holding each reference test in the adult South Sudan RDS
SSD_OUTCOMES <- c(igm = "baseline_igm_result", pcr = "baseline_pcr_result")

#' Load South Sudan data (patients >= 14 years, or < 14 with `children = TRUE`).
#'
#' @param outcome reference test defining `infect_bin`: "igm" (anti-HEV IgM by
#'   ELISA, primary analysis) or "pcr" (HEV RNA, adults only).
load_ssd <- function(children = FALSE, binary = FALSE, outcome = c("igm", "pcr")) {
  outcome <- match.arg(outcome)
  if (children) {
    if (outcome != "igm") stop("PCR results are only available for patients aged 14 and over")
    dat <- readRDS(here("dat_ssd_adjusted_children.rds"))
  } else {
    dat <- readRDS(here("dat_ssd_adjusted.rds")) %>%
      mutate(infect_bin = as.numeric(.data[[SSD_OUTCOMES[[outcome]]]])) %>%
      select(infect_bin, everything(), -all_of(SSD_OUTCOMES))
  }
  stopifnot(all(dat$infect_bin %in% c(0, 1)), !any(SSD_OUTCOMES %in% names(dat)))
  if (binary) dat <- to_binary(dat, SSD_SYMPTOMS)
  dat
}


#' Load Bangladesh data.
load_bd <- function(binary = FALSE) {
  dat <- readRDS(here("dat_bd_adjusted.rds"))
  if (binary) dat <- to_binary(dat, BD_SYMPTOMS)
  dat
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

#' Add one 0/1 column per international definition (`intl_MSF`, ...), using
#' the definition as operationalised in that population.
add_international_flags <- function(dat, defs) {
  for (nm in names(defs)) dat[[paste0("intl_", nm)]] <- as.integer(rowSums(dat[defs[[nm]]]) >= 1)
  dat
}

#' Combined SSD + BD data on shared symptoms (0/1), with post-stratification
#' weights so each study population contributes 50% regardless of size. Weights
#' sum to the total number of patients, so weighted 2x2 counts are in patients.
#' IgM is the reference test in both populations.
#'
#' The `intl_*` columns flag the international definitions as evaluated in each
#' country (Table 4), computed *before* restricting to the shared symptoms, so
#' they keep epigastric pain (MSF) and joint pain (IMC) for South Sudan.
load_joint <- function() {
  intl <- paste0("intl_", names(INTERNATIONAL_DEFS$ssd))
  ssd <- load_ssd(binary = TRUE, outcome = "igm") %>%
    add_international_flags(INTERNATIONAL_DEFS$ssd) %>%
    select(infect_bin, all_of(SSD_TO_JOINT), all_of(intl)) %>%
    mutate(dataset = "South Sudan")
  bd <- load_bd(binary = TRUE) %>%
    add_international_flags(INTERNATIONAL_DEFS$bd) %>%
    select(infect_bin, all_of(JOINT_SYMPTOMS), all_of(intl)) %>%
    mutate(dataset = "Bangladesh")
  bind_rows(ssd, bd) %>%
    group_by(dataset) %>%
    mutate(weight = 0.5 / n()) %>%
    ungroup() %>%
    mutate(weight = weight * n())
}
