# Data loading and symptom definitions -------------------------------------
#
# All analyses use the three analytic datasets in the repo root. Each row is
# one patient with acute jaundice (onset <= 30 days before the visit);
# `infect_bin` is 1 if anti-HEV IgM was detected by ELISA.
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

#' Load South Sudan data (patients >= 14 years, or < 14 with `children = TRUE`).
load_ssd <- function(children = FALSE, binary = FALSE) {
  file <- if (children) "dat_ssd_adjusted_children.rds" else "dat_ssd_adjusted.rds"
  dat <- readRDS(here(file))
  if (binary) dat <- to_binary(dat, SSD_SYMPTOMS)
  dat
}

#' Load Bangladesh data.
load_bd <- function(binary = FALSE) {
  dat <- readRDS(here("dat_bd_adjusted.rds"))
  if (binary) dat <- to_binary(dat, BD_SYMPTOMS)
  dat
}

#' Combined SSD + BD data on shared symptoms (0/1), with post-stratification
#' weights so each study population contributes 50% regardless of size. Weights
#' sum to the total number of patients, so weighted 2x2 counts are in patients.
load_joint <- function() {
  ssd <- load_ssd(binary = TRUE) %>%
    select(infect_bin, all_of(SSD_TO_JOINT)) %>%
    mutate(dataset = "South Sudan")
  bd <- load_bd(binary = TRUE) %>%
    select(infect_bin, all_of(JOINT_SYMPTOMS)) %>%
    mutate(dataset = "Bangladesh")
  bind_rows(ssd, bd) %>%
    group_by(dataset) %>%
    mutate(weight = 0.5 / n()) %>%
    ungroup() %>%
    mutate(weight = weight * n())
}
