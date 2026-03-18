#!/usr/bin/env Rscript
# =========================================================
# 06_predictors_zbi.R
# ---------------------------------------------------------
# Purpose:
#   Univariate predictor analysis for:
#     a) Baseline caregiver burden (ZBI_pre)
#     b) Change in burden (ΔZBI = post – pre)
#   Includes clinical scales, socio-demographic factors,
#   and caregiving-related variables.
#
# Usage:
#   source("06_predictors_zbi.R")
#
# Inputs:
#   - Objects from 02_longitudinal_anova.R:
#       zbi_long, who5_long, bdi_long, df_pat, crf_raw
#
# Outputs:
#   - Output/zbi_predictors_univariate.csv
#   - Output/zbi_moderator_ancova.csv
# =========================================================

source("00_setup.R")


# =========================================================
# 1) Load longitudinal analysis results
# =========================================================
#
# Provides: zbi_long, who5_long, bdi_long, df_pat, crf_raw
# Each long-format object contains: patient, time, score

source("02_longitudinal_anova.R")


# =========================================================
# 2) Build baseline dataset (pre-values only)
# =========================================================

baseline <- zbi_long %>%
  filter(time == "pre") %>%
  select(patient, zbi_pre = score) %>%
  left_join(
    zbi_long %>%
      filter(time == "post") %>%
      select(patient, zbi_post = score),
    by = "patient"
  ) %>%
  left_join(
    who5_long %>%
      filter(time == "pre") %>%
      select(patient, who5_pre = score),
    by = "patient"
  ) %>%
  left_join(
    bdi_long %>%
      filter(time == "pre") %>%
      select(patient, bdi_pre = score),
    by = "patient"
  ) %>%
  left_join(
    df_pat %>% select(patient_id, group),
    by = c("patient" = "patient_id")
  )


# =========================================================
# 3) Add socio-demographic and caregiving variables from CRF
# =========================================================

# Detect patient ID column in CRF automatically
idcol_crf <- detect_patient_id_col(crf_raw)

crf_for_join <- crf_raw %>%
  mutate(
    patient = as.character(.data[[idcol_crf]]),

    # --- Demographics ---
    age            = suppressWarnings(as.numeric(id1_1_ihr_alter)),
    sex            = as.factor(id2_2_ihr_geschlecht),
    marital_status = as.factor(id3_3_ihr_familienstand),
    education      = as.factor(id4_4_ihr_hochster_schulabschluss),

    # --- Relationship to patient ---
    rel_to_patient = as.factor(
      id9_8_welche_beziehung_haben_sie_zur_patientin_zum_patienten_sind_sie
    ),
    years_known = suppressWarnings(
      as.numeric(
        id85_11_wie_lange_kennen_sie_schon_die_patientin_den_patienten_in_jahren
      )
    ),

    # --- Household and additional care recipients ---
    same_household = as.factor(
      id16_12_leben_sie_mit_der_patientin_dem_patienten_in_einem_haushalt_zusammen
    ),
    n_children = suppressWarnings(
      as.numeric(id24_15_1_fur_wie_viele_im_haushalt_lebende_kinder)
    ),
    n_children_u21 = suppressWarnings(
      as.numeric(id25_15_1_1_fur_wie_viele_kinder_unter_21_jahren)
    ),
    n_adult_care = suppressWarnings(
      as.numeric(
        id26_15_2_fur_wie_viele_im_haushalt_lebende_pflegebedurftige_erwachsene
      )
    ),

    # --- Care level and duration ---
    care_level = suppressWarnings(
      as.numeric(id22_14_welchen_pflegegrad_hat_die_patientin_der_patient)
    ),
    years_care = suppressWarnings(
      as.numeric(
        id89_13_seit_wie_vielen_jahren_pflegen_sie_die_patientin_den_patienten
      )
    ),

    # --- Employment ---
    employed  = as.factor(id52_23_sind_sie_selbst_derzeit_erwerbstatig),
    work_hours = suppressWarnings(
      as.numeric(
        id117_wie_viele_stunden_sind_s_ie_in_der_woche_erwerbstatig_angabe_als_dezimalzahl_z_b_1_5
      )
    ),
    worktime_reduced = as.factor(
      id60_haben_sie_im_letzten_monat_ihre_wochentliche_arbeitszeit_wegen_ihrer_pflegerischen_aufgaben_verkurzt
    )
  ) %>%
  select(
    patient,
    age, sex, marital_status, education,
    rel_to_patient, years_known,
    same_household, n_children, n_children_u21, n_adult_care,
    care_level, years_care,
    employed, work_hours, worktime_reduced
  )

# Join CRF variables into baseline dataset
baseline <- baseline %>%
  left_join(crf_for_join, by = "patient")


# =========================================================
# 4) Compute ΔZBI (post – pre)
# =========================================================

zbi_change <- zbi_long %>%
  select(patient, time, score) %>%
  pivot_wider(names_from = time, values_from = score) %>%
  mutate(delta_zbi = post - pre)

baseline <- baseline %>%
  left_join(
    zbi_change %>% select(patient, delta_zbi),
    by = "patient"
  ) %>%
  filter(!is.na(zbi_pre))   # retain only cases with ZBI baseline


# =========================================================
# 5) Derive additional predictors
# =========================================================

baseline <- baseline %>%
  mutate(
    # numeric group coding: 0 = control, 1 = intervention
    group_num = if_else(group == "Intervention", 1, 0, missing = 0)
  )


# =========================================================
# 6) Define predictor lists
# =========================================================

# Full predictor list (clinical + demographic + caregiving)
predictors <- c(

  # Clinical baseline scales
  "who5_pre",
  "bdi_pre",

  # Group assignment
  "group_num",

  # Demographics
  "age",
  "sex",
  "marital_status",
  "education",

  # Relationship to patient
  "rel_to_patient",
  "years_known",

  # Household and caregiving context
  "same_household",
  "n_children",
  "n_children_u21",
  "n_adult_care",

  # Caregiving intensity
  "care_level",
  "years_care",

  # Employment
  "employed",
  "work_hours",
  "worktime_reduced"
)

# Moderator list: excludes group_num (already in model as main effect)
moderator_predictors <- c(
  "who5_pre",
  "bdi_pre",
  "age",
  "sex",
  "marital_status",
  "education",
  "rel_to_patient",
  "years_known",
  "same_household",
  "n_children",
  "n_children_u21",
  "n_adult_care",
  "care_level",
  "years_care",
  "employed",
  "work_hours",
  "worktime_reduced"
)


# =========================================================
# 7) Analysis functions
# =========================================================

run_uni_lm <- function(data, outcome, predictor) {
  # Fit a univariate linear model and return tidy results.
  # Returns NULL silently if the model fails (e.g. zero variance).
  form <- reformulate(predictor, response = outcome)
  fit  <- try(lm(form, data = data), silent = TRUE)

  if (inherits(fit, "try-error")) return(NULL)

  broom::tidy(fit, conf.int = TRUE) %>%
    filter(term != "(Intercept)") %>%
    mutate(
      outcome   = outcome,
      predictor = predictor
    )
}

run_moderator_ancova <- function(data, predictor) {
  # Fit an ANCOVA model with group × predictor interaction.
  # Model: zbi_post ~ zbi_pre + group * predictor
  # Returns NULL silently if the model fails.
  form <- as.formula(paste("zbi_post ~ zbi_pre + group *", predictor))
  fit  <- try(lm(form, data = data), silent = TRUE)

  if (inherits(fit, "try-error")) return(NULL)

  broom::tidy(fit, conf.int = TRUE) %>%
    mutate(
      outcome   = "zbi_post",
      predictor = predictor,
      model     = "Moderator ANCOVA"
    )
}


# =========================================================
# 8) Run analyses
# =========================================================

res_baseline <- map_dfr(
  predictors,
  ~ run_uni_lm(baseline, outcome = "zbi_pre", predictor = .x)
) %>%
  mutate(model = "Baseline ZBI (univariate)")

res_delta <- map_dfr(
  predictors,
  ~ run_uni_lm(baseline, outcome = "delta_zbi", predictor = .x)
) %>%
  mutate(model = "Delta ZBI (univariate)")

res_moderator <- map_dfr(
  moderator_predictors,
  ~ run_moderator_ancova(baseline, .x)
)

# Extract interaction terms only
res_moderator_interaction <- res_moderator %>%
  filter(grepl("^group.*:", term) | grepl(":group", term))

# Combine all results
res_all <- bind_rows(
  res_baseline,
  res_delta,
  res_moderator_interaction
) %>%
  relocate(outcome, predictor, model, .before = term)


# =========================================================
# 9) Save results
# =========================================================

if (!dir.exists("Output")) {
  dir.create("Output", recursive = TRUE, showWarnings = FALSE)
}

write_csv(res_all,                    "Output/zbi_predictors_univariate.csv")
write_csv(res_moderator_interaction,  "Output/zbi_moderator_ancova.csv")

message("Saved:")
message(" - Output/zbi_predictors_univariate.csv")
message(" - Output/zbi_moderator_ancova.csv")

