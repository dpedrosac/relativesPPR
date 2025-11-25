#!/usr/bin/env Rscript
# Script: 05_outlier_sensitivity.R
# Purpose: sensitivity analysis and clinical outlier review
# Inputs: zbi_long (from longitudinal analysis)
# Outputs: Output/zbi_top_changes_for_clinical_review.csv,
#          Output/anova_sensitivity_summary_brief.csv
# Version: 1.0
# Date: 2025-11-24

library(dplyr)
library(tidyr)
library(readr)
library(afex)

# 1) Ensure ZBI data exists
if (!exists("zbi_long")) stop("zbi_long fehlt – vorher 02_longitudinal_anova.R ausführen.")

# 2) Keep only participants with pre AND post
zbi_pairs_ids <- zbi_long %>%
  count(patient) %>%
  filter(n >= 2) %>%
  pull(patient)

zbi_long_pairs <- zbi_long %>%
  filter(patient %in% zbi_pairs_ids)

# 3) Wide form / compute change
zbi_changes <- zbi_long_pairs %>%
  select(patient, group, time, score) %>%
  pivot_wider(names_from = time, values_from = score) %>%
  mutate(change = post - pre) %>%
  arrange(desc(abs(change)))

# 4) Export top 10 changes (with review note column)
zbi_review <- zbi_changes %>%
  mutate(review_note = NA_character_) %>%
  select(patient, group, pre, post, change, review_note)

if (!dir.exists("Output")) dir.create("Output")
write_csv(zbi_review, "Output/zbi_top_changes_for_clinical_review.csv")

cat("\n✅ File 'zbi_top_changes_for_clinical_review.csv' saved – please review clinically.\n")

# 5) ANOVA with and without the top 2 changes
if (nrow(zbi_long_pairs) >= 3) {
  fit_all <- afex::aov_car(
    score ~ group * time + Error(patient/time),
    data = zbi_long_pairs,
    factorize = FALSE
  )
  nice_all <- afex::nice(fit_all, es = "pes")

  # Top-2 outliers (can be adjusted after clinical review)
  outlier_ids <- head(zbi_changes$patient, 2)
  zbi_no_outliers <- zbi_long_pairs %>%
    filter(!patient %in% outlier_ids)

  fit_no <- afex::aov_car(
    score ~ group * time + Error(patient/time),
    data = zbi_no_outliers,
    factorize = FALSE
  )
  nice_no <- afex::nice(fit_no, es = "pes")

  # 6) Save short report with/without outliers
  report <- tibble::tibble(
    model = c("mit Ausreißern", "ohne Ausreißer"),
    p_value_group_time = c(
      nice_all$p.value[nice_all$Effect == "group:time"],
      nice_no$p.value[nice_no$Effect == "group:time"]
    ),
    eta_sq_partial = c(
      nice_all$pes[nice_all$Effect == "group:time"],
      nice_no$pes[nice_no$Effect == "group:time"]
    )
  )

  write_csv(report, "Output/anova_sensitivity_summary_brief.csv")
  cat("\n✅ Sensitivity comparison saved to 'anova_sensitivity_summary_brief.csv'\n")

} else {
  cat("\n⚠️ Not enough complete pre/post pairs for sensitivity ANOVA.\n")
}