#!/usr/bin/env Rscript
# =========================================================
# 05_outlier_sensitivity.R
# ---------------------------------------------------------
# Purpose:
#   Run a sensitivity analysis and clinical outlier review
#   for ZBI change scores.
#
# Provides:
#   - Selection of patients with complete pre/post ZBI data
#   - Wide-format change-score table
#   - Export of cases for clinical review
#   - Repeated-measures ANOVA with and without the largest
#     changes
#   - Brief CSV summary of sensitivity results
#
# Input:
#   - zbi_long (created in 02_longitudinal_anova.R)
#
# Output:
#   - Output/zbi_top_changes_for_clinical_review.csv
#   - Output/anova_sensitivity_summary_brief.csv
#
# Dependencies:
#   00_setup.R
#   zbi_long must already exist in the environment
# =========================================================


# =========================================================
# 1) Ensure ZBI data are available
# =========================================================

if (!exists("zbi_long")) {
  stop("zbi_long is missing — run 02_longitudinal_anova.R first.")
}


# =========================================================
# 2) Keep only patients with both pre and post values
# =========================================================

zbi_pairs_ids <- zbi_long %>%
  dplyr::count(patient) %>%
  dplyr::filter(n >= 2) %>%
  dplyr::pull(patient)

zbi_long_pairs <- zbi_long %>%
  dplyr::filter(patient %in% zbi_pairs_ids)


# =========================================================
# 3) Convert to wide format and compute change scores
# =========================================================

zbi_changes <- zbi_long_pairs %>%
  dplyr::select(patient, group, time, score) %>%
  tidyr::pivot_wider(names_from = time, values_from = score) %>%
  dplyr::mutate(change = post - pre) %>%
  dplyr::arrange(dplyr::desc(abs(change)))


# =========================================================
# 4) Export cases for clinical review
# =========================================================

zbi_review <- zbi_changes %>%
  dplyr::mutate(review_note = NA_character_) %>%
  dplyr::select(patient, group, pre, post, change, review_note)

if (!dir.exists("Output")) {
  dir.create("Output")
}

readr::write_csv(zbi_review, "Output/zbi_top_changes_for_clinical_review.csv")

cat("\nFile 'zbi_top_changes_for_clinical_review.csv' saved — please review clinically.\n")


# =========================================================
# 5) Run ANOVA with and without the two largest changes
# =========================================================

if (nrow(zbi_long_pairs) >= 3) {
  fit_all <- afex::aov_car(
    score ~ group * time + Error(patient / time),
    data = zbi_long_pairs,
    factorize = FALSE
  )

  nice_all <- afex::nice(fit_all, es = "pes")

  # Top-2 potential outliers
  # Adjust after clinical review if needed
  outlier_ids <- head(zbi_changes$patient, 2)

  zbi_no_outliers <- zbi_long_pairs %>%
    dplyr::filter(!patient %in% outlier_ids)

  fit_no <- afex::aov_car(
    score ~ group * time + Error(patient / time),
    data = zbi_no_outliers,
    factorize = FALSE
  )

  nice_no <- afex::nice(fit_no, es = "pes")


  # =========================================================
  # 6) Save brief report with and without outliers
  # =========================================================

  report <- tibble::tibble(
    model = c("with outliers", "without outliers"),
    p_value_group_time = c(
      nice_all$p.value[nice_all$Effect == "group:time"],
      nice_no$p.value[nice_no$Effect == "group:time"]
    ),
    eta_sq_partial = c(
      nice_all$pes[nice_all$Effect == "group:time"],
      nice_no$pes[nice_no$Effect == "group:time"]
    )
  )

  readr::write_csv(report, "Output/anova_sensitivity_summary_brief.csv")

  cat("\nSensitivity summary saved to 'anova_sensitivity_summary_brief.csv'\n")

} else {
  cat("\nNot enough complete pre/post pairs for sensitivity ANOVA.\n")
}
