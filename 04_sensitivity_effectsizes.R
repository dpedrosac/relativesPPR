#!/usr/bin/env Rscript
# =========================================================
# 03_sensitivity_effectsizes.R
# ---------------------------------------------------------
# Purpose:
#   Run sensitivity analyses for ZBI including imputation,
#   outlier inspection, repeated-measures ANOVA, and effect
#   size estimation.
#
# Provides:
#   - Numeric recoding of ZBI items
#   - Raw and imputed ZBI total scores
#   - Pre/post long-format ZBI data
#   - Identification of largest individual changes
#   - Outlier overview
#   - Repeated-measures ANOVA with and without outliers
#   - Within-group and between-group effect sizes
#
# Notes:
#   - No plots are generated in this script.
#   - Output is written to the Output/ directory.
#
# Input:
#   - Objects created upstream, including:
#       zbi_raw
#       df_pat
#   - 00_setup.R
#
# Output:
#   - Output/zbi_imputed_overview.csv
#   - Output/zbi_top_changes_for_clinical_review.csv
#   - Output/zbi_outliers_top2.csv
#   - Output/anova_sensitivity_summary.rds
#   - Output/zbi_effectsizes.html
#
# Dependencies:
#   00_setup.R
# =========================================================


# =========================================================
# 0) Load setup and configure afex
# =========================================================
source("00_setup.R")
afex::afex_options(type = 3)

if (!dir.exists("Output")) {
  dir.create("Output", recursive = TRUE, showWarnings = FALSE)
}

# =========================================================
# 1) Recode ZBI items to numeric values
# =========================================================

zbi_items <- names(zbi_raw)[stringr::str_detect(names(zbi_raw), "^id\\d+_\\d+")]
stopifnot(length(zbi_items) > 0)

# Keep mapping consistent with the other scripts
map_levels <- c(
  "zu keiner zeit" = 0,
  "nie" = 0,
  "selten" = 1,
  "manchmal" = 2,
  "häufig" = 3,
  "haeufig" = 3,
  "oft" = 3,
  "ziemlich oft" = 3,
  "sehr oft" = 4,
  "fast immer" = 4,
  "immer" = 4
)

zbi_num <- zbi_raw |>
  dplyr::mutate(
    dplyr::across(
      dplyr::all_of(zbi_items),
      ~ {
        v <- stringr::str_to_lower(stringr::str_squish(as.character(.x)))
        num <- suppressWarnings(readr::parse_number(v))
        ifelse(!is.na(num), num, unname(map_levels[v]))
      }
    )
  )


# =========================================================
# 2) Compute raw sum score (rule: <= 4 missing items)
# =========================================================

zbi_scored <- zbi_num |>
  dplyr::rowwise() |>
  dplyr::mutate(
    na_count = sum(is.na(dplyr::c_across(dplyr::all_of(zbi_items)))),
    zbi_total_raw = ifelse(
      na_count <= 4,
      sum(dplyr::c_across(dplyr::all_of(zbi_items)), na.rm = TRUE),
      NA_real_
    )
  ) |>
  dplyr::ungroup()


# =========================================================
# 3) Merge data and perform regression-based imputation
# =========================================================
# Imputation model uses group only.

zbi_merged <- zbi_scored |>
  dplyr::left_join(df_pat, by = c("patient" = "patient_id")) |>
  dplyr::mutate(
    group = factor(group, levels = c("Control", "Intervention"))
  )

mod_zbi <- stats::lm(zbi_total_raw ~ group, data = zbi_merged)
zbi_merged$zbi_pred <- predict(mod_zbi, newdata = zbi_merged)

zbi_merged <- zbi_merged |>
  dplyr::mutate(
    zbi_total_imp = dplyr::if_else(is.na(zbi_total_raw), zbi_pred, zbi_total_raw),
    zbi_total_use = zbi_total_imp
  )

cat("\nMissing ZBI values after imputation (zbi_total_imp):\n")
print(table(is.na(zbi_merged$zbi_total_imp)))

readr::write_csv(
  zbi_merged |>
    dplyr::select(patient, zbi_total_raw, zbi_total_imp, zbi_total_use, group),
  "Output/zbi_imputed_overview.csv"
)

cat("zbi_imputed_overview.csv saved.\n")


# =========================================================
# 4) Create time variable, assign pre/post, and deduplicate
# =========================================================
zbi_merged <- zbi_scored |>
  dplyr::left_join(df_pat, by = c("patient" = "patient_id")) |>
  dplyr::mutate(
    group = factor(group, levels = c("Control", "Intervention"))
  )

mod_zbi <- stats::lm(zbi_total_raw ~ group, data = zbi_merged)
zbi_merged$zbi_pred <- predict(mod_zbi, newdata = zbi_merged)

zbi_merged <- zbi_merged |>
  dplyr::mutate(
    zbi_total_imp = dplyr::if_else(is.na(zbi_total_raw), zbi_pred, zbi_total_raw),
    zbi_total_use = zbi_total_imp
  )

zbi_time <- zbi_merged |>
  get_when() |>
  dplyr::arrange(patient, when_dt)

zbi_long <- zbi_time |>
  add_pre_post() |>
  dplyr::select(patient, time, score = zbi_total_use, when_dt) |>
  dplyr::left_join(df_pat, by = c("patient" = "patient_id")) |>
  dplyr::mutate(
    group = factor(group, levels = c("Control", "Intervention")),
    time = factor(time, levels = c("pre", "post")),
    score = as.numeric(score)
  )
  
# Keep the most recent entry per patient and time point
zbi_long_clean <- zbi_long |>
  dplyr::group_by(patient, group, time) |>
  dplyr::arrange(dplyr::desc(when_dt), .by_group = TRUE) |>
  dplyr::slice(1) |>
  dplyr::ungroup()

# Keep complete pre/post pairs only
zbi_pairs_ids <- zbi_long_clean |>
  dplyr::group_by(patient) |>
  dplyr::filter(dplyr::n_distinct(time) == 2) |>
  dplyr::ungroup() |>
  dplyr::pull(patient) |>
  unique()

zbi_long_pairs <- zbi_long_clean |>
  dplyr::filter(patient %in% zbi_pairs_ids)


# =========================================================
# 5) Create change table for clinical review
# =========================================================

zbi_changes <- zbi_long_clean |>
  tidyr::pivot_wider(
    id_cols = c(patient, group),
    names_from = time,
    values_from = score
  ) |>
  dplyr::mutate(change = post - pre) |>
  dplyr::arrange(dplyr::desc(abs(change)))
  
readr::write_csv(
  zbi_changes |>
    dplyr::select(patient, group, pre, post, change),
  "Output/zbi_top_changes_for_clinical_review.csv"
)

cat("zbi_top_changes_for_clinical_review.csv saved.\n")


# =========================================================
# 6) Identify outliers based on largest changes
# =========================================================

outlier_tbl <- zbi_changes |>
  dplyr::mutate(abs_change = abs(change)) |>
  dplyr::slice_max(order_by = abs_change, n = 2, with_ties = TRUE) |>
  dplyr::arrange(dplyr::desc(abs_change))

cat("\nTop changes (potential outliers):\n")
print(
  outlier_tbl |>
    dplyr::select(patient, group, pre, post, change) |>
    dplyr::mutate(change = round(change, 2))
)

readr::write_csv(outlier_tbl, "Output/zbi_outliers_top2.csv")
cat("File saved: Output/zbi_outliers_top2.csv\n")

outlier_ids <- unique(outlier_tbl$patient)


# =========================================================
# 7) Run ANOVA and compute effect sizes
# =========================================================

cat("\nContingency table (zbi_long_pairs):\n")
print(table(zbi_long_pairs$group, zbi_long_pairs$time, useNA = "ifany"))

# Condition: enough data and both groups represented at pre and post
if (nrow(zbi_long_pairs) >= 3) {
  tab_gt <- table(zbi_long_pairs$group, zbi_long_pairs$time)
  valid_groups <- sum(rowSums(tab_gt) > 0)
} else {
  tab_gt <- matrix(0, nrow = 0, ncol = 0)
  valid_groups <- 0
}

if (
  nrow(zbi_long_pairs) >= 3 &&
  valid_groups >= 2 &&
  all(tab_gt > 0)
) {
  # Full ANOVA with all cases
  fit_all <- afex::aov_car(
    score ~ group * time + Error(patient / time),
    data = zbi_long_pairs,
    factorize = FALSE
  )

  nice_all <- afex::nice(fit_all, es = "pes")

  # Remove outliers
  zbi_no_out <- zbi_long_pairs |>
    dplyr::filter(!patient %in% outlier_ids)

  if (nrow(zbi_no_out) >= 3) {
    tab_gt_no <- table(zbi_no_out$group, zbi_no_out$time)
    valid_groups_no <- sum(rowSums(tab_gt_no) > 0)
  } else {
    tab_gt_no <- matrix(0, nrow = 0, ncol = 0)
    valid_groups_no <- 0
  }

  if (
    nrow(zbi_no_out) >= 3 &&
    valid_groups_no >= 2 &&
    all(tab_gt_no > 0)
  ) {
    fit_no <- afex::aov_car(
      score ~ group * time + Error(patient / time),
      data = zbi_no_out,
      factorize = FALSE
    )

    nice_no <- afex::nice(fit_no, es = "pes")

  } else {
    cat(
      "\nNo meaningful group structure remains after removing outliers; ",
      "second ANOVA is skipped.\n",
      sep = ""
    )
    nice_no <- NULL
  }

  # Effect sizes: within-group dz
  eff_within <- function(df, grp) {
    w <- df |>
      dplyr::filter(group == grp) |>
      dplyr::select(patient, time, score) |>
      tidyr::pivot_wider(names_from = time, values_from = score) |>
      dplyr::mutate(delta = post - pre) |>
      dplyr::pull(delta)

    n_valid <- sum(is.finite(w))
    s <- stats::sd(w, na.rm = TRUE)

    if (n_valid < 2 || is.na(s) || s == 0) {
      return(c(dz = NA_real_, n = n_valid, ci_low = NA_real_, ci_high = NA_real_))
    }

    m <- mean(w, na.rm = TRUE)
    dz <- m / s
    se <- sqrt((1 / n_valid) + (dz^2 / (2 * (n_valid - 1))))
    ci <- dz + c(-1, 1) * 1.96 * se

    c(dz = dz, n = n_valid, ci_low = ci[1], ci_high = ci[2])
  }

  wz_c <- eff_within(zbi_long_pairs, "Control")
  wz_i <- eff_within(zbi_long_pairs, "Intervention")

  # Between-group: Hedges g on change scores
  ch <- zbi_long_pairs |>
    dplyr::select(group, patient, time, score) |>
    tidyr::pivot_wider(names_from = time, values_from = score) |>
    dplyr::mutate(delta = post - pre)

  dC <- ch |>
    dplyr::filter(group == "Control") |>
    dplyr::pull(delta)

  dI <- ch |>
    dplyr::filter(group == "Intervention") |>
    dplyr::pull(delta)

  if (
    length(dC) < 2 ||
    length(dI) < 2 ||
    is.na(stats::var(dC, na.rm = TRUE)) ||
    is.na(stats::var(dI, na.rm = TRUE)) ||
    stats::var(dC, na.rm = TRUE) == 0 ||
    stats::var(dI, na.rm = TRUE) == 0
  ) {
    cat(
      "\nToo few change scores or zero variance in one group; ",
      "Hedges g is not computed.\n",
      sep = ""
    )
    g <- NA_real_
    ci_g <- c(NA_real_, NA_real_)

  } else {
    m_diff <- mean(dI, na.rm = TRUE) - mean(dC, na.rm = TRUE)

    s_pooled <- sqrt(
      ((length(dI) - 1) * stats::var(dI, na.rm = TRUE) +
         (length(dC) - 1) * stats::var(dC, na.rm = TRUE)) /
        (length(dI) + length(dC) - 2)
    )

    d <- m_diff / s_pooled
    J <- 1 - (3 / (4 * (length(dI) + length(dC)) - 9))
    g <- d * J

    se_g <- sqrt(
      (length(dI) + length(dC)) / (length(dI) * length(dC)) +
        (g^2 / (2 * (length(dI) + length(dC) - 2)))
    )

    ci_g <- g + c(-1, 1) * 1.96 * se_g
  }

  # Save results object
  sens_list <- list(
    with_outliers = nice_all,
    without_outliers = nice_no,
    effect_sizes = list(
      within_control = wz_c,
      within_intervention = wz_i,
      between_delta_g = c(
        g = g,
        ci_low = ci_g[1],
        ci_high = ci_g[2],
        n_I = length(dI),
        n_C = length(dC)
      )
    )
  )

  saveRDS(sens_list, "Output/anova_sensitivity_summary.rds")
  cat("anova_sensitivity_summary.rds saved.\n")

  # Compact HTML summary
  es_tbl <- dplyr::tibble(
    Metric = c("dz (within) Control", "dz (within) Intervention", "g (between Δ)"),
    Effect = c(round(wz_c["dz"], 3), round(wz_i["dz"], 3), round(g, 3)),
    CI_low = c(round(wz_c["ci_low"], 3), round(wz_i["ci_low"], 3), round(ci_g[1], 3)),
    CI_high = c(round(wz_c["ci_high"], 3), round(wz_i["ci_high"], 3), round(ci_g[2], 3)),
    n = c(
      wz_c["n"],
      wz_i["n"],
      paste0("I=", length(dI), ", C=", length(dC))
    )
  )

  gt::gtsave(
    gt::gt(es_tbl) |>
      gt::tab_header(
        title = gt::md("**ZBI effect sizes (pre→post)**"),
        subtitle = "dz: within-group, g: between-group (change)"
      ),
    "Output/zbi_effectsizes.html"
  )

  cat("zbi_effectsizes.html saved.\n")

} else {
  cat(
    "\nNot enough meaningful data for ANOVA/effect sizes ",
    "(minimum: 2 groups with pre+post).\n",
    sep = ""
  )
}

cat("\n--- Sensitivity analysis complete ---\n")
