#!/usr/bin/env Rscript
# Script: 04_sensitivity_effectsizes.R
# Purpose: sensitivity analyses, imputation, outlier checks and effect sizes (no plots)
# Version: 1.0
# Date: 2025-11-24
############################################################################################
# Setup
############################################################################################

source("00_setup.R")
need_extra <- c("afex", "gt")
miss <- setdiff(need_extra, rownames(installed.packages()))
if (length(miss)) install.packages(miss, dependencies = TRUE)
invisible(lapply(need_extra, library, character.only = TRUE))
afex::afex_options(type = 3)

# 1) Read data
data_dir <- "Data/exports_demapped_251013"
read_latest <- function(filename) {
  readr::read_csv(file.path(data_dir, filename), show_col_types = FALSE) |>
    janitor::clean_names()
}

zbi_raw <- read_latest("cleaned_export_qform_5_2025-10-13-09-54_16_demapped.csv")
df_pat  <- read_latest("export_patients_5_2025-10-13-09-52_demapped.csv") |>
  dplyr::mutate(
    group = dplyr::if_else(active, "Intervention", "Control", missing = "Control"),
    group = factor(group, levels = c("Control", "Intervention"))
  ) |>
  dplyr::select(patient_id, group)

if (!dir.exists("Output")) dir.create("Output", recursive = TRUE, showWarnings = FALSE)

# 2) ZBI items -> numeric
zbi_items <- names(zbi_raw)[stringr::str_detect(names(zbi_raw), "^id\\d+_\\d+")]
stopifnot(length(zbi_items) > 0)

map_levels <- c(
  "nie" = 0, "selten" = 1, "manchmal" = 2,
  "häufig" = 3, "haeufig" = 3, "oft" = 3, "ziemlich oft" = 3,
  "sehr oft" = 4, "fast immer" = 4, "immer" = 4
)

zbi_num <- zbi_raw |>
  dplyr::mutate(dplyr::across(dplyr::all_of(zbi_items), ~{
    v <- stringr::str_to_lower(stringr::str_squish(as.character(.x)))
    num <- suppressWarnings(readr::parse_number(v))
    ifelse(!is.na(num), num, unname(map_levels[v]))
  }))

# 3) Unimputed sum score (rule: ≤4 missing)
zbi_scored <- zbi_num |>
  dplyr::rowwise() |>
  dplyr::mutate(
    na_count = sum(is.na(dplyr::c_across(dplyr::all_of(zbi_items)))),
    zbi_total_raw = ifelse(na_count <= 4, sum(dplyr::c_across(dplyr::all_of(zbi_items)), na.rm = TRUE), NA_real_)
  ) |>
  dplyr::ungroup()

# 4) Merge + regression imputation (by group)
zbi_merged <- zbi_scored |>
  dplyr::left_join(df_pat, by = c("patient" = "patient_id")) |>
  dplyr::mutate(group = factor(group, levels = c("Control", "Intervention")))

mod_zbi <- stats::lm(zbi_total_raw ~ group, data = zbi_merged)
zbi_merged$zbi_pred <- predict(mod_zbi, newdata = zbi_merged)

zbi_merged <- zbi_merged |>
  dplyr::mutate(
    zbi_total_imp = dplyr::if_else(is.na(zbi_total_raw), zbi_pred, zbi_total_raw),
    zbi_total_use = zbi_total_imp
  )

cat("\nFehlende ZBI nach Imputation (zbi_total_imp):\n")
print(table(is.na(zbi_merged$zbi_total_imp)))

readr::write_csv(
  zbi_merged |> dplyr::select(patient, zbi_total_raw, zbi_total_imp, zbi_total_use, group),
  "Output/zbi_imputed_overview.csv"
)
cat("✅ zbi_imputed_overview.csv gespeichert.\n")

# 5) Time variable + pre/post + deduplication
zbi_time <- zbi_scored |> get_when() |> dplyr::arrange(patient, when_dt)

zbi_long <- zbi_time |>
  dplyr::left_join(zbi_merged |> dplyr::select(patient, zbi_total_use), by = "patient") |>
  add_pre_post() |>
  dplyr::select(patient, time, score = zbi_total_use, when_dt) |>
  dplyr::left_join(df_pat, by = c("patient" = "patient_id")) |>
  dplyr::mutate(
    group = factor(group, levels = c("Control", "Intervention")),
    time  = factor(time,  levels = c("pre", "post")),
    score = as.numeric(score)
  )

# for each patient & time: keep the latest entry
zbi_long_clean <- zbi_long |>
  dplyr::group_by(patient, group, time) |>
  dplyr::arrange(dplyr::desc(when_dt), .by_group = TRUE) |>
  dplyr::slice(1) |>
  dplyr::ungroup()

# only complete pairs
zbi_pairs_ids <- zbi_long_clean |>
  dplyr::group_by(patient) |>
  dplyr::filter(dplyr::n_distinct(time) == 2) |>
  dplyr::ungroup() |>
  dplyr::pull(patient) |>
  unique()

zbi_long_pairs <- zbi_long_clean |> dplyr::filter(patient %in% zbi_pairs_ids)

# 6) Change list (for clinical review)
zbi_changes <- zbi_long_clean |>
  dplyr::select(patient, when_dt, group, time, score) |>
  tidyr::pivot_wider(names_from = time, values_from = score, values_fn = dplyr::last) |>
  dplyr::mutate(change = post - pre) |>
  dplyr::arrange(dplyr::desc(abs(change)))

readr::write_csv(
  zbi_changes |> dplyr::select(patient, when_dt, group, pre, post, change),
  "Output/zbi_top_changes_for_clinical_review.csv"
)
cat("✅ zbi_top_changes_for_clinical_review.csv gespeichert.\n")

# === Outliers (top changes) identify, show, save ===
outlier_tbl <- zbi_changes %>%
  dplyr::mutate(abs_change = abs(change)) %>%
  dplyr::slice_max(order_by = abs_change, n = 2, with_ties = TRUE) %>%
  dplyr::arrange(dplyr::desc(abs_change))

cat("\n===== Top-Changes (potential outliers) =====\n")
print(
  outlier_tbl %>%
    dplyr::select(patient, group, pre, post, change) %>%
    dplyr::mutate(change = round(change, 2))
)

readr::write_csv(outlier_tbl, "Output/zbi_outliers_top2.csv")
cat("✅ Datei gespeichert: Output/zbi_outliers_top2.csv\n")

outlier_ids <- unique(outlier_tbl$patient)

# 7) ANOVA + effect sizes
if (nrow(zbi_long_pairs) >= 3) {
  fit_all <- afex::aov_car(score ~ group * time + Error(patient/time),
                           data = zbi_long_pairs, factorize = FALSE)
  nice_all <- afex::nice(fit_all, es = "pes")

  # remove strongest 2 changes
  zbi_no_out <- zbi_long_pairs |> dplyr::filter(!patient %in% outlier_ids)
  fit_no <- afex::aov_car(score ~ group * time + Error(patient/time),
                          data = zbi_no_out, factorize = FALSE)
  nice_no <- afex::nice(fit_no, es = "pes")

  eff_within <- function(df, grp) {
    w <- df |> dplyr::filter(group == grp) |>
      dplyr::select(patient, time, score) |>
      tidyr::pivot_wider(names_from = time, values_from = score) |>
      dplyr::mutate(delta = post - pre) |> dplyr::pull(delta)
    m <- mean(w, na.rm = TRUE); s <- stats::sd(w, na.rm = TRUE)
    dz <- m / s
    n  <- sum(is.finite(w))
    # 95%-CI (Morris & DeShon approximation)
    se <- sqrt((1/n) + (dz^2 / (2*(n-1))))
    ci <- dz + c(-1,1) * 1.96 * se
    c(dz = dz, n = n, ci_low = ci[1], ci_high = ci[2])
  }
  wz_c <- eff_within(zbi_long_pairs, "Control")
  wz_i <- eff_within(zbi_long_pairs, "Intervention")

  ch <- zbi_long_pairs |>
    dplyr::select(group, patient, time, score) |>
    tidyr::pivot_wider(names_from = time, values_from = score) |>
    dplyr::mutate(delta = post - pre)

  dC <- ch |> dplyr::filter(group == "Control") |> dplyr::pull(delta)
  dI <- ch |> dplyr::filter(group == "Intervention") |> dplyr::pull(delta)
  m_diff <- mean(dI, na.rm = TRUE) - mean(dC, na.rm = TRUE)
  s_pooled <- sqrt(((length(dI)-1)*stats::var(dI, na.rm = TRUE) +
                     (length(dC)-1)*stats::var(dC, na.rm = TRUE)) /
                    (length(dI) + length(dC) - 2))
  d <- m_diff / s_pooled
  J <- 1 - (3 / (4*(length(dI) + length(dC)) - 9))  # Hedges correction
  g <- d * J
  se_g <- sqrt((length(dI) + length(dC)) / (length(dI) * length(dC)) + (g^2 / (2*(length(dI) + length(dC) - 2))))
  ci_g <- g + c(-1,1) * 1.96 * se_g

  sens_list <- list(with_outliers = nice_all, without_outliers = nice_no,
                    effect_sizes = list(
                      within_control = wz_c,
                      within_intervention = wz_i,
                      between_delta_g = c(g = g, ci_low = ci_g[1], ci_high = ci_g[2], n_I = length(dI), n_C = length(dC))
                    ))
  saveRDS(sens_list, "Output/anova_sensitivity_summary.rds")
  cat("✅ anova_sensitivity_summary.rds gespeichert.\n")

  es_tbl <- dplyr::tibble(
    Kennzahl = c("dz (within) Control", "dz (within) Intervention", "g (between Δ)"),
    Effekt  = c(round(wz_c["dz"], 3), round(wz_i["dz"], 3), round(g, 3)),
    CI_low  = c(round(wz_c["ci_low"], 3), round(wz_i["ci_low"], 3), round(ci_g[1], 3)),
    CI_high = c(round(wz_c["ci_high"], 3), round(wz_i["ci_high"], 3), round(ci_g[2], 3)),
    n       = c(wz_c["n"], wz_i["n"], paste0("I=", length(dI), ", C=", length(dC)))
  )

  gt::gtsave(
    gt::gt(es_tbl) |>
      gt::tab_header(title = gt::md("**ZBI effect sizes (pre→post)**"),
                     subtitle = "dz: within-group, g: between-group (Change)"),
    "Output/zbi_effectsizes.html"
  )
  cat("✅ zbi_effectsizes.html gespeichert.\n")

} else {
  cat("\n⚠ Nicht genug vollständige pre/post-Paare für ANOVA/Effektstärken.\n")
}

cat("\n--- Sensitivity analysis END ---\n")