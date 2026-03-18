#!/usr/bin/env Rscript
# =========================================================
# 02_longitudinal_anova.R
# ---------------------------------------------------------
# Purpose:
#   Prepare longitudinal scale data for repeated-measures
#   analyses and sensitivity analyses.
#   Provides:
#     - Scale preparation for ZBI, WHO-5, BDI-II, and CRF
#     - Item recoding to numeric values
#     - Scale-level imputation
#     - Raw, imputed, and preferred ZBI scores
#     - Pre/post assignment per patient
#     - Mixed ANOVA (group × time)
#     - Identification of largest individual changes
#     - GLMM sensitivity analyses
#
# Notes:
#   - Table 1 is created first in 01_clean_baseline.R.
#   - Imputation is performed here, after baseline reporting.
#   - ZBI outputs include:
#       zbi_total_raw = raw score, up to 4 missing items allowed
#       zbi_total_imp = imputed score, missing items replaced
#                       by person mean
#       zbi_total_use = preferred score for ANOVA; imputed if
#                       available, otherwise raw
#
# Input:
#   - cleaned_export_qform_5_2025-10-13-09-54_12_demapped.csv
#   - cleaned_export_qform_5_2025-10-13-09-54_15_demapped.csv
#   - cleaned_export_qform_5_2025-10-13-09-54_16_demapped.csv
#   - cleaned_export_qform_5_2025-10-13-09-54_27_demapped.csv
#   - export_patients_5_2025-10-13-09-52_demapped.csv
#
# Output:
#   - Output/anova_ergebnisse.html
#   - Output/zbi_changes.csv
#   - Output/crf_changes.csv (if available)
#   - Output/residuals_shapiro_per_scale.csv
#   - Output/glmmTMB_models_all_scales.rds (if available)
#
# Dependencies:
#   00_setup.R
# =========================================================


# =========================================================
# 0) Load setup and configure afex
# =========================================================
source("00_setup.R")
afex::afex_options(type = 3)


# =========================================================
# 1) Prepare ZBI and perform imputation
# =========================================================

zbi_items <- names(zbi_raw)[stringr::str_detect(names(zbi_raw), "^id\\d+_\\d+($|_)")]
if (length(zbi_items) == 0) {
  stop("No ZBI item columns found.")
}

# Small fix: include "zu keiner zeit" in addition to "nie"
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

zbi_items <- names(zbi_raw) |>
  grep(pattern = "^id\\d{1,2}_\\d{1,2}_", x = _, value = TRUE)

zbi_dat <- zbi_raw |>
  mutate(
    zbi_miss_n = rowSums(is.na(across(all_of(zbi_items))))
  ) |>
  arrange(patient, authored) |>
  group_by(patient) |>
  mutate(
    visit_n = row_number()
  ) |>
  ungroup()

zbi_imp <- impute_scale(
  df = zbi_num,
  items = zbi_items,
  score_name = "zbi_total_imp",
  max_missing_count = 4
) |>
  dplyr::rowwise() |>
  dplyr::mutate(
    na_count = sum(is.na(dplyr::c_across(dplyr::all_of(zbi_items)))),
    zbi_total_raw = ifelse(
      na_count <= 4,
      sum(dplyr::c_across(dplyr::all_of(zbi_items)), na.rm = TRUE),
      NA_real_
    ),
    zbi_total_use = dplyr::coalesce(zbi_total_imp, zbi_total_raw)
  ) |>
  dplyr::ungroup() |>
  get_when()

summary_check <- zbi_imp |>
  dplyr::summarise(
    n_rows = dplyr::n(),
    na_raw = sum(is.na(zbi_total_raw)),
    na_imp = sum(is.na(zbi_total_imp)),
    na_use = sum(is.na(zbi_total_use)),
    na_raw_imp_ok = sum(is.na(zbi_total_raw) & !is.na(zbi_total_imp)),
    too_many_missing = sum(is.na(zbi_total_imp))
  )

print(summary_check)


# =========================================================
# 2) Prepare WHO-5 and perform imputation
# =========================================================

who5_items <- names(who5_raw)[stringr::str_detect(names(who5_raw), "^id[1-5]_")]

who5_map <- c(
  "die ganze zeit" = 5,
  "meistens" = 4,
  "etwas mehr als die hälfte der zeit" = 3,
  "etwas weniger als die hälfte der zeit" = 2,
  "ab und zu" = 1,
  "nie" = 0
)

who5_num <- who5_raw |>
  dplyr::mutate(
    dplyr::across(
      dplyr::all_of(who5_items),
      ~ {
        v <- stringr::str_to_lower(stringr::str_squish(as.character(.x)))
        unname(who5_map[v])
      }
    )
  )

who5_imp <- impute_scale(
  df = who5_num,
  items = who5_items,
  score_name = "who5_total_imp",
  max_missing_prop = 0.2
)

who5_scored <- who5_imp |>
  dplyr::rowwise() |>
  dplyr::mutate(
    n_ans = sum(!is.na(dplyr::c_across(dplyr::all_of(who5_items)))),
    sum_i = sum(dplyr::c_across(dplyr::all_of(who5_items)), na.rm = TRUE),
    who5_raw = ifelse(n_ans >= 4, (sum_i / n_ans) * 5, NA_real_),
    who5_pct = who5_raw * 4
  ) |>
  dplyr::ungroup() |>
  get_when()


# =========================================================
# 3) Prepare BDI-II and perform imputation
# =========================================================

# Small fix: exclude id25 (comment column)
bdi_items <- names(bdi_raw)[
  stringr::str_detect(names(bdi_raw), "^id\\d{1,2}_\\d+") &
    !stringr::str_detect(names(bdi_raw), "^id25_")
]

if (length(bdi_items) < 10) {
  stop("BDI item columns appear incomplete.")
}

bdi_num <- bdi_raw |>
  dplyr::mutate(
    dplyr::across(
      dplyr::all_of(bdi_items),
      ~ {
        suppressWarnings(readr::parse_number(as.character(.x)))
      }
    )
  )

k_bdi <- length(bdi_items)

bdi_imp <- impute_scale(
  df = bdi_num,
  items = bdi_items,
  score_name = "bdi2_total_imp",
  max_missing_count = 4
)

bdi_scored <- bdi_imp |>
  dplyr::rowwise() |>
  dplyr::mutate(
    n_ans = sum(!is.na(dplyr::c_across(dplyr::all_of(bdi_items)))),
    sum_i = sum(dplyr::c_across(dplyr::all_of(bdi_items)), na.rm = TRUE),
    bdi2_total = ifelse(
      n_ans >= ceiling(0.8 * k_bdi),
      (sum_i / n_ans) * k_bdi,
      NA_real_
    )
  ) |>
  dplyr::ungroup() |>
  get_when()


# =========================================================
# 4) Prepare CRF / fatigue and perform imputation
# =========================================================

crf_item_names <- c(
  "id48_1_sie_haben_aufgrund_ihrer_pflegeaufgaben_ein_gefuehl_von_muede_sein",
  "id49_1_sie_haben_sich_haeufig_muede_gefuehlt",
  "id50_1_sie_fuehlen_sich_koerperlich_erschoepft",
  "id51_1_sie_haben_haeufig_das_gefuehl_geistig_ausgelaugt_zu_sein",
  "id52_1_sie_empfinden_ihre_pflegeaufgabe_als_koerperlich_anstrengend",
  "id53_1_sie_empfinden_ihre_pflegeaufgabe_als_geistig_anstrengend"
)

crf_item_names <- intersect(crf_item_names, names(crf_raw))

if (length(crf_item_names) > 0) {
  message("CRF items found: ", paste(crf_item_names, collapse = ", "))

  crf_num <- crf_raw |>
    dplyr::mutate(
      dplyr::across(
        dplyr::all_of(crf_item_names),
        ~ {
          suppressWarnings(readr::parse_number(as.character(.x)))
        }
      )
    )

  crf_imp <- impute_scale(
    df = crf_num,
    items = crf_item_names,
    score_name = "crf_total_imp",
    max_missing_prop = 0.2
  )

  crf_scored <- crf_imp |>
    dplyr::rowwise() |>
    dplyr::mutate(
      n_ans = sum(!is.na(dplyr::c_across(dplyr::all_of(crf_item_names)))),
      sum_i = sum(dplyr::c_across(dplyr::all_of(crf_item_names)), na.rm = TRUE),
      crf_total = ifelse(
        n_ans >= ceiling(0.8 * length(crf_item_names)),
        (sum_i / n_ans) * length(crf_item_names),
        NA_real_
      )
    ) |>
    dplyr::ungroup() |>
    get_when()

} else {
  warning("No matching CRF items found; CRF analysis skipped.")
  crf_scored <- NULL
}


# =========================================================
# 5) Create long-format data and run ANOVA
# =========================================================

make_long_with_group <- function(df_scored, score_var, instrument_label) {
  if (is.null(df_scored)) return(NULL)

  df_long <- df_scored |>
    add_pre_post() |>
    dplyr::select(patient, time, score = {{ score_var }}, when_dt) |>
    dplyr::left_join(df_pat, by = c("patient" = "patient_id")) |>
    dplyr::mutate(
      instrument = instrument_label,
      group = factor(group, levels = c("Control", "Intervention")),
      time = factor(time, levels = c("pre", "post")),
      score = as.numeric(score)
    )

  df_long
}

zbi_long  <- make_long_with_group(zbi_imp, zbi_total_use, "ZBI")
who5_long <- make_long_with_group(who5_scored, who5_pct, "WHO-5 (0–100 %)")
bdi_long  <- make_long_with_group(bdi_scored, bdi2_total, "BDI-II (0–63)")
crf_long  <- make_long_with_group(crf_scored, crf_total, "CRF")

long_list <- list(
  "WHO-5 (0–100 %)" = who5_long,
  "BDI-II (0–63)" = bdi_long,
  "ZBI" = zbi_long,
  "CRF" = crf_long
)

long_list <- long_list[!vapply(long_list, is.null, logical(1))]

run_anova_summary <- function(dlong, instrument_label) {
  ids_both <- dlong |>
    dplyr::count(patient) |>
    dplyr::filter(n >= 2) |>
    dplyr::pull(patient)

  d2 <- dlong |>
    dplyr::filter(patient %in% ids_both) |>
    dplyr::mutate(
      group = factor(group, levels = c("Control", "Intervention")),
      time = factor(time, levels = c("pre", "post")),
      score = as.numeric(score)
    )

  wide_tmp <- d2 |>
    dplyr::select(patient, time, score) |>
    tidyr::pivot_wider(names_from = time, values_from = score)

  if (nrow(dplyr::filter(wide_tmp, !is.na(pre) & !is.na(post))) < 3) {
    return(NULL)
  }

  fit <- afex::aov_car(
    score ~ group * time + Error(patient / time),
    data = d2,
    factorize = FALSE
  )

  res <- afex::nice(fit, es = "pes") |>
    dplyr::as_tibble()

  if (!"Effect" %in% names(res)) {
    if ("term" %in% names(res)) res <- dplyr::rename(res, Effect = term)
    if ("Term" %in% names(res)) res <- dplyr::rename(res, Effect = Term)
  }

  if (!all(c("df1", "df2") %in% names(res))) {
    if (all(c("num.Df", "den.Df") %in% names(res))) {
      res <- dplyr::rename(res, df1 = num.Df, df2 = den.Df)
    } else if ("df" %in% names(res)) {
      res <- tidyr::separate(
        res,
        df,
        into = c("df1", "df2"),
        sep = ",\\s*",
        convert = TRUE,
        fill = "right"
      )
    }
  }

  if (!"pes" %in% names(res) && "ges" %in% names(res)) {
    res <- dplyr::rename(res, pes = ges)
  }

  keep <- intersect(c("Effect", "F", "df1", "df2", "p.value", "pes"), names(res))

  res |>
    dplyr::mutate(instrument = instrument_label) |>
    dplyr::select(instrument, dplyr::all_of(keep))
}

anova_list <- purrr::imap(long_list, ~ run_anova_summary(.x, .y))
anova_results <- purrr::compact(anova_list)

anova_table <- if (length(anova_results)) {
  dplyr::bind_rows(anova_results)
} else {
  dplyr::tibble(
    instrument = character(),
    Effect = character(),
    F = numeric(),
    df1 = numeric(),
    df2 = numeric(),
    p.value = numeric(),
    pes = numeric()
  )
}

anova_table <- anova_table |>
  dplyr::mutate(
    F = round(as.numeric(F), 2),
    df1 = suppressWarnings(round(as.numeric(df1), 0)),
    df2 = suppressWarnings(round(as.numeric(df2), 0)),
    p.value = signif(as.numeric(p.value), 3),
    pes = round(as.numeric(pes), 3)
  )

anova_gt <- anova_table |>
  gt::gt() |>
  gt::tab_header(
    title = gt::md("**ANOVA results (group × time)**"),
    subtitle = "Only patients with complete pre- and post-values"
  ) |>
  gt::sub_missing(columns = dplyr::everything(), missing_text = "-")


# =========================================================
# 6) Identify largest individual ZBI and CRF changes
# =========================================================

zbi_pairs <- zbi_long |>
  dplyr::count(patient) |>
  dplyr::filter(n >= 2) |>
  dplyr::pull(patient)

zbi_wide <- zbi_long |>
  dplyr::filter(patient %in% zbi_pairs) |>
  dplyr::select(patient, group, time, score) |>
  tidyr::pivot_wider(names_from = time, values_from = score) |>
  dplyr::mutate(change = post - pre) |>
  dplyr::arrange(dplyr::desc(abs(change)))

if (!is.null(crf_long)) {
  crf_pairs <- crf_long |>
    dplyr::count(patient) |>
    dplyr::filter(n >= 2) |>
    dplyr::pull(patient)

  crf_wide <- crf_long |>
    dplyr::filter(patient %in% crf_pairs) |>
    dplyr::select(patient, group, time, score) |>
    tidyr::pivot_wider(names_from = time, values_from = score) |>
    dplyr::mutate(change = post - pre) |>
    dplyr::arrange(dplyr::desc(abs(change)))
}


# =========================================================
# 7) Save ANOVA outputs and change files
# =========================================================

if (!dir.exists("Output")) {
  dir.create("Output", recursive = TRUE, showWarnings = FALSE)
}

gt::gtsave(anova_gt, "Output/anova_ergebnisse.html")
utils::write.csv(zbi_wide, "Output/zbi_changes.csv", row.names = FALSE)

if (exists("crf_wide")) {
  utils::write.csv(crf_wide, "Output/crf_changes.csv", row.names = FALSE)
}

message(
  "Files saved: Output/anova_ergebnisse.html, zbi_changes.csv",
  if (exists("crf_wide")) ", crf_changes.csv" else ""
)


# =========================================================
# 8) Residual normality and GLMM models
# =========================================================
resid_diag_list <- purrr::imap(long_list, function(dlong, instrument_label) {
  ids_both <- dlong |>
    dplyr::count(patient) |>
    dplyr::filter(n >= 2) |>
    dplyr::pull(patient)

  d2 <- dlong |>
    dplyr::filter(patient %in% ids_both) |>
    dplyr::mutate(
      group = factor(group, levels = c("Control", "Intervention")),
      time = factor(time, levels = c("pre", "post")),
      score = as.numeric(score)
    )

  if (nrow(d2) < 3) return(NULL)

  lm_fit <- stats::lm(score ~ group * time, data = d2)
  res <- stats::residuals(lm_fit)

  if (length(res) < 3 || length(res) > 5000) {
    shapiro_p <- NA_real_
  } else {
    shapiro_p <- tryCatch(
      stats::shapiro.test(res)$p.value,
      error = function(e) NA_real_
    )
  }

  dplyr::tibble(
    instrument = instrument_label,
    n_obs = nrow(d2),
    shapiro_p = shapiro_p
  )
})

resid_diag_table <- resid_diag_list |>
  purrr::compact() |>
  dplyr::bind_rows()

if (!dir.exists("Output")) {
  dir.create("Output", recursive = TRUE, showWarnings = FALSE)
}

utils::write.csv(
  resid_diag_table,
  "Output/residuals_shapiro_per_scale.csv",
  row.names = FALSE
)

message(
  "Shapiro-Wilk residual summary saved to ",
  "Output/residuals_shapiro_per_scale.csv"
)


# =========================================================
# 8b) GLMM with glmmTMB (nbinom2) as a more robust alternative
# =========================================================

if ("glmmTMB" %in% rownames(installed.packages())) {
  glmm_models <- purrr::imap(long_list, function(dlong, instrument_label) {
    ids_both <- dlong |>
      dplyr::count(patient) |>
      dplyr::filter(n >= 2) |>
      dplyr::pull(patient)

    d2 <- dlong |>
      dplyr::filter(patient %in% ids_both) |>
      dplyr::mutate(
        group = factor(group, levels = c("Control", "Intervention")),
        time = factor(time, levels = c("pre", "post")),
        score = as.numeric(score)
      )

    if (nrow(d2) < 3) return(NULL)

    # GLMM: score ~ group * time + random intercept for patient
    fit_glmm <- tryCatch(
      glmmTMB::glmmTMB(
        score ~ group * time + (1 | patient),
        data = d2,
        family = glmmTMB::nbinom2()
      ),
      error = function(e) {
        message(
          "GLMM for ", instrument_label,
          " could not be estimated: ", e$message
        )
        return(NULL)
      }
    )

    if (is.null(fit_glmm)) return(NULL)

    list(
      instrument = instrument_label,
      model = fit_glmm,
      summary = summary(fit_glmm)
    )
  })

  glmm_models <- purrr::compact(glmm_models)
  saveRDS(glmm_models, file = "Output/glmmTMB_models_all_scales.rds")

  message("GLMM models saved to Output/glmmTMB_models_all_scales.rds")

} else {
  message("Package 'glmmTMB' not installed; GLMM sensitivity analysis skipped.")
}
