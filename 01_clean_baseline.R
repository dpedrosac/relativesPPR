#!/usr/bin/env Rscript
# =========================================================
# 01_clean_baseline.R
# ---------------------------------------------------------
# Purpose:
#   Create the baseline dataset without imputation.
#   Table 1 is based exclusively on observed, non-imputed
#   raw data. Imputation is performed later in
#   02_longitudinal_anova.R.
#
# Notes:
#   - This script performs no imputation.
#   - Baseline values are derived from the earliest available
#     time point per patient.
#
# Changes:
#   2026-03-12:
#     ZBI had multiple combinations, so data were previously
#     not uniquely identifiable.
#
# Input:
#   - cleaned_export_qform_5_2025-10-13-09-54_16_demapped.csv
#   - export_patients_5_2025-10-13-09-52_demapped.csv
#   - cleaned_export_qform_5_2025-10-13-09-54_27_demapped.csv
#   - cleaned_export_qform_5_2025-10-13-09-54_15_demapped.csv
#
# Output:
#   - Output/table1_baseline_<timestamp>.html
#
# Dependencies:
#   00_setup.R
# =========================================================


# =========================================================
# 0) Load setup and helpers
# =========================================================
source("00_setup.R")

# =========================================================
# 1) Harmonize IDs
# =========================================================

idcol_zbi <- detect_patient_id_col(zbi_raw)
if (idcol_zbi != "patient_id") {
  zbi_raw <- zbi_raw |>
    dplyr::rename(patient_id = !!rlang::sym(idcol_zbi))
}

idcol_crf <- detect_patient_id_col(crf_raw)
if (idcol_crf != "patient_id") {
  crf_raw <- crf_raw |>
    dplyr::rename(patient_id = !!rlang::sym(idcol_crf))
}

standardize_id <- function(x) {
  trimws(as.character(x))
}

zbi_raw$patient_id <- standardize_id(zbi_raw$patient_id)
crf_raw$patient_id <- standardize_id(crf_raw$patient_id)
df_pat$patient_id  <- standardize_id(df_pat$patient_id)


# =========================================================
# 2) CRF: age / sex / marital status (robust extraction)
# =========================================================

# 2a) Select age column heuristically
# Exclude questionnaire-like metadata fields where possible
age_candidates <- grep(
  "(alter|geburtsjahr|geburtsdatum|age|birth)",
  names(crf_raw),
  ignore.case = TRUE,
  value = TRUE
)

pick_age <- function(df, cands) {
  if (!length(cands)) return(NA_character_)

  # Restrict to candidates that can plausibly be parsed as numeric
  cands_num <- cands[
    vapply(
      df[cands],
      function(x) {
        is.numeric(x) || is.integer(x) || all(suppressWarnings(!is.na(as.numeric(x))))
      },
      logical(1)
    )
  ]

  if (!length(cands_num)) return(NA_character_)

  # Prefer columns with plausible median (10..110) and non-zero variance
  score <- sapply(cands_num, function(nm) {
    x <- suppressWarnings(as.numeric(df[[nm]]))
    x <- x[is.finite(x)]

    if (!length(x)) return(-Inf)

    med <- stats::median(x, na.rm = TRUE)
    sdv <- stats::sd(x, na.rm = TRUE)

    as.numeric((med >= 10 && med <= 110)) +
      as.numeric(is.finite(sdv) && sdv > 0.5)
  })

  best <- cands_num[which.max(score)]

  if (!length(best) || is.infinite(score[which.max(score)]) || max(score) == 0) {
    NA_character_
  } else {
    best
  }
}

age_col2 <- pick_age(crf_raw, age_candidates)

crf_raw$alter <- if (is.na(age_col2)) {
  NA_real_
} else {
  suppressWarnings(as.numeric(crf_raw[[age_col2]]))
}


# 2b) Detect sex
sex_candidates <- c(
  "id2_2_ihr_geschlecht",
  "geschlecht",
  "sex",
  "gender",
  "biologisches_geschlecht",
  "sexus"
)

sex_col <- intersect(sex_candidates, names(crf_raw))[1]

if (is.na(sex_col) || is.null(sex_col)) {
  # Heuristic token search across character/factor columns
  chr_cols <- names(crf_raw)[
    vapply(crf_raw, function(x) is.character(x) || is.factor(x), logical(1))
  ]

  hit <- NA_character_

  for (nm in chr_cols) {
    v <- tolower(trimws(as.character(crf_raw[[nm]])))
    if (any(grepl("\\b(w|f|m|weiblich|männlich|maennlich|female|male|divers|non[- ]binary|nb)\\b", v))) {
      hit <- nm
      break
    }
  }

  sex_col <- hit
}

if (is.na(sex_col) || is.null(sex_col)) {
  crf_raw$geschlecht <- factor(
    NA_character_,
    levels = c("weiblich", "männlich", "divers")
  )
} else {
  sx <- tolower(trimws(as.character(crf_raw[[sex_col]])))

  sx <- dplyr::case_when(
    sx %in% c("w", "f", "weiblich", "female", "frau") ~ "weiblich",
    sx %in% c("m", "male", "männlich", "maennlich", "mann") ~ "männlich",
    sx %in% c("divers", "d", "non-binary", "non binary", "nb") ~ "divers",
    TRUE ~ NA_character_
  )

  crf_raw$geschlecht <- factor(
    sx,
    levels = c("weiblich", "männlich", "divers")
  )
}


# 2c) Detect marital status
fs_candidates <- c(
  "id3_3_ihr_familienstand",
  "familienstand",
  "family_status",
  "marital_status"
)

fs_col <- intersect(fs_candidates, names(crf_raw))[1]

default_fs_levels <- c(
  "Ledig",
  "Verheiratet",
  "Geschieden",
  "Verwitwet",
  "Sonstiges"
)

if (is.na(fs_col) || is.null(fs_col)) {
  crf_raw$familienstand <- factor(
    NA_character_,
    levels = default_fs_levels
  )
} else {
  fs <- tolower(trimws(as.character(crf_raw[[fs_col]])))

  fs <- dplyr::case_when(
    grepl("ledig", fs) ~ "Ledig",
    grepl("verheiratet|lebenspartnerschaft|eingetragene", fs) ~ "Verheiratet",
    grepl("geschieden", fs) ~ "Geschieden",
    grepl("verwitwet", fs) ~ "Verwitwet",
    TRUE ~ "Sonstiges"
  )

  crf_raw$familienstand <- factor(
    fs,
    levels = default_fs_levels
  )
}


# 2d) Create CRF subset
crf_clean <- crf_raw |>
  dplyr::select(patient_id, alter, geschlecht, familienstand) |>
  dplyr::distinct()


# =========================================================
# 3) ZBI: identify items and map response values robustly
# =========================================================

zbi_items <- names(zbi_raw)[stringr::str_detect(names(zbi_raw), "^id\\d+_\\d+")]
if (length(zbi_items) == 0) {
  zbi_items <- names(zbi_raw)[stringr::str_detect(names(zbi_raw), "^id\\d+")]
}
if (length(zbi_items) == 0) {
  stop("No ZBI item columns found.")
}

map_levels <- c(
  "zu keiner zeit" = 0,   # corresponds to "never"
  "nie" = 0,
  "selten" = 1,
  "manchmal" = 2,
  "häufig" = 3,
  "haeufig" = 3,
  "oft" = 3,
  "ziemlich oft" = 3,
  "sehr oft" = 4,
  "fast immer" = 4,
  "immer" = 4,
  "trifft nicht zu" = 0,
  "nie/überhaupt nicht" = 0,
  "gelegentlich" = 2,
  "hin und wieder" = 2,
  "häufig/oft" = 3,
  "ziemlich häufig" = 3,
  "sehr häufig" = 4,
  "meistens" = 4
)

parse_likert <- function(x) {
  v <- tolower(trimws(as.character(x)))
  num <- suppressWarnings(readr::parse_number(v))  # e.g. "3 - häufig"
  out <- ifelse(!is.na(num), num, unname(map_levels[v]))
  suppressWarnings(as.numeric(out))
}

zbi_num <- zbi_raw |>
  dplyr::mutate(
    dplyr::across(dplyr::all_of(zbi_items), parse_likert)
  )


# =========================================================
# 4) Compute ZBI score (no imputation; max 4 missing items)
# =========================================================

zbi_scored <- zbi_num |>
  dplyr::rowwise() |>
  dplyr::mutate(
    na_count = sum(is.na(dplyr::c_across(dplyr::all_of(zbi_items)))),
    zbi_total_raw = ifelse(
      na_count <= 4,
      sum(dplyr::c_across(dplyr::all_of(zbi_items)), na.rm = TRUE),
      NA_real_
    ),
    zbi_kat_raw = cut(
      zbi_total_raw,
      breaks = c(-Inf, 20, 40, 60, Inf),
      labels = c("keine/geringe", "leicht", "moderat", "schwer"),
      right = TRUE
    )
  ) |>
  dplyr::ungroup()


# =========================================================
# 5a) Read and score BDI (BDI-II, 21 items, 0-3)
# =========================================================

bdi_raw <- read_latest("cleaned_export_qform_5_2025-10-13-09-54_15_demapped.csv")

# Detect and standardize patient_id
idcol_bdi <- detect_patient_id_col(bdi_raw)
if (idcol_bdi != "patient_id") {
  bdi_raw <- bdi_raw |>
    dplyr::rename(patient_id = !!rlang::sym(idcol_bdi))
}

bdi_raw$patient_id <- trimws(as.character(bdi_raw$patient_id))

# Item columns (id<number> or id<number>_<number>)
# Exclude id25_ (comments)
bdi_items <- names(bdi_raw)[
  stringr::str_detect(names(bdi_raw), "^id\\d{1,2}_\\d+") &
    !stringr::str_detect(names(bdi_raw), "^id25_")
]

if (length(bdi_items) == 0) {
  stop("No BDI item columns found (file 15).")
}

# Robust text-to-number mapping
bdi_map <- c(
  "0" = 0,
  "1" = 1,
  "2" = 2,
  "3" = 3,
  "gar nicht" = 0,
  "überhaupt nicht" = 0,
  "ueberhaupt nicht" = 0,
  "trifft nicht zu" = 0,
  "nie" = 0,
  "leicht" = 1,
  "ein wenig" = 1,
  "ein bisschen" = 1,
  "manchmal" = 1,
  "mäßig" = 2,
  "maessig" = 2,
  "ziemlich" = 2,
  "oft" = 2,
  "stark" = 3,
  "sehr stark" = 3,
  "fast immer" = 3,
  "immer" = 3
)

parse_bdi <- function(x) {
  v <- tolower(trimws(as.character(x)))
  num <- suppressWarnings(readr::parse_number(v))  # captures "2 - mäßig" etc.
  out <- ifelse(!is.na(num), num, unname(bdi_map[v]))
  suppressWarnings(as.numeric(out))
}

bdi_num <- bdi_raw |>
  dplyr::mutate(
    dplyr::across(dplyr::all_of(bdi_items), parse_bdi)
  )

# Compute BDI score without imputation:
# allow at most 2 missing items, otherwise NA
# if <= 2 items are missing, sum observed items only
bdi_scored <- bdi_num |>
  dplyr::rowwise() |>
  dplyr::mutate(
    bdi_na_count = sum(is.na(dplyr::c_across(dplyr::all_of(bdi_items)))),
    bdi_total_raw = ifelse(
      bdi_na_count <= 2,
      sum(dplyr::c_across(dplyr::all_of(bdi_items)), na.rm = TRUE),
      NA_real_
    ),
    bdi_kat_raw = dplyr::case_when(
      !is.na(bdi_total_raw) & bdi_total_raw <= 13 ~ "minimal",
      !is.na(bdi_total_raw) & bdi_total_raw <= 19 ~ "leicht",
      !is.na(bdi_total_raw) & bdi_total_raw <= 28 ~ "moderat",
      !is.na(bdi_total_raw) & bdi_total_raw >= 29 ~ "schwer",
      TRUE ~ NA_character_
    )
  ) |>
  dplyr::ungroup() |>
  dplyr::select(patient_id, authored, bdi_total_raw, bdi_kat_raw)


# =========================================================
# 6) Merge data and define baseline
# =========================================================

dat <- zbi_scored |>
  dplyr::left_join(bdi_scored, by = c("patient_id", "authored")) |>
  dplyr::left_join(crf_clean, by = "patient_id") |>
  dplyr::left_join(df_pat, by = "patient_id") |>
  dplyr::mutate(
    group = factor(group, levels = c("Control", "Intervention"))
  )

dat_baseline <- dat |>
  get_when() |>
  dplyr::arrange(patient_id, when_dt) |>
  dplyr::group_by(patient_id) |>
  dplyr::slice_head(n = 1) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    geschlecht = if (!is.factor(geschlecht)) {
      factor(geschlecht, levels = c("weiblich", "männlich", "divers"))
    } else {
      geschlecht
    },
    familienstand = if (!is.factor(familienstand)) {
      factor(familienstand, levels = default_fs_levels)
    } else {
      familienstand
    }
  )


# =========================================================
# 7) Create Table 1 (gtsummary)
# =========================================================

tbl1_baseline <- dat_baseline |>
  dplyr::select(
    zbi_total_raw,
    zbi_kat_raw,
    bdi_total_raw,
    bdi_kat_raw,
    alter,
    geschlecht,
    familienstand,
    group
  ) |>
  gtsummary::tbl_summary(
    by = group,
    type = list(
      gtsummary::all_continuous() ~ "continuous",
      geschlecht ~ "categorical",
      familienstand ~ "categorical"
    ),
    statistic = list(
      gtsummary::all_continuous() ~ "{mean} ± {sd}",
      gtsummary::all_categorical() ~ "{n} ({p}%)"
    ),
    percent = "column",
    missing = "ifany",
    label = list(
      zbi_total_raw ~ "ZBI total score (raw)",
      zbi_kat_raw   ~ "ZBI category (raw)",
      bdi_total_raw ~ "BDI total score (raw)",
      bdi_kat_raw   ~ "BDI category (raw)",
      alter         ~ "Age (years)",
      geschlecht    ~ "Sex",
      familienstand ~ "Marital status"
    )
  ) |>
  gtsummary::add_n() |>
  gtsummary::bold_labels() |>
  gtsummary::modify_caption(
    "**Table 1 (baseline, no imputation).** Mean ± SD; n (%)."
  )

# Optional: add p-values
# tbl1_baseline <- tbl1_baseline |>
#   gtsummary::add_p(test = everything() ~ "wilcox.test")


gt_tbl <- gtsummary::as_gt(tbl1_baseline)
print(gt_tbl)

if (!dir.exists("Output")) {
  dir.create("Output", recursive = TRUE, showWarnings = FALSE)
}

stamp <- format(Sys.time(), "%Y%m%d-%H%M%S")
out_file <- normalizePath(
  file.path("Output", paste0("table1_baseline_", stamp, ".html")),
  winslash = "/",
  mustWork = FALSE
)

gt::gtsave(gt_tbl, out_file, inline_css = TRUE)
message("Table saved to: ", out_file)

