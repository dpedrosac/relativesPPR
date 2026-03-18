#!/usr/bin/env Rscript
# =========================================================
# 00_setup.R
# ---------------------------------------------------------
# Purpose:
#   Central setup script loaded by all analysis scripts.
#   Provides:
#     - Package installation and loading
#     - Global options
#     - Data reading helpers
#     - Date/time parsing utilities
#     - Pre/post time-point assignment
#     - Scale imputation helpers
#     - Instrument-specific helpers (BDI, ZBI)
#     - read raw data into workspace
#
# Input:
#   - Data/exports_demapped_251013/ (CSV files)
#
# =========================================================


# =========================================================
# 0) Install missing packages and load all dependencies
# =========================================================

pkgs <- c(
  "readr", "dplyr", "lubridate", "janitor", "broom",
  "gtsummary", "gt", "stringr", "tidyr",
  "survival", "tableone", "rlang", "purrr",
  "glmmTMB", "ggplot2", "scales", "afex"
)

install_missing <- function(p) {
  missing_pkgs <- setdiff(p, rownames(installed.packages()))
  if (length(missing_pkgs) > 0) install.packages(missing_pkgs)
}

install_missing(pkgs)
invisible(lapply(pkgs, library, character.only = TRUE))

options(dplyr.summarise.inform = FALSE)


# =========================================================
# 1) Paths and data reading
# =========================================================

# Directory containing the latest CSV exports (2025-10-13)
data_dir <- normalizePath(
  "Data/exports_demapped_251013",
  winslash  = "/",
  mustWork  = FALSE
)

read_latest <- function(filename) {
  # Read a CSV from the data directory with robust encoding and NA handling
  readr::read_csv(
    file.path(data_dir, filename),
    locale         = readr::locale(encoding = "UTF-8"),  # fallback: "Latin1"
    na             = c("", "NA", "NaN", "NULL"),
    guess_max      = 100000,
    show_col_types = FALSE
  ) |>
    janitor::clean_names()
}


# =========================================================
# 2) Helper: detect patient ID column
# =========================================================

detect_patient_id_col <- function(df) {
  # Returns the first matching candidate column name for patient ID
  candidates <- c("patient", "patient_id", "pat_id", "pid", "id_patient", "id")
  hit <- intersect(candidates, names(df))

  if (length(hit) == 0) {
    stop(
      "No patient ID column found. Available columns: ",
      paste(names(df), collapse = ", ")
    )
  }

  hit[1]
}


# =========================================================
# 3) Date/time parsing: unified POSIXct output
# =========================================================
#
# Supports: ISO, German d.m.Y, m/d/Y, and Excel serial numbers.

parse_any_datetime <- function(x, tz = "UTC") {

  if (inherits(x, "POSIXct")) return(x)
  if (inherits(x, "Date"))    return(as.POSIXct(x, tz = tz))
  if (is.factor(x))           x <- as.character(x)

  # Excel serial numbers (plausible range ~1997–2120)
  if (is.numeric(x)) {
    ok  <- !is.na(x) & x > 10000 & x < 80000
    out <- rep(NA_real_, length(x))
    out[ok] <- as.numeric(
      as.POSIXct(as.Date(x[ok], origin = "1899-12-30"), tz = tz)
    )
    return(as.POSIXct(out, origin = "1970-01-01", tz = tz))
  }

  x <- as.character(x)

  # Broad parser covering German and ISO formats, with and without time
  p <- suppressWarnings(lubridate::parse_date_time(
    x,
    orders = c(
      "Ymd HMS", "Ymd HM", "Ymd",
      "Y-m-d H:M:S", "Y-m-d H:M", "Y-m-d",
      "dmY HMS", "dmY HM", "dmY",
      "d.m.Y H:M:S", "d.m.Y H:M", "d.m.Y",
      "mdY HMS", "mdY HM", "mdY",
      "m/d/Y H:M:S", "m/d/Y H:M", "m/d/Y"
    ),
    tz    = tz,
    quiet = TRUE
  ))

  # Fallback: try standard single-format parsers if all NA
  if (all(is.na(p))) {
    p <- suppressWarnings(lubridate::ymd(x, tz = tz, quiet = TRUE))
    if (all(is.na(p))) p <- suppressWarnings(lubridate::dmy(x, tz = tz, quiet = TRUE))
    if (all(is.na(p))) p <- suppressWarnings(lubridate::mdy(x, tz = tz, quiet = TRUE))
    if (!all(is.na(p))) return(as.POSIXct(p, tz = tz))
  }

  as.POSIXct(p, tz = tz)
}


# =========================================================
# 4) Create "when_dt" column (POSIXct) from multiple candidates
# =========================================================
#
# Falls back to row-number ordering per patient if no date found.

get_when <- function(df, tz = "UTC") {

  candidates <- intersect(
    c("when_dt", "authored", "created_at", "created",
      "timestamp", "time", "date", "datum", "submitted_at"),
    names(df)
  )

  # No date column available: use row order as fallback
  if (length(candidates) == 0) {
    id_col <- detect_patient_id_col(df)
    return(
      df |>
        dplyr::group_by(dplyr::across(dplyr::all_of(id_col))) |>
        dplyr::mutate(
          when_dt = as.POSIXct(dplyr::row_number(), origin = "1970-01-01", tz = tz)
        ) |>
        dplyr::ungroup()
    )
  }

  # Parse each candidate column and coalesce results
  dd <- df
  for (v in candidates) {
    dd[[paste0(v, "_dt")]] <- parse_any_datetime(dd[[v]], tz = tz)
  }

  dt_cols      <- paste0(candidates, "_dt")
  dd$when_dt   <- purrr::reduce(dd[dt_cols], dplyr::coalesce)

  # Second fallback: if all parsed values are NA, use row order
  if (all(is.na(dd$when_dt))) {
    id_col <- detect_patient_id_col(dd)
    dd <- dd |>
      dplyr::group_by(dplyr::across(dplyr::all_of(id_col))) |>
      dplyr::mutate(
        when_dt = as.POSIXct(dplyr::row_number(), origin = "1970-01-01", tz = tz)
      ) |>
      dplyr::ungroup()
  }

  dd
}


# =========================================================
# 5) Assign pre/post labels based on earliest/latest row
# =========================================================

add_pre_post <- function(df_long, time_col = "when_dt", id_col = NULL) {
  # Labels first row per patient as "pre" and last as "post".
  # Rows in between are dropped. Single-row patients are labelled "single".

  stopifnot(time_col %in% names(df_long))

  if (is.null(id_col)) id_col <- detect_patient_id_col(df_long)

  df_long |>
    dplyr::filter(!is.na(.data[[time_col]])) |>
    dplyr::arrange(.data[[id_col]], .data[[time_col]]) |>
    dplyr::group_by(.data[[id_col]]) |>
    dplyr::mutate(
      rk   = dplyr::row_number(),
      mrk  = max(rk),
      time = dplyr::case_when(
        mrk == 1  ~ "single",
        rk  == 1  ~ "pre",
        rk  == mrk ~ "post",
        TRUE       ~ NA_character_
      )
    ) |>
    dplyr::ungroup() |>
    dplyr::filter(time %in% c("pre", "post")) |>
    dplyr::mutate(time = factor(time, levels = c("pre", "post")))
}


# =========================================================
# 6) Scale imputation (mean substitution within tolerance)
# =========================================================
#
# Arguments:
#   df                 - data frame
#   items              - character vector of item column names
#   score_name         - name for the resulting sum score column
#   max_missing_count  - max number of missing items allowed
#   max_missing_prop   - max proportion of missing items allowed

impute_scale <- function(df, items, score_name,
                         max_missing_count = NULL,
                         max_missing_prop  = NULL) {
  stopifnot(length(items) > 0)

  score_sym <- rlang::sym(score_name)
  n_items   <- length(items)

  thr <- if (!is.null(max_missing_count)) {
    max_missing_count / n_items
  } else if (!is.null(max_missing_prop)) {
    max_missing_prop
  } else {
    0
  }

  df |>
    dplyr::rowwise() |>
    dplyr::mutate(
      n_obs       = sum(!is.na(dplyr::c_across(dplyr::all_of(items)))),
      prop_miss   = 1 - n_obs / n_items,
      person_mean = ifelse(
        n_obs > 0,
        mean(dplyr::c_across(dplyr::all_of(items)), na.rm = TRUE),
        NA_real_
      )
    ) |>
    dplyr::mutate(
      dplyr::across(
        dplyr::all_of(items),
        ~ ifelse(is.na(.x) & prop_miss <= thr, person_mean, .x)
      ),
      !!score_sym := dplyr::if_else(
        prop_miss <= thr,
        sum(dplyr::c_across(dplyr::all_of(items)), na.rm = TRUE),
        NA_real_
      )
    ) |>
    dplyr::ungroup() |>
    dplyr::select(-n_obs, -prop_miss, -person_mean)
}


# =========================================================
# 7) BDI helper: select item columns (exclude id25*)
# =========================================================

get_bdi_items <- function(df) {
  # Returns column names matching BDI item pattern, excluding id25_* range
  names(df)[
    stringr::str_detect(names(df), "^id\\d{1,2}_\\d+") &
      !stringr::str_detect(names(df), "^id25_")
  ]
}


# =========================================================
# 8) ZBI helper: recode response categories to 0–4
# =========================================================
#
# Maps both "zu keiner zeit" and "nie" to 0
# to handle minor label variations across data exports.

recode_zbi_item <- function(x) {
  dplyr::recode(
    x,
    "zu keiner zeit" = 0L,
    "nie"            = 0L,
    "selten"         = 1L,
    "manchmal"       = 2L,
    "ziemlich oft"   = 3L,
    "fast immer"     = 4L,
    .default         = NA_integer_
  )
}


# =========================================================
# Optional: quick sanity check (uncomment to run)
# =========================================================

# df_test <- read_latest("fragebogen_scores.csv") |>
#   get_when() |>
#   add_pre_post()
# dplyr::glimpse(df_test)
# summary(df_test$when_dt)
# table(df_test$time, useNA = "ifany")



# =========================================================
# 9) Read data
# =========================================================

who5_raw <- read_latest("cleaned_export_qform_5_2025-10-13-09-54_12_demapped.csv")
bdi_raw  <- read_latest("cleaned_export_qform_5_2025-10-13-09-54_15_demapped.csv")
zbi_raw  <- read_latest("cleaned_export_qform_5_2025-10-13-09-54_16_demapped.csv")
crf_raw  <- read_latest("cleaned_export_qform_5_2025-10-13-09-54_27_demapped.csv")

df_pat <- read_latest("export_patients_5_2025-10-13-09-52_demapped.csv") |>
  dplyr::mutate(
    group = dplyr::if_else(active, "Intervention", "Control", missing = "Control"),
    group = factor(group, levels = c("Control", "Intervention"))
  ) |>
  dplyr::select(patient_id, group)

