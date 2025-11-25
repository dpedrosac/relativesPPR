#!/usr/bin/env Rscript
# Script for analyzing ParkProReakt data (project 2022–2025)
# Focus: changes in relatives' burden during the intervention
# Code: Liz Wappler & David Pedrosa
#
# Version: 1.2
# Date: 2025-11-24
# Notes: structural adjustments and additional plausibility checks
############################################################################################
# Setup
############################################################################################

pkgs <- c( # necessary packages to load or install
  "readr", "dplyr", "lubridate", "janitor", "gtsummary", "gt",
  "stringr", "tidyr", "survival", "tableone", "rlang", "purrr"
)

install_missing <- function(p) {
  missing <- setdiff(p, rownames(installed.packages()))
  if (length(missing) > 0) install.packages(missing)
}
install_missing(pkgs)
invisible(lapply(pkgs, library, character.only = TRUE))

options(
  dplyr.summarise.inform = FALSE
)

# 1) Define the folder that contains all raw data (exports from 2025-10-13)
data_dir <- normalizePath(
  "Data/exports_demapped_251013",
  winslash = "/", mustWork = FALSE
)

read_latest <- function(filename) {
  readr::read_csv(
    file.path(data_dir, filename),
    locale = readr::locale(encoding = "UTF-8"),  # if necessary, consider "Latin1"
    na = c("", "NA", "NaN", "NULL"),
    guess_max = 100000,
    show_col_types = FALSE
  ) |>
    janitor::clean_names()
}

# 2) Helper: robustly detect the patient-ID column
detect_patient_id_col <- function(df) {
  cand <- c("patient", "patient_id", "pat_id", "pid", "id_patient", "id")
  hit  <- intersect(cand, names(df))
  if (length(hit) == 0) {
    stop("Keine Patienten-ID-Spalte gefunden. Verfügbare Spalten: ",
         paste(names(df), collapse = ", "))
  }
  hit[1]
}

# 3) Date/time parsing: unify to POSIXct
#    Supports ISO, d.m.Y[ H:M(:S)], m/d/Y and Excel serials
parse_any_datetime <- function(x, tz = "UTC") {
  if (inherits(x, "POSIXct")) return(x)
  if (inherits(x, "Date"))    return(as.POSIXct(x, tz = tz))
  if (is.factor(x)) x <- as.character(x)

  # Excel-serial numbers (plausible dates ~ 1997..2120)
  if (is.numeric(x)) {
    ok  <- !is.na(x) & x > 10000 & x < 80000
    out <- rep(NA_real_, length(x))
    out[ok] <- as.numeric(
      as.POSIXct(as.Date(x[ok], origin = "1899-12-30"), tz = tz)
    )
    return(as.POSIXct(out, origin = "1970-01-01", tz = tz))
  }

  x <- as.character(x)

  # wide parser incl. German and ISO formats with/without time
  p <- suppressWarnings(
    lubridate::parse_date_time(
      x,
      orders = c(
        "Ymd HMS", "Ymd HM", "Ymd",
        "Y-m-d H:M:S", "Y-m-d H:M", "Y-m-d",
        "dmY HMS", "dmY HM", "dmY", "d.m.Y H:M:S", "d.m.Y H:M", "d.m.Y",
        "mdY HMS", "mdY HM", "mdY", "m/d/Y H:M:S", "m/d/Y H:M", "m/d/Y"
      ),
      tz = tz, quiet = TRUE
    )
  )

  if (all(is.na(p))) {
    p <- suppressWarnings(lubridate::ymd(x, tz = tz, quiet = TRUE))
    if (all(is.na(p))) p <- suppressWarnings(lubridate::dmy(x, tz = tz, quiet = TRUE))
    if (all(is.na(p))) p <- suppressWarnings(lubridate::mdy(x, tz = tz, quiet = TRUE))
    if (!all(is.na(p))) return(as.POSIXct(p, tz = tz))
  }

  as.POSIXct(p, tz = tz)
}

# 4) Safely create a 'when' column (POSIXct) from multiple candidates
#    Fallback: use row order per patient if no date is available
get_when <- function(df, tz = "UTC") {
  cand <- intersect(
    c(
      "when_dt", "authored", "created_at", "created", "timestamp", "time",
      "date", "datum", "submitted_at"
    ),
    names(df)
  )

  if (length(cand) == 0) {
    idcol <- detect_patient_id_col(df)
    return(
      df |>
        dplyr::group_by(dplyr::across(dplyr::all_of(idcol))) |>
        dplyr::mutate(
          when_dt = as.POSIXct(dplyr::row_number(), origin = "1970-01-01", tz = tz)
        ) |>
        dplyr::ungroup()
    )
  }

  dd <- df
  for (v in cand) {
    dd[[paste0(v, "_dt")]] <- parse_any_datetime(dd[[v]], tz = tz)
  }
  dt_cols <- paste0(cand, "_dt")

  # Coalesce across similar (POSIXct) vectors
  dd$when_dt <- purrr::reduce(dd[dt_cols], dplyr::coalesce)

  if (all(is.na(dd$when_dt))) {
    idcol <- detect_patient_id_col(dd)
    dd <- dd |>
      dplyr::group_by(dplyr::across(dplyr::all_of(idcol))) |>
      dplyr::mutate(
        when_dt = as.POSIXct(dplyr::row_number(), origin = "1970-01-01", tz = tz)
      ) |>
      dplyr::ungroup()
  }
  dd
}

# 5) Determine pre/post based on earliest/latest row per patient
add_pre_post <- function(df_long, time_col = "when_dt", id_col = NULL) {
  stopifnot(time_col %in% names(df_long))
  if (is.null(id_col)) id_col <- detect_patient_id_col(df_long)

  df_long |>
    dplyr::filter(!is.na(.data[[time_col]])) |>
    dplyr::arrange(.data[[id_col]], .data[[time_col]]) |>
    dplyr::group_by(.data[[id_col]]) |>
    dplyr::mutate(
      rk = dplyr::row_number(),
      mrk = max(rk),
      time = dplyr::case_when(
        mrk == 1  ~ "single",
        rk == 1   ~ "pre",
        rk == mrk ~ "post",
        TRUE      ~ NA_character_
      )
    ) |>
    dplyr::ungroup() |>
    dplyr::filter(time %in% c("pre", "post")) |>
    dplyr::mutate(time = factor(time, levels = c("pre", "post")))
}

# 6) Scale imputation (score_name as string; safe NSE)
#    - items: vector of item column names
#    - score_name: name of the summed score column (string)
#    - max_missing_count / max_missing_prop: tolerance for imputation
impute_scale <- function(df, items, score_name,
                         max_missing_count = NULL, max_missing_prop = NULL) {
  stopifnot(length(items) > 0)
  score_sym <- rlang::sym(score_name)

  n_items <- length(items)
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
      n_obs = sum(!is.na(dplyr::c_across(dplyr::all_of(items)))),
      prop_miss = 1 - n_obs / n_items,
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
        as.numeric(NA)
      )
    ) |>
    dplyr::ungroup() |>
    dplyr::select(-n_obs, -prop_miss, -person_mean)
}

# End setup

# Optional quick test (commented out)
# df_test <- read_latest("fragebogen_scores.csv") |>
#   get_when() |>
#   add_pre_post()
# dplyr::glimpse(df_test)
# summary(df_test$when_dt)
# table(df_test$time, useNA = "ifany")