#!/usr/bin/env Rscript
# Script: 03_plots_main_analysis.R
# Purpose: prepare long-format data and create plots for the main analysis
# Version: 1.2
# Date: 2025-11-24
############################################################################################
# Setup
############################################################################################

source("00_setup.R")
need_extra <- c("afex", "ggplot2", "purrr")
missing_extra <- setdiff(need_extra, rownames(installed.packages()))
if (length(missing_extra)) install.packages(missing_extra, dependencies = TRUE)
invisible(lapply(need_extra, library, character.only = TRUE))
afex::afex_options(type = 3)

# 1) Read data locally
data_dir <- "Data/exports_demapped_251013"
read_latest <- function(filename) {
  readr::read_csv(file.path(data_dir, filename), show_col_types = FALSE) |>
    janitor::clean_names()
}

# WHO-5 relatives (questionnaire 12)
who5_raw <- read_latest("cleaned_export_qform_5_2025-10-13-09-54_12_demapped.csv")

# BDI-II relatives (questionnaire 15)
bdi_raw  <- read_latest("cleaned_export_qform_5_2025-10-13-09-54_15_demapped.csv")

# ZBI relatives (questionnaire 16)
zbi_raw  <- read_latest("cleaned_export_qform_5_2025-10-13-09-54_16_demapped.csv")

# CRF relatives (questionnaire 27)
crf_raw  <- read_latest("cleaned_export_qform_5_2025-10-13-09-54_27_demapped.csv")

# Patient table with group assignment
df_pat <- read_latest("export_patients_5_2025-10-13-09-52_demapped.csv") |>
  dplyr::mutate(
    group = dplyr::if_else(active, "Intervention", "Control", missing = "Control"),
    group = factor(group, levels = c("Control", "Intervention"))
  ) |>
  dplyr::select(patient_id, group)

# 2) Prepare scales (same as in 02_longitudinal_anova.R)
# -> ZBI
zbi_items <- names(zbi_raw)[stringr::str_detect(names(zbi_raw), "^id\\d+_\\d+_")]
if (length(zbi_items) == 0) {
  zbi_items <- names(zbi_raw)[stringr::str_detect(names(zbi_raw), "^id\\d+_\\d+($|_)")]
}
map_levels <- c(
  "nie" = 0, "selten" = 1, "manchmal" = 2,
  "häufig" = 3, "haeufig" = 3, "oft" = 3,
  "ziemlich oft" = 3, "sehr oft" = 4,
  "fast immer" = 4, "immer" = 4
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

# -> WHO-5
who5_items <- names(who5_raw)[stringr::str_detect(names(who5_raw), "^id[1-5]_")]
who5_map <- c(
  "die ganze zeit" = 5,
  "meistens" = 4,
  "etwas mehr als die hälfte der zeit" = 3,
  "etwas weniger als die hälfte der zeit" = 2,
  "ab und zu" = 1,
  "zu keinem zeitpunkt" = 0
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
    who5_raw_mean5 = ifelse(n_ans >= 4, (sum_i / n_ans) * 5, NA_real_),
    who5_pct = who5_raw_mean5 * 4
  ) |>
  dplyr::ungroup() |>
  get_when()

# -> BDI-II
bdi_items <- names(bdi_raw)[
  stringr::str_detect(names(bdi_raw), "^id\\d{1,2}_\\d+") &
    !stringr::str_detect(names(bdi_raw), "^id25_")]
k_bdi <- length(bdi_items)

bdi_num <- bdi_raw |>
  dplyr::mutate(
    dplyr::across(
      dplyr::all_of(bdi_items),
      ~ suppressWarnings(readr::parse_number(as.character(.x)))
    )
  )

bdi_imp <- impute_scale(
  df = bdi_num,
  items = bdi_items,
  score_name = "bdi2_total_imp",
  max_missing_count = 4
)

bdi_scored <- bdi_imp |>
  dplyr::select(-starts_with("id25")) |>
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

# -> CRF (Fatigue)
crf_item_names <- c(
  "id48_1_sie_haben_aufgrund_ihrer_pflegeaufgaben_ein_gefuehl_von_muede_sein",
  "id49_1_sie_haben_sich_haeufig_muede_gefuehlt",
  "id50_1_sie_fuehlen_sich_koerperlich_erschoepft",
  "id51_1_sie_haben_haeufig_das_gefuehl_geistig_ausgelaugt_zu_sein",
  "id52_1_sie_empfinden_ihre_pflegeaufgabe_als_koerperlich_anstrengend",
  "id53_1_sie_empfinden_ihre_pflegeaufgabe_als_geistig_anstrengend"
)

crf_item_names <- intersect(crf_item_names, names(crf_raw))

crf_scored <- NULL
if (length(crf_item_names) > 0) {
  crf_num <- crf_raw |>
    dplyr::mutate(
      dplyr::across(
        dplyr::all_of(crf_item_names),
        ~ suppressWarnings(readr::parse_number(as.character(.x)))
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
}

# 3) Attach pre/post + group (long datasets)
make_long_with_group <- function(df_scored, score_var, instrument_label) {
  if (is.null(df_scored)) return(NULL)
  df_scored |>
    add_pre_post() |>
    dplyr::select(patient, time, score = {{ score_var }}, when_dt) |>
    dplyr::left_join(
      df_pat,
      by = c("patient" = "patient_id")
    ) |>
    dplyr::mutate(
      instrument = instrument_label,
      group = factor(group, levels = c("Control", "Intervention")),
      time = factor(time, levels = c("pre", "post")),
      score = as.numeric(score)
    )
}

zbi_long  <- make_long_with_group(zbi_imp,     zbi_total_use, "ZBI")
who5_long <- make_long_with_group(who5_scored, who5_pct,      "WHO-5 (0–100 %)")
bdi_long  <- make_long_with_group(bdi_scored,  bdi2_total,    "BDI-II (0–63)")
crf_long  <- make_long_with_group(crf_scored,  crf_total,     "CRF")

all_long <- list(zbi_long, who5_long, bdi_long, crf_long) |>
  purrr::discard(is.null) |>
  dplyr::bind_rows() |>
  dplyr::filter(!is.na(score), !is.na(time), !is.na(group))

# 4) Plot: mean ± SE per group, pre vs post
p <- ggplot2::ggplot(
  all_long,
  ggplot2::aes(x = time, y = score, color = group, group = group)
) +
  ggplot2::stat_summary(fun = mean, geom = "line", linewidth = 1) +
  ggplot2::stat_summary(fun = mean, geom = "point", size = 3) +
  ggplot2::stat_summary(fun.data = ggplot2::mean_se,
                        geom = "errorbar",
                        width = 0.1) +
  ggplot2::facet_wrap(~ instrument, scales = "free_y") +
  ggplot2::labs(
    title = "Course of scores (pre vs. post) by group",
    x = "Time point",
    y = "Score (mean ± SE)",
    color = "Group"
  ) +
  ggplot2::theme_minimal(base_size = 14) +
  ggplot2::theme(
    strip.text = ggplot2::element_text(face = "bold"),
    legend.position = "bottom"
  )
print(p)

# 5) Save plot for manuscript
if (!dir.exists("Output")) dir.create("Output")
ggplot2::ggsave(
  filename = "Output/verlauf_scores.png",
  plot = p,
  width = 9,
  height = 6,
  dpi = 300
)
message("✅ Plot saved: Output/verlauf_scores.png")