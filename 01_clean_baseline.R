# 01_clean_baseline.R — Baseline (ohne Imputation)

# 0) Setup/Helfer laden
source("00_setup.R")

# 1) Daten einlesen

zbi_raw <- read_latest("cleaned_export_qform_5_2025-10-13-09-54_16_demapped.csv")
df_pat  <- read_latest("export_patients_5_2025-10-13-09-52_demapped.csv") |>
  dplyr::mutate(
    group = dplyr::if_else(active, "Intervention", "Control", missing = "Control"),
    group = factor(group, levels = c("Control","Intervention"))
  ) |>
  dplyr::select(patient_id, group)

# CRF (Demografie)
crf_raw <- read_latest("cleaned_export_qform_5_2025-10-13-09-54_27_demapped.csv")


# 2) IDs vereinheitlichen

idcol_zbi <- detect_patient_id_col(zbi_raw)
if (idcol_zbi != "patient_id") {
  zbi_raw <- zbi_raw |> dplyr::rename(patient_id = !!rlang::sym(idcol_zbi))
}
idcol_crf <- detect_patient_id_col(crf_raw)
if (idcol_crf != "patient_id") {
  crf_raw <- crf_raw |> dplyr::rename(patient_id = !!rlang::sym(idcol_crf))
}

standardize_id <- function(x) trimws(as.character(x))
zbi_raw$patient_id <- standardize_id(zbi_raw$patient_id)
crf_raw$patient_id <- standardize_id(crf_raw$patient_id)
df_pat$patient_id  <- standardize_id(df_pat$patient_id)


# 3) CRF: Alter / Geschlecht / Familienstand (robust)

# 3a) Alter-Spalte per Heuristik wählen (keine "questionnaire"-Felder o.ä.)
age_candidates <- grep("(alter|geburtsjahr|geburtsdatum|age|birth)", names(crf_raw),
                       ignore.case = TRUE, value = TRUE)

pick_age <- function(df, cands) {
  if (!length(cands)) return(NA_character_)
  # auf numerisch prüfbare Kandidaten einschränken
  cands_num <- cands[vapply(df[cands], function(x) {
    is.numeric(x) || is.integer(x) || all(suppressWarnings(!is.na(as.numeric(x))))
  }, logical(1))]
  if (!length(cands_num)) return(NA_character_)
  # plausibler Median (10..110) & Varianz > 0.5
  score <- sapply(cands_num, function(nm) {
    x <- suppressWarnings(as.numeric(df[[nm]]))
    x <- x[is.finite(x)]
    if (!length(x)) return(-Inf)
    med <- stats::median(x, na.rm = TRUE)
    sdv <- stats::sd(x, na.rm = TRUE)
    as.numeric((med >= 10 && med <= 110)) + as.numeric(is.finite(sdv) && sdv > 0.5)
  })
  best <- cands_num[which.max(score)]
  if (!length(best) || is.infinite(score[which.max(score)]) || max(score) == 0) NA_character_ else best
}

age_col2 <- pick_age(crf_raw, age_candidates)
crf_raw$alter <- if (is.na(age_col2)) NA_real_ else suppressWarnings(as.numeric(crf_raw[[age_col2]]))

# 3b) Geschlecht erkennen 
sex_candidates <- c("id2_2_ihr_geschlecht","geschlecht","sex","gender","biologisches_geschlecht","sexus")
sex_col <- intersect(sex_candidates, names(crf_raw))[1]
if (is.na(sex_col) || is.null(sex_col)) {
  # heuristisch nach Tokens suchen
  chr_cols <- names(crf_raw)[vapply(crf_raw, function(x) is.character(x) || is.factor(x), logical(1))]
  hit <- NA_character_
  for (nm in chr_cols) {
    v <- tolower(trimws(as.character(crf_raw[[nm]])))
    if (any(grepl("\\b(w|f|m|weiblich|männlich|maennlich|female|male|divers|non[- ]binary|nb)\\b", v))) {
      hit <- nm; break
    }
  }
  sex_col <- hit
}
if (is.na(sex_col) || is.null(sex_col)) {
  crf_raw$geschlecht <- factor(NA_character_, levels = c("weiblich","männlich","divers"))
} else {
  sx <- tolower(trimws(as.character(crf_raw[[sex_col]])))
  sx <- dplyr::case_when(
    sx %in% c("w","f","weiblich","female","frau") ~ "weiblich",
    sx %in% c("m","male","männlich","maennlich","mann") ~ "männlich",
    sx %in% c("divers","d","non-binary","non binary","nb") ~ "divers",
    TRUE ~ NA_character_
  )
  crf_raw$geschlecht <- factor(sx, levels = c("weiblich","männlich","divers"))
}

# 3c) Familienstand 
fs_candidates <- c("id3_3_ihr_familienstand","familienstand","family_status","marital_status")
fs_col <- intersect(fs_candidates, names(crf_raw))[1]
default_fs_levels <- c("Ledig","Verheiratet","Geschieden","Verwitwet","Sonstiges")

if (is.na(fs_col) || is.null(fs_col)) {
  crf_raw$familienstand <- factor(NA_character_, levels = default_fs_levels)
} else {
  fs <- tolower(trimws(as.character(crf_raw[[fs_col]])))
  fs <- dplyr::case_when(
    grepl("ledig", fs) ~ "Ledig",
    grepl("verheiratet|lebenspartnerschaft|eingetragene", fs) ~ "Verheiratet",
    grepl("geschieden", fs) ~ "Geschieden",
    grepl("verwitwet", fs) ~ "Verwitwet",
    TRUE ~ "Sonstiges"
  )
  crf_raw$familienstand <- factor(fs, levels = default_fs_levels)
}

# 3d) CRF-Subset
crf_clean <- crf_raw |>
  dplyr::select(patient_id, alter, geschlecht, familienstand) |>
  dplyr::distinct()


# 4) ZBI: Items & Mapping robust

zbi_items <- names(zbi_raw)[stringr::str_detect(names(zbi_raw), "^id\\d+_\\d+")]
if (length(zbi_items) == 0) zbi_items <- names(zbi_raw)[stringr::str_detect(names(zbi_raw), "^id\\d+")]
if (length(zbi_items) == 0) stop("Keine ZBI-Itemspalten gefunden.")

map_levels <- c(
  "nie"=0,"selten"=1,"manchmal"=2,"häufig"=3,"haeufig"=3,"oft"=3,"ziemlich oft"=3,
  "sehr oft"=4,"fast immer"=4,"immer"=4,
  # Zusatz-Varianten
  "trifft nicht zu"=0,"nie/überhaupt nicht"=0,
  "gelegentlich"=2,"hin und wieder"=2,
  "häufig/oft"=3,"ziemlich häufig"=3,
  "sehr häufig"=4,"meistens"=4
)
parse_likert <- function(x) {
  v <- tolower(trimws(as.character(x)))
  num <- suppressWarnings(readr::parse_number(v))  # z.B. "3 - häufig"
  out <- ifelse(!is.na(num), num, unname(map_levels[v]))
  suppressWarnings(as.numeric(out))
}
zbi_num <- zbi_raw |> dplyr::mutate(dplyr::across(dplyr::all_of(zbi_items), parse_likert))

# 5) ZBI-Score (ohne Imputation; ≤4 fehlende Items)

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
      breaks = c(-Inf,20,40,60,Inf),
      labels = c("keine/geringe","leicht","moderat","schwer"),
      right = TRUE
    )
  ) |>
  dplyr::ungroup()


# 5a) BDI einlesen & scorieren (BDI-II, 21 Items, 0-3)

bdi_raw <- read_latest("cleaned_export_qform_5_2025-10-13-09-54_15_demapped.csv")

# patient_id erkennen & vereinheitlichen
idcol_bdi <- detect_patient_id_col(bdi_raw)
if (idcol_bdi != "patient_id") {
  bdi_raw <- bdi_raw |> dplyr::rename(patient_id = !!rlang::sym(idcol_bdi))
}
bdi_raw$patient_id <- trimws(as.character(bdi_raw$patient_id))

# Item-Spalten (id<Zahl> oder id<Zahl>_<Zahl>)
bdi_items <- names(bdi_raw)[stringr::str_detect(names(bdi_raw), "^id\\d+_\\d+")]
if (length(bdi_items) == 0) {
  stop("Keine BDI-Itemspalten gefunden (Datei 15).")
}

# Text→Zahl Mapping (robust)
bdi_map <- c(
  "0"=0,"1"=1,"2"=2,"3"=3,
  "gar nicht"=0,"überhaupt nicht"=0,"ueberhaupt nicht"=0,"trifft nicht zu"=0,"nie"=0,
  "leicht"=1,"ein wenig"=1,"ein bisschen"=1,"manchmal"=1,
  "mäßig"=2,"maessig"=2,"ziemlich"=2,"oft"=2,
  "stark"=3,"sehr stark"=3,"fast immer"=3,"immer"=3
)

parse_bdi <- function(x) {
  v <- tolower(trimws(as.character(x)))
  num <- suppressWarnings(readr::parse_number(v)) # fängt "2 - mäßig" etc. ab
  out <- ifelse(!is.na(num), num, unname(bdi_map[v]))
  suppressWarnings(as.numeric(out))
}

bdi_num <- bdi_raw |>
  dplyr::mutate(
    dplyr::across(dplyr::all_of(bdi_items), parse_bdi)
  )

# BDI-Score (ohne Imputation). Regel: max 2 fehlende Items, sonst NA.
# (Wenn <=2 fehlen, summieren nur die vorhandenen – konservativ ohne Imputation.)
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
  dplyr::select(patient_id, bdi_total_raw, bdi_kat_raw)

# 6) Merge & Baseline bestimmen

dat <- zbi_scored |>
  dplyr::left_join(bdi_scored, by = "patient_id") |>   #  <<< NEU: BDI dazu
  dplyr::left_join(crf_clean, by = "patient_id") |>
  dplyr::left_join(df_pat,  by = "patient_id") |>
  dplyr::mutate(group = factor(group, levels = c("Control","Intervention")))

dat_baseline <- dat |>
  get_when() |>
  dplyr::arrange(patient_id, when_dt) |>
  dplyr::group_by(patient_id) |>
  dplyr::slice_head(n = 1) |>
  dplyr::ungroup() |>
  # Sicherheits-Cast: Kategorien als Faktor halten (auch bei all-NA)
  dplyr::mutate(
    geschlecht    = if (!is.factor(geschlecht)) factor(geschlecht, levels = c("weiblich","männlich","divers")) else geschlecht,
    familienstand = if (!is.factor(familienstand)) factor(familienstand, levels = default_fs_levels) else familienstand
  )


# 7) Tabelle 1 (gtsummary)

tbl1_baseline <- dat_baseline |>
  dplyr::select(
    zbi_total_raw, zbi_kat_raw,
    bdi_total_raw, bdi_kat_raw,          #  <<< NEU
    alter, geschlecht, familienstand,
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
      zbi_total_raw ~ "ZBI-Gesamtscore (roh)",
      zbi_kat_raw   ~ "ZBI-Kategorie (roh)",
      bdi_total_raw ~ "BDI-Gesamtscore (roh)",        #  <<< NEU
      bdi_kat_raw   ~ "BDI-Kategorie (roh)",          #  <<< NEU
      alter         ~ "Alter (Jahre)",
      geschlecht    ~ "Geschlecht",
      familienstand ~ "Familienstand"
    )
  ) |>
  gtsummary::add_n() |>
  gtsummary::bold_labels() |>
  gtsummary::modify_caption("**Tabelle 1 (Baseline, ohne Imputation).** Mittelwert ± SD; n (%).")

# Optional: p-Werte
# tbl1_baseline <- tbl1_baseline |> gtsummary::add_p(test = everything() ~ "wilcox.test")


# 8) Anzeigen & robust speichern

gt_tbl <- gtsummary::as_gt(tbl1_baseline)
print(gt_tbl)

if (!dir.exists("Output")) dir.create("Output", recursive = TRUE, showWarnings = FALSE)
stamp    <- format(Sys.time(), "%Y%m%d-%H%M%S")
out_file <- normalizePath(file.path("Output", paste0("table1_baseline_", stamp, ".html")),
                          winslash = "/", mustWork = FALSE)
gt::gtsave(gt_tbl, out_file, inline_css = TRUE)
message("✅ Tabelle gespeichert unter: ", out_file)

if (interactive()) {
  if (requireNamespace("rstudioapi", quietly = TRUE)) rstudioapi::viewer(out_file)
  if (file.exists(out_file)) utils::browseURL(out_file)
}
