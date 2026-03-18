# 06_mice_sensitivity.R
# Sensitivitätsanalyse mit Mehrfach-Imputation (MICE)
# - Imputation auf Ebene der Skalen (pre/post)
# - Modell: ΔScore = post - pre ~ group
#   für ZBI, WHO-5, BDI-II

# 0) Setup 

source("00_setup.R")

need_extra <- c("mice", "dplyr", "tidyr", "broom", "readr")
miss <- setdiff(need_extra, rownames(installed.packages()))
if (length(miss)) install.packages(miss, dependencies = TRUE)
invisible(lapply(need_extra, library, character.only = TRUE))

# Hauptanalyse-Skript laden, damit zbi_long / who5_long / bdi_long / crf_long / df_pat existieren
source("02_longitudinal_anova.R")

# Sicherheit:
stopifnot(exists("zbi_long"), exists("who5_long"), exists("bdi_long"), exists("df_pat"))

# 1) Wide-Datensatz (pre/post je Skala) 

make_wide <- function(dlong, prefix) {
  dlong %>%
    dplyr::select(patient, group, time, score) %>%
    tidyr::pivot_wider(
      names_from  = time,
      values_from = score,
      names_prefix = paste0(prefix, "_")
    ) %>%
    # beachte: group ist gleich pro Patient, daher first()
    dplyr::group_by(patient) %>%
    dplyr::summarise(
      group = dplyr::first(group),
      dplyr::across(dplyr::starts_with(prefix), ~ dplyr::first(.x)),
      .groups = "drop"
    )
}

zbi_wide  <- make_wide(zbi_long,  "zbi")
who5_wide <- make_wide(who5_long, "who5")
bdi_wide  <- make_wide(bdi_long,  "bdi")

# CRF nur, falls vorhanden
crf_wide <- NULL
if (exists("crf_long") && !is.null(crf_long)) {
  crf_wide <- make_wide(crf_long, "crf")
}

# alles zusammenführen
dat_mice <- zbi_wide %>%
  dplyr::full_join(who5_wide %>% dplyr::select(-group), by = "patient") %>%
  dplyr::full_join(bdi_wide  %>% dplyr::select(-group), by = "patient")

if (!is.null(crf_wide)) {
  dat_mice <- dat_mice %>%
    dplyr::full_join(crf_wide %>% dplyr::select(-group), by = "patient")
}

# Gruppenzugehörigkeit von df_pat sicherstellen (falls irgendwo NA)
dat_mice <- dat_mice %>%
  dplyr::left_join(df_pat, by = c("patient" = "patient_id")) %>%
  dplyr::mutate(
    group = dplyr::coalesce(group.x, group.y),
    group = factor(group, levels = c("Control","Intervention"))
  ) %>%
  dplyr::select(-group.x, -group.y)

# Fälle ohne Gruppe raus
dat_mice <- dat_mice %>%
  dplyr::filter(!is.na(group))

# optional: Datensatz einmal speichern, um zu sehen, was imputiert wird
if (!dir.exists("Output")) dir.create("Output", recursive = TRUE, showWarnings = FALSE)
readr::write_csv(dat_mice, "Output/mice_input_wide_scales.csv")

cat("Basisdatensatz für MICE gespeichert: Output/mice_input_wide_scales.csv\n")

# 2) Mehrfach-Imputation mit mice 

# nur die relevanten Variablen für MICE auswählen
vars_for_mice <- c(
  "group",
  "zbi_pre","zbi_post",
  "who5_pre","who5_post",
  "bdi_pre","bdi_post"
)

# falls CRF existiert, dazunehmen
vars_for_mice <- intersect(vars_for_mice, names(dat_mice))  # robust
if ("crf_pre" %in% names(dat_mice) && "crf_post" %in% names(dat_mice)) {
  vars_for_mice <- c(vars_for_mice, "crf_pre","crf_post")
}

mice_data <- dat_mice %>%
  dplyr::select(patient, dplyr::all_of(vars_for_mice))

# Patienten behalten, die überhaupt irgendwas haben
mice_data <- mice_data %>%
  dplyr::filter(rowSums(!is.na(dplyr::across(-patient))) > 0)

# group als Faktor
mice_data <- mice_data %>%
  dplyr::mutate(group = factor(group, levels = c("Control","Intervention")))

# predictorMatrix anpassen: patient nicht imputieren / nicht als Prädiktor
ini <- mice::mice(mice_data[,-1], maxit = 0, printFlag = FALSE)
predMat <- ini$predictorMatrix
predMat[,] <- 1
diag(predMat) <- 0  # keine Selbstvorhersage

# MICE laufen lassen (z.B. m = 20)
mids_obj <- mice::mice(
  data   = mice_data[,-1],   # ohne patient
  m      = 20,
  method = "pmm",
  predictorMatrix = predMat,
  seed   = 20251130,
  printFlag = FALSE
)

cat("MICE-Imputation abgeschlossen (20 Datensätze).\n")

# 3) Modelle: ΔScore ~ group  (gepoolt) 

run_delta_model <- function(mids, pre_name, post_name, label) {
  form_call <- bquote({
    delta <- .(as.name(post_name)) - .(as.name(pre_name))
    stats::lm(delta ~ group)
  })
  fit <- with(mids, eval(form_call))
  pooled <- pool(fit)
  sum_tab <- summary(pooled, conf.int = TRUE)
  
  # nur Gruppeneffekt
  res <- sum_tab %>%
    dplyr::filter(term == "groupIntervention") %>%
    dplyr::mutate(instrument = label) %>%
    dplyr::select(
      instrument, term, estimate, std.error, statistic, p.value,
      conf.low, conf.high
    )
  res
}

res_list <- list()

if (all(c("zbi_pre","zbi_post") %in% colnames(mids_obj$data))) {
  res_list[["ZBI"]] <- run_delta_model(mids_obj, "zbi_pre","zbi_post","ZBI")
}
if (all(c("who5_pre","who5_post") %in% colnames(mids_obj$data))) {
  res_list[["WHO-5"]] <- run_delta_model(mids_obj, "who5_pre","who5_post","WHO-5")
}
if (all(c("bdi_pre","bdi_post") %in% colnames(mids_obj$data))) {
  res_list[["BDI-II"]] <- run_delta_model(mids_obj, "bdi_pre","bdi_post","BDI-II")
}
if (all(c("crf_pre","crf_post") %in% colnames(mids_obj$data))) {
  res_list[["CRF"]] <- run_delta_model(mids_obj, "crf_pre","crf_post","CRF")
}

mice_results <- dplyr::bind_rows(res_list)

readr::write_csv(mice_results, "Output/mice_delta_group_effects.csv")
print(mice_results)

cat("Ergebnistabelle gespeichert unter Output/mice_delta_group_effects.csv\n")
cat("    (Effekt = mittlere Veränderung Intervention vs. Control, nach MICE)\n")