#!/usr/bin/env Rscript
# =========================================================
# 07_plots_predictors.R
# ---------------------------------------------------------
# Purpose:
#   Create forest plots for:
#     A) Baseline predictors of ZBI_pre
#     B) Moderator effects from ANCOVA (group × predictor)
#
# Input:
#   - Output/zbi_predictors_univariate.csv
#   - res_moderator (from analysis script)
#
# Output:
#   - Output/forest_baseline_zbi_A.png
#   - Output/forest_moderators_zbi_B.png
# =========================================================

source("00_setup.R")

if (!dir.exists("Output")) {
  dir.create("Output", recursive = TRUE, showWarnings = FALSE)
}

# =========================================================
# 1) Read baseline predictor results
# =========================================================

uni <- read_csv("Output/zbi_predictors_univariate.csv", show_col_types = FALSE)

stopifnot(all(c("estimate", "conf.low", "conf.high",
                "p.value", "outcome", "predictor") %in% names(uni)))

# =========================================================
# 2) Plot A: Baseline predictors
# =========================================================

uni_base_plot <- uni %>%
  filter(outcome == "zbi_pre") %>%
  mutate(
    sig = case_when(
      is.na(p.value) ~ "n.s.",
      p.value < 0.01 ~ "p < .01",
      p.value < 0.05 ~ "p < .05",
      TRUE           ~ "n.s."
    ),
    label = case_when(
      term == predictor ~ predictor,  # numeric predictors
      TRUE ~ paste0(predictor, ": ", str_remove(term, paste0("^", predictor)))
    ),
    label = str_replace_all(label, "_", " "),
    label = forcats::fct_reorder(label, estimate)
  )

# ---------- manual renaming ----------
manual_renaming <- TRUE
if (manual_renaming) {
  label_map <- c(
    "who5_pre" = "WHO-5 well-being (baseline)",
    "bdi_pre" = "Depression (BDI, baseline)",
    "group_num" = "Intervention group",
    "age" = "Age (years)",
    "sexWeiblich" = "Female sex",
    "marital_statusVerheiratet/ eingetragene Lebenspartnerschaft" = "Married / civil partnership",
    "marital_statusVerwitwet" = "Widowed",
    "educationHauptschulabschluss" = "Primary education",
    "educationRealschulabschluss" = "Secondary education",
    "rel_to_patientLebenspartner / Ehemann" = "Partner / husband",
    "rel_to_patientLebenspartnerin / Ehefrau" = "Partner / wife",
    "rel_to_patientTochter / Sohn" = "Daughter / son",
    "years_known" = "Years known to patient",
    "same_householdTRUE" = "Same household",
    "n_children" = "Number of children",
    "n_children_u21" = "Children under 21",
    "n_adult_care" = "Adults cared for",
    "years_care" = "Years of caregiving",
    "employedTRUE" = "Currently employed",
    "work_hours" = "Weekly work hours"
  )

  uni_base_plot <- uni_base_plot %>%
    mutate(label = dplyr::recode(term, !!!label_map))
}

uni_base_plot <- uni_base_plot %>%
  mutate(
    label       = forcats::fct_reorder(label, estimate),
    stats_label = sprintf("%.2f [%.2f, %.2f]", estimate, conf.low, conf.high)
  )

# ---------- build plot ----------
xlim_vals  <- range(c(uni_base_plot$conf.low, uni_base_plot$conf.high), na.rm = TRUE)
x_text_pos <- xlim_vals[2] + 0.05 * diff(xlim_vals)   # 5% beyond right edge

p_A <- ggplot(uni_base_plot, aes(x = estimate, y = label, color = sig)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.8, color = "grey40") +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high),
                 width = 0.2, linewidth = 0.9, color = "grey30") +
  geom_point(size = 3.0, aes(fill = sig), shape = 21, color = "black") +
  geom_text(aes(y = label, x = x_text_pos, label = stats_label),
            hjust = 0, color = "black", size = 4) +
  scale_color_manual(values = c(
    "n.s."    = "grey60",
    "p < .05" = "darkorange3",
    "p < .01" = "firebrick3"
  )) +
  scale_x_continuous(
    name   = "Coefficient estimate (95% CI)",
    breaks = scales::breaks_pretty(n = 5)
  ) +
  coord_cartesian(
    xlim = c(xlim_vals[1], xlim_vals[2] + 0.15 * diff(xlim_vals)),
    clip = "off"
  ) +
  labs(
    title = "A. Predictors of baseline burden",
    x     = "Coefficient estimate (95% CI)",
    y     = NULL,
    color = "Significance"
  ) +
  theme_minimal(base_size = 20) +
  theme(
    axis.line.x        = element_line(color = "black"),
    axis.ticks.x       = element_line(color = "black"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    legend.position    = "bottom",
    plot.margin        = margin(t = 20, r = 130, b = 20, l = 20)
  )

ggsave("Output/forest_baseline_zbi_A.png",
       p_A, width = 9, height = 6.5, dpi = 300)

# =========================================================
# 3) Read moderator results
# =========================================================

mod_plot <- res_moderator %>%
  filter(grepl("group", term) & grepl(":", term)) %>%
  filter(!(is.na(estimate) & is.na(conf.low) & is.na(conf.high))) %>%
  mutate(
    sig = case_when(
      is.na(p.value) ~ "n.s.",
      p.value < 0.01 ~ "p < .01",
      p.value < 0.05 ~ "p < .05",
      TRUE           ~ "n.s."
    ),
    other_var = ifelse(str_starts(term, "group"),
                       str_remove(term, "^group[^:]*:"),
                       str_remove(term, ":group[^:]*$")),
    other_var = str_trim(other_var),
    other_label = recode(other_var,
      "who5_pre" = "Well‑being (WHO‑5)",
      "bdi_pre" = "Depression (BDI)",
      "age" = "Age (years)",
      "sexWeiblich" = "Sex (female)",
      "marital_statusVerheiratet/ eingetragene Lebenspartnerschaft" = "Married / partnered",
      "marital_statusVerwitwet" = "Widowed",
      "educationHauptschulabschluss" = "Primary education",
      "educationRealschulabschluss" = "Secondary education",
      "rel_to_patientLebenspartner / Ehemann" = "Partner / husband",
      "rel_to_patientLebenspartnerin / Ehefrau" = "Partner / wife",
      "rel_to_patientTochter / Sohn" = "Daughter / son",
      "years_known" = "Years known to patient",
      "same_householdTRUE" = "Same household",
      "n_children" = "Number of children",
      "n_children_u21" = "Children under 21",
      "n_adult_care" = "Adults cared for",
      "care_level" = "Patient care‑level",
      "years_care" = "Years of caregiving",
      "employedTRUE" = "Employed",
      "work_hours" = "Weekly work hours",
      "worktime_reducedTRUE" = "Reduced work hours",
      .default = str_replace_all(other_var, "_", " ")
    ),
    short_label = paste0("Group × ", other_label),
    stats_label = sprintf("%.2f [%.2f, %.2f]", estimate, conf.low, conf.high)
  ) %>%
  arrange(estimate) %>%
  mutate(short_label = factor(short_label, levels = unique(short_label)))

# =========================================================
# 4) Plot B: Moderator effects
# =========================================================

xlim_vals  <- range(c(mod_plot$conf.low, mod_plot$conf.high), na.rm = TRUE)
x_text_pos <- xlim_vals[2] + 0.05 * diff(xlim_vals)   # 5% beyond right limit

p_B <- ggplot(mod_plot, aes(x = estimate, y = short_label, color = sig)) +
  geom_vline(xintercept = 0, linetype = "dashed",
             linewidth = 0.8, color = "grey40") +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high),
                 width = 0.15, linewidth = 0.9, color = "grey30") +
  geom_point(size = 3.2, aes(fill = sig), shape = 21, color = "black") +
  geom_text(aes(y = short_label, x = x_text_pos, label = stats_label),
            hjust = 0, color = "black", size = 4) +
  scale_color_manual(values = c(
    "n.s."    = "grey60",
    "p < .05" = "darkorange3",
    "p < .01" = "firebrick3"
  )) +
  scale_x_continuous(
    name = "Interaction coefficient (95% CI)",
    breaks = scales::breaks_pretty(n = 5)
  ) +
  coord_cartesian(
    xlim = c(xlim_vals[1], xlim_vals[2] + 0.20 * diff(xlim_vals)),
    clip = "off"
  ) +
  labs(
    title = "B. Moderator effects",
    x     = "Interaction coefficient (95% CI)",
    y     = NULL,
    color = "Significance"
  ) +
  theme_minimal(base_size = 20) +
  theme(
    axis.line.x        = element_line(color = "black"),
    axis.ticks.x       = element_line(color = "black"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    legend.position    = "bottom",
    plot.margin        = margin(t = 20, r = 140, b = 20, l = 20)
  )

ggsave("Output/forest_moderators_zbi_B.png",
       p_B, width = 10, height = 7.5, dpi = 300)

# =========================================================
# 5) Log message
# =========================================================

message("Saved:")
message(" - Output/forest_baseline_zbi_A.png")
message(" - Output/forest_moderators_zbi_B.png")

