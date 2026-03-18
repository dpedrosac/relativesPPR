# 08_mcid_responder.R
# MCID-Analyse + Responder + logistische Regression + Grafiken

source("00_setup.R")
source("02_longitudinal_anova.R")  # liefert zbi_long, who5_long, bdi_long

library(dplyr)
library(tidyr)
library(broom)
library(ggplot2)
library(forcats)

if (!dir.exists("Output")) dir.create("Output")


# 1) Veränderungen berechnen


make_delta <- function(df_long, label) {
  df_long %>%
    select(patient, group, time, score) %>%
    pivot_wider(names_from = time, values_from = score) %>%
    mutate(delta = post - pre,
           instrument = label)
}

zbi_delta  <- make_delta(zbi_long,  "ZBI")
who5_delta <- make_delta(who5_long, "WHO-5")
bdi_delta  <- make_delta(bdi_long,  "BDI-II")


# 2) MCIDs

sd_zbi <- sd(zbi_delta$pre, na.rm = TRUE)
MCID_ZBI  <- 0.5 * sd_zbi       # distribution-based

MCID_WHO5 <- 10                 # WHO guideline
MCID_BDI  <- 5                  # validated cutoff


# 3) Responder definieren


zbi_delta <- zbi_delta %>% mutate(
  responder = delta <= -MCID_ZBI,
  worsener  = delta >= +MCID_ZBI
)

who5_delta <- who5_delta %>% mutate(
  responder = delta >= MCID_WHO5,
  worsener  = delta <= -MCID_WHO5
)

bdi_delta <- bdi_delta %>% mutate(
  responder = delta <= -MCID_BDI,
  worsener  = delta >= +MCID_BDI
)


# 4) Responder-Raten pro Gruppe


count_responder <- function(df) {
  df %>%
    group_by(group) %>%
    summarise(
      N = n(),
      responders = sum(responder, na.rm = TRUE),
      worseners  = sum(worsener,  na.rm = TRUE),
      resp_rate = responders / N
    )
}

resp_zbi  <- count_responder(zbi_delta)
resp_who5 <- count_responder(who5_delta)
resp_bdi  <- count_responder(bdi_delta)

write.csv(resp_zbi,  "Output/mcid_responder_zbi.csv", row.names = FALSE)
write.csv(resp_who5, "Output/mcid_responder_who5.csv", row.names = FALSE)
write.csv(resp_bdi,  "Output/mcid_responder_bdi.csv", row.names = FALSE)


# 5) Logistische Regression (OR)

or_model <- function(df) {
  glm(responder ~ group, data = df, family = binomial()) %>% tidy()
}

or_zbi  <- or_model(zbi_delta) %>% mutate(instrument = "ZBI")
or_who5 <- or_model(who5_delta) %>% mutate(instrument = "WHO-5")
or_bdi  <- or_model(bdi_delta) %>% mutate(instrument = "BDI-II")

or_all <- bind_rows(or_zbi, or_who5, or_bdi)
write.csv(or_all, "Output/mcid_or_all.csv", row.names = FALSE)


# 6) Grafiken


plot_responder <- function(df, title) {
  df %>%
    ggplot(aes(x = group, y = resp_rate, fill = group)) +
    geom_col(width = 0.6) +
    geom_text(aes(label = scales::percent(resp_rate, accuracy = 1)),
              vjust = -0.4, size = 5) +
    scale_y_continuous(labels = scales::percent_format()) +
    theme_minimal(base_size = 16) +
    guides(fill = "none") +
    ggtitle(title)
}

# Barplots
p_zbi  <- plot_responder(resp_zbi,  "Clinically Relevant Improvement (ZBI)")
p_who5 <- plot_responder(resp_who5, "Clinically Relevant Improvement (WHO-5)")
p_bdi  <- plot_responder(resp_bdi,  "Clinically Relevant Improvement (BDI-II)")

ggsave("Output/mcid_plot_zbi.png",  p_zbi,  width = 7, height = 5)
ggsave("Output/mcid_plot_who5.png", p_who5, width = 7, height = 5)
ggsave("Output/mcid_plot_bdi.png",  p_bdi,  width = 7, height = 5)


# 7) Forest Plot der ORs

forest_data <- or_all %>%
  filter(term == "groupIntervention") %>%
  mutate(
    OR = exp(estimate),
    CI_low = exp(estimate - 1.96 * std.error),
    CI_high = exp(estimate + 1.96 * std.error),
    instrument = fct_inorder(instrument)
  )

forest_plot <- forest_data %>%
  ggplot(aes(y = instrument, x = OR)) +
  geom_point(size = 4, color = "darkblue") +
  geom_errorbarh(aes(xmin = CI_low, xmax = CI_high), height = 0.2) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  theme_minimal(base_size = 16) +
  xlab("Odds Ratio (Responder)") +
  ylab("") +
  ggtitle("Odds Ratios for Clinically Relevant Improvement")

ggsave("Output/mcid_forest_or.png", forest_plot, width = 7, height = 5)

cat("MCID-Skript vollständig abgeschlossen.\n")