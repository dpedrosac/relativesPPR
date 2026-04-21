# Caregiver Burden in Relatives of Patients with Parkinson’s Disease

## Overview
This repository contains the full statistical analysis pipeline for a longitudinal intervention study on caregiver burden, wellbeing, and depressive symptoms in relatives of patients with Parkinson’s disease.

The project emphasizes **reproducibility**, **transparent reporting**, and **robust sensitivity analyses**.

---

## Outcomes

- **Caregiver burden** — Zarit Burden Interview (ZBI)
- **Psychological wellbeing** — WHO-5 Well-Being Index
- **Depressive symptoms** — Beck Depression Inventory II (BDI-II)

---

## Repository Structure

```bash
├── Data/
│   └── Raw questionnaire exports
│
├── Output/
│   └── Tables, figures, statistical results
│
├── 00_setup.R
├── 01_clean_baseline.R
├── 02_longitudinal_anova.R
├── 03_plots_main_analysis.R
├── 04_sensitivity_effectsizes.R
├── 05_outlier_sensitivity.R
├── 06_predictors_zbi.R
├── 07_plots_predictors.R
├── 08_mcid_responder.R
├── 09_mice_sensitivity.R
└── 10_elasticnet_moderators.R
```

---

## Quick Start

Run the main workflow in sequence:

```r
source("00_setup.R")
source("01_clean_baseline.R")
source("02_longitudinal_anova.R")
source("03_plots_main_analysis.R")
```

Optional analyses:

```r
source("04_sensitivity_effectsizes.R")
source("05_outlier_sensitivity.R")
source("09_mice_sensitivity.R")
source("06_predictors_zbi.R")
source("07_plots_predictors.R")
source("10_elasticnet_moderators.R")
source("08_mcid_responder.R")
```

---

## Analytical Workflow

### Setup & Preprocessing
- **00_setup.R**  
  Environment setup, package loading, helper functions, and data import

- **01_clean_baseline.R**  
  Baseline characteristics (Table 1), scale scoring (ZBI, BDI-II), descriptive statistics

### Main Analysis
- **02_longitudinal_anova.R**  
  Data reshaping, pre/post identification, and mixed ANOVA models

  ```r
  score ~ group * time + Error(patient/time)
  ```

- **03_plots_main_analysis.R**  
  Score trajectories and pre–post visualizations

### Sensitivity Analyses
- **04_sensitivity_effectsizes.R**  
  - Regression-based imputation of missing values
  - Detection of extreme change scores
  - Effect sizes: within-group (*dz*), between-group (*Hedges g*)

- **05_outlier_sensitivity.R**  
  Influence of extreme observations on model estimates

- **09_mice_sensitivity.R**  
  Multiple imputation (MICE) with pooled regression results

### Predictor & Exploratory Analyses
- **06_predictors_zbi.R**  

  ```r
  ZBI_pre ~ predictors
  ΔZBI ~ predictors
  ```

- **07_plots_predictors.R**  
  Forest plots and regression visualizations

- **10_elasticnet_moderators.R**  
  Elastic Net models for identifying moderators of treatment response

### Clinical Interpretation
- **08_mcid_responder.R**  
  Responder classification and logistic regression for clinically meaningful change

---

## Software Requirements

- R ≥ 4.2

**Core packages:**
- tidyverse
- afex
- glmmTMB
- mice
- glmnet
- ggplot2
- gtsummary

All dependencies are loaded within the scripts.

---

## Data Availability

Raw data cannot be shared due to ethical and privacy constraints.  
This repository provides analysis code to ensure reproducibility.

---

## Use of Large Language Models

Large language models (LLMs) were used during script preparation for **code polishing, formatting, and readability improvements**.  
All analyses, study design decisions, variable definitions, and statistical reasoning were developed independently by the authors based on their own scientific judgment and the original study data.

---

## Author

Clinical research project on caregiver burden in Parkinson’s disease.

**Supervisor:** David Pedrosa  
Department of Neurology  
Philipps-University Marburg / UKGM

---

## License

Code is provided for scientific transparency and reproducibility.
