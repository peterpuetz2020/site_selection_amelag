## Compute quality criteria for wastewater treatment plants -----------------

# -------------------------------------------------------------------------
# This script creates site-level metrics that reflect data quality.
# The resulting data is saved for later selection.
# -------------------------------------------------------------------------

# Load shared helpers and paths. Only pacman is loaded explicitly here because
# setup.R manages the remaining dependencies for all scripts.
if (!require("pacman")) install.packages("pacman")
pacman::p_load(here)
source(here("Scripts", "setup.R"), encoding = "UTF-8")

# Metrics used for ranking. Each metric is also ranked later; the mean rank
# becomes the composite quality score.
rel_vars <- c("riip_share", "mae", "bloq_share")

# Read and sort input data.
df <- readRDS(here(read_data_here, "data_example.RDS")) %>%
  arrange(site_name, date, lab)

## Share of values below limit of quantification ----------------------------

# bloq_share captures the proportion of measurements falling below the limit of
# quantification. High shares can indicate weak signal quality.
df <- df %>%
  group_by(site_name) %>%
  mutate(bloq_share = mean(below_loq)) %>%
  ungroup()

## RIIP share calculation ----------------------------------------------------

# RIIP share detects implausible jumps between consecutive samples. Both IQR
# outliers and point-validation outliers (extreme swings followed by reversals)
# are counted to penalise unstable series.
temp <- df %>%
  filter(below_loq == 0) %>%
  group_by(site_name, lab) %>%
  mutate(
    date_lag = as.numeric(date - lag(date, 1)),
    R = (viral_load / lag(viral_load))^(1 / date_lag),
    R_outlier = identify_outliers_IQR_factor(R),
    PV_outlier = ifelse(
      (R > quantile(R, .75, na.rm = TRUE) & lead(R) < quantile(R, .25, na.rm = TRUE)) |
        (R < quantile(R, .25, na.rm = TRUE) & lead(R) > quantile(R, .75, na.rm = TRUE)),
      1, 0
    ),
    PV_outlier = ifelse(is.na(PV_outlier), 0, PV_outlier),
    PV_outlier = accumulate(PV_outlier, ~ if (.x == 1 && .y == 1) 0 else .y),
    is_outlier = ifelse(R_outlier == 1 | PV_outlier == 1, 1, 0),
    is_outlier = ifelse(is.na(is_outlier), 0, is_outlier)
  ) %>%
  filter(!is.na(R)) %>%
  group_by(site_name) %>%
  summarise(riip_share = mean(is_outlier, na.rm = TRUE))

## MAE calculation -----------------------------------------------------------

# Compute median absolute error (MAE) relative to loess-smoothed values.
# Steps:
# 1) Remove values below LOQ and log-transform measurements.
# 2) Drop sparse site/lab combinations where smoothing is unreliable.
# 3) Fit adaptive loess curves and measure absolute deviations.
df <- df %>%
  filter(below_loq == 0) %>%
  mutate(logvalue = log10(viral_load)) %>%
  mutate(date = as.Date(format(
    as.POSIXct(date, format = "%Y-%m-%d %H:%M:%S"),
    format = "%Y-%m-%d"
  ))) %>%
  arrange(site_name, date, lab)

# Drop small site/lab combinations where loess is unreliable due to too few
# observations.
df <- df %>%
  group_by(site_name, lab) %>%
  mutate(n = n()) %>%
  ungroup() %>%
  filter(n >= min_obs) %>%
  select(-n)

# Count non-NAs per site/lab.
df <- df %>%
  group_by(site_name, lab) %>%
  mutate(across(
    logvalue,
    .fns = list(non_nas = ~ sum(!is.na(.))),
    .names = "{fn}_{col}"
  )) %>%
  ungroup()

# Compute loess-smoothed values.
df <- loess_fun("logvalue", "non_nas_logvalue")

# Compute median absolute deviations between observed and smoothed values
# on the log scale for each site.
df_var <- df %>%
  group_by(site_name, lab) %>%
  mutate(dev_abs_15 = abs(loess_logvalue - logvalue)) %>%
  group_by(site_name) %>%
  summarise(mae = median(dev_abs_15, na.rm = TRUE), .groups = "drop")

## Combine metrics and rank --------------------------------------------------

# Merge metrics, rank them so higher rank values indicate better quality, and
# average the ranks into the composite indicator used for site selection.
df <- df %>%
  left_join(df_var, by = "site_name") %>%
  left_join(temp, by = "site_name") %>%
  select(site_name, lab, state, pop_share, mae, riip_share, bloq_share) %>%
  distinct(site_name, .keep_all = TRUE) %>%
  mutate(across(
    all_of(rel_vars),
    .fns = list(rank = ~ rank(-.)),
    .names = "{col}_{fn}"
  )) %>%
  rowwise() %>%
  mutate(rank_mean = mean(c_across(all_of(paste0(rel_vars, "_rank"))))) %>%
  ungroup()

# Persist results for site selection script.
saveRDS(df, here(read_data_here, "quality_data.RDS"))
