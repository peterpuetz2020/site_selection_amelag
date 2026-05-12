## Plot aggregated viral load curves ---------------------------------------

# -------------------------------------------------------------------------
# This script aggregates viral load trajectories across all wastewater sites
# and produces a comparison plot between all sites and the selected subset.
# -------------------------------------------------------------------------

# Load shared helpers and paths. setup.R handles dependencies and helper functions.
if (!require("pacman")) install.packages("pacman")
pacman::p_load(here)
source(here("Scripts", "setup.R"), encoding = "UTF-8")

# Read and prepare the full dataset.
df <- readRDS(here(read_data_here, "data.RDS")) %>%
  mutate(log_viral_load = log10(viral_load + 1))
  
# Create weekly grouping based on Thursdays. th_week is a running index that
# increments each Thursday and is later used to align weekly aggregates.
thursday_data <- df %>%
  distinct(date) %>%
  mutate(
    Tag = lubridate::wday(date, week_start = 1),
    is_thursday = ifelse(Tag == 4, 1, 0),
    th_week = cumsum(is_thursday)
  ) %>%
  dplyr::select(date, th_week)

# Join weekly grouping metadata.
df_agg <- df %>% left_join(thursday_data, by = "date")

# Add warm-up days before sampling to harmonise rolling means.
new_rows <- df_agg %>%
  group_by(site_name) %>%
  summarise(min_date = min(date, na.rm = TRUE), .groups = "drop") %>%
  rowwise() %>%
  do(
    data.frame(
      site_name = .$site_name,
      date = seq(.$min_date - days(7), .$min_date - days(1), by = "day")
    )
  ) %>%
  ungroup()

# Bind artificial warm-up rows, add calendar metadata, and average measurements
# within each site/week cell. Keep only rows with population coverage and
# measured viral load to avoid propagating NA series.
df_agg <- bind_rows(df_agg, new_rows) %>%
  arrange(site_name, date) %>%
  mutate(
    KW = lubridate::week(date),
    Jahr = lubridate::year(date)
  ) %>%
  group_by(site_name, th_week) %>%
  mutate_at(vars(contains("viral_load")), ~ mean(., na.rm = TRUE)) %>%
  filter(!is.na(pop_covered), !is.na(log_viral_load)) %>%
  filter(date == max(date, na.rm = TRUE)) %>%
  ungroup()

# Center viral loads per site before aggregation so that site-specific level
# differences do not dominate the combined curve. First subtract the weekly mean
# per th_week, then remove residual site-level offsets, and finally add weekday
# information used for filtering.
df_agg <- df_agg %>%
  group_by(th_week) %>%
  mutate(mean_log_viral_load = mean(log_viral_load)) %>%
  ungroup() %>%
  mutate(log_viral_load_dev = log_viral_load - mean_log_viral_load) %>%
  group_by(site_name, lab) %>%
  mutate(log_viral_load_dev = mean(log_viral_load_dev)) %>%
  ungroup() %>%
  mutate(log_viral_load = log_viral_load - log_viral_load_dev) %>%
  select(-mean_log_viral_load, -log_viral_load_dev) %>%
  mutate(Tag = lubridate::wday(date, week_start = 1))

# Aggregate all sites and selected sites using population coverage weights.
df_agg <- bind_rows(
  aggregation(df = df_agg, weighting = TRUE) %>%
    mutate(all_sites = "ja"),
  aggregation(df = df_agg %>% filter(kept_sites == "ja"), weighting = TRUE) %>%
    mutate(all_sites = "nein")
) %>%
  arrange(all_sites, date)

# Calculate loess estimates with standard errors. The weights
# emphasize observations with lower variance.
pred <- df_agg %>%
  group_by(all_sites) %>%
  nest() %>%
  mutate(pred = purrr::map(data, ~ predict(
    loess.as(
      .x$obs[!is.na(.x$log_viral_load)],
      .x$log_viral_load[!is.na(.x$log_viral_load)],
      degree = 2,
      weights = sqrt(.x$weights[!is.na(.x$log_viral_load)]),
      user.span = .15
    ),
    newdata = data.frame(x = .x$obs),
    se = TRUE
  ))) %>%
  dplyr::select(pred) %>%
  unnest(cols = c(pred))

# Convert nested predictions into plain columns aligned with df_agg rows.
pred_list <- pred[, "pred"]$pred

# Store number of observations per group.
reps <- df_agg %>%
  group_by(all_sites) %>%
  summarise(n = n()) %>%
  pull(n)

# Combine predictions, compute confidence intervals, and back-transform values
# to the original scale.
df_agg <- df_agg %>%
  add_column(
    loess_optimized = extract_prediction(lis = pred_list, extract = "fit"),
    loess_optimized_se = extract_prediction(lis = pred_list, extract = "se.fit"),
    loess_optimized_df = extract_prediction(pred_list, "df") %>%
      map2(., reps, ~ rep(.x, .y)) %>%
      unlist()
  ) %>%
  mutate(
    loess_optimized_pw_lb = loess_optimized - qt(0.975, loess_optimized_df) *
      loess_optimized_se,
    loess_optimized_pw_ub = loess_optimized + qt(0.975, loess_optimized_df) *
      loess_optimized_se
  ) %>%
  add_column(site_name = "Aggregiert") %>%
  dplyr::select(
    site_name,
    all_sites,
    n_non_na,
    date,
    log_viral_load,
    contains("loess_optimized"),
    -contains("d_df")
  ) %>%
  mutate_at(vars(contains("loess")), list(orig = ~ 10 ^ . - 1)) %>%
  mutate(viral_load = 10 ^ log_viral_load - 1)

# Compare average width of confidence bands for all sites vs. selected sites.
df_agg %>%
  mutate(width = loess_optimized_pw_ub_orig - loess_optimized_pw_lb_orig) %>%
  group_by(all_sites) %>%
  summarise(m = mean(width)) %>%
  mutate(ratio =  m[2] /  m[1])

# Concordance between all-site and selected-site curves.
df_ccc <- df_agg %>%
  select(date, all_sites, loess_optimized_orig) %>%
  pivot_wider(
    names_from = all_sites,
    values_from = loess_optimized_orig
  )

# Lin's Concordance Correlation Coefficient (CCC).
ccc <- CCC(df_ccc$ja, df_ccc$nein)
ccc$rho.c
# Pearson correlation coefficient.
cor(df_ccc$ja, df_ccc$nein)

# Peak and trough timing. Reshape to wide format first.
df_wide <- df_agg %>%
  select(date, all_sites, loess_optimized_orig) %>%
  pivot_wider(
    names_from = all_sites,
    values_from = loess_optimized_orig
  ) %>%
  arrange(date)

# Identify turning points for both series.
tp_ja <- find_turning_points(df_wide$ja, df_wide$date) %>%
  mutate(series = "ja")

tp_nein <- find_turning_points(df_wide$nein, df_wide$date) %>%
  mutate(series = "nein")

# Combine turning points and calculate timing differences.
turning_points <- bind_rows(tp_ja, tp_nein)
timing_diff <- turning_points %>% 
  filter(
    # Peaks.
    (type == "peak" &
       (
         (year(date) == 2022 & month(date) == 12) |
           (year(date) == 2023 & month(date) == 3)  |
           (year(date) == 2023 & month(date) == 2)  |
           (year(date) == 2023 & month(date) == 12)
       )) |
      
      # Troughs.
      (type == "trough" &
         (
           (year(date) == 2023 & month(date) == 7) |
             (year(date) == 2023 & month(date) == 6) |
             (year(date) == 2024 & month(date) == 4)
         ))
  ) %>% 
  arrange(type, date) %>%
  group_by(type, series) %>%
  mutate(id = row_number()) %>%
  ungroup() %>%
  tidyr::pivot_wider(
    names_from = series,
    values_from = date
  ) %>%
  mutate(
    diff_days = as.numeric(nein - ja),
    abs_diff_days = abs(diff_days)
  )
timing_diff %>% 
  mutate(mean_abs_diff_days = mean(abs_diff_days))

# Convert viral load units to thousands for presentation.
df_agg <- df_agg %>%
  mutate_at(vars(viral_load, contains("orig")), ~ ./1000)

# Prepare facet labels and locale for month abbreviations.
custom_labels <- c("ja" = "All sites", "nein" = "Sites selected for 2025")
Sys.setlocale("LC_TIME", "C")

# Build plot with ribbons for confidence intervals, points for raw data, and
# lines for loess fits.
p_trans <- df_agg %>%
  ggplot(aes(x = date, y = viral_load)) +
  geom_ribbon(
    aes(ymin = loess_optimized_pw_lb_orig, ymax = loess_optimized_pw_ub_orig),
    fill = "lightblue",
    linemitre = 200
  ) +
  geom_point(colour = "grey") +
  geom_line(aes(date, y = loess_optimized_orig)) +
  facet_wrap(~all_sites, ncol = 1, labeller = labeller(all_sites = custom_labels)) +
  scale_x_date(date_breaks = "2 month", date_labels = "%b \n%y", expand = c(0.05, 0.05)) +
  labs(
    y = expression(atop("Viral load in wastewater", atop(paste("in gene copies / liter (in thousand)")))),
    x = "Date"
  ) +
  theme_palatino

# Display plot.
print(p_trans)

# Save plot.
ggsave(here(results_here, paste0("curves_compared.png")), p_trans, width = 5.25, height = 3.5, dpi = 300)
