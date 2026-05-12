## Select sites based on quality metrics and population coverage ------------

# -------------------------------------------------------------------------
# This script applies deterministic rules to choose wastewater sites per state
# such that population coverage thresholds are met while prioritising high data
# quality.
# -------------------------------------------------------------------------

# Load shared helpers and paths. Only pacman/here are loaded directly; setup.R
# provides the rest of the shared context.
if (!require("pacman")) install.packages("pacman")
pacman::p_load(here)
source(here("Scripts", "setup.R"), encoding = "UTF-8")

# Read computed quality data from previous step.
df <- readRDS(here(read_data_here, "quality_data.RDS"))

# Focus sites that must always remain in the selection (placeholder list).
focus_sites <- c("XXX")

# Thresholds controlling the optimisation.
threshold    <- 0.15  # 15% population coverage per state
small_cutoff <- 0.02  # 2% minimum population share

# Find the smallest feasible subset that reaches the population threshold,
# breaking ties by the highest average quality rank.
best_subset <- function(idx_set, threshold, pop, rank_mean) {
  idx_set <- unique(idx_set)
  if (length(idx_set) == 0L) return(NULL)
  if (sum(pop[idx_set]) < threshold) return(NULL)

  for (k in seq_along(idx_set)) {
    combs <- combn(idx_set, k)
    if (k == 1L) combs <- matrix(combs, nrow = 1)

    pop_sums   <- colSums(matrix(pop[combs], nrow = k))
    valid_cols <- which(pop_sums >= threshold)

    if (length(valid_cols) > 0L) {
      valid_combs <- combs[, valid_cols, drop = FALSE]
      rm_sums     <- colSums(matrix(rank_mean[valid_combs], nrow = k))
      best_col    <- which.max(rm_sums)
      return(as.integer(valid_combs[, best_col]))
    }
  }
  NULL
}

# Process each state independently and recombine results with a logical keep flag.
df_selected <- df %>%
  group_by(state) %>%
  group_modify(~ {
    state_df <- .x

    pop   <- state_df$pop_share
    rmean <- state_df$rank_mean
    n     <- nrow(state_df)
    all_idx <- seq_len(n)

    focus_idx    <- which(state_df$site_name %in% focus_sites)
    nonfocus_idx <- setdiff(all_idx, focus_idx)

    big_idx   <- nonfocus_idx[pop[nonfocus_idx] >= small_cutoff]
    small_idx <- nonfocus_idx[pop[nonfocus_idx] <  small_cutoff]

    # If total coverage already fails the threshold, keep everything to avoid
    # data loss and flag for manual review.
    if (sum(pop) < threshold) {
      state_df$keep <- TRUE
      return(state_df)
    }

    chosen  <- NULL
    sum_big <- sum(pop[big_idx])

    # Case 1: big non-focus sites alone can reach threshold. Choose the smallest
    # combination that meets coverage while maximising rank_mean.
    if (sum_big >= threshold) {
      chosen <- best_subset(big_idx, threshold, pop, rmean)
      if (is.null(chosen)) {
        state_df$keep <- TRUE
        return(state_df)
      }

    } else {
      # Case 2: big non-focus sites cannot reach threshold. Start with all big
      # sites, then top up with focus and small sites as needed.
      chosen  <- big_idx
      cum_pop <- sum_big
      remaining_threshold <- threshold - cum_pop

      # First try to satisfy threshold with focus sites.
      add_focus <- best_subset(focus_idx, remaining_threshold, pop, rmean)
      if (!is.null(add_focus)) {
        chosen  <- sort(c(chosen, add_focus))
        cum_pop <- sum(pop[chosen])
      }

      # If still below threshold, augment with high-ranked small sites.
      if (cum_pop < threshold) {
        remaining_threshold <- threshold - cum_pop
        add_small <- best_subset(small_idx, remaining_threshold, pop, rmean)
        if (!is.null(add_small)) {
          chosen  <- sort(c(chosen, add_small))
          cum_pop <- sum(pop[chosen])
        }
      }

      # If no combination hits the target, retain all sites for review.
      if (cum_pop < threshold) {
        state_df$keep <- TRUE
        return(state_df)
      }
    }

    # Special rule: if a single site already crosses 15%, also include the
    # second-largest site (where available) to improve robustness.
    if (length(chosen) == 1 && pop[chosen] >= threshold && n >= 2) {
      remaining   <- setdiff(all_idx, chosen)
      second_site <- remaining[order(pop[remaining], decreasing = TRUE)][1]
      chosen      <- sort(c(chosen, second_site))
    }

    # Ensure at least 2 sites per state if possible to provide redundancy.
    if (n >= 2 && length(chosen) < 2) {
      remaining <- setdiff(all_idx, chosen)
      needed    <- 2 - length(chosen)
      add_idx   <- remaining[order(rmean[remaining], decreasing = TRUE)][seq_len(needed)]
      chosen    <- sort(c(chosen, add_idx))
    }

    # Always keep focus sites and flag selected ones.
    chosen <- sort(unique(c(chosen, focus_idx)))
    state_df$keep <- all_idx %in% chosen
    state_df
  }) %>%
  ungroup()

# Show final site selection; keep indicates retained sites.
df_selected %>% 
  select(state, site_name, keep)
