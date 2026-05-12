# Overview

This repository contains an R-based pipeline for evaluating wastewater treatment plants, comparing aggregated viral load trends, and selecting sites that meet population coverage thresholds while maintaining measurement quality.

## General information and advice
The analyses were conducted using R 4.5.1 (64-bit, Windows). Recreate the project environment with `renv` by installing the same R version, then running `install.packages("renv")` and `renv::restore()` from the project root.

If you use RStudio, open `site_selection_amelag.Rproj` in the main directory first. If you do not use RStudio, run commands from the project root so the `here` package can resolve the `Data/`, `Results/`, and `Scripts/` directories correctly.

## Repository structure
- **Scripts/** – modular R scripts that implement each stage of the analysis:
  - `setup.R` – shared configuration, dependency loading, helper utilities (adaptive loess smoothing, aggregation, plotting theme, outlier detection).
  - `compute_quality_criteria.R` – computes site-level quality metrics (below-LOQ share, RIIP outlier share, median absolute error), ranks sites, and stores results.
  - `plot_curves.R` – aggregates viral load trajectories across all sites and the selected subset, then saves comparison plots.
  - `select_sites.R` – applies deterministic selection rules per state to meet population coverage targets while prioritizing high-quality sites.
  - `main.R` – orchestrates the full workflow by sourcing the other scripts in sequence.
- **Data/** – input/output data files. The pipeline expects `data.RDS` for full time-series inputs, `data_example.RDS` for exemplary sites, and `quality_data.RDS` for exemplary sites with calculated quality metrics.
- **Results/** – generated output figure (`curves_compared.png`).
- **renv/**, `renv.lock` – dependency lock files for reproducible environments.

## Requirements
- R with the following packages, installed automatically via `pacman::p_load` in `setup.R`: `tidyverse`, `here`, `fANCOVA`, `padr`, `extrafont`, `DescTools`, and `pacman`.
- Fonts: Palatino ("Palatino Linotype") is referenced for plotting; run `extrafont::font_import()` once per machine if needed.

## Running the pipeline
From the project root, run the orchestrator:

```r
Rscript Scripts/main.R
```

This loads shared helpers, generates aggregated curve plots, computes quality criteria, and selects sites.
