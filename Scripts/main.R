# Main runner script to execute project steps in order
# 1. Load required packages and functions
# 2. Generate aggregated curves
# 3. Compute quality criteria for exemplary treatment plants
# 4. Select treatment plants according to several criteria and conditions

if (!require("pacman")) install.packages("pacman")
pacman::p_load(here)

message("1) Install and load packages, load functions...")
source(here("Scripts", "setup.R"), encoding = "UTF-8")

message("2) Plotting aggregated curves...")
source(here("Scripts", "plot_curves.R"), encoding = "UTF-8")

message("3) Computing quality criteria...")
source(here("Scripts", "compute_quality_criteria.R"), encoding = "UTF-8")

message("4) Selecting sites...")
source(here("Scripts", "select_sites.R"), encoding = "UTF-8")