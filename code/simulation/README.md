# Simulation Study Reproduction Guide

This directory contains the necessary R scripts to fully reproduce the simulation study results in section 3 of the manuscript.

## System Requirements & Dependencies

**Hardware Considerations:**
Running the full simulation (`1_simulations.R`) is computationally intensive because it fits multiple models (IMR, Soft-Impute, MCCI) across 500 replications for 4 dimension settings.

**Software & Packages:**
The codebase is written in R. Prior to execution, ensure the following packages are installed:
*   `tidyverse` 
*   `magrittr`
*   `scales`
*   `kableExtra` 
*   `RSSthemes` 
*   `IMR` (authors' package)

## Execution Workflow


### Step 1: Run the Simulations
```r
source("R/simulation/1_simulations.R")
```
*   **Output:** Generates three `.rds` result files (dataframes) saved to `data/simulation/`:
    *   `results_scenario_1.rds` # Setting 1
    *   `results_scenario_2_part_1.rds` #  Setting 2 (figure 1)
    *   `results_scenario_2_part_2.rds` # Setting 2 (figure 2)

### Step 2: Generate Tables and Figures
```r
source("R/simulation/2_tables_plots.R")
```
*   **Output:** 
    *   Prints the LaTeX code for the table (Table 1).
    *   Saves the two figures to `data/simulation/` as `figure_1_scenario_2.pdf` and `figure_2_scenario_2.pdf`.

---
