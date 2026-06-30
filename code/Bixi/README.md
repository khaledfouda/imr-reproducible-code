# Bixi Application Reproduction Guide

This directory contains the necessary R scripts to fully reproduce the Bixi bike-sharing data application presented in the manuscript.

##  Dependencies

 The BKTR model fitting step (`2_fit_BKTR.R`) is computationally intensive.
It requires fitting 250 model variations, with each fit taking approximately 22 minutes on an m4 cpu core.
Total estimated time is 4 days.

The number of replications is controlled by `NUM_SPLITS`, defined in `config_default.R`
(full setting: 50) and loaded by `helper.R`. To reduce computation time, lower this
value in `config_default.R` (or create a `config.R` in this folder that sets
`NUM_SPLITS`, which overrides the default). The repository-level `run_all_lite.sh`
script does this automatically for a fast test.


**Packages:**
*   `IMR`
*   `BKTR` (for the benchmark method)
*   `tidyverse` (includes `dplyr`, `tidyr`)
*   `data.table`
*   `bench`
*   `lubridate`
*   `kableExtra` (for rendering the results table)

## Workflow


### Step 1: Generate Train/Test Splits
```r
source("code/Bixi/1_generate_train_test_splits.R")
```
*  Creates missing data scenarios across 50 replications at various training proportions (55% to 75%).
*    Generates 300 `.rds` data split files saved in `data/Bixi/train_test_splits/`.

### Step 2: Fit the BKTR Model 
```r
source("code/Bixi/2_fit_BKTR.R")
```
*    Individual model fits are saved in `output/Bixi/model_fits/BKTR_all_fits/`.
The aggregated results dataframe is saved to `output/Bixi/results_bktr.rds`.

### Step 3: Fit the IMR Models
```r
source("code/Bixi/3_fit_IMR.R")
```
*  Saves the aggregated results table to `output/Bixi/results_imr.rds`.

### Step 4: Generate the Results Table
```r
source("code/Bixi/4_generate_results_table.R")
```
* reads the saved results and generates (prints) the latex code for Table 2. 


---
