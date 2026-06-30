# Output Directory

This directory contains the model outputs, fitted objects, and generated figures for the three analyses in the paper.

## `movielens/`

This directory stores the results and plots for the MovieLens application.

The following files and folders are included in this directory:

| File / Folder | Contents |
|------|----------|
| `plot_full_model.pdf` | Box plot of the fitted ratings across demographic groups |
| `model_fits/` | Directory containing all saved model objects |

## `Bixi/`

This directory stores the results for the Bixi bike-sharing data application. 

The following files and folders are included in this directory:

| File / Folder | Contents |
|------|----------|
| `results_bktr.rds` | Aggregated results table for the BKTR model |
| `results_imr.rds` | Aggregated results table for the IMR models |
| `model_fits/BKTR_all_fits/` | Directory containing BKTR model fits for all train/test splits |

## `simulation/`

This directory stores the generated data, tables, and plots from the simulation study. 

The following files are included in this directory:

| File | Contents |
|------|----------|
| `results_scenario_1.rds` | Results data for Setting 1 |
| `results_scenario_2_part_1.rds` | Results data for Setting 2 (Figure 1) |
| `results_scenario_2_part_2.rds` | Results data for Setting 2 (Figure 2) |
| `figure_1_scenario_2.pdf` | Generated plot for Setting 2 (Figure 1) |
| `figure_2_scenario_2.pdf` | Generated plot for Setting 2 (Figure 2) |
