# Simulation Study


## Directory Structure

- **Code (`code/simulation/`):** Contains all R scripts required to execute the simulations and generate the tables and figures.
- **Output (`output/simulation/`):** Contains the model outputs (fit files) and generated figures.

## Workflow
- Requirements:
    - R packages: `tidyverse`, `magrittr`, `scales`, `kableExtra`, `RSSthemes`, `IMR`

To fully reproduce the simulation results, execute the scripts in the following order.


The full simulation fits IMR, Soft-Impute, and MCCI across 500 replications for 4
dimension settings. To run fewer, lower `NUM_REPLICATIONS` in `config_default.R`.

### 1. Run the Simulations
- **`1_simulations.R`**: Runs the simulations for all settings. Generates three `.rds` result files (dataframes) saved to `output/simulation/`:
    - `results_scenario_1.rds` # Setting 1
    - `results_scenario_2_part_1.rds` # Setting 2 (figure 1)
    - `results_scenario_2_part_2.rds` # Setting 2 (figure 2)

### 2. Generate Tables and Figures
- **`2_tables_plots.R`**: Reads the simulation output files and generates the tables and plots.
- Prints the LaTeX code for the table (Table 1).
- Saves the two figures to `output/simulation/` as `figure_1_scenario_2.pdf` and `figure_2_scenario_2.pdf`.

## Running the analysis:
Make sure you are in the project root directory and that all requirements are installed.
```bash
Rscript code/simulation/1_simulations.R
Rscript code/simulation/2_tables_plots.R
```
