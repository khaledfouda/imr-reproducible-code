#!/bin/bash
# ============================================================================
# Reproduction Script for:
#   "Incomplete Matrix Regression"
#
# Usage:
#   chmod +x run_all.sh
#   ./run_all.sh
#
# ============================================================================

set -e  # Exit immediately on any error

echo "============================================"
echo " IMR Reproduction Pipeline"
echo " Started at: $(date)"
echo "============================================"

# run configuration functions

source ./helper.sh

# Define the variables
set_simulation_reps 500 1000 1e-6 1e6 5
set_bixi_splits 50
set_movielens_cv FALSE
set_glocalk_epochs 20 30 500 1000


# ---------- 1. Simulation Study ----------
echo ""
echo "========== 1. Simulation Study =========="

echo "  [1/2] Running simulations (1_simulations.R)..."
Rscript ./code/simulation/1_simulations.R

echo "  [2/2] Generating tables and figures (2_tables_plots.R)..."
Rscript ./code/simulation/2_tables_plots.R

echo "  -> Simulation complete. Outputs in output/simulation/"

# ---------- 2. Bixi Application ----------
echo ""
echo "========== 2. Bixi Application =========="

echo "  [1/4] Generating train/test splits (1_generate_train_test_splits.R)..."
Rscript ./code/Bixi/1_generate_train_test_splits.R

echo "  [2/4] Fitting BKTR models (2_fit_BKTR.R)..."
Rscript ./code/Bixi/2_fit_BKTR.R

echo "  [3/4] Fitting IMR models (3_fit_IMR.R)..."
Rscript ./code/Bixi/3_fit_IMR.R

echo "  [4/4] Generating results table (4_generate_results_table.R)..."
Rscript ./code/Bixi/4_generate_results_table.R

echo "  -> Bixi complete. Outputs in output/Bixi/"

# ---------- 3. MovieLens Application ----------
echo ""
echo "========== 3. MovieLens Application =========="

echo "  [1/4] Fitting R-based models (1_fit_MovieLens.R)..."
Rscript ./code/movielens/1_fit_MovieLens.R

echo "  [2/4] Preparing data for GLocal-K (2_1_prepare_python_data.R)..."
Rscript ./code/movielens/2_1_prepare_python_data.R

echo "  [3/4] Fitting GLocal-K PyTorch model (2_2_GlocalK_torch.py)..."
python ./code/movielens/2_2_GlocalK_torch.py

echo "  [4/4] Generating results table and plot (3_generate_results_table.R)..."
Rscript ./code/movielens/3_generate_results_table.R

echo "  -> MovieLens complete. Outputs in output/movielens/"

# ---------- Done ----------
echo ""
echo "============================================"
echo " Finished successfully!"
echo " Completed at: $(date)"
echo ""
echo " All outputs are in the output/ directory."
echo "============================================"
