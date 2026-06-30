#!/bin/bash
# Global configuration variables for R and Python scripts
set_simulation_reps() {
  local num_rep=$1
  local maxit=$2
  local thresh=$3
  local max_grid_size=$4
  local mcci_nfolds=$5
  
  cat <<EOF > ./code/simulation/config.R
NUM_REPLICATIONS <- ${num_rep}
CONVERGENCE <- IMR::imr_convergence(maxit=${maxit}, thresh=${thresh})
MAX_GRID_SIZE <- ${max_grid_size}
MCCI_NFOLDS <- ${mcci_nfolds}

EOF
  echo "  Generated ./code/simulation/config.R with NUM_REPLICATIONS=${num_rep}, maxit=${maxit}, and thresh=${thresh}"
}

set_bixi_splits() {
  local num_splits=$1
  
  cat <<EOF > ./code/Bixi/config.R
NUM_SPLITS <- ${num_splits}
EOF
  echo "  Generated ./code/Bixi/config.R with NUM_SPLITS=${num_splits}"
}

set_movielens_cv() {
  local run_cv=$1
  
  cat <<EOF > ./code/movielens/config.R
run_cross_validation <- ${run_cv}
EOF
  echo "  Generated ./code/movielens/config.R with run_cross_validation=${run_cv}"
}

set_glocalk_epochs() {
  local epoch_p=$1
  local epoch_f=$2
  local max_epoch_p=$3
  local max_epoch_f=$4
  
  cat <<EOF > ./code/movielens/config.py
epoch_p = ${epoch_p}
epoch_f = ${epoch_f}
max_epoch_p = ${max_epoch_p}
max_epoch_f = ${max_epoch_f}
EOF
  echo "  Generated ./code/movielens/config.py with epoch_p=${epoch_p}, epoch_f=${epoch_f}, max_epoch_p=${max_epoch_p}, max_epoch_f=${max_epoch_f}"
}