# Parts 2, 3 and 4: tuned fits. Each part is a table of cells (one row per
# fit) and they all use the same worker, tuned_block().
#
# Columns of a table of cells:
#   task      : the cells of one task share their random numbers and are run
#               by one worker; one result file per task.
#   cell_id   : unique name of the cell.
#   seed      : seed of the data set. Within a task the cells have the same
#               seed, so they share every component that their factor does
#               not change (common random numbers).
#   seed_fit  : seed of the validation split and of the fits.
#   the other columns are arguments of generate_simulated_data().

# Part 2, variability by source. For each data set and each source, `draws`
# tuning runs in which only that source is redrawn:
#   split: a new validation split;  noise: a new noise matrix;
#   mask : a new mask at the same rate.
design_sources <- function(cfg) {
  cells <- expand.grid(
    draw = seq_len(cfg$sources$draws), source = cfg$sources$sources,
    seed = cfg$sources$seeds, stringsAsFactors = FALSE
  )
  new_seed <- 10000 * cells$seed + cells$draw
  cells$seed_fit <- ifelse(cells$source == "split", new_seed, cells$seed)
  cells$seed_noise <- ifelse(cells$source == "noise", new_seed, NA)
  cells$seed_mask <- ifelse(cells$source == "mask", new_seed, NA)
  cells$task <- sprintf("seed_%02d_%s", cells$seed, cells$source)
  cells$cell_id <- sprintf("%s_%03d", cells$task, cells$draw)
  return(cells)
}

# Part 3, identification. In each replication: every n, with the generating
# components identified or not.
design_ident <- function(cfg) {
  cells <- expand.grid(
    n = cfg$ident$n, identified = c(TRUE, FALSE),
    seed = seq_len(cfg$ident$reps)
  )
  cells$seed_fit <- cells$seed
  cells$task <- sprintf("rep_%04d", cells$seed)
  cells$cell_id <- sprintf(
    "%s_n%d_%s", cells$task, cells$n, ifelse(cells$identified, "id", "gen")
  )
  return(cells)
}

# Part 4, missingness. In each replication: MCAR (mar_strength = 0) and the
# MAR strengths, on the same theta, noise and uniforms.
design_missing <- function(cfg) {
  cells <- expand.grid(
    mar_strength = cfg$missing$mar_strength,
    seed = seq_len(cfg$missing$reps)
  )
  cells$mechanism <- ifelse(cells$mar_strength == 0, "mcar", "mar")
  cells$seed_fit <- cells$seed
  cells$task <- sprintf("rep_%04d", cells$seed)
  cells$cell_id <- sprintf("%s_mar_%s", cells$task, cells$mar_strength)
  return(cells)
}

# Worker: run the cells of one task, one after the other.
# deciles = TRUE adds the errors by decile of the MAR score (part 4).
tuned_block <- function(cells, cfg, deciles = FALSE) {
  groups <- if (deciles) sens_score_groups(cfg, cells[1, ])
  results <- lapply(seq_len(nrow(cells)), function(i) {
    sens_cell(cells[i, ], cfg, groups)
  })
  return(list(
    rows = dplyr::bind_rows(lapply(results, `[[`, "row")),
    history = dplyr::bind_rows(lapply(results, `[[`, "history"))
  ))
}
