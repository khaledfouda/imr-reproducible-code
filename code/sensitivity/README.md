# Sensitivity analysis

Simulations of the sensitivity section: starting values and the monotone
objective (part 1), variability by source (part 2), identification (part 3)
and the missingness mechanism (part 4).

## Run

From the root of this repository, in a terminal:

```sh
OPENBLAS_NUM_THREADS=1 Rscript code/sensitivity/1_run.R lite     # test, minutes
OPENBLAS_NUM_THREADS=1 Rscript code/sensitivity/1_run.R draft    # fewer replications
OPENBLAS_NUM_THREADS=1 Rscript code/sensitivity/1_run.R full     # the whole design
Rscript code/sensitivity/2_summaries.R
```

To read the results, run `code/sensitivity/3_overview.R` in RStudio, section
by section (after `2_summaries.R`, which the figure of part 1 needs).

- Parts can be named after the size: `1_run.R draft missing ident`. The parts
  are `missing`, `ident`, `sources` and `init`.
- `OPENBLAS_NUM_THREADS=1` gives each worker one BLAS thread. Without it, each
  of the 8 workers may start one thread per core. It must be set before R
  starts.
- A run that is stopped continues where it stopped: every task has its own
  result file and finished tasks are skipped. A task that failed leaves a
  `.err` file with the error message and is tried again by the next run.
- `draft` and `full` write to the same folder, and the draft tasks are a
  subset of the full ones: a full run after a draft run only adds the missing
  tasks. Delete `output/sensitivity/` to start again, for example after a
  change to the package.

## Files

| File | Content |
|---|---|
| `settings.R` | `sens_config()`: every number of the study, for the three sizes |
| `runner.R` | `run_tasks()`, `collect_tasks()`: parallel runs with one result file per task |
| `helpers.R` | the data set of a cell, the tuned fit, the metrics |
| `tuned.R` | parts 2, 3 and 4: their tables of cells and their worker |
| `init.R` | part 1: its two stages |
| `1_run.R` | runs the parts and saves the results |
| `2_summaries.R` | the numbers behind each table and figure |
| `3_overview.R` | prints readable tables and draws four figures, to read the results; saves nothing |
| `4_plots.R` | the four figures of the paper (one per case), saved to `output/sensitivity/paper/`; prints the numbers quoted in the text |

A cell is one fit. A task is a group of cells that share their random numbers
(one replication with all its factor levels) and is run by one worker. Every
seed is written in the table of cells, so the results do not depend on the
number of workers.

## Output (`output/sensitivity/`)

- `<part>/<task>.rds`: the result of one task; `<part>/progress.log`.
- `results_<part>.rds`: a list with
  - `rows`: one row per fit (seeds, factor levels, data facts, tuning outputs,
    errors);
  - `history`: the tuning paths;
  - `extra`: for part 1, the fits (coefficients and objective after each
    iteration);
  - `config`, `info`: the settings and the versions of the run.
- `summaries.rds`: the output of `2_summaries.R`.

Main columns of `rows`:

- `status`: `ok`, or `skipped` when a row or column of the draw had fewer
  than two observed entries.
- `rrmse_all`, `rrmse_obs`, `rrmse_unobs`: relative error of the fitted matrix
  on all, observed and unobserved entries.
- `<C>_<estimate>_<target>` for `C` in `beta`, `gamma`, `M`: relative error of
  the `raw` or identified (`id`) estimate from the generating component
  (`gen`) or from its part of the identified truth (`T`).
- `gap_<C>`: distance from the generating component to the identified truth.
- `lambda_m`, `rank_in`, `rank_fit`, `n_iter`, `converged`: tuning outputs.
- Part 4: `rrmse_g01` to `rrmse_g10`, the error by decile of the MAR score
  (1 = least often observed); `pi_min`, `pi_max`.
- Part 1: `regime`, `start`, `level`, `objective`, `certificate`,
  `max_rel_increase`, `n_increase`, and the differences `d_*` from the
  default-start fit of the same data set, regime and threshold.
