# A small parallel engine. It knows nothing about IMR.
#
# A task is one unit of work for one worker process. Its result is saved in
# its own file, <dir>/<task name>.rds, so an interrupted run resumes: tasks
# that already have a file are skipped.

# Run worker(task, ...) for every task that has no result file yet.
#   tasks    : named list; the names are the file names.
#   worker   : function(task, ...) returning
#              list(rows = <data frame>, history = <data frame>, extra = <list>)
#              (history and extra are optional).
#   n_workers: number of forked R processes.
run_tasks <- function(tasks, worker, dir, n_workers, ...) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  result_file <- function(name) file.path(dir, paste0(name, ".rds"))
  error_file <- function(name) file.path(dir, paste0(name, ".err"))

  pending <- names(tasks)[!file.exists(result_file(names(tasks)))]
  message(sprintf(
    "%s: %d tasks, %d already done, %d to run on %d workers",
    basename(dir), length(tasks), length(tasks) - length(pending),
    length(pending), n_workers
  ))

  run_one <- function(name) {
    start <- Sys.time()
    result <- tryCatch(worker(tasks[[name]], ...), error = function(e) e)
    if (inherits(result, "error")) {
      # no result file is written, so the next run tries this task again
      writeLines(conditionMessage(result), error_file(name))
      return("failed")
    }
    # write under a temporary name first: an interrupted run then leaves no
    # half-written result file
    saveRDS(result, paste0(result_file(name), ".tmp"))
    file.rename(paste0(result_file(name), ".tmp"), result_file(name))
    unlink(error_file(name))
    cat(
      sprintf("%s  %s  %.0f s\n", format(Sys.time(), "%H:%M:%S"), name,
              as.numeric(Sys.time() - start, units = "secs")),
      file = file.path(dir, "progress.log"), append = TRUE
    )
    return("done")
  }

  # mc.preschedule = FALSE: a worker takes the next task when it is free, so
  # tasks of unequal length are balanced
  status <- parallel::mclapply(pending, run_one, mc.cores = n_workers,
                               mc.preschedule = FALSE)
  failed <- pending[!vapply(status, identical, logical(1), "done")]
  if (length(failed) > 0) {
    message(sprintf("%d tasks failed (see the .err files): %s",
                    length(failed), paste(failed, collapse = ", ")))
  }
  return(invisible(failed))
}

# Read every result file of a folder and stack the results.
collect_tasks <- function(dir) {
  files <- list.files(dir, pattern = "\\.rds$", full.names = TRUE)
  results <- lapply(files, readRDS)
  names(results) <- sub("\\.rds$", "", basename(files))
  return(list(
    rows = dplyr::bind_rows(lapply(results, `[[`, "rows")),
    history = dplyr::bind_rows(lapply(results, `[[`, "history")),
    extra = lapply(results, `[[`, "extra")
  ))
}
