#' @export
run_method <- function(task, method, parameters, suppress_output = TRUE) {
  summary <- data.frame(row.names = 1)

  # Add the counts to the parameters
  arglist <- c(list(counts = task$counts), parameters)

  # Include start cell if method requires it
  if ("start_cell_id" %in% formalArgs(method$run_fun)) {
    arglist$start_cell_id <- task$special_cells$start_cell_id
  }

  # Run model on task with given parameters. Suppress output if need be.
  time0 <- Sys.time()
  if (suppress_output) {
    capture.output({
      model <- do.call(method$run_fun, arglist)
    })
  } else {
    model <- do.call(method$run_fun, arglist)
  }
  time1 <- Sys.time()
  summary$time_method <- as.numeric(difftime(time1, time0, units = "sec"))

  rownames(summary) <- NULL

  lst(model, summary)
}
