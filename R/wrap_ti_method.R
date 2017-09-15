#' Return all TI descriptions
#'
#' @param as_tibble Whether or not to return the descriptions as a tibble
#'
#' @importFrom utils lsf.str
#' @importFrom dynutils list_as_tibble
#' @export
get_descriptions <- function(as_tibble = TRUE) {
  requireNamespace("dyneval")
  functions <- lsf.str(asNamespace("dyneval"))
  description_functions <- functions[grep("description_", functions)]
  descriptions <- lapply(description_functions, function(fun_name) {
    do.call(fun_name, args = list(), envir = asNamespace("dyneval"))
  })
  if (as_tibble) {
    list_as_tibble(descriptions)
  } else {
    descriptions
  }
}

#' Check which packages need to be installed for all TI methods to function correctly
#'
#' @export
check_dependencies <- function() {
  for (descr in get_descriptions(as_tibble = FALSE)) {
    required_packages <- c(descr$package_loaded, descr$package_required)
    installed <- required_packages %in% rownames(installed.packages())
    if (any(!installed)) {
      message(sQuote(descr$name), " requires the following packages still to be installed: ", paste(sQuote(required_packages[!installed]), collapse = ", "))
    }
  }
}

create_description <- function(
  name,
  short_name,
  package_loaded,
  package_required,
  par_set,
  properties,
  run_fun,
  plot_fun
) {
  lst(
    name,
    short_name,
    package_loaded,
    package_required,
    par_set,
    properties,
    run_fun,
    plot_fun
  )
}

#' Run a method on a set of tasks with a set of parameters
#' @param tasks The tasks on which to evaluate.
#' @param method The method to evaluate.
#' @param parameters The parameters to evaluate with.
#' @param give_start_cell Whether a start cell should be provided even though a method doesn't require it.
#' @param give_end_cells Whether end cells should be provided even though a method doesn't require it.
#' @param give_cell_grouping Whether a cell grouping should be provided even though a method doesn't require it.
#' @param timeout Kill execution after a given amount of time.
#'
#' @importFrom utils capture.output
#' @importFrom dynutils wait_or_kill
#' @export
execute_method <- function(
  tasks,
  method,
  parameters,
  give_start_cell = FALSE,
  give_end_cells = FALSE,
  give_cell_grouping = FALSE,
  timeout = Inf
) {
  # dry run of the wait_or_kill method
  dry_run <- wait_or_kill({1}, wait_time = 5, function(x) x, .1)
  # Run the method on each of the tasks
  dynutils::wait_or_kill(
    wait_time = timeout,
    cancel_output_fun = function(time) stop("Timeout after ", time, " seconds"),
    check_interval = 1,
    verbose = FALSE,
    globals = c("tasks", "method", "parameters", "give_start_cell", "give_end_cells", "give_cell_grouping"),
    packages = c("dyneval", "dynutils", method$package_loaded),
    expr = {
      # Load required namespaces
      for (pack in method$package_required) {
        suppressMessages(do.call(requireNamespace, list(pack)))
      }

      # Disable seed setting
      orig_setseed <- base::set.seed
      override_setseed(function(i) {})

      # Run method on each task
      outputs <- lapply(seq_len(nrow(tasks)), function(i) {
        task <- extract_row_to_list(tasks, i)

        summary <- data.frame(row.names = 1)

        # Add the counts to the parameters
        arglist <- c(list(counts = task$counts), parameters)

        # Add prior information
        # Including the task is only used for perturbing the gold standard, not used by any real methods
        param_names <- c("start_cell_id", "end_cell_id", "cell_grouping", "task")
        param_bools <- c(give_start_cell, give_end_cells, give_cell_grouping, FALSE)

        for (i in seq_along(param_names)) {
          param_name <- param_names[[i]]
          param_bool <- param_bools[[i]]

          if (param_name %in% formalArgs(method$run_fun) | param_bool) {
            arglist[[param_name]] <-
              if (param_name == "task") {
                task
              } else if (param_name == "cell_grouping") {
                task$cell_grouping
              } else {
                task$special_cells[[param_name]]
              }
          }
        }

        # Run model on task with given parameters. Suppress output if need be.
        time0 <- Sys.time()

        # set working directory to a temporary directory, to avoid polluting the working directory
        # if a method starts writing output to the working directory
        oldwd <- getwd()
        setwd(tempdir())

        tryCatch({
          model <- do.call(method$run_fun, arglist)
        }, finally = {
          # make sure working directory is always set to original even when errored
          setwd(oldwd)
        })
        time1 <- Sys.time()
        summary$time_method <- as.numeric(difftime(time1, time0, units = "sec"))

        rownames(summary) <- NULL

        list(model = model, summary = summary)
      })

      override_setseed(orig_setseed)

      outputs
    }
  )
}
