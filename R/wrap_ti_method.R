#' Return all TI descriptions
#'
#' @param as_tibble Whether or not to return the descriptions as a tibble
#'
#' @importFrom utils lsf.str
#' @importFrom dynutils list_as_tibble
#' @export
get_descriptions <- function(as_tibble = T) {
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
  for (descr in get_descriptions(as_tibble = F)) {
    required_packages <- c(descr$package_loaded, descr$package_required)
    installed <- required_packages %in% rownames(installed.packages())
    if (any(!installed)) {
      warning(sQuote(descr$name), " requires the following packages still to be installed: ", paste(sQuote(required_packages[!installed]), collapse = ", "))
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
#' @param tasks the tasks on which to evaluate
#' @param method the method to evaluate
#' @param parameters the parameters to evaluate with
#' @param give_start_cell whether a start cell should be provided even though a method doesn't require it
#' @param give_endt_cells whether end cells should be provided even though a method doesn't require it
#' @param give_cell_grouping whether a cell grouping should be provided even though a method doesn't require it
#' @param suppress_output whether or not to suppress the outputted messages
#'
#' @importFrom utils capture.output
#' @export
execute_method <- function(tasks, method, parameters, give_start_cell = FALSE, give_end_cells = FALSE, give_cell_grouping = FALSE, suppress_output = TRUE) {
  # Run the method on each of the tasks
  method_futures <- future(
    {
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
        # Include start cell if method requires it
        if ("start_cell_id" %in% formalArgs(method$run_fun)) {
          # Some methods do not require a start cell, but can use it, determined by whether the default start_cell_id is NULL
          # Give the start cell if the start cell is needed, OR if give_start_cell is TRUE
          if(
            (!is.null(as.list(args(method$run_fun))$start_cell_id)) |
            give_start_cell
          ) {
            arglist$start_cell_id <- task$special_cells$start_cell_id
          }
        }

        # Include end cells if method requires it
        if ("end_cell_ids" %in% formalArgs(method$run_fun)) {
          # Some methods do not require end cells, but can use it, determined by whether the default end_cell_ids is NULL
          # Give the end cells if the end cells are needed, OR if give_end_cells is TRUE
          if(
            (!is.null(as.list(args(method$run_fun))$end_cell_ids)) |
            give_end_cells
          ) {
            arglist$end_cell_ids <- task$special_cells$end_cell_ids
          }
        }

        # Include cell_grouping if method requires it
        if ("cell_grouping" %in% formalArgs(method$run_fun)) {
          # Some methods do not require a grouping, but can use it, determined by whether the default grouping is NULL
          # Give the grouping if the grouping is needed, OR if give_cell_grouping is TRUE
          if(
            (!is.null(as.list(args(method$run_fun))$cell_grouping)) |
            give_cell_grouping
          ) {
            arglist$cell_grouping <- task$cell_grouping
          }
        }

        # Include task (only used for perturbing the gold standard, not used by any real methods)
        if ("task" %in% formalArgs(method$run_fun)) {
          arglist$task <- task
        }

        # Run model on task with given parameters. Suppress output if need be.
        time0 <- Sys.time()

        oldwd <- getwd() # set working directory, to avoid polluting the working directory
        setwd(tempdir())
        tryCatch({
          if (suppress_output) {
            capture.output({
              model <- do.call(method$run_fun, arglist)
            })
          } else {
            model <- do.call(method$run_fun, arglist)
          }
        }, finally=setwd(oldwd)) # make sure working directory is always set to original even when errored
        time1 <- Sys.time()
        summary$time_method <- as.numeric(difftime(time1, time0, units = "sec"))

        rownames(summary) <- NULL

        list(model = model, summary = summary)
      })

      override_setseed(orig_setseed)

      outputs
    },
    globals = c("tasks", "method", "parameters", "suppress_output", "give_start_cell", "give_cell_grouping"),
    packages = c("dyneval", "dynutils", method$package_loaded),
    evaluator = plan("multisession")
  )
  value(method_futures)
}
