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
    if (!is.null(descr$make_command)) {
      pkg_dir <- find.package("dyneval")
      if (pkg_dir == getwd()) {
        pkg_dir <- paste0(pkg_dir, "/inst")
      }
      system(paste0("bash ", pkg_dir, "/", descr$make_command))
    }
  }
}

get_dyneval_install_path <- function() {
  p <- "~/.dyneval"
  if (!dir.exists(p)) dir.create(p, recursive = T, showWarnings = F)
  p
}

create_description <- function(
  name,
  short_name,
  package_loaded,
  package_required,
  par_set,
  properties,
  run_fun,
  plot_fun,
  make_command = NULL
) {
  lst(
    name,
    short_name,
    package_loaded,
    package_required,
    par_set,
    properties,
    run_fun,
    plot_fun,
    make_command
  )
}

#' Run a method on a set of tasks with a set of parameters
#' @param tasks the tasks on which to evaluate
#' @param method the method to evaluate
#' @param parameters the parameters to evaluate with
#' @param suppress_output whether or not to suppress the outputted messages
#'
#' @importFrom utils capture.output
#' @export
execute_method <- function(tasks, method, parameters, suppress_output = TRUE) {
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
          arglist$start_cell_id <- task$special_cells$start_cell_id
        }

        # Include cell_grouping if method requires it
        if ("cell_grouping" %in% formalArgs(method$run_fun)) {
          arglist$cell_grouping <- task$cell_grouping
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

        list(model = model, summary = summary)
      })

      override_setseed(orig_setseed)

      outputs
    },
    globals = c("tasks", "method", "parameters", "suppress_output"),
    packages = c("dyneval", "dynutils", method$package_loaded),
    evaluator = plan("multisession")
  )
  value(method_futures)
}
