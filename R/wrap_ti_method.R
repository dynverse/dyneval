#' Return all TI descriptions
#'
#' @param as_tibble Whether or not to return the descriptions as a tibble
#'
#' @importFrom utils lsf.str
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
    descriptions %>% setNames(descriptions %>% map_chr(~.$short_name))
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

#' Create a TI method description
#'
#' @param name The name of the TI method
#' @param short_name A short name for the method, max 8 characters
#' @param package_loaded The packages that need to be loaded before executing the method
#' @param package_required The packages that need to be installed before executing the method
#' @param par_set A bunch of parameters created by \code{\link[ParamHelpers]{makeParamSet}}
#' @param properties Several descriptive properties of the method. WIP.
#' @param run_fun A function to run the TI, needs to have 'counts' as its first param.
#' @param plot_fun A function to plot the results of a TI, needs to have 'prediction' as its first param.
#' @param override_runfun_params Whether or not to override the default parameters
#' of \code{run_fun} with those described in \code{par_set}.
create_description <- function(
  name,
  short_name,
  package_loaded,
  package_required,
  par_set,
  properties,
  run_fun,
  plot_fun,
  override_runfun_params = TRUE
) {
  if (override_runfun_params) {
    default_params <- par_set %>%
      generateDesignOfDefaults(trafo = TRUE) %>%
      ParamHelpers::dfRowToList(par_set, 1)

    # TODO: check params does not contain 'y' or 'y_1', 'y_2', ... etc

    if(!all(names(default_params) %in% formalArgs(run_fun))) {
      stop("Not all default params described in par_set are listed in the run_fun.")
    }

    formals(run_fun)[names(default_params)] <- default_params
  }
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
#' @param debug Setting debug to \code{TRUE} will avoid running the method in a separate R session
#'   using \code{\link[dynutils]{wait_or_kill}} and run the method directly. Note that the timeout functionality
#'   will not work when \code{debug} is \code{TRUE}.
#'
#' @importFrom utils capture.output
#' @importFrom readr read_file
#' @importFrom stringr str_length
#' @export
execute_method <- function(
  tasks,
  method,
  parameters,
  give_start_cell = FALSE,
  give_end_cells = FALSE,
  give_cell_grouping = FALSE,
  timeout = Inf,
  debug = FALSE
) {
  # Run method on each task
  lapply(seq_len(nrow(tasks)), function(i) {
    # start the timer
    time0 <- Sys.time()

    # get the task
    task <- dynutils::extract_row_to_list(tasks, i)

    # Add the counts to the list parameters
    arglist <- c(list(counts = task$counts), parameters)

    # Add prior information
    # Including the task is only used for perturbing the gold standard, not used by any real methods
    # TODO: This should be uniformised to allow many different outputs
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

    # create a temporary directory to set as working directory,
    # to avoid polluting the working directory if a method starts
    # producing files :angry_face:
    tmp_dir <- tempfile(pattern = method$short_name)
    dir.create(tmp_dir)
    old_wd <- getwd()
    setwd(tmp_dir)

    # disable seed setting
    # a method shouldn't set seeds during regular execution,
    # it should be left up to the user instead
    orig_setseed <- base::set.seed
    setseed_detection_file <- tempfile(pattern = "seedsetcheck")

    # run the method and catch the error, if necessary
    out <-
      tryCatch({
        if (debug) {
          cat("Running ", method$name, " on ", task$name, " in debug mode!\n", sep = "")

          # Load packages into global environment
          for (pack in method$package_loaded) {
            suppressMessages(do.call(require, list(pack)))
          }

          # Execute method
          model <- execute_method_internal(method, arglist, setseed_detection_file)
        } else {
          # dry run of the wait_or_kill method
          dry_run <- wait_or_kill({1}, wait_time = 5, function(x) x, .1)

          # Run the method on each of the tasks
          model <- dynutils::wait_or_kill(
            wait_time = timeout,
            cancel_output_fun = function(time) stop("Timeout after ", time, " seconds"),
            check_interval = 1,
            verbose = FALSE,
            globals = c("method", "arglist", "tmp_dir", "old_wd", "setseed_detection_file"),
            packages = c("dyneval", "dynutils", method$package_loaded),
            expr = {
              tryCatch({
                # create a temporary directory to set as working directory,
                # to avoid polluting the working directory if a method starts
                # producing files :angry_face:
                setwd(tmp_dir)

                # run method
                execute_method_internal(method, arglist, setseed_detection_file)
              }, finally = {
                # return to old_wd
                # (in the likely event that this is a different thread)
                setwd(old_wd)
              })
            }
          )
        }

        c(model, list(error = NULL))
      }, error = function(e) {
        list(model = NULL, time1 = time0, time2 = Sys.time(), error = e)
      })

    # retrieve the model and error message
    model <- out$model
    error <- out$error
    time1 <- out$time1
    time2 <- out$time2

    # check whether the method produced output files and
    # wd to previous state
    num_files_created <- length(list.files(tmp_dir, recursive = TRUE))
    setwd(old_wd)

    # Temporary fix: do not remove the whole tmp wd, as futures might still be in this wd.
    # TODO: solve it
    # unlink(tmp_dir, recursive = TRUE, force = TRUE)
    for (x in list.dirs(tmp_dir, recursive = FALSE)) {
      unlink(x, recursive = TRUE, force = TRUE)
    }
    file.remove(list.files(tmp_dir))

    # read how many seeds were set and
    # restore environment to previous state
    num_setseed_calls <-
      if (file.exists(setseed_detection_file)) {
        stringr::str_length(readr::read_file(setseed_detection_file))
      } else {
        0
      }
    if (file.exists(setseed_detection_file)) {
      file.remove(setseed_detection_file)
    }
    dynutils::override_setseed(orig_setseed)

    # stop the timer
    time3 <- Sys.time()

    # create a summary tibble
    summary <- tibble(
      method_name = method$name,
      method_short_name = method$short_name,
      task_id = task$id,
      time_setup = as.numeric(difftime(time1, time0, units = "sec")),
      time_method = as.numeric(difftime(time2, time1, units = "sec")),
      time_cleanup = as.numeric(difftime(time3, time2, units = "sec")),
      error = list(error),
      num_files_created = num_files_created,
      num_setseed_calls = num_setseed_calls
    )

    # return output
    list(model = model, summary = summary)
  })
}

#' Internal method for executing a method
#'
#' If you're reading this, you're supposed to be using \code{execute_method} instead.
#'
#' @param method A TI method wrapper
#' @param arglist The arguments to apply to the method
#' @param setseed_detection_file A file to which will be written if a method
#'   uses the set.seed function.
#'
#' @export
#' @importFrom readr write_file
execute_method_internal <- function(method, arglist,setseed_detection_file) {
  # disable seed setting
  # a method shouldn't set seeds during regular execution,
  # it should be left up to the user instead
  new_setseed <- function(i) {
    readr::write_file("1", setseed_detection_file, append = TRUE)
  }
  dynutils::override_setseed(new_setseed)

  # Load required namespaces
  for (pack in method$package_required) {
    suppressMessages(do.call(requireNamespace, list(pack)))
  }

  # measure second time point
  time1 <- Sys.time()

  # execute method and return model
  out <- do.call(method$run_fun, arglist)

  # measure third time point
  time2 <- Sys.time()

  list(time1 = time1, time2 = time2, out = out)
}
