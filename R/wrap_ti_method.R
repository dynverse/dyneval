#' Return all TI descriptions
#'
#' @param as_tibble Whether or not to return the descriptions as a tibble
#'
#' @importFrom utils lsf.str
#' @export
get_descriptions <- function(as_tibble = T) {
  functions <- lsf.str("package:dyneval")
  description_functions <- functions[grep("description_", functions)]
  descriptions <- lapply(description_functions, do.call, list())
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
  descrs <- get_descriptions()
  for (descr_fun in description_functions) {
    descr <- do.call(descr_fun, list())
    required_packages <- c(descr$package_load, descr$package_installed)
    installed <- required_packages %in% rownames(installed.packages())
    if (any(!installed)) {
      warning(sQuote(descr$name), " requires the following packages still to be installed: ", paste(sQuote(required_packages[!installed]), collapse = ", "))
    }
  }
}

#' @importFrom utils capture.output
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

