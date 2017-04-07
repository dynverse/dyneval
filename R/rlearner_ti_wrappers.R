#' TI R Learner
#'
#' @import mlr
#' @export
#' @rdname RLearner
RLearnerTI <- function(cl, package, par.set, par.vals = list(), properties = character(0L), name = cl, short.name = cl, note = "", callees = character(0L)) {
  addClasses(
    mlr:::makeRLearnerInternal(cl, "ti", package, par.set, par.vals, properties, name, short.name, note, callees),
    c(cl, "RLearnerTI")
  )
}

#' @export
read_task_data <- function(type, ti_type, name, data_directory = "data") {
  readRDS(paste0(data_directory, "/", type, "/", ti_type, "/", name, ".rds"))
}

#' @export
get_task_identifier <- function(task) {
  task[c("type", "ti_type", "name")]
}

#' @export
get_task_data <- function(.task, .subset) {
  # does not work as it should yet
  # can we use mlr::getTaskData for this?
  .task
}

#' @export
wrap_ti_task_data <- function(ti_type, name, expression, state_names, state_network, state_percentages, sample_info = NULL, feature_info = NULL, ...) {
  abstract_wrapper(
    "ti",
    ti_type,
    name,
    state_names,
    state_network,
    state_percentages,
    expression = expression,
    sample_info = sample_info,
    feature_info = feature_info,
    ...
  )
}

#' @export
wrap_ti_prediction <- function(ti_type, name, state_names, state_network, state_percentages, task_id, ...) {
  abstract_wrapper(
    "ti_pred",
    ti_type,
    name,
    state_names,
    state_network,
    state_percentages,
    task_id = task_id,
    ...
  )
}

abstract_wrapper <- function(type, ti_type, name, state_names, state_network, state_percentages, ...) {
  if (!is.data.frame(state_network) || any(colnames(state_network) != c("from", "to", "length"))) {
    stop(sQuote("state_network"), " should be a data frame with exactly three columns named ", sQuote("from"),
         ", ", sQuote("to"), " and ", sQuote("length"), ".")
  }
  if (any(!state_network$from %in% state_names) || any(!state_network$to %in% state_names)) {
    stop("Not all states in ", sQuote("state_network"), " are in ", sQuote("state_names"), ".")
  }
  if (!is.data.frame(state_percentages) || ncol(state_percentages) != length(state_names)+1) {
    stop(sQuote("state_network"), " should be a data frame with exactly N+1 columns, where N is the number of states as defined by ", sQuote("state_names"), ".")
  }
  if (any(!state_names %in% colnames(state_percentages))) {
    stop("Not all states in ", sQuote("state_percentages"), " are in ", sQuote("state_names"), ".")
  }

  l <- list(
    type = type,
    ti_type = ti_type,
    name = name,
    state_names = state_names,
    state_network = state_network,
    state_percentages = state_percentages,
    ...
  )
  class(l) <- paste0("dyneval::ti_wrapper")
  l
}

is_ti_wrapper <- function(object) {
  class(object) == "dyneval::ti_wrapper"
}


