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

read_ti_task_data <- function(dataset_name) {
  readRDS(paste0("data/", dataset_name, ".rds"))
}

make_ti_task_data <- function(ti_type, name, expression, state_network, state_percentages, sample_info = NULL, feature_info = NULL) {
  l <- list(
    type = "ti",
    ti_type = ti_type,
    name = name,
    expression = expression,
    state_network = state_network,
    state_percentages = state_percentages,
    sample_info = sample_info,
    feature_info = feature_info
  )
  l
}


get_task_data <- function(.task, .subset) {
  # does not work as it should
  # can we use mlr::getTaskData for this?
  .task
}

wrap_ti_output <- function(state_network, state_percentages) {
  l <- list(state_network = state_network, state_percentages = state_percentages)
  l
}
