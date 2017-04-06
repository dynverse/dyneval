#' @import mlr
#' @export
#' @rdname RLearner
RLearnerTI <- function(cl, package, par.set, par.vals = list(), properties = character(0L), name = cl, short.name = cl, note = "", callees = character(0L)) {
  addClasses(
    mlr:::makeRLearnerInternal(cl, "ti", package, par.set, par.vals, properties, name, short.name, note, callees),
    c(cl, "RLearnerTI")
  )
}

make_task_data <- function(expression, gold_structure, gold_cells) {
  list(
    expression = expression,
    structure = gold_structure,
    cells = gold_cells
  )
}


get_task_data <- function(.task, .subset) {
  # does not work as it should
  # can we use mlr::getTaskData for this?
  .task
}

wrap_ti_output <- function(structure, cells) {
  l <- list(structure = structure, cells = cells)
  l
}
