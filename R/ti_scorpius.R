makeRLearner.ti.scorpius <- function() {
  makeRLearnerTI(
    cl = "ti.scorpius",
    package = c("SCORPIUS"),
    par.set = makeParamSet(
      makeIntegerLearningParam(id = "num_dimensions", lower = 1L, default = 3L),
      makeIntegerLearningParam(id = "num_clusters", lower = 2L, default = 4L)
    )
  )
}

trainLearner.ti.scorpius = function(.task, .subset, num_dimensions, num_clusters) {
  data <- get_task_data(.task, .subset) # this will not work yet, but is based on mlr

  expression <- data[["expression"]]
  dist <- SCORPIUS::correlation.distance(expression)
  space <- SCORPIUS::reduce.dimensionality(dist, ndim = num_dimensions)
  traj <- SCORPIUS::infer.trajectory(space, k = num_clusters)

  state_network <- tibble::data_frame(from = "A", to = "B", length = 1)
  state_percentages <- tibble::data_frame(cellid = rownames(expression), A = 1 - traj$time, B = 1 - A)

  wrap_ti_output(state_network, state_percentages)
}
