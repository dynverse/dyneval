#' @export
makeRLearner.ti.scorpius <- function() {
  makeRLearnerTI(
    cl = "ti.scorpius",
    package = c("SCORPIUS"),
    par.set = makeParamSet(
      makeIntegerLearningParam(id = "num_dimensions", lower = 1L, default = 3L),
      makeIntegerLearningParam(id = "num_clusters", lower = 2L, default = 4L)
    ),
    properties = c("linear", "dimred_samples", "dimred_traj", "pseudotime"), # What to add?
    name = "SCORPIUS",
    short.name = "SCORPIUS"
  )
}

#' @importFrom SCORPIUS correlation.distance reduce.dimensionality infer.trajectory
#' @importFrom tibble data_frame
#' @export
trainLearner.ti.scorpius <- function(.task, .subset, num_dimensions, num_clusters) {
  # subsetting will not work yet, but the function is already provided
  data <- get_task_data(.task, .subset)

  expression <- data[["expression"]]
  dist <- SCORPIUS::correlation.distance(expression)
  space <- SCORPIUS::reduce.dimensionality(dist, ndim = num_dimensions)
  traj <- SCORPIUS::infer.trajectory(space, k = num_clusters)

  state_names <- c("A", "B")
  state_network <- tibble::data_frame(from = "A", to = "B", length = 1)
  state_percentages <- tibble::data_frame(id = rownames(expression), A = 1 - traj$time, B = 1 - A)

  wrap_ti_prediction(
    ti_type = "linear",
    name = "SCORPIUS",
    state_names,
    state_network,
    state_percentages,
    task_id = get_task_identifier(.task),
    dimred_samples = space,
    dimred_traj = traj$path,
    pseudotime = traj$time
  )
}

#' @import ggplot2
#' @importFrom viridis scale_color_viridis
#' @export
plotLearner.ti.scorpius <- function(.task, .subset, ti_predictions) {
  sample_df <- data.frame(
    ti_predictions$dimred_samples,
    pseudotime = ti_predictions$pseudotime
  )
  traj_df <- data.frame(
    ti_predictions$dimred_traj
  )
  ggplot() +
    geom_point(aes(Comp1, Comp2, colour = pseudotime), sample_df) +
    geom_path(aes(Comp1, Comp2), traj_df) +
    coord_equal() +
    scale_color_viridis()
}
