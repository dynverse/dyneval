#' @import ParamHelpers
#' @import mlr
#' @export
makeRLearner.ti.scorpius <- function() {
  makeRLearnerTI(
    cl = "ti.scorpius",
    package = c("SCORPIUS"),
    par.set = makeParamSet(
      makeIntegerLearnerParam(id = "num_dimensions", lower = 2L, default = 3L, upper = 20L),
      makeIntegerLearnerParam(id = "num_clusters", lower = 2L, default = 4L, upper = 20L),
      makeDiscreteLearnerParam(id = "distance_method", default = "spearman", values = c("spearman", "pearson", "kendall"))
    ),
    properties = c("tibble", "dimred", "dimred_traj", "pseudotime"),
    name = "SCORPIUS",
    short.name = "SCORPIUS"
  )
}

#' @export
trainLearner.ti.scorpius <- function(.learner, .task, .subset, ...) {
  NULL
}

#' @export
predictLearner.ti.scorpius <- function(.learner, .model, .newdata, ...) {
  outs <- lapply(.newdata$counts, run.ti.scorpius, ...)
  outs_tib <- to_tibble(outs)
  outs_tib
}

#' @importFrom SCORPIUS correlation.distance reduce.dimensionality infer.trajectory
#' @importFrom tibble data_frame
run.ti.scorpius <- function(counts, num_dimensions = 3, num_clusters = 4, distance_method = "spearman") {
  expression <- log2(as.matrix(counts)+1)

  dist <- SCORPIUS::correlation.distance(expression, method = distance_method)
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
    dimred_samples = space,
    dimred_traj = traj$path,
    pseudotime = traj$time
  )
}

#' @import ggplot2
#' @importFrom viridis scale_color_viridis
#' @export
plotLearner.ti.scorpius <- function(ti_predictions) {
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
