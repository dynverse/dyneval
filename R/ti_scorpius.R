#' @import ParamHelpers
#' @import mlr
#' @export
description_scorpius <- function() {
  list(
    name = "SCORPIUS",
    short_name = "SCORPIUS",
    package_loaded = c("tidyverse", "ggplot2"),
    package_installed = c("SCORPIUS"),
    par_set = makeParamSet(
      makeDiscreteParam(id = "distance_method", default = "spearman", values = c("spearman", "pearson", "kendall")),
      makeIntegerParam(id = "num_dimensions", lower = 2L, default = 3L, upper = 20L),
      makeIntegerParam(id = "num_clusters", lower = 2L, default = 4L, upper = 20L, special.vals = list(NULL)),
      makeNumericParam(id = "thresh", lower = -5L, upper = 5L, default = -3L, trafo = function(x) 10^x),
      makeIntegerParam(id = "maxit", lower = 0, upper = 50, default = 10),
      makeNumericParam(id = "stretch", lower = 0, upper = 5, default = 0),
      makeDiscreteParam(id = "smoother", default = "smooth.spline", values = c("smooth.spline", "lowess", "periodic.lowess"))

    ),
    properties = c("tibble", "dimred", "dimred_traj", "pseudotime"),
    run_fun = run_scorpius,
    plot_fun = plot_scorpius
  )
}

#' @importFrom SCORPIUS correlation.distance reduce.dimensionality infer.trajectory
#' @importFrom tibble data_frame
run_scorpius <- function(counts,
                         num_dimensions = 3, num_clusters = 4, distance_method = "spearman",
                         thresh = .001, maxit = 10, stretch = 0, smoother = "smooth.spline") {
  expression <- log2(as.matrix(counts)+1)

  dist <- SCORPIUS::correlation.distance(expression, method = distance_method)
  space <- SCORPIUS::reduce.dimensionality(dist, ndim = num_dimensions)
  traj <- SCORPIUS::infer.trajectory(space, k = num_clusters, thresh = thresh, maxit = maxit, stretch = stretch, smoother = smoother)

  milestone_ids <- c("milestone_A", "milestone_B")
  milestone_network <- tibble::data_frame(from = milestone_ids[[1]], to = milestone_ids[[2]], length = 1)
  milestone_percentages <- bind_rows(
    tibble::data_frame(id = rownames(expression), state = milestone_ids[[1]], percentage = 1 - traj$time),
    tibble::data_frame(id = rownames(expression), state = milestone_ids[[2]], percentage = traj$time)
  )

  wrap_ti_prediction(
    ti_type = "linear",
    name = "SCORPIUS",
    cell_ids = rownames(counts),
    milestone_ids = milestone_ids,
    milestone_network = milestone_network,
    milestone_percentages = milestone_percentages,
    dimred_samples = space,
    dimred_traj = traj$path,
    pseudotime = traj$time
  )
}

#' @import ggplot2
#' @importFrom viridis scale_color_viridis
#' @export
plot_scorpius <- function(ti_predictions) {
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
