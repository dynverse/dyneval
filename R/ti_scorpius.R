#' Description for SCORPIUS
#' @export
description_scorpius <- function() create_description(
  name = "SCORPIUS",
  short_name = "SCORPIUS",
  package_loaded = c(),
  package_required = c("SCORPIUS"),
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

run_scorpius <- function(counts,
                         num_dimensions = 3, num_clusters = 4, distance_method = "spearman",
                         thresh = .001, maxit = 10, stretch = 0, smoother = "smooth.spline") {
  requireNamespace("SCORPIUS")

  expression <- log2(as.matrix(counts)+1)

  dist <- SCORPIUS::correlation_distance(expression, method = distance_method)
  space <- SCORPIUS::reduce_dimensionality(dist, ndim = num_dimensions)
  traj <- SCORPIUS::infer_trajectory(space, k = num_clusters, thresh = thresh, maxit = maxit, stretch = stretch, smoother = smoother)

  milestone_ids <- c("milestone_A", "milestone_B")
  milestone_network <- data_frame(from = milestone_ids[[1]], to = milestone_ids[[2]], length = 1, directed=TRUE)
  milestone_percentages <- bind_rows(
    data_frame(cell_id = rownames(expression), milestone_id = milestone_ids[[1]], percentage = 1 - traj$time),
    data_frame(cell_id = rownames(expression), milestone_id = milestone_ids[[2]], percentage = traj$time)
  )

  wrap_ti_prediction(
    ti_type = "linear",
    id = "SCORPIUS",
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
