#' @import ParamHelpers
#' @import mlr
#' @export
description_embeddr <- function() {
  list(
    name = "embeddr",
    short_name = "embeddr",
    package_load = c("ggplot2"),
    package_installed = c("scater", "embeddr"),
    par_set = makeParamSet(
      makeDiscreteParam(id = "kernel", default = "nn", values = c("nn", "dist", "heat")),
      makeDiscreteParam(id = "metric", default = "correlation", values = c("correlation", "euclidean", "cosine")),
      makeNumericParam(id = "nn_pct", lower = -2, upper = log10(10), default = 0, trafo = function(x) 10^x),
      makeNumericParam(id = "eps", lower = -5L, upper = 5L, default = 0, trafo = function(x) 10^x),
      makeNumericParam(id = "t", lower = -5L, upper = 5L, default = 0, trafo = function(x) 10^x),
      makeDiscreteParam(id = "symmetrize", default = "mean", values = c("mean", "ceil", "floor")),
      makeDiscreteParam(id = "measure_type", default = "unorm", values = c("unorm", "norm")),
      makeIntegerParam(id = "p", lower = 2, upper = 10, default = 2),
      makeNumericParam(id = "thresh", lower = -5L, upper = 5L, default = -3L, trafo = function(x) 10^x),
      makeIntegerParam(id = "maxit", lower = 0, upper = 50, default = 10),
      makeNumericParam(id = "stretch", lower = 0, upper = 5, default = 2),
      makeDiscreteParam(id = "smoother", default = "smooth.spline", values = c("smooth.spline", "lowess", "periodic.lowess"))
    ),
    properties = c(),
    run_fun = run_embeddr,
    plot_fun = plot_embeddr
  )
}

#' @export
run_embeddr <- function(counts,
                        kernel = "nn", metric = "correlation",
                        nn_pct = 1, eps = 1, t = 1,
                        symmetrize = "mean", measure_type = "unorm", p = 2,
                        thresh = .001, maxit = 10, stretch = 2, smoother = "smooth.spline") {
  nn = round(log(nrow(counts)) * nn_pct)

  sce <- scater::newSCESet(countData = t(counts))
  sce <- embeddr::embeddr(sce, kernel = kernel, metric = metric, nn = nn, eps = eps, t = t, symmetrize = symmetrize, measure_type = measure_type, p = p)
  sce <- embeddr::fit_pseudotime(sce, thresh = thresh, maxit = maxit, stretch = stretch, smoother = smoother)

  results <- as(sce@phenoData, "data.frame")
  pseudotime <- results$pseudotime

  dimred_samples <- sce@reducedDimension
  dimred_traj <- results %>% dplyr::arrange(pseudotime) %>% select(starts_with("trajectory"))

  milestone_ids <- c("milestone_A", "milestone_B")
  milestone_network <- tibble::data_frame(from = milestone_ids[[1]], to = milestone_ids[[2]], length = 1)
  milestone_percentages <- bind_rows(
    tibble(
      cell_id = rownames(counts), milestone_id=milestone_ids[[1]], percentage=1-pseudotime
    ),
    tibble(
      cell_id = rownames(counts), milestone_id=milestone_ids[[2]], percentage=pseudotime
    )
  )

  wrap_ti_prediction(
    ti_type = "linear",
    id = "embeddr",
    cell_ids = rownames(counts),
    milestone_ids = milestone_ids,
    milestone_network = milestone_network,
    milestone_percentages = milestone_percentages,
    dimred_samples = dimred_samples,
    dimred_traj = dimred_traj
  )
}

#' @import ggplot2
#' @importFrom viridis scale_color_viridis
#' @export
plot_embeddr <- function(ti_predictions) {
  dimred_samples <- ti_predictions$dimred_samples
  dimred_traj <- ti_predictions$dimred_traj
  ggplot(data.frame(dimred_samples), aes(component_1, component_2)) +
    geom_point(alpha = .65, size = 3.5) +
    geom_path(aes(trajectory_1, trajectory_2), dimred_traj, size = 1.5, alpha = 0.8, linetype = 2)
}
