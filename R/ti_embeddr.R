#' Description for Embeddr
#' @export
description_embeddr <- function() create_description(
  name = "embeddr",
  short_name = "embeddr",
  package_loaded = c(),
  package_required = c("scater", "embeddr"),
  par_set = makeParamSet(
    makeDiscreteParam(id = "kernel", default = "nn", values = c("nn", "dist", "heat")),
    makeDiscreteParam(id = "metric", default = "correlation", values = c("correlation", "euclidean", "cosine")),
    makeNumericParam(id = "nn_pct", lower = -2, upper = log10(10), default = 0, trafo = function(x) 10^x),
    makeNumericParam(id = "eps", lower = -5L, upper = 5L, default = 0, trafo = function(x) 10^x),
    makeNumericParam(id = "t", lower = -5L, upper = 5L, default = 0, trafo = function(x) 10^x),
    makeDiscreteParam(id = "symmetrize", default = "mean", values = c("mean", "ceil", "floor")),
    makeDiscreteParam(id = "measure_type", default = "unorm", values = c("unorm", "norm")),
    makeIntegerParam(id = "p", lower = 2L, upper = 10L, default = 2L),
    makeNumericParam(id = "thresh", lower = -5L, upper = 5L, default = -3L, trafo = function(x) 10^x),
    makeIntegerParam(id = "maxit", lower = 0L, upper = 50L, default = 10L),
    makeNumericParam(id = "stretch", lower = 0, upper = 5, default = 2),
    makeDiscreteParam(id = "smoother", default = "smooth.spline", values = c("smooth.spline", "lowess", "periodic.lowess"))
  ),
  properties = c(),
  run_fun = run_embeddr,
  plot_fun = plot_embeddr
)

run_embeddr <- function(counts,
                        kernel = "nn",
                        metric = "correlation",
                        nn_pct = 1,
                        eps = 1,
                        t = 1,
                        symmetrize = "mean",
                        measure_type = "unorm",
                        p = 2,
                        thresh = .001,
                        maxit = 10,
                        stretch = 2,
                        smoother = "smooth.spline") {
  requireNamespace("scater")
  requireNamespace("embeddr")

  # calculate nn param
  nn <- round(log(nrow(counts)) * nn_pct)

  # load data in scater
  sce <- scater::newSCESet(countData = t(counts))

  # run embeddr
  sce <- embeddr::embeddr(
    sce,
    kernel = kernel,
    metric = metric,
    nn = nn,
    eps = eps,
    t = t,
    symmetrize = symmetrize,
    measure_type = measure_type,
    p = p
  )

  # fit pseudotime
  sce <- embeddr::fit_pseudotime(
    sce,
    thresh = thresh,
    maxit = maxit,
    stretch = stretch,
    smoother = smoother
  )

  # construct milestone network
  milestone_ids <- c("milestone_A", "milestone_B")

  milestone_network <- data_frame(
    from = milestone_ids[[1]],
    to = milestone_ids[[2]],
    length = 1,
    directed = FALSE
  )

  # construct progressions
  progressions <- as(sce@phenoData, "data.frame") %>%
    rownames_to_column(var = "cell_id") %>%
    mutate(from = milestone_ids[[1]], to = milestone_ids[[2]]) %>%
    select(cell_id, from, to, percentage = pseudotime)

  # creating extra output for visualisation purposes
  dimred_samples <- sce@reducedDimension %>%
    as.data.frame() %>%
    rownames_to_column("cell_id")
  dimred_traj <- as(sce@phenoData, "data.frame") %>%
    arrange(pseudotime) %>%
    select(pseudotime, starts_with("trajectory_"))

  # return output
  wrap_ti_prediction(
    ti_type = "linear",
    id = "embeddr",
    cell_ids = rownames(counts),
    milestone_ids = milestone_ids,
    milestone_network = milestone_network,
    progressions = progressions,
    dimred_samples = dimred_samples,
    dimred_traj = dimred_traj
  )
}

plot_embeddr <- function(predictions) {
  ggplot() +
    geom_point(aes(component_1, component_2), predictions$dimred_samples, alpha = .65, size = 3.5) +
    geom_path(aes(trajectory_1, trajectory_2), predictions$dimred_traj, size = 1.5, alpha = 0.8, linetype = 2) +
    cowplot::theme_cowplot()
}
