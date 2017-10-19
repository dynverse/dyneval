#' Description for SCORPIUS
#' @export
description_scorpius <- function() create_description(
  name = "SCORPIUS",
  short_name = "SCORPIUS",
  package_loaded = c(),
  package_required = c("SCORPIUS"),
  par_set = makeParamSet(
    makeDiscreteParam(id = "distance_method", default = "spearman", values = c("spearman", "pearson", "kendall")),
    makeIntegerParam(id = "ndim", lower = 2L, default = 3L, upper = 20L),
    makeIntegerParam(id = "k", lower = 0L, default = 4L, upper = 20L),
    makeNumericParam(id = "thresh", lower = -5, upper = 5, default = -3, trafo = function(x) 10^x),
    makeIntegerParam(id = "maxit", lower = 0L, upper = 50L, default = 10L),
    makeNumericParam(id = "stretch", lower = 0, upper = 5, default = 0),
    makeDiscreteParam(id = "smoother", default = "smooth.spline", values = c("smooth.spline", "lowess", "periodic.lowess"))
  ),
  properties = c(),
  run_fun = run_scorpius,
  plot_fun = plot_scorpius
)

run_scorpius <- function(counts,
                         ndim = 3,
                         k = 4,
                         distance_method = "spearman",
                         thresh = .001,
                         maxit = 10,
                         stretch = 0,
                         smoother = "smooth.spline") {
  requireNamespace("SCORPIUS")

  # if k is too low, turn off clustering
  if (k < 2) {
    k <- NULL
  }

  # transform counts
  expr <- log2(as.matrix(counts)+1)

  # calculate distances between cells
  dist <- SCORPIUS::correlation_distance(expr, method = distance_method)

  # perform dimensionality reduction
  space <- SCORPIUS::reduce_dimensionality(dist, ndim = ndim)

  # infer a trajectory through the data
  traj <- SCORPIUS::infer_trajectory(
    space,
    k = k,
    thresh = thresh,
    maxit = maxit,
    stretch = stretch,
    smoother = smoother
  )

  # return output
  wrap_linear_ti_prediction(
    id = "SCORPIUS",
    cell_ids = rownames(counts),
    pseudotimes = traj$time,
    space = space,
    traj = traj
  )
}

#' @importFrom viridis scale_colour_viridis
#' @importFrom reshape2 melt
plot_scorpius <- function(prediction) {
  requireNamespace("SCORPIUS")
  requireNamespace("MASS")

  space <- prediction$space[,1:2]
  ranges <- apply(space, 2, range)
  maxrange <- apply(ranges, 2, diff) %>% max

  limits <- ranges %>% sweep(1, c(-.5, .5) * maxrange, "+")

  space_df <- space %>%
    as.data.frame() %>%
    set_colnames(paste0("Comp", seq_len(ncol(.)))) %>%
    rownames_to_column("cell_id") %>%
    mutate(time = prediction$pseudotimes)

  kde_out <- MASS::kde2d(space[, 1], space[, 2], lims = as.vector(limits))
  z_melt <- reshape2::melt(kde_out$z)
  kde_df <- data.frame(
    Comp1 = kde_out$x[z_melt$Var1],
    Comp2 = kde_out$y[z_melt$Var2],
    density = z_melt$value
  )

  traj_df <- prediction$traj$path[,1:2] %>%
    as.data.frame() %>%
    set_colnames(paste0("Comp", seq_len(ncol(.))))

  g <- ggplot() +
    geom_path(aes(Comp1, Comp2), alpha = 0, data.frame(ranges %>% sweep(1, c(-.2, .2) * maxrange, "+"))) +
    stat_contour(aes(Comp1, Comp2, z = density), geom = "polygon", kde_df, alpha = 0.02, breaks = seq(.2, 1, by = .1)) +
    geom_point(aes(Comp1, Comp2, colour = time), space_df) +
    geom_path(aes(Comp1, Comp2), traj_df) +
    viridis::scale_colour_viridis() +
    labs(colour = "Pseudotime") +
    theme(legend.position = c(.92, .12))

  process_dyneval_plot(g, prediction$id, expand = FALSE)
}
