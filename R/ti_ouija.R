#' Description for Ouijia
#' @export
description_ouija <- function() create_description(
  name = "ouija",
  short_name = "ouija",
  package_required = c("ouija", "rstan"),
  package_loaded = c("coda"),
  par_set = makeParamSet(
    makeNumericParam(id = "iter", lower = log(2), default = log(10000), upper = log(50000), trafo = function(x) round(exp(x))),
    makeDiscreteParam(id = "response_type", default = "switch", values = c("switch", "transient")),
    makeDiscreteParam(id = "inference_type", default = "hmc", values = c("hmc", "vb")),
    makeLogicalParam(id = "normalise_expression", default = TRUE)
  ),
  properties = c(),
  run_fun = run_ouija,
  plot_fun = plot_ouija
)

run_ouija <- function(
    counts,
    iter = 1000, # default is actually 10'000.
    response_type = "switch",
    inference_type = "hmc",
    normalise_expression = TRUE
  ) {
  requireNamespace("ouija")
  requireNamespace("rstan")
  requireNamespace("coda")

  # TODO: ouija assumes a small panel of marker genes chosen a priori!
  # marker genes should also be a possible form of prior information
  expr <- log2(counts + 1)

  # write compiled instance of the stanmodel to HDD
  rstan::rstan_options(auto_write = TRUE)

  # run ouija
  oui <- ouija::ouija(
    x = expr,
    iter = iter,
    response_type = response_type,
    inference_type = inference_type,
    normalise_expression = normalise_expression
  )

  # obtain the pseudotimes
  pseudotimes <- ouija::map_pseudotime(oui) %>% dynutils::scale_minmax()

  # run pca for visualisation purposes
  space <- stats::prcomp(expr)$x[,1:2] %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    as_data_frame() %>%
    mutate(pseudotime = pseudotimes)

  # produce output objects
  milestone_ids <- c("milestone_A", "milestone_B")
  milestone_network <- data_frame(
    from = milestone_ids[[1]],
    to = milestone_ids[[2]],
    length = 1,
    directed = FALSE
  )
  progressions <- data_frame(
    cell_id = rownames(counts),
    from = milestone_ids[[1]],
    to = milestone_ids[[2]],
    percentage = pseudotimes
  )

  # extract data for visualisation
  # adapted from ouija::plot_switch_times(oui)
  # to avoid saving the whole oui object
  k_trace <- rstan::extract(oui$fit, "k")$k
  kmean <- colMeans(k_trace)
  t0 <- rstan::extract(oui$fit, "t0")$t0
  t0_means <- colMeans(t0)
  t0_interval <- coda::HPDinterval(coda::mcmc(t0))
  t0_df <- data_frame(t0_mean = t0_means, lower = t0_interval[, 1], upper = t0_interval[, 2], kmean = kmean)
  t0_df$Gene <- colnames(oui$Y[, oui$response_type == "switch"])
  t0_df$Gene <- factor(t0_df$Gene, t0_df$Gene[order(t0_means)])

  # return output
  wrap_ti_prediction(
    ti_type = "linear",
    id = "ouija",
    cell_ids = rownames(counts),
    milestone_ids = milestone_ids,
    milestone_network = milestone_network,
    progressions = progressions,
    t0_df = t0_df,
    space = space
  )
}

#' @importFrom cowplot theme_cowplot
#' @importFrom viridis viridis_pal scale_colour_viridis
plot_ouija <- function(prediction, type = c("dimred", "switch_times")) {
  type <- match.arg(type)

  switch(
    type,
    dimred = {
      ggplot(prediction$space) +
        geom_point(aes(PC1, PC2, colour = pseudotime)) +
        viridis::scale_colour_viridis() +
        cowplot::theme_cowplot()
    },
    switch_times = {
      requireNamespace("ouija")

      t0_df <- prediction$t0_df

      # adapted from ouija::plot_switch_times(prediction$oui)
      # to avoid saving the whole oui file
      vpal <- viridis::viridis_pal()(8)

      ggplot(t0_df, aes(x = Gene, y = t0_mean, fill = kmean)) +
        geom_errorbar(aes(ymin = lower, ymax = upper), color = "grey60", width = 0.5, alpha = 0.5) +
        coord_flip() +
        geom_point(color = "grey50", shape = 21, size = 3) +
        ylab("Switch point") +
        scale_fill_gradient2(name = "Regulation", low = vpal[1], high = vpal[5]) +
        scale_color_gradient2(name = "Regulation", low = vpal[1], high = vpal[5]) +
        theme(legend.position = "top") +
        cowplot::theme_cowplot()
    }
  )
}
