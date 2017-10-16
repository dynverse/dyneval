#' Description for phenopath
#' @export
description_phenopath <- function() create_description(
  name = "phenopath",
  short_name = "phenopat",
  package_required = c("phenopath"),
  package_loaded = c(),
  par_set = makeParamSet(
    makeIntegerParam(id = "thin", lower = 2L, upper = 500L, default = 40L),
    makeDiscreteParam(id = "z_init", default = 1, values = list(1, 2, 3, 4, 5, "random")),
    makeLogicalParam(id = "model_mu", default = FALSE),
    makeLogicalParam(id = "scale_y", default = TRUE)
  ),
  properties = c(),
  run_fun = run_phenopath,
  plot_fun = plot_phenopath
)

run_phenopath <- function(counts,
                          thin = 40,
                          z_init = 1,
                          model_mu = FALSE,
                          scale_y = TRUE
) {
  requireNamespace("phenopath")

  expr <- log2(counts + 1)

  # run phenopath
  fit <- phenopath::phenopath(
    exprs_obj = expr,
    x = rep(1, nrow(counts)),
    elbo_tol = 1e-6,
    thin = thin,
    z_init = z_init,
    model_mu = model_mu,
    scale_y = scale_y
  )
  pseudotimes <- phenopath::trajectory(fit) %>% dynutils::scale_minmax()

  # run pca for visualisation purposes
  space <- stats::prcomp(expr)$x[,1:2] %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    as_data_frame() %>%
    mutate(pseudotime = pseudotimes)

  # wrap output
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

  # return output
  wrap_ti_prediction(
    ti_type = "linear",
    id = "phenopath",
    cell_ids = rownames(counts),
    milestone_ids = milestone_ids,
    milestone_network = milestone_network,
    progressions = progressions,
    space = space
  )
}

#' @importFrom viridis scale_colour_viridis
#' @importFrom cowplot theme_cowplot
plot_phenopath <- function(prediction) {
  ggplot(prediction$space) +
    geom_point(aes(PC1, PC2, colour = pseudotime)) +
    viridis::scale_colour_viridis() +
    cowplot::theme_cowplot()
}
