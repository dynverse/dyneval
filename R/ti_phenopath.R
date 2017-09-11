#' Description for phenopath
#' @export
description_phenopath <- function() create_description(
  name = "phenopath",
  short_name = "phenopath",
  package_required = c("phenopath"),
  package_loaded = c(),
  par_set = makeParamSet(
    makeIntegerParam(id = "thin", lower = 2, upper = 500, default = 40),
    makeDiscreteParam(id = "z_init", default = 1, values = list(1, 2, 3, 4, 5, "random")),
    makeLogicalParam(id="model_mu", default=FALSE),
    makeLogicalParam(id="scale_y", default=TRUE)
  ),
  properties = c(),
  run_fun = run_phenopath,
  plot_fun = plot_phenopath
)

run_phenopath <- function(
    counts,
    thin = 40,
    z_init = 1,
    model_mu = FALSE,
    scale_y = TRUE
  ) {
  requireNamespace("phenopath")

  fit <- phenopath::phenopath(counts, rep(1, nrow(counts)), elbo_tol = 1e-6, thin = thin, z_init=z_init, model_mu=model_mu, scale_y=scale_y)
  pseudotimes <- phenopath::trajectory(fit)
  pseudotimes <- (pseudotimes - min(pseudotimes))/(max(pseudotimes) - min(pseudotimes))

  milestone_ids <- c("milestone_A", "milestone_B")
  milestone_network <- tibble::data_frame(from = milestone_ids[[1]], to = milestone_ids[[2]], length = 1)
  progressions <- tibble(
    cell_id = rownames(counts),
    from = milestone_ids[[1]],
    to = milestone_ids[[2]],
    percentage=pseudotimes
  )

  wrap_ti_prediction(
    ti_type = "linear",
    id = "phenopath",
    cell_ids = rownames(counts),
    milestone_ids = milestone_ids,
    milestone_network = milestone_network,
    progressions = progressions,
    fit = fit
  )
}

#' @import ggplot2
plot_phenopath <- plot_default
