#' Description for pseudogp
#' @export
description_pseudogp <- function() create_description(
  name = "pseudogp",
  short_name = "pseudogp",
  package_loaded = c("pseudogp"),
  package_required = c("rstan", "coda", "MCMCglmm"),
  par_set = makeParamSet(
    makeNumericParam(id = "smoothing_alpha", lower = 1, upper = 20, default = 10),
    makeNumericParam(id = "smoothing_beta", lower = 1, upper = 20, default = 3),
    makeNumericParam(id = "pseudotime_mean", lower = 0, upper = 1, default = 0.5),
    makeNumericParam(id = "pseudotime_var", lower = 0.01, upper = 1, default = 1),
    makeIntegerParam(id = "chains", lower = 1L, default = 1L, upper = 20L),
    makeNumericParam(id = "iter", lower = log(5), default = log(50), upper = log(1000), trafo = function(x) round(exp(x))), # default is 1000
    makeLogicalVectorParam(id = "dimreds", len = length(list_dimred_methods()), default = names(list_dimred_methods()) == "pca"),
    makeDiscreteParam(id = "initialise_from", values=c("random", "principal_curve", "pca"), default="random")
  ),
  properties = c(),
  run_fun = run_pseudogp,
  plot_fun = plot_pseudogp
)

run_pseudogp <- function(counts,
                         dimreds = names(list_dimred_methods()) == "pca",
                         chains = 1,
                         iter = 1000,
                         smoothing_alpha = 10,
                         smoothing_beta = 3,
                         pseudotime_mean = 0.5,
                         pseudotime_var = 1,
                         initialise_from = "random") {
  requireNamespace("pseudogp")
  requireNamespace("rstan")
  requireNamespace("coda")
  requireNamespace("MCMCglmm")

  # log transform counts
  expr <- log2(counts + 1)

  # perform dimreds
  spaces <- list_dimred_methods()[dimreds] %>%
    map(~.(expr, 2)) # only 2 dimensions per dimred are allowed

  # fit probabilistic pseudotime model
  fit <- pseudogp::fitPseudotime(
    X = spaces,
    smoothing_alpha = smoothing_alpha,
    smoothing_beta = smoothing_beta,
    iter = iter,
    chains = chains,
    initialise_from = initialise_from,
    pseudotime_var = pseudotime_var,
    pseudotime_mean = pseudotime_mean
  )

  # extract pseudotime
  pst <- rstan::extract(fit, pars = "t")$t
  tmcmc <- coda::mcmc(pst)
  pseudotimes <- MCMCglmm::posterior.mode(tmcmc)

  # collect data for visualisation purposes
  # code is adapted from pseudogp::posteriorCurvePlot
  pst <- rstan::extract(fit, pars = "t", permute = FALSE)
  lambda <- rstan::extract(fit, pars = "lambda", permute = FALSE)
  sigma <- rstan::extract(fit, pars = "sigma", permute = FALSE)

  # return output
  wrap_linear_ti_prediction(
    id = "pseudogp",
    cell_ids = rownames(counts),
    pseudotimes = pseudotimes,
    spaces = spaces,
    chains = chains,
    pst = pst,
    lambda = lambda,
    sigma = sigma
  )
}

#' @importFrom cowplot plot_grid
plot_pseudogp <- function(prediction) {
  # code is adapted from pseudogp::posteriorCurvePlot
  requireNamespace("pseudogp")

  spaces <- prediction$spaces
  chains <- prediction$chains
  pst <- prediction$pst
  lambda <- prediction$lambda
  sigma <- prediction$sigma

  Ns <- length(spaces)

  plots <- lapply(seq_len(Ns), function(i) {
    l <- lambda[, , (2 * i - 1):(2 * i), drop = FALSE]
    s <- sigma[, , (2 * i - 1):(2 * i), drop = FALSE]
    sp <- spaces[[i]]
    plt <- pseudogp:::makeEnvelopePlot(
      pst,
      l,
      s,
      sp,
      chains,
      posterior_mean = TRUE,
      ncurves = min(50, nrow(sp)),
      nnt = 80,
      point_colour = "darkred",
      curve_colour = "black",
      point_alpha = 1,
      curve_alpha = .05,
      use_cowplot = TRUE,
      standardize_ranges = FALSE
    )
  })
  cowplot::plot_grid(
    plotlist = plots,
    labels = names(spaces)
  )
}

