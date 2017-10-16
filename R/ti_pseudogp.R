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
    makeIntegerParam(id = "chains", lower = 1L, upper = 20L, default = 1L),
    makeNumericParam(id = "iter", lower = 1, upper = 3, default = 3, trafo = function(x) floor(10^x)),
    makeLogicalVectorParam(id = "dimreds", len = length(list_dimred_methods()), default = names(list_dimred_methods()) == "pca"),
    makeDiscreteParam(id = "initialise_from", values=c("random", "principal_curve", "pca"), default="random")
  ),
  properties = c(),
  run_fun = run_pseudogp,
  plot_fun = plot_pseudogp
)

run_pseudogp <- function(
  counts,
  dimreds = names(list_dimred_methods()) == "pca",
  chains = 1,
  iter = 1000,
  smoothing_alpha = 10,
  smoothing_beta = 3,
  pseudotime_mean = 0.5,
  pseudotime_var = 1,
  initialise_from = "random"
) {
  requireNamespace("pseudogp")
  requireNamespace("rstan")
  requireNamespace("coda")
  requireNamespace("MCMCglmm")

  spaces <- list_dimred_methods()[dimreds] %>% map(~.(counts, 2)) # only 2 dimensions are allowed

  le_fit <- pseudogp::fitPseudotime(
    X = spaces,
    smoothing_alpha = smoothing_alpha,
    smoothing_beta = smoothing_beta,
    iter = iter,
    chains = chains,
    initialise_from = initialise_from,
    pseudotime_var = pseudotime_var,
    pseudotime_mean = pseudotime_mean)

  num_samples <- iter / 2

  pst <- rstan::extract(le_fit, pars = "t", permute = FALSE)
  lambda <- rstan::extract(le_fit, pars = "lambda", permute = FALSE)
  sigma <- rstan::extract(le_fit, pars = "sigma", permute = FALSE)

  # not necessary, but can be used for plotting: the times of every sample across chains
  sample_posterior_times <-
    map_df(seq_len(chains), function(chain_id) {
      mat <- pst[,chain_id,] %>% coda::mcmc() %>% as.matrix()
      colnames(mat) <- rownames(counts)
      mat %>%
        reshape2::melt(varnames = c("sample_id", "cell_id"), value.name = "time") %>%
        as_data_frame() %>%
        mutate(cell_id = as.character(cell_id), chain_id = chain_id)
    })

  # calculate the final pseudotime by averaging over the modes of every chain
  chain_posterior_times <-
    map_df(seq_len(chains), function(chain_id) {
      vec <- pst[,chain_id,] %>% coda::mcmc() %>% MCMCglmm::posterior.mode()
      data_frame(chain_id, cell_id = rownames(counts), time = vec)
    })

  # calculate progressions and milestone network
  progressions <- chain_posterior_times %>%
    group_by(cell_id) %>%
    summarise(percentage = mean(time)) %>%
    mutate(from = "M1", to = "M2") %>%
    select(cell_id, from, to, percentage)

  milestone_network <- tibble(from = "M1", to = "M2", length = 1, directed = TRUE)

  wrap_ti_prediction(
    ti_type = "linear",
    id = "pseudogp",
    cell_ids = rownames(counts),
    milestone_ids = c("M1", "M2"),
    milestone_network = milestone_network,
    progressions = progressions,
    dimreds_samples = spaces,
    sample_posterior_times = sample_posterior_times,
    chain_posterior_times = chain_posterior_times,
    num_samples = num_samples,
    le_fit = le_fit
  )
}

plot_pseudogp <- function(prediction) {
  requireNamespace("pseudogp")

  pseudogp::posteriorCurvePlot(
    X = prediction$dimreds_samples,
    fit = prediction$le_fit,
    nsamples = min(prediction$num_samples, 50),
    posterior_mean = TRUE)
}

