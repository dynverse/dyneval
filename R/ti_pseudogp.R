#' Description for pseudogp
#' @export
description_pseudogp <- function() {
  list(
    name = "pseudogp",
    short_name = "pseudogp",
    package_load = c(),
    package_installed = c("pseudogp", "rstan"),
    par_set = makeParamSet(
      makeNumericParam(id = "smoothing_alpha", lower = 1, upper = 20, default = 10),
      makeNumericParam(id = "smoothing_beta", lower = 1, upper = 20, default = 3),
      makeNumericParam(id = "pseudotime_mean", lower = 0, upper = 1, default = 0.5),
      makeNumericParam(id = "pseudotime_var", lower = 0.01, upper = 1, default = 1),
      makeIntegerParam(id = "chains", lower = 1L, upper = 20L, default = 1L),
      makeIntegerParam(id = "iter", lower = 10L, upper = 1000L, default = 1000L),
      makeDiscreteParam(id = "dimred_names", values=c("pca", "mds", "tsne"), default="pca"),
      makeDiscreteParam(id = "intialise_from", values=c("random", "principal_curve", "pca"), default="random")
    ),
    properties = c(),
    run_fun = run_pseudogp,
    plot_fun = plot_pseudogp
  )
}

run_pseudogp <- function(
  counts,
  dimred_names="pca",
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
  dimred_funcs <- map(dimred_names, ~getFromNamespace(paste0("dimred_", .), "dyneval"))

  spaces <- map(dimred_funcs, ~.(counts, 2)) # only 2 dimensions are allowed
  #ggplot(spaces[[1]] %>% as.data.frame) + geom_point(aes(Comp1, Comp2))

  le_fit <- pseudogp::fitPseudotime(spaces, smoothing_alpha, smoothing_beta, iter = iter, chains = chains, initialise_from = initialise_from, pseudotime_var=pseudotime_var, pseudotime_mean=pseudotime_mean)
  nsamples <- iter/2
  #posteriorCurvePlot(spaces, le_fit, nsamples=nsamples, posterior_mean = TRUE)

  pst <- rstan::extract(le_fit, pars = "t", permute = FALSE)
  lambda <- rstan::extract(le_fit, pars = "lambda", permute = FALSE)
  sigma <- rstan::extract(le_fit, pars = "sigma", permute = FALSE)

  # not necessary, but can be used for plotting: the times of every sample across chains
  sample_posterior_times <- map(seq_len(chains), ~mcmc(pst[, ., ])) %>%
    invoke(rbind, .)
  colnames(sample_posterior_times) <- rownames(counts)
  sample_posterior_times <- sample_posterior_times %>%
    reshape2::melt(varnames=c("sample_id", "cell_id"), value.name="time") %>%
    mutate(chain_id = floor((sample_id-1)/(max(sample_id) / chains)) + 1)

  # calculate the final pseudotime by averaging over the modes of every chain
  chain_posterior_times <- map(seq_len(chains), ~posterior.mode(mcmc(pst[, ., ]))) %>%
    invoke(rbind, .)
  colnames(chain_posterior_times) <- rownames(counts)
  chain_posterior_times <- chain_posterior_times %>% reshape2::melt(varnames = c("chain_id", "cell_id"), value.name="time")

  # calculate progressions and milestone network
  progressions <- chain_posterior_times %>%
    group_by(cell_id) %>%
    summarise(percentage=mean(time)) %>%
    mutate(from="M1", to="M2") %>%
    mutate(cell_id = as.character(cell_id))

  milestone_network <- tibble(from="M1", to="M2", length=1)

  wrap_ti_prediction(
    ti_type = "linear",
    id = "pseudogp",
    cell_ids = rownames(counts),
    milestone_ids = c("M1", "M2"),
    milestone_network = milestone_network %>% select(from, to, length),
    progressions = progressions %>% select(cell_id, from, to, percentage),
    dimreds_samples = spaces,
    sample_posterior_times = sample_posterior_times,
    chain_posterior_times = chain_posterior_times,
    le_fit = le_fit
  )
}

plot_pseudogp <- function(prediction) {
  # variability between samples and chains
  prediction$sample_posterior_times %>%
    group_by(cell_id) %>%
    mutate(mean_time=mean(time)) %>%
    ggplot() +
    geom_boxplot(aes(mean_time, time, group=cell_id, color=factor(chain_id)), width=0.01) +
    facet_wrap(~chain_id)

  # variability between chains, eg. to check whether there is a structure present in the data
  prediction$chain_posterior_times %>%
    group_by(cell_id) %>%
    mutate(mean_time=mean(time)) %>%
    ggplot() +
    geom_boxplot(aes(mean_time, time, group=cell_id))

  posteriorCurvePlot(prediction$dimreds_samples, prediction$le_fit, nsamples=50, posterior_mean = TRUE)
}
