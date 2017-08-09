#' @import tidyverse
generate_toy_linear <- function(ncells = 100) {
  milestone_network <- tibble(from="M1", to="M2", length=1)
  milestone_ids <- unique(c(milestone_network$from, milestone_network$to))
  cell_ids <- paste0("C", seq_len(ncells))

  progressions <- tibble(cell_id = cell_ids) %>%
    mutate(
      from="M1",
      to="M2",
      percentage=runif(length(.))
    )

  dyneval::wrap_ti_prediction(
    "linear",
    "toy",
    cell_ids,
    milestone_ids,
    milestone_network,
    progressions = progressions
  )
}


change_one_cell <- function(task) {
  the_chosen_one <- sample(task$state_percentages$id, 1)


}

#' @import tidyverse


library(tidyverse)
gs <- generate_toy_linear()
dyneval::plot_default(gs)

datasets <- list(gs=gs)


scores <- dyneval::compute_coranking(gs$geodesic_dist, gs$geodesic_dist)$summary

dummy_method <- function(task) list(
  name = "dummy",
  short_name = "dummy",
  package_load = c(),
  package_installed = c(),
  par_set = ParamHelpers::makeParamSet(
    # makeDiscreteParam(id = "kernel", default = "nn", values = c("nn", "dist", "heat")),
    # makeDiscreteParam(id = "metric", default = "correlation", values = c("correlation", "euclidean", "cosine")),
    # makeNumericParam(id = "nn_pct", lower = -2, upper = log10(10), default = 0, trafo = function(x) 10^x),
    # makeNumericParam(id = "eps", lower = -5L, upper = 5L, default = 0, trafo = function(x) 10^x),
    # makeNumericParam(id = "t", lower = -5L, upper = 5L, default = 0, trafo = function(x) 10^x),
    # makeDiscreteParam(id = "symmetrize", default = "mean", values = c("mean", "ceil", "floor")),
    # makeDiscreteParam(id = "measure_type", default = "unorm", values = c("unorm", "norm")),
    # makeIntegerParam(id = "p", lower = 2, upper = 10, default = 2),
    # makeNumericParam(id = "thresh", lower = -5L, upper = 5L, default = -3L, trafo = function(x) 10^x),
    # makeIntegerParam(id = "maxit", lower = 0, upper = 50, default = 10),
    # makeNumericParam(id = "stretch", lower = 0, upper = 5, default = 2),
    # makeDiscreteParam(id = "smoother", default = "smooth.spline", values = c("smooth.spline", "lowess", "periodic.lowess"))
  ),
  properties = c(),
  run_fun = function(counts) task,
  plot_fun = function(task) plot(1:10)
)


fun <- make_obj_fun(dummy_method(gs))
result <- fun(list(), to_tibble(list(gs)))
scores <- attr(result, "extras")$.summary



