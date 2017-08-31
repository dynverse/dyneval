#' Description for random linear
#' @export
description_random_linear <- function() {
  list(
    name = "random_linear",
    short_name = "random_linear",
    package_load = c(),
    package_installed = c(),
    par_set = makeParamSet(
    ),
    properties = c(),
    run_fun = run_random_linear,
    plot_fun = plot_default
  )
}

run_random_linear <- function(counts) {
  milestone_network <- tibble(from=1, to=2, length=1)
  progressions <- tibble(cell_id=rownames(counts), from=1, to=2) %>%
    mutate(percentage=runif(n()))

  wrap_ti_prediction(
    ti_type = "linear",
    id = "random_linear",
    cell_ids = rownames(counts),
    milestone_ids = c(1,2),
    milestone_network = milestone_network,
    progressions = progressions
  )
}

