#' Description for Ouijia
#' @export
description_ouija <- function() create_description(
  name = "ouija",
  short_name = "ouija",
  package_required = c("ouija"),
  package_loaded = c(),
  par_set = makeParamSet(
    makeIntegerParam(id = "iter", lower = 2, upper = 500, default = 20), # default 10000!
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
    iter = 20,
    response_type = "switch",
    inference_type = "hmc",
    normalise_expression = TRUE
  ) {
  requireNamespace("ouija")

  oui <- ouija::ouija(
    counts,
    iter=iter,
    response_type=response_type,
    inference_type=inference_type,
    normalise_expression=normalise_expression
  )
  pseudotimes <- ouija::map_pseudotime(oui)

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
    id = "ouija",
    cell_ids = rownames(counts),
    milestone_ids = milestone_ids,
    milestone_network = milestone_network,
    progressions = progressions,
    oui = oui
  )
}

#' @import ggplot2
plot_ouija <- function(prediction) {
  ouija::plot_switch_times(prediction$oui)
}
