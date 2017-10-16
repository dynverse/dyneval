#' Description for Waterfall
#' @export
description_waterfall <- function() create_description(
  name = "Waterfall",
  short_name = "waterfll", # max 8 chars
  package_loaded = c(),
  package_required = c("Waterfall"),
  par_set = makeParamSet(
    makeIntegerParam(id = "num_clusters", lower = 2L, default = 10L, upper = 20L)
  ),
  properties = c("pseudotime"),
  run_fun = run_waterfall,
  plot_fun = plot_waterfall
)

run_waterfall <- function(counts, num_clusters = 10) {
  requireNamespace("Waterfall")

  # Run waterfall
  expr <- log2(counts+1)
  ps <- Waterfall::pseudotimeprog.foo(t(expr), k = num_clusters)

  # create output
  cell_ids <- rownames(counts)
  wrap_linear_ti_prediction(
    id = "Waterfall",
    cell_ids = cell_ids,
    pseudotimes = ps$pseudotime,
    ps = ps
  )
}

plot_waterfall <- function(ti_predictions) {
  requireNamespace("Waterfall")

  ps <- ti_predictions$ps
  Waterfall::plot_waterfall(ps)
}
