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

  # log transform
  expr <- log2(counts+1)

  # run waterfall
  ps <- Waterfall::pseudotimeprog.foo(t(expr), k = num_clusters)

  # return output
  wrap_linear_ti_prediction(
    id = "Waterfall",
    cell_ids = rownames(counts),
    pseudotimes = ps$pseudotime,
    ps = ps
  )
}

plot_waterfall <- function(prediction) {
  requireNamespace("Waterfall")

  Waterfall::plot_waterfall(prediction$ps)
}
