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
  milestone_ids <- paste0("milestone_", 1:2)
  milestone_network <- data_frame(from = milestone_ids[[1]], to = milestone_ids[[2]], length = 1, directed=TRUE)
  progressions <- data_frame(cell_id = cell_ids, from = milestone_ids[[1]], to = milestone_ids[[2]], percentage = ps$pseudotime)

  # wrap and return
  wrap_ti_prediction(
    ti_type = "linear",
    id = "Waterfall",
    cell_ids = cell_ids,
    milestone_ids = milestone_ids,
    milestone_network = milestone_network,
    progressions = progressions,
    pseudotime = setNames(cell_ids, ps$pseudotime),
    ps = ps
  )
}

plot_waterfall <- function(ti_predictions) {
  requireNamespace("Waterfall")

  ps <- ti_predictions$ps
  Waterfall::plot_waterfall(ps)
}
