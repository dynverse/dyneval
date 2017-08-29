#' @export
description_waterfall <- function() {
  list(
    name = "Waterfall",
    short_name = "Waterf", # max 8 chars
    package_load = c("matrixStats", "limma", "MASS", "ape", "RColorBrewer", "RHmm"),
    package_installed = c(),
    par_set = makeParamSet(
      makeIntegerParam(id = "num_clusters", lower = 2, default = 10, upper = 20)
    ),
    properties = c("pseudotime"),
    run_fun = run_waterfall,
    plot_fun = plot_waterfall
  )
}

#' @export
run_waterfall <- function(counts, k = 10) {
  ## load class definition and functions
  code_path <- paste0(path.package("dyneval"), "/extra_code/Waterfall")
  source(paste0(code_path, "/waterfall.R"))

  # Run waterfall
  expr <- log2(counts+1)
  ps <- pseudotimeprog.foo(t(expr), k, color = "red", plot = F)

  # create output
  cell_ids <- rownames(counts)
  milestone_ids <- paste0("milestone_", 1:2)
  milestone_network <- data_frame(from = milestone_ids[[1]], to = milestone_ids[[2]], length = 1)
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

#' @importFrom ggplot2 ggplot geom_point aes scale_colour_distiller theme
plot_waterfall <- function(ti_predictions) {
  ps <- ti_predictions$ps
  ggplot() +
    geom_point(aes(pseudotime, pseudotime.y, colour = pseudotime), ps, size = 5) +
    scale_colour_distiller(palette = "RdBu") +
    theme(legend.position = "none")
}





