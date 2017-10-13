#' Description for monocle DDRTree
#' @export
description_monocle2_ddrtree <- function() abstract_monocle_description("DDRTree")

#' Description for monocle ICA
#' @export
description_monocle1_ica <- function() abstract_monocle_description("ICA")

# These don't work yet.
#
# #' Description for monocle SimplePPT
# #' @export
# description_monocle2_simpleppt <- function() abstract_monocle_description("SimplePPT")
#
# #' Description for monocle L1-graph
# #' @export
# description_monocle2_l1graph <- function() abstract_monocle_description("L1-graph")
#
# #' Description for monocle SGL-tree
# #' @export
# description_monocle2_sgltree <- function() abstract_monocle_description("SGL-tree")

abstract_monocle_description <- function(reduction_method) {
  if(reduction_method == "DDRTree") {
    par_set <- makeParamSet(
      makeDiscreteParam(id = "reduction_method", values = "DDRTree", default = "DDRTree"),
      makeIntegerParam(id = "max_components", lower = 2L, default = 2L, upper = 20L),
      makeDiscreteParam(id = "norm_method", default = "vstExprs", values = c("vstExprs", "log", "none")),
      makeLogicalParam(id = "auto_param_selection", default = TRUE)
    )
  } else {
    par_set <- makeParamSet(
      makeDiscreteParam(id = "reduction_method", values = reduction_method, default = reduction_method),
      makeIntegerParam(id = "max_components", lower = 2L, default = 2L, upper = 20L),
      makeDiscreteParam(id = "norm_method", default = "vstExprs", values = c("vstExprs", "log", "none"))
    )
  }

  run_fun <- run_monocle
  formals(run_fun)$reduction_method <- reduction_method

  short_name <- c(
    "DDRTree" = "mnclDDR",
    "ICA" = "mnclICA",
    "tSNE" = "mncltSNE",
    "SimplePPT" = "mnclSPPT",
    "L1-graph" = "mnclL1gr",
    "SGL-tree" = "mnclSGLT"
  )
  create_description(
    name = glue::glue("monocle with {reduction_method}"),
    short_name = short_name[reduction_method],
    package_loaded = c("monocle"),
    package_required = c("BiocGenerics"),
    par_set = par_set,
    properties = c(),
    run_fun = run_fun,
    plot_fun = plot_monocle
  )
}

#' @importFrom igraph as_data_frame
run_monocle <- function(counts,
                        reduction_method,
                        start_cell_id = NULL,
                        max_components = 2,
                        norm_method = "vstExprs",
                        auto_param_selection = TRUE,
                        num_paths = NULL) {
  requireNamespace("monocle")
  requireNamespace("BiocGenerics")

  # TODO: implement num_paths prior

  # just in case
  if (is.factor(norm_method)) norm_method <- as.character(norm_method)

  # load in the new dataset
  pd <- new("AnnotatedDataFrame", data.frame(row.names = rownames(counts)))
  fd <- new("AnnotatedDataFrame", data.frame(row.names = colnames(counts), gene_short_name = colnames(counts)))
  cds <- monocle::newCellDataSet(t(counts), pd, fd)

  # estimate size factors and dispersions
  cds <- BiocGenerics::estimateSizeFactors(cds)
  cds <- BiocGenerics::estimateDispersions(cds)

  # reduce the dimensionality
  cds <- monocle::reduceDimension(
    cds,
    max_components = max_components,
    reduction_method = reduction_method,
    norm_method = norm_method,
    auto_param_selection = auto_param_selection
  )

  # order the cells
  cds <- monocle::orderCells(cds, num_paths = num_paths)

  # extract the igraph and which cells are on the trajectory
  gr <- monocle::minSpanningTree(cds)
  is_trajectory <- setNames(rep(TRUE, nrow(counts)), rownames(counts))

  # convert to milestone representation
  out <- dynutils::simplify_sample_graph(
    edges = igraph::as_data_frame(gr, "edges") %>%
      rename(length = weight) %>%
      mutate(directed = FALSE),
    is_trajectory = is_trajectory,
    is_directed = FALSE
  )

  # wrap output
  wrap_ti_prediction(
    ti_type = "tree",
    id = paste0("monocle_", reduction_method),
    cell_ids = rownames(counts),
    milestone_ids = out$milestone_ids,
    milestone_network = out$milestone_network,
    progressions = out$progressions,
    cds = cds
  )
}

plot_monocle <- function(predictions) {
  requireNamespace("monocle")
  monocle::plot_cell_trajectory(predictions$cds)
}
