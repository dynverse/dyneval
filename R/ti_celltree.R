#' Description for celltree maptpx
#' @export
description_celltree_maptpx <- function() abstract_celltree_description("maptpx")

#' Description for celltree gibbs
#' @export
description_celltree_gibbs <- function() abstract_celltree_description("Gibbs")

#' Description for celltree vem
#' @export
description_celltree_vem <- function() abstract_celltree_description("VEM")

abstract_celltree_description <- function(method) {
  if (method == "maptpx") {
    par_set <- makeParamSet(
      makeDiscreteParam(id = "method", values = "maptpx", default = "maptpx"),
      makeIntegerParam(id = "num_topics_lower", lower = 2L, upper = 15L, default = 2L),
      makeIntegerParam(id = "num_topics_upper", lower = 2L, upper = 15L, default = 15L),
      makeNumericParam(id = "sd_filter", lower = log(.01), upper = log(5.0), default = log(.5), special.vals = list(FALSE), trafo = exp),
      makeNumericParam(id = "tot_iter", lower = log(10^4), upper = log(10^7), default = log(10^6), trafo = function(x) as.integer(round(exp(x)))),
      makeNumericParam(id = "tolerance", lower = log(.001), upper = log(.5), default = log(.05), trafo = exp),
      makeNumericParam(id = "width_scale_factor", lower = 1.01, default = 1.2, upper = 2),
      forbidden = quote(num_topics_lower > num_topics_upper)
    )
  } else if (method == "Gibbs") {
    par_set <- makeParamSet(
      makeDiscreteParam(id = "method", values = "Gibbs", default = "Gibbs"),
      makeIntegerParam(id = "num_topics", lower = 2L, default = 4L, upper = 15L),
      makeNumericParam(id = "sd_filter", lower = log(.01), upper = log(5.0), default = log(.5), special.vals = list(FALSE), trafo = exp),
      makeNumericParam(id = "tot_iter", lower = log(50), upper = log(500), default = log(200), trafo = function(x) round(exp(x))),
      makeNumericParam(id = "tolerance", lower = log(10^-7), upper = log(10^-3), default = log(10^-5), trafo = exp),
      makeNumericParam(id = "width_scale_factor", lower = log(.1), default = log(1.2), upper = log(100), trafo = exp)
    )
  } else if (method == "VEM") {
    par_set <- makeParamSet(
      makeDiscreteParam(id = "method", values = "VEM", default = "VEM"),
      makeIntegerParam(id = "num_topics", lower = 2L, default = 4L, upper = 15L),
      makeNumericParam(id = "sd_filter", lower = log(.01), upper = log(5.0), default = log(.5), special.vals = list(FALSE), trafo = exp),
      makeNumericParam(id = "tot_iter", lower = log(10^4), upper = log(10^7), default = log(10^6), trafo = function(x) round(exp(x))),
      makeNumericParam(id = "tolerance", lower = log(10^-7), upper = log(10^-3), default = log(10^-5), trafo = exp),
      makeNumericParam(id = "width_scale_factor", lower = log(.1), default = log(1.5), upper = log(100), trafo = exp)
    )
  }

  run_fun <- run_celltree
  formals(run_fun)$method <- method

  create_description(
    name = glue::glue("cellTree with {method}"),
    short_name = glue::glue("CT{method}"),
    package_loaded = c(),
    package_required = c("cellTree"),
    par_set = par_set,
    properties = c(),
    run_fun = run_fun,
    plot_fun = plot_celltree
  )
}

#' @importFrom igraph degree distances get.vertex.attribute induced_subgraph
#' @importFrom reshape2 melt
run_celltree <- function(counts,
                         start_cell_id = NULL,
                         cell_grouping = NULL,
                         method = "maptpx",
                         num_topics_lower = 2,
                         num_topics_upper = 15,
                         num_topics = num_topics_lower:num_topics_upper,
                         sd_filter = .5,
                         tot_iter = 1e6,
                         tolerance = .05,
                         width_scale_factor = 1.5
) {
  requireNamespace("cellTree")

  expression <- log2(counts+1)

  # infer the LDA model
  lda_out <- cellTree::compute.lda(
    t(expression) + min(expression) + 1,
    k.topics = num_topics,
    method = method,
    log.scale = FALSE,
    sd.filter = sd_filter,
    tot.iter = tot_iter,
    tol = tolerance)

  # put the parameters for the backbones in a list,
  # for adding optional cell_grouping and (if grouping is given) start group
  backbone_params <- list(
    lda.results = lda_out,
    width.scale.factor = width_scale_factor,
    only.mst = FALSE,
    merge.sequential.backbone = FALSE
  )

  # if these parameters are available, add them to the list
  if(!is.null(cell_grouping)) {
    backbone_params$grouping <- cell_grouping %>% slice(match(cell_id, rownames(counts))) %>% pull(group_id)
    if(!is.null(start_cell_id)) {
      backbone_params$start.group.label <- cell_grouping %>% filter(cell_id == start_cell_id) %>% pull(group_id)
    }
  }

  # construct the backbone tree
  mst_tree <- do.call(cellTree::compute.backbone.tree, backbone_params)

  # simplify sample graph to just its backbone
  edges <- igraph::as_data_frame(mst_tree, "edges") %>% select(from, to, length = weight) %>% mutate(directed = FALSE)
  is_trajectory <- igraph::V(mst_tree)$is.backbone %>% setNames(names(igraph::V(mst_tree)))
  out <- dynutils::simplify_sample_graph(edges, is_trajectory, is_directed = FALSE)

  # wrap output
  wrap_ti_prediction(
    ti_type = "tree",
    id = "cellTree",
    cell_ids = rownames(counts),
    milestone_ids = out$milestone_ids,
    milestone_network = out$milestone_network,
    progressions = out$progressions,
    mst_tree = mst_tree
  )
}

plot_celltree <- function(prediction) {
  requireNamespace("cellTree")
  requireNamespace("igraph")
  igraph::plot.igraph(prediction$mst_tree)
}
