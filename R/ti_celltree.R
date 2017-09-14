#' Description for celltree maptpx
#' @export
description_celltree_maptpx <- function() create_description(
  name = "cellTree with maptpx",
  short_name = "CTmaptpx",
  package_loaded = c("cellTree"),
  package_required = c(),
  par_set = makeParamSet(
    makeDiscreteParam(id = "method", values = "maptpx", default = "maptpx"),
    makeIntegerParam(id = "num_topics_lower", lower = 2L, upper = 15L, default = 2),
    makeIntegerParam(id = "num_topics_upper", lower = 2L, upper = 15L, default = 15),
    makeNumericParam(id = "sd_filter", lower = log(.01), upper = log(5.0), default = log(.5), special.vals = list(FALSE), trafo = exp),
    makeNumericParam(id = "tot_iter", lower = log(10^4), upper = log(10^7), default = log(10^6), trafo = function(x) round(exp(x))),
    makeNumericParam(id = "tolerance", lower = log(.001), upper = log(.5), default = log(.05), trafo = exp),
    makeNumericParam(id = "width_scale_factor", lower = 1.01, default = 1.2, upper = 2),
    forbidden = quote(num_topics_lower > num_topics_upper)
  ),
  properties = c(),
  run_fun = run_celltree,
  plot_fun = plot_celltree
)

#' Description for celltree gibbs
#' @export
description_celltree_gibbs <- function() create_description(
  name = "cellTree with Gibbs",
  short_name = "CTGibbs",
  package_loaded = c("cellTree"),
  package_required = c(),
  par_set = makeParamSet(
    makeDiscreteParam(id = "method", values = "Gibbs", default = "Gibbs"),
    makeIntegerParam(id = "num_topics", lower = 2L, default = 4L, upper = 15L),
    makeNumericParam(id = "sd_filter", lower = log(.01), upper = log(5.0), default = log(.5), special.vals = list(FALSE), trafo = exp),
    makeNumericParam(id = "tot_iter", lower = log(50), upper = log(500), default = log(200), trafo = function(x) round(exp(x))),
    makeNumericParam(id = "tolerance", lower = log(10^-7), upper = log(10^-3), default = log(10^-5), trafo = exp),
    makeNumericParam(id = "width_scale_factor", lower = log(.1), default = log(1.2), upper = log(100), trafo = exp)
  ),
  properties = c(),
  run_fun = run_celltree,
  plot_fun = plot_celltree
)

#' Description for celltree VEM
#' @export
description_celltree_vem <- function() create_description(
  name = "cellTree with VEM",
  short_name = "CTVEM",
  package_loaded = c("cellTree"),
  package_required = c(),
  par_set = makeParamSet(
    makeDiscreteParam(id = "method", values = "VEM", default = "VEM"),
    makeIntegerParam(id = "num_topics", lower = 2L, default = 4L, upper = 15L),
    makeNumericParam(id = "sd_filter", lower = log(.01), upper = log(5.0), default = log(.5), special.vals = list(FALSE), trafo = exp),
    makeNumericParam(id = "tot_iter", lower = log(10^4), upper = log(10^7), default = log(10^6), trafo = function(x) round(exp(x))),
    makeNumericParam(id = "tolerance", lower = log(10^-7), upper = log(10^-3), default = log(10^-5), trafo = exp),
    makeNumericParam(id = "width_scale_factor", lower = log(.1), default = log(1.5), upper = log(100), trafo = exp)
  ),
  properties = c(),
  run_fun = run_celltree,
  plot_fun = plot_celltree
)

#' @importFrom igraph degree distances get.vertex.attribute induced_subgraph
#' @importFrom reshape2 melt
run_celltree <- function(counts,
                         start_cell_id=NULL,
                         cell_grouping=NULL,
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

  lda_out <- cellTree::compute.lda(
    t(expression) + min(expression) + 1,
    k.topics = num_topics,
    method = method,
    log.scale = FALSE,
    sd.filter = sd_filter,
    tot.iter = tot_iter,
    tol = tolerance)

  # put the parameters for the backbones in separate list, for adding optional cell_grouping and (if grouping is given) start group
  backbone_params <- list(lda_out, width.scale.factor = width_scale_factor, only.mst = FALSE, merge.sequential.backbone = FALSE)

  if(!is.null(cell_grouping)) {
    backbone_params$grouping <- cell_grouping %>% slice(match(cell_id, rownames(counts))) %>% pull(group_id)
    if(!is.null(start_cell_id)) {
      backbone_params$start.group.label <- cell_grouping %>% filter(cell_id == start_cell_id) %>% pull(group_id)
    }
  }

  mst_tree <- do.call(cellTree::compute.backbone.tree, backbone_params)

  backbone_gr <- igraph::induced_subgraph(mst_tree, igraph::get.vertex.attribute(mst_tree, "is.backbone"))
  tomerge <- names(igraph::V(backbone_gr))[igraph::degree(backbone_gr) == 2]
  backbone <- igraph::as_long_data_frame(backbone_gr) %>% select(from = from_name, to = to_name, weight, arrow.mode)

  for (node in tomerge) {
    subgraph <- backbone %>% filter((from == node) | (to == node))
    includeds <- unlist(subgraph$included)
    newnodes <- subgraph %>% {c(.$from, .$to)} %>% keep(~.!=node)

    backbone <- backbone %>% filter((from != node) & (to != node)) %>% bind_rows(list(from=newnodes[[1]], to=newnodes[[2]], included=list(c(includeds, node))))
  }
  backbone <- backbone %>% tidyr::unnest()

  backbonenodes <- names(igraph::V(mst_tree))[igraph::get.vertex.attribute(mst_tree, "is.backbone")]
  sidenodes <- names(igraph::V(mst_tree))[!igraph::get.vertex.attribute(mst_tree, "is.backbone")]
  centralnodes <- unique(c(backbone$from, backbone$to))
  names(centralnodes) <- seq_along(centralnodes)
  sidenodes2backbone <-  igraph::as_long_data_frame(mst_tree) %>% select(from = from_name, to = to_name, weight, arrow.mode) %>% filter(to %in% sidenodes) %$% set_names(from, to)

  percentages <- tibble()
  for (node in names(igraph::V(mst_tree))) {
    if (node %in% names(sidenodes2backbone)) {
      realnode <- as.character(sidenodes2backbone[[node]])
    } else {
      realnode <- node
    }

    if(realnode %in% centralnodes) {
      percentages <- percentages %>% bind_rows(tibble(milestone=as.character(which(centralnodes == realnode)), cell=node, percentage=1))
    } else {
      centralnodesoi <- backbone %>% filter(included == realnode) %>% {c(.$from, .$to)}
      distances <- igraph::distances(mst_tree, realnode, centralnodesoi)
      percentages <- percentages %>% bind_rows(tibble(milestone=as.character(match(centralnodesoi, centralnodes)), cell=node, percentage=1-distances[1, ]/sum(distances)))
    }
  }

  milestone_percentages <- percentages %>% mutate(cell_id = rownames(expression)[as.integer(cell)], milestone_id = paste0("milestone_", milestone)) %>% select(cell_id, milestone_id, percentage)
  milestone_ids <- paste0("milestone_", names(centralnodes))

  # rename milestones
  milestone_network <- backbone %>% tidyr::nest(included) %>% dplyr::select(from, to) %>% mutate(
    from = milestone_ids[match(from, centralnodes)],
    to = milestone_ids[match(to, centralnodes)],
    length = 1,
    directed=TRUE
  )

  # wrap output
  wrap_ti_prediction(
    ti_type = "tree",
    id = "cellTree",
    cell_ids = rownames(counts),
    milestone_ids = milestone_ids,
    milestone_network = milestone_network,
    milestone_percentages = milestone_percentages,
    mst_tree = mst_tree
  )
}

plot_celltree <- function(ti_predictions) {
  requireNamespace("cellTree")
  cellTree::ct.plot.topics(ti_predictions$mst_tree)
}
