#' @export
makeRLearner.ti.celltree <- function() {
  makeRLearnerTI(
    cl = "ti.celltree",
    package = c("cellTree"),
    par.set = makeParamSet(
      makeIntegerLearningParam(id = "num_topics", lower = 3L, default = 10L)
      # makeIntegerLearningParam(id = "num_clusters", lower = 2L, default = 4L)
    ),
    # properties = c("linear", "dimred_samples", "dimred_traj", "pseudotime"), # What to add?
    properties = c(),
    name = "celltree",
    short.name = "celltree"
  )
}


#' @importFrom igraph degree distances get.vertex.attribute induced_subgraph
#' @importFrom reshape2 melt
#' @import dplyr
#' @import cellTree
#' @import magrittr
#'
#' @export
trainLearner.ti.celltree <- function(.task, .subset, num_topics, width_scale_factor) {
  # subsetting will not work yet, but the function is already provided
  data <- get_task_data(.task, .subset)
  expression = data$expression

  #expression <- t(SCORPIUS::quant.scale(t(data$expression), 0))

  #num_topics = 4;width_scale_factor=1.5

  lda.results <- compute.lda(t(expression) + min(expression) + 1, k.topics=num_topics, method="maptpx", log.scale=F, sd.filter=0)
  dists <- get.cell.dists(lda.results)
  #mst.tree <- compute.backbone.tree(lda.results, grouping = goldstandard$cellinfo$piecestateid, width.scale.factor=1.5, start.group.label = "1")
  #mst.tree <- compute.backbone.tree(lda.results, only.mst = T)
  mst.tree <- compute.backbone.tree(lda.results, width.scale.factor=width_scale_factor)
  #ct.plot.topics(mst.tree)
  #ct.plot.grouping(mst.tree)


  backbone <- induced_subgraph(mst.tree, get.vertex.attribute(mst.tree, "is.backbone"))
  tomerge <- names(V(backbone))[igraph::degree(backbone) == 2]
  backbone <- backbone %>% as_data_frame
  for(node in tomerge) {
    subgraph <- backbone %>% filter((from == node) | (to == node))
    includeds <- unlist(subgraph$included)
    newnodes <- subgraph %>% {c(.$from, .$to)} %>% keep(~.!=node)

    backbone <- backbone %>% filter((from != node) & (to != node)) %>% bind_rows(list(from=newnodes[[1]], to=newnodes[[2]], included=list(c(includeds, node))))
  }
  backbone <- backbone %>% tidyr::unnest()

  backbonenodes <- names(V(mst.tree))[get.vertex.attribute(mst.tree, "is.backbone")]
  sidenodes <- names(V(mst.tree))[!get.vertex.attribute(mst.tree, "is.backbone")]
  centralnodes <- unique(c(backbone$from, backbone$to))
  names(centralnodes) <- seq_along(centralnodes)
  sidenodes2backbone <- mst.tree %>% as_data_frame %>% filter(to %in% sidenodes) %$% set_names(from, to)

  percentages <- tibble()
  for(node in names(V(mst.tree))) {
    if (node %in% names(sidenodes2backbone)) {
      realnode <- as.character(sidenodes2backbone[[node]])
    } else {
      realnode <- node
    }

    if(realnode %in% centralnodes) {
      percentages <- percentages %>% bind_rows(tibble(milestone=as.character(which(centralnodes == realnode)), cell=node, percentage=1))
    } else {
      centralnodesoi <- backbone %>% filter(included == realnode) %>% {c(.$from, .$to)}
      distances <- igraph::distances(mst.tree, realnode, centralnodesoi)
      percentages <- percentages %>% bind_rows(tibble(milestone=as.character(match(centralnodesoi, centralnodes)), cell=node, percentage=1-distances[1, ]/sum(distances)))
    }
  }
  state_percentages <- reshape2::acast(percentages %>% mutate(cell=as.numeric(cell)) %>% arrange(cell), cell~milestone, value.var="percentage", fill=0) %>%
    as.data.frame() %>% mutate(id=rownames(expression))
  state_names <- names(centralnodes)

  state_percentages <- state_percentages[,c("id", state_names)]

  # rename states
  state_network <- backbone %>% tidyr::nest(included) %>% dplyr::select(from, to) %>% mutate(
    from = names(centralnodes)[match(from, centralnodes)],
    to = names(centralnodes)[match(to, centralnodes)]
  )
  state_network <- state_network %>% mutate(from = paste0("state_", from), to = paste0("state_", to), length=1)
  state_names <- paste0("state_", state_names)
  colnames(state_percentages) <- c("id", state_names)

  # wrap output
  wrap_ti_prediction(
    ti_type = "tree",
    name = "celltree",
    state_names,
    state_network,
    state_percentages,
    task_id = get_task_identifier(.task),
    mst.tree = mst.tree
  )
}

#' @export
plotLearner.ti.celltree <- function(ti_predictions) {
  cellTree::ct.plot.grouping(ti_predictions$mst.tree)
  #monocle::plot_cell_trajectory(ti_predictions$cds)
}
