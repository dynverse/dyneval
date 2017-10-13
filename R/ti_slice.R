#' Description for SLICE
#' @export
description_slice <- function() create_description(
  name = "SLICE",
  short_name = "SLICE",
  package_loaded = c("SLICE"),
  package_required = c("igraph"),
  par_set = makeParamSet(
    makeDiscreteParam(id = "lm.method", default = "clustering", values = c("clustering", "graph")),
    makeDiscreteParam(id = "model.type", default = "tree", values = c("tree", "graph")),
    makeDiscreteParam(id = "ss.method", default = "all", values = c("all", "top", "pcst")),
    makeNumericParam(id = "ss.threshold", default=0.25, lower=0, upper=1),
    makeDiscreteParam(id = "community.method", default = "louvain", values = c("fast_greedy", "edge_betweenness", "label_prop", "leading_eigen","louvain","spinglass", "walktrap", "auto")),
    makeDiscreteParam(id = "cluster.method", default = "kmeans", values = c("kmeans", "pam")),
    makeDiscreteParam(id = "k", default = 0, values = c(0, 3:20)),
    makeIntegerParam(id = "k.max", lower = 3L, upper = 20L, default = 10L),
    makeIntegerParam(id = "B", lower = 3L, upper = 500L, default = 100L),
    makeDiscreteParam(id = "k.opt.method", default = "firstmax", values = c("firstmax", "globalmax", "Tibs2001SEmax", "firstSEmax", "globalSEmax"))
    # makeIntegerParam(id = "B.size", lower = 10L, upper = 10000L, default = 1000L),
    # makeIntegerParam(id = "B.num", lower = 1L, upper = 100L, default = 1L),
    # makeIntegerParam(id = "clustering.k", lower = 1L, upper = 100L, default = 1L),
    # makeDiscreteParam(id = "calculation", default = "bootstrap", values = c("bootstrap", "deterministic"))
  ),
  properties = c(),
  run_fun = run_slice,
  plot_fun = plot_slice
)

run_slice <- function(
  counts,
  cell_grouping = NULL,
  lm.method = "clustering",
  model.type = "tree",
  ss.method = "all",
  ss.threshold = 0.25,
  community.method = "louvain",
  cluster.method = "kmeans",
  k = 0,
  k.max = 10,
  B = 100,
  k.opt.method = "firstmax"
) {
  requireNamespace("igraph")
  requireNamespace("SLICE")

  if (k == 0) k = NULL

  # if cell_grouping is not given, fill it with 1's
  if(!is.null(cell_grouping)) {
    cellidentity <- cell_grouping %>%
      slice(match(rownames(counts), cell_id)) %>%
      pull(group_id) %>%
      factor()
  } else {
    cellidentity <- factor(rep(1, nrow(counts)))
  }

  # wrap data
  sc <- SLICE::construct(
    exprmatrix = as.data.frame(t(counts)),
    cellidentity = cellidentity
  )

  # Should we provide a better km?
  # According to the documentation, km should be:
  # A symmetric matrix encoding the functional similarity of genes;
  # the row names and column names must be official NCBI gene symbols.
  num_genes <- ncol(counts)
  km <- matrix(
    runif(num_genes * num_genes),
    ncol = num_genes,
    dimnames = list(colnames(counts), colnames(counts))
  )

  # calculate the entropy of individual cells
  sc <- SLICE::getEntropy(sc, km = km)

  # reduce expression space
  sc <- SLICE::getRDS(
    sc,
    method = "pca",
    num_dim = 2,
    log.base = 2,
    do.center = TRUE,
    do.scale = FALSE,
    use.cor = TRUE,
    min.var = 0,
    min.cells = 0
  )

  # infer entropy-directed cell lineage model
  sc <- SLICE::getLineageModel(
    sc,
    lm.method = lm.method,
    model.type = model.type,
    ss.method = ss.method,
    ss.threshold = ss.threshold,
    community.method = community.method,
    cluster.method = cluster.method,
    k = k,
    k.max = k.max,
    B = B,
    k.opt.method = k.opt.method,
    do.plot = F
  )

  # extract the stable state to which each cell belongs
  states <- sc@model$cells.df %>%
    rownames_to_column("cell_id") %>%
    slice(match(rownames(counts), cell_id)) %>%
    mutate(state = paste0("slice.ss.", slice.state)) %>%
    select(cell_id, state)

  # each stable state is a milestone
  lin_model <- sc@model$lineageModel
  milestone_ids <- names(igraph::V(lin_model))
  milestone_network <- lin_model %>%
    igraph::as_data_frame() %>%
    rename(length = weight) %>%
    mutate(directed = TRUE)

  # now extract the pseudotimes
  # this is not directly available for us,
  # we will use the method's functions to constuct small trajectories
  # between every stable state (SLICE uses principal curves)
  pseudotimes <- map_df(seq_len(nrow(milestone_network)), function(i) {
    from <- milestone_network[i, 1]
    to <- milestone_network[i, 2]

    fromid <- match(from, milestone_ids)
    toid <- match(to, milestone_ids)

    # has to use pc here, because the other method does not return all cells
    sc_tmp <- SLICE::getTrajectories(
      sc,
      method = "pc",
      start = fromid,
      end = toid,
      do.plot = FALSE,
      do.trim = FALSE
    )
    sc_tmp@transitions[[1]]$i.pseudotime %>%
      rownames_to_column("cell_id") %>%
      mutate(from = from, to = to) %>%
      select(cell_id, from, to, percentage = ptime)
  })

  # some cells can have strange percentages, late in one trajectory, early in the other,
  # even if the state of the cell is neither within the from or to of an edge
  # we therefore filter everything here:
  #  - check whether the state of a cell is in the from or to
  #  - get the earliest timepoint
  progressions <- pseudotimes %>%
    left_join(states, by = "cell_id") %>%
    filter((state == from) | (state == to)) %>%
    group_by(cell_id) %>%
    arrange(percentage) %>%
    filter(row_number() == 1) %>%
    select(-state) %>%
    ungroup()

  wrap_ti_prediction(
    ti_type = "tree",
    id = "SLICE",
    cell_ids = rownames(counts),
    milestone_ids = milestone_ids,
    milestone_network = milestone_network,
    progressions = progressions,
    sc = sc
  )
}

plot_slice <- function(prediction) {
  requireNamespace("igraph")

  sc <- prediction$sc
  list2env(sc@model, environment())

  # Copied from the code of SLICE itself
  # This code is painful to watch
  ss.cells.df <- cells.df

  g <- ggplot() + ggtitle("Inferred Lineage Model") + labs(x="PC1", y="PC2")
  g <- g + geom_point(data=subset(cells.df, slice.realcell==1 & slice.stablestate != "NA" ), aes(x=x, y=y, col=slice.state, size=entropy))
  g <- g + geom_point(data=cells.df[which(cells.df$slice.realcell==0), ], aes(x=x, y=y, size=entropy), col="black")
  g <- g + geom_point(data=subset(cells.df, slice.realcell==1 & slice.stablestate == "NA" ), aes(x=x, y=y, col=slice.state, size=entropy))

  edge.df <- as.data.frame(igraph::get.edgelist(lineageModel))
  edge.df$src.x <- edge.df$src.y <- edge.df$dst.x <- edge.df$dst.y <- 0
  for (ei in 1:dim(edge.df)[1]) {
    src.id <- which(rownames(ss.cells.df) == as.character(edge.df$V1[ei]))
    dst.id <- which(rownames(ss.cells.df) == as.character(edge.df$V2[ei]))
    edge.df$src.x[ei] <- ss.cells.df$x[src.id]
    edge.df$src.y[ei] <- ss.cells.df$y[src.id]
    edge.df$dst.x[ei] <- ss.cells.df$x[dst.id]
    edge.df$dst.y[ei] <- ss.cells.df$y[dst.id]
  }
  g <- g + geom_segment(data=edge.df, aes(x=src.x, y=src.y, xend=dst.x, yend=dst.y), size=I(2), linetype="solid", col=I("black"), alpha=0.6, arrow=arrow(), na.rm=TRUE)

  g <- g +
    theme(strip.background = element_rect(colour = 'white', fill = 'white')) +
    theme(panel.border = element_blank(), axis.line = element_line()) +
    theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) +
    theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) +
    theme(panel.background = element_rect(fill='white')) +
    theme(legend.key = element_blank()) +
    theme(axis.line.x = element_line(size=0.5), axis.line.y=element_line(size=0.5))

  g
}
