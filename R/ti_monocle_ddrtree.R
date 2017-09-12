#' Description for Monocle DDRtree
#' @export
description_monocle_ddrtree <- function() create_description(
  name = "monocle with DDRtree",
  short_name = "monocDDR",
  package_loaded = c("monocle", "igraph", "reshape2"),
  package_required = c(),
  par_set = makeParamSet(
    makeIntegerParam(id = "num_dimensions", lower = 2L, default = 2L, upper = 20L),
    makeDiscreteParam(id = "norm_method", default = "vstExprs", values = c("vstExprs", "log", "none")),
    makeIntegerParam(id = "maxIter", lower = 1L, default = 20L, upper = 100L),
    makeNumericParam(id = "sigma", lower = 0, default = .001, upper = 100),
    makeLogicalParam(id = "lambda_null", default = T),
    makeNumericParam(id = "lambda", lower = 0, default = 5, upper = 100),
    makeLogicalParam(id = "ncenter_null", default = T),
    makeIntegerParam(id = "ncenter", lower = 3, default = 5, upper = 20),
    makeNumericParam(id = "param.gamma", lower = 0, default = 20, upper = 1e5),
    makeNumericParam(id = "tol", lower = 0, default = .001, upper = 10),
    makeLogicalParam(id = "auto_param_selection", default = T)
  ),
  properties = c("tibble"),#, "dimred", "dimred_traj", "pseudotime"), # todo: implement other outputs
  run_fun = run_monocle_ddrtree,
  plot_fun = plot_monocle_ddrtree
)

#' @importFrom igraph degree all_shortest_paths distances
#' @importFrom reshape2 melt
run_monocle_ddrtree <- function(counts,
                                start_cell_id = NULL,
                                num_dimensions = 2,
                                norm_method = "vstExprs",
                                maxIter = 20,
                                sigma = 0.001,
                                lambda_null = T,
                                lambda = NULL,
                                ncenter_null = T,
                                ncenter = NULL,
                                param.gamma = 20,
                                tol = 0.001,
                                auto_param_selection = T) {
  requireNamespace("monocle")

  if (lambda_null) lambda <- NULL
  if (ncenter_null) ncenter <- NULL
  if (is.factor(norm_method)) norm_method <- as.character(norm_method)

  # load in the new dataset
  expr <- counts
  #expr <- log2(as.matrix(counts)+1)
  featureData <- new("AnnotatedDataFrame", data.frame(row.names = colnames(expr), gene_short_name = colnames(expr)))
  cds_1 <- newCellDataSet(t(expr), featureData = featureData)

  # estimate sparameters
  cds_1 <- estimateSizeFactors(cds_1)
  cds_1 <- estimateDispersions(cds_1)

  # reduce dimension
  cds_2 <- reduceDimension(cds_1,
                           max_components = num_dimensions, norm_method = norm_method,
                           maxIter = maxIter, sigma = sigma, lambda = lambda, ncenter = ncenter,
                           param.gamma = param.gamma, tol = tol, auto_param_selection = auto_param_selection)

  # order the cells
  cds_3 <- orderCells(cds_2, reverse = FALSE)

  # retrieve the graph and the root cell
  gr <- cds_3@auxOrderingData$DDRTree$pr_graph_cell_proj_tree

  if(!is.null(start_cell_id)) {
    root <- cds_3@auxOrderingData$DDRTree$root_cell
  } else {
    root <- start_cell_id
  }

  # find the branching cells and the terminal cells using the degree
  deg <- igraph::degree(gr, mode = c("all"))
  milestone_ids <- names(deg)[deg != 2]
  branching <- names(deg)[deg > 2]
  terminal <- names(deg)[deg == 1]


  asp <- igraph::all_shortest_paths(gr, from = root, to = milestone_ids, mode = c("all"))
  asp2 <- lapply(asp$res, function(path) {
    last_bit <- tail(which(path$name %in% milestone_ids), 2)
    if (length(last_bit) == 1) {
      path[last_bit]
    } else {
      path[last_bit[[1]]:last_bit[[2]]]
    }
  })
  asp2 <- asp2[sapply(asp2, length) > 1]
  dist_m <- igraph::distances(gr, v = milestone_ids, to = milestone_ids)

  milestone_network <- bind_rows(lapply(asp2, function(p) {
    from_ <- head(p, 1)$name
    to_ <- tail(p, 1)$name
    data_frame(from = paste0("milestone_", from_), to = paste0("milestone_", to_), length = dist_m[[from_, to_]])
  })) %>% mutate(length = length / max(length))

  milestone_percentages <- bind_rows(lapply(asp2, function(path) {
    dists <- t(igraph::distances(gr, v = path[c(1, length(path))], to = path))
    pct <- 1 - t(apply(dists, 1, function(x) x / sum(x)))
    pct %>%
      reshape2::melt(varnames = c("cell_id", "milestone_id"), value.name = "percentage") %>%
      mutate(cell_id = as.character(cell_id), milestone_id = paste0("milestone_", as.character(milestone_id)))
  })) %>% as_tibble() %>% filter(percentage > 0) %>% distinct(cell_id, milestone_id, percentage)
  milestone_ids <- paste0("milestone_", milestone_ids)

  # wrap output
  wrap_ti_prediction(
    ti_type = "tree",
    id = "monocleDDRtree",
    cell_ids = rownames(counts),
    milestone_ids = milestone_ids,
    milestone_network = milestone_network,
    milestone_percentages = milestone_percentages,
    cds = cds_3
  )
}

plot_monocle_ddrtree <- function(ti_predictions) {
  requireNamespace("monocle")
  monocle::plot_cell_trajectory(ti_predictions$cds)
}
