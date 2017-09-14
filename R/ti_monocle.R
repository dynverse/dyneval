#' Description for monocle DDRTree
#' @export
description_monocle_ddrtree <- function() generic_monocle_description("DDRTree")

#' Description for monocle PQTree
#' @export
description_monocle_pqtree <- function() generic_monocle_description("ICA")

generic_monocle_description <- function(reduction_method) {
  if(reduction_method == "DDRTree") {
    par_set = makeParamSet(
      makeIntegerParam(id = "num_dimensions", lower = 2L, default = 2L, upper = 20L),
      makeDiscreteParam(id = "norm_method", default = "vstExprs", values = c("vstExprs", "log", "none")),
      makeIntegerParam(id = "maxIter", lower = 1L, default = 20L, upper = 100L),
      makeNumericParam(id = "sigma", lower = 0, default = .001, upper = 100),
      makeLogicalParam(id = "lambda_null", default = TRUE),
      makeNumericParam(id = "lambda", lower = 0, default = 5, upper = 100),
      makeLogicalParam(id = "ncenter_null", default = TRUE),
      makeIntegerParam(id = "ncenter", lower = 3, default = 5, upper = 20),
      makeNumericParam(id = "param.gamma", lower = 0, default = 20, upper = 1e5),
      makeNumericParam(id = "tol", lower = 0, default = .001, upper = 10),
      makeLogicalParam(id = "auto_param_selection", default = TRUE)
    )
  } else if(reduction_method == "ICA"){
    par_set = makeParamSet(
      makeIntegerParam(id = "num_dimensions", lower = 2L, default = 2L, upper = 20L),
      makeDiscreteParam(id = "norm_method", default = "vstExprs", values = c("vstExprs", "log", "none")),
      makeNumericParam(id = "lambda", lower = 0, default = 5, upper = 100),
      makeLogicalParam(id = "ncenter_null", default = TRUE),
      makeNumericParam(id = "tol", lower = 0, default = .001, upper = 10),
      makeLogicalParam(id = "auto_param_selection", default = TRUE)
    )
  }

  create_description(
    name = glue::glue("monocle with {ifelse(reduction_method == 'DDRTree', 'DDRTree', 'ICA and PQTree')}"),
    short_name = ifelse(reduction_method == "DDRTree", "monocDDR", "monocPQ"),
    package_loaded = c("monocle", "igraph", "reshape2"),
    package_required = c(),
    par_set = par_set,
    properties = c("tibble"),#, "dimred", "dimred_traj", "pseudotime"), # todo: implement other outputs
    run_fun = ifelse(reduction_method == "DDRTree", run_monocle_ddrtree, run_monocle_pqtree),
    plot_fun = plot_monocle
  )
}

run_monocle_ddrtree <- function(...) {
  args <- list(...)
  print(args)
  args$reduction_method <- "DDRTree"

  do.call(run_monocle, args)
}

run_monocle_pqtree <- function(...) {
  args <- list(...)
  args$reduction_method <- "ICA"

  do.call(run_monocle, args)
}

#' @importFrom igraph degree all_shortest_paths distances
#' @importFrom reshape2 melt
run_monocle <- function(counts,
                        reduction_method,
                        start_cell_id = NULL,
                        num_dimensions = 2,
                        norm_method = "vstExprs",
                        maxIter = 20,
                        sigma = 0.001,
                        lambda_null = TRUE,
                        lambda = NULL,
                        ncenter_null = TRUE,
                        ncenter = NULL,
                        param.gamma = 20,
                        tol = 0.001,
                        auto_param_selection = TRUE) {
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
  if(reduction_method == "DDRTree") {
    cds_2 <- reduceDimension(cds_1,
                             reduction_method = reduction_method,
                             max_components = num_dimensions, norm_method = norm_method,
                             maxIter = maxIter, sigma = sigma, lambda = lambda, ncenter = ncenter,
                             param.gamma = param.gamma, tol = tol, auto_param_selection = auto_param_selection)
  } else if (reduction_method == "ICA") {
    cds_2 <- reduceDimension(cds_1,
                             reduction_method = reduction_method,
                             max_components = num_dimensions, norm_method = norm_method,
                             tol = tol, auto_param_selection = auto_param_selection)
  }


  # order the cells
  cds_3 <- orderCells(cds_2, reverse = FALSE)

  orderingData <- cds_3@auxOrderingData[[reduction_method]] # location of ordering data depends on reduction method
  # retrieve the graph and the root cell, also depends on reduction method (-_-)
  if (reduction_method == "DDRTree") {
    gr <- orderingData$pr_graph_cell_proj_tree
  } else if (reduction_method == "ICA") {
    gr <- orderingData$cell_ordering_tree
  }

  if(is.null(start_cell_id)) {
    root <- orderingData$root_cell
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
    print(last_bit)
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
  })) %>% mutate(length = length / max(length), directed=TRUE)

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
    id = ifelse(reduction_method == "DDRTree", "monocDDR", "monocPQ"),
    cell_ids = rownames(counts),
    milestone_ids = milestone_ids,
    milestone_network = milestone_network,
    milestone_percentages = milestone_percentages,
    cds = cds_3
  )
}

plot_monocle <- function(ti_predictions) {
  requireNamespace("monocle")
  monocle::plot_cell_trajectory(ti_predictions$cds)
}
