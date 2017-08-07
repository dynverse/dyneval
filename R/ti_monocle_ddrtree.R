#' @import ParamHelpers
#' @export
description_monocle_ddrtree <- function() {
  list(
    name = "monocle_DDRtree",
    short_name = "monocle with DDRtree",
    package_load = c("monocle", "igraph", "tidyverse", "reshape2"),
    package_installed = c(),
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
}

#' @importFrom igraph degree all_shortest_paths distances
#' @importFrom reshape2 melt
#' @import dplyr
#' @import monocle
#'
#' @export
run_monocle_ddrtree <- function(counts,
                                num_dimensions = 2, norm_method = "vstExprs",
                                maxIter = 20, sigma = 0.001, lambda_null = T, lambda = NULL,
                                ncenter_null = T, ncenter = NULL, param.gamma = 20, tol = 0.001,
                                auto_param_selection = T) {
  if (lambda_null) lambda <- NULL
  if (ncenter_null) ncenter <- NULL
  if (is.factor(norm_method)) norm_method <- as.character(norm_method)

  # load in the new dataset
  cds_1 <- newCellDataSet(t(as.matrix(counts)))

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
  root <- cds_3@auxOrderingData$DDRTree$root_cell

  # find the branching cells and the terminal cells using the degree
  deg <- igraph::degree(gr, mode = c("all"))
  state_names <- names(deg)[deg != 2]
  branching <- names(deg)[deg > 2]
  terminal <- names(deg)[deg == 1]


  asp <- igraph::all_shortest_paths(gr, from = root, to = state_names, mode = c("all"))
  asp2 <- lapply(asp$res, function(path) {
    last_bit <- tail(which(path$name %in% state_names), 2)
    if (length(last_bit) == 1) {
      path[last_bit]
    } else {
      path[last_bit[[1]]:last_bit[[2]]]
    }
  })
  asp2 <- asp2[sapply(asp2, length) > 1]
  dist_m <- igraph::distances(gr, v = state_names, to = state_names)

  ids <- rownames(counts)
  state_network <- bind_rows(lapply(asp2, function(p) {
    from_ <- head(p, 1)$name
    to_ <- tail(p, 1)$name
    data_frame(from = paste0("state_", from_), to = paste0("state_", to_), length = dist_m[[from_, to_]])
  })) %>% mutate(length = length / max(length))

  state_percentages <- bind_rows(lapply(asp2, function(path) {
    dists <- t(igraph::distances(gr, v = path[c(1, length(path))], to = path))
    pct <- 1 - t(apply(dists, 1, function(x) x / sum(x)))
    pct %>%
      reshape2::melt(varnames = c("id", "state"), value.name = "percentage") %>%
      mutate(id = as.character(id), state = paste0("state_", as.character(state)))
  })) %>% as_data_frame() %>% filter(percentage > 0)
  state_names <- paste0("state_", state_names)

  # wrap output
  wrap_ti_prediction(
    ti_type = "tree",
    name = "monocleDDRtree",
    ids = ids,
    state_names = state_names,
    state_network = state_network,
    state_percentages = state_percentages,
    cds = cds_3
  )
}

#' @importFrom monocle plot_cell_trajectory
#' @export
plot_monocle_ddrtree <- function(ti_predictions) {
  monocle::plot_cell_trajectory(ti_predictions$cds)
}
