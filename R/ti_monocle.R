#' @export
makeRLearner.ti.monocle <- function() {
  makeRLearnerTI(
    cl = "ti.mococle",
    package = c("monocle"),
    par.set = makeParamSet(
      makeIntegerLearningParam(id = "num_dimensions", lower = 1L, default = 2L)
      # makeIntegerLearningParam(id = "num_clusters", lower = 2L, default = 4L)
    ),
    # properties = c("linear", "dimred_samples", "dimred_traj", "pseudotime"), # What to add?
    properties = c(),
    name = "monocle",
    short.name = "monocle"
  )
}

#' @import monocle
#' @importFrom igraph degree all_shortest_paths distances
#' @importFrom reshape2 melt
#' @import dplyr
#'
#' @export
trainLearner.ti.monocle <- function(.task, .subset, num_dimensions) {
  # subsetting will not work yet, but the function is already provided
  data <- get_task_data(.task, .subset)

  expression <- t(SCORPIUS::quant.scale(t(data$expression), 0))

  cds_1 <- monocle::newCellDataSet(t(as.matrix(expression)))
  cds_2 <- monocle::reduceDimension(cds_1, max_components = num_dimensions)
  cds_3 <- monocle::orderCells(cds_2, reverse = FALSE)

  gr <- cds_3@auxOrderingData$DDRTree$pr_graph_cell_proj_tree
  root <- cds_3@auxOrderingData$DDRTree$root_cell

  deg <- igraph::degree(gr, mode = c("all"))
  state_names <- names(deg)[deg != 2]

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

  state_network <- bind_rows(lapply(asp2, function(p) {
    from <- p %>% head(1) %>% .$name
    to <- p %>% tail(1) %>% .$name
    data_frame(from, to, length = dist_m[[from, to]])
  })) %>% mutate(length = length / max(length))

  pct_melted <- bind_rows(lapply(asp2, function(path) {
    dists <- t(igraph::distances(gr, v = path[c(1, length(path))], to = path))
    pct <- 1 - t(apply(dists, 1, function(x) x / sum(x)))
    pct %>%
      reshape2::melt(varnames = c("id", "waypoint")) %>%
      mutate(id = as.character(id), waypoint = as.character(waypoint))
  }))

  state_percentages <- pct_melted %>%
    unique %>%
    spread(waypoint, value, fill = 0) %>%
    slice(match(rownames(expression), id))

  state_percentages <- state_percentages[,c("id", state_names)]

  #' rename states
  state_network <- state_network %>% mutate(from = paste0("state_", from), to = paste0("state_", to))
  state_names <- paste0("state_", state_names)
  colnames(state_percentages) <- c("id", state_names)

  #' wrap output
  wrap_ti_prediction(
    ti_type = "tree",
    name = "monocle",
    state_names,
    state_network,
    state_percentages,
    task_id = get_task_identifier(.task),
    cds = cds_3
  )
}

#' @importFrom monocle plot_cell_trajectory
#' @export
plotLearner.ti.monocle <- function(ti_predictions) {
  monocle::plot_cell_trajectory(ti_predictions$cds)
}
