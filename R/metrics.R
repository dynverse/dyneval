
#' Calculate Earth Mover's Distance between cells in a trajectory
#'
#' @param traj the trajectory
#'
#' @export
#'
#' @importFrom igraph graph_from_data_frame E distances
#' @importFrom transport transport
#' @import dplyr
emdist <- function(traj) {
  state_network <- traj$state_network
  state_names <- traj$state_names
  state_percentages <- traj$state_percentages

  # calculate the shortest path distances between milestones
  gr <- igraph::graph_from_data_frame(state_network, directed = T, vertices = state_names)
  milestone_distances <- igraph::distances(gr, weights = igraph::E(gr)$length)

  # transport percentages data
  pct <- as.matrix(state_percentages[,-1])
  rownames(pct) <- state_percentages$id

  # calculate the distances between cells
  # could also implement myself https://people.cs.umass.edu/~mcgregor/papers/13-approx1.pdf
  cell_dists <- expand.grid(from = rownames(pct), to = rownames(pct))
  cell_dists$dist <- pbapply::pbsapply(seq_len(nrow(cell_dists)), function(i) {
    a <- pct[cell_dists$from[[i]],]
    b <- pct[cell_dists$to[[i]],]
    transport::transport(a, b, costm = milestone_distances, method = "revsimplex") %>%
      mutate(dist = milestone_distances[cbind(from,to)], mult = mass * dist) %>%
      .$mult %>%
      sum
  })

  cell_dists
}

#' Plot the Earth Mover's distances in a heatmap
#'
#' @param traj the trajectory (less than 500 cells is recommended)
#' @param emdist the Earth Mover's distances as calculated by \code{\link{emdist}}

#' @export
#'
#' @importFrom reshape2 acast
#' @importFrom pheatmap pheatmap
plot_emdist <- function(traj, emdist) {
  state_percentages <- traj$state_percentages
  pct <- as.data.frame(state_percentages[,-1])
  rownames(pct) <- state_percentages$id
  dist_mat <- reshape2::acast(emdist, from~to, value.var = "dist")
  pheatmap::pheatmap(dist_mat, cluster_rows = F, cluster_cols = F, annotation_col = pct, annotation_row = pct)
}
