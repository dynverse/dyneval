
#' Calculate Earth Mover's Distance between cells in a trajectory
#'
#' @param traj the trajectory
#' @param parallel compute distances in parallel
#'
#' @export
#'
#' @importFrom igraph graph_from_data_frame E distances
#' @importFrom transport transport
#' @import dplyr
compute_em_dist <- function(traj, parallel = T) {
  state_network <- traj$state_network
  state_names <- traj$state_names
  state_percentages <- traj$state_percentages

  # calculate the shortest path distances between milestones
  gr <- igraph::graph_from_data_frame(state_network, directed = T, vertices = state_names)
  milestone_distances <- igraph::distances(gr, weights = igraph::E(gr)$length)

  # transport percentages data
  pct <- as.matrix(state_percentages[,-1])
  ids <- state_percentages$id
  rownames(pct) <- ids

  # calculate the distances between cells
  # could also implement myself https://people.cs.umass.edu/~mcgregor/papers/13-approx1.pdf
  sapply_fun <-
    if (parallel) {
      function(X, FUN, ...) {
        unlist(parallel::mclapply(X, FUN, mc.cores = 8, ...))
      }
    } else {
      pbapply::pbsapply
    }

  cell_dists <- expand.grid(from = rownames(pct), to = rownames(pct)) %>% filter(as.integer(from) < as.integer(to))
  cell_dists$dist <- sapply_fun(seq_len(nrow(cell_dists)), function(i) {
    a <- pct[cell_dists$from[[i]],]
    b <- pct[cell_dists$to[[i]],]
    transport::transport(a, b, costm = milestone_distances, method = "shortsimplex") %>%
      mutate(dist = milestone_distances[cbind(from,to)], mult = mass * dist) %>%
      .$mult %>%
      sum
  })

  dist_m <- cell_dists %>% reshape2::acast(from~to, value.var = "dist", fill = 0)
  dist_m2 <- matrix(0, nrow = length(ids), ncol = length(ids), dimnames = list(ids, ids))
  dist_m2[rownames(dist_m), colnames(dist_m)] <- dist_m
  (dist_m2 + t(dist_m2))
}

#' Plot the Earth Mover's distances in a heatmap
#'
#' @param traj the trajectory (less than 500 cells is recommended)
#' @param emdist the Earth Mover's distances as calculated by \code{\link{emdist}}

#' @export
#'
#' @importFrom reshape2 acast
#' @importFrom pheatmap pheatmap
plot_emdist <- function(traj, dist) {
  state_percentages <- traj$state_percentages
  pct <- as.data.frame(state_percentages[,-1])
  rownames(pct) <- state_percentages$id
  pheatmap::pheatmap(
    dist,
    cluster_rows = F,
    cluster_cols = F,
    annotation_col = pct,
    annotation_row = pct,
    legend = F,
    legend_labels = F,
    legend_breaks = F,
    border_color = NA,
    show_rownames = F,
    show_colnames = F)
}

#' Compute the coranking matrix and
#'
#' @param gold_dist A data frame containing the pairwise distances in the original space
#' @param pred_dist A data frame containing the pairwise distances in the new space
#'
#' @export
#'
#' @importFrom coRanking coranking LCMC
compute_coranking <- function(gold_dist, pred_dist) {
  gold_dist <- gold_dist + runif(length(gold_dist), 0, 1e-20)
  pred_dist <- pred_dist + runif(length(pred_dist), 0, 1e-20)
  diag(gold_dist) <- 0
  diag(pred_dist) <- 0

  corank <- coRanking::coranking(gold_dist, pred_dist, input = "dist")

  # coRanking::imageplot(corank)
  # pheatmap::pheatmap(corank, cluster_rows = F, cluster_cols = F, color = colorRampPalette(c("white", "black"))(100), show_rownames = F, show_colnames = F)

  lcmc <- coRanking::LCMC(corank)

  auc_lcmc <- mean(lcmc)

  list(
    corank = corank,
    lcmc = lcmc,
    auc_lcmc = auc_lcmc
  )
}
