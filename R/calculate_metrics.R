#' Calculate the performance of a model with respect to a task
#'
#' @param task the original task
#' @param model the predicted model
#' @param metrics which metrics to evaluate:
#' \enumerate{
#'   \item Coranking of geodesic distances: \code{"mean_R_nx"}, \code{"auc_R_nx"}, \code{"Q_local"}, \code{"Q_global"}
#'   \item Spearman correlation of geodesic distances: \code{"correlation"}
#'   \item Isomorphism of the two networks: \code{"isomorphic"}
#'   \item GEDEVO Graph Edit Distance: \code{"ged"}
#'   \item Earth Mover's Distance on orbit counts: \code{"net_emd"}
#'   \item Genetic Algorithm for aligning small graphs: \code{"robbie_network_score"}
#'   \item Mantel test p-value: \code{"mantel_pval"}
#' }
#'
#' @importFrom igraph is_isomorphic_to graph_from_data_frame
#'
#' @export
calculate_metrics <- function(task, model, metrics) {
  summary_list <- list()

  if (any(c("mean_R_nx", "auc_R_nx", "Q_local", "Q_global", "correlation", "mantel_pval") %in% metrics)) {
    # compute coranking
    time0 <- Sys.time()
    coranking <- compute_coranking(task$geodesic_dist, model$geodesic_dist)
    summary_list <- c(summary_list, coranking$summary)
    time1 <- Sys.time()
    summary_list$time_coranking <- as.numeric(difftime(time1, time0, units = "sec"))

    # compute corrrelation
    time0 <- Sys.time()
    summary_list$correlation <- cor(task$geodesic_dist %>% as.vector, model$geodesic_dist %>% as.vector, method="spearman")
    time1 <- Sys.time()
    summary_list$time_correlation <- as.numeric(difftime(time1, time0, units = "sec"))

    # compute mantel pval
    time0 <- Sys.time()
    mantel <- vegan::mantel(task$geodesic_dist, model$geodesic_dist, permutations = 100)
    summary_list$mantel_pval <- -log10(mantel$signif)
    time1 <- Sys.time()
    summary_list$time_mantel <- as.numeric(difftime(time1, time0, units = "sec"))
  }

  net1 <- dynutils::simplify_milestone_network(model$milestone_network)
  net2 <- dynutils::simplify_milestone_network(task$milestone_network) %>% filter(to != "FILTERED_CELLS")

  # Compute the milestone network isomorphic
  if ("isomorphic" %in% metrics) {
    time0 <- Sys.time()

    summary_list$isomorphic <- (is_isomorphic_to(
      graph_from_data_frame(net1),
      graph_from_data_frame(net2)
    ))+0
    time1 <- Sys.time()
    summary_list$time_isomorphic <- as.numeric(difftime(time1, time0, units = "sec"))
  }

  # Compute the milestone network GED
  if ("ged" %in% metrics) {
    time0 <- Sys.time()
    summary_list$ged <- calculate_ged(net1, net2)
    time1 <- Sys.time()
    summary_list$time_ged <- as.numeric(difftime(time1, time0, units = "sec"))
  }

  # Compute Netdist EMD (see scripts/wouter/network_scores_tests.R)
  if ("net_emd" %in% metrics) {
    time0 <- Sys.time()
    gdd1 <- netdist::gdd(net1 %>% igraph::graph_from_data_frame())
    gdd2 <- netdist::gdd(net2 %>% igraph::graph_from_data_frame())
    summary_list$net_emd <- netdist::net_emd(gdd1, gdd2)
    time1 <- Sys.time()
    summary_list$time_net_emd <- 1-as.numeric(difftime(time1, time0, units = "sec"))
  }

  if ("robbie_network_score" %in% metrics) {
    time0 <- Sys.time()
    summary_list$robbie_network_score <- calculate_robbie_network_score(net1, net2)
    time1 <- Sys.time()
    summary_list$time_robbie_network_score <- 1-as.numeric(difftime(time1, time0, units = "sec"))
  }

  summary <- as_tibble(summary_list)

  list(coranking = coranking, summary = summary)
}