#' Calculate the performance of a model with respect to a task
#'
#' @param task the original task
#' @param model the predicted model
#' @param metrics which metrics to evaluate:
#' \enumerate{
#'   \item Coranking of geodesic distances: \code{"mean_R_nx"}, \code{"auc_R_nx"}, \code{"Q_local"}, \code{"Q_global"}
#'   \item Spearman correlation of geodesic distances: \code{"correlation"}
#'   \item Mantel test p-value: \code{"mantel_pval"}
#'   \item Edge flip score: \code{"edge_flip"}
#'   \item RF MSE: \code{"rf_mse"}, \code{"rf_rsq"}
#' }
#'
#' @importFrom igraph is_isomorphic_to graph_from_data_frame
#' @importFrom testthat expect_equal
#'
#' @export
calculate_metrics <- function(
  task,
  model,
  metrics = c("mean_R_nx", "auc_R_nx", "Q_local", "Q_global", "correlation", "mantel_pval", "edge_flip", "rf_mse", "rf_rsq")
) {

  testthat::expect_equal(task$cell_ids, rownames(task$geodesic_dist))
  testthat::expect_equal(task$cell_ids, colnames(task$geodesic_dist))

  if (!is.null(model)) {
    testthat::expect_equal(task$cell_ids, model$cell_ids)
    testthat::expect_equal(task$cell_ids, rownames(model$geodesic_dist))
    testthat::expect_equal(task$cell_ids, colnames(model$geodesic_dist))
  }

  summary_list <- list()

  if (any(c("mean_R_nx", "auc_R_nx", "Q_local", "Q_global", "correlation", "mantel_pval") %in% metrics)) {

    if (!is.null(model)) {
      task$geodesic_dist[is.infinite(task$geodesic_dist)] <- .Machine$double.xmax
      model$geodesic_dist[is.infinite(model$geodesic_dist)] <- .Machine$double.xmax

      # compute coranking
      time0 <- Sys.time()
      coranking <- compute_coranking(task$geodesic_dist, model$geodesic_dist)
      summary_list <- c(summary_list, coranking$summary)
      time1 <- Sys.time()
      summary_list$time_coranking <- as.numeric(difftime(time1, time0, units = "sec"))

      # compute corrrelation
      time0 <- Sys.time()
      summary_list$correlation <- cor(task$geodesic_dist %>% as.vector, model$geodesic_dist %>% as.vector, method = "spearman")
      time1 <- Sys.time()
      summary_list$time_correlation <- as.numeric(difftime(time1, time0, units = "sec"))

      # compute mantel pval
      time0 <- Sys.time()
      mantel <- vegan::mantel(task$geodesic_dist, model$geodesic_dist, permutations = 100)
      summary_list$mantel_pval <- -log10(mantel$signif)
      time1 <- Sys.time()
      summary_list$time_mantel <- as.numeric(difftime(time1, time0, units = "sec"))
    } else {
      summary_list <- c(summary_list, list(mean_R_nx = 0, auc_R_nx = 0, Q_local = 0, Q_global = 0, correlation = 0, mantel_pval = 0))
    }
  }

  if ("edge_flip" %in% metrics) {
    if (!is.null(model)) {
      net1 <- model$milestone_network %>% filter(to != "FILTERED_CELLS")
      net2 <- task$milestone_network %>% filter(to != "FILTERED_CELLS")

      time0 <- Sys.time()
      summary_list$edge_flip <- calculate_edge_flip(net1, net2)
      time1 <- Sys.time()
      summary_list$time_edge_flip <- as.numeric(difftime(time1, time0, units = "sec"))
    } else {
      summary_list$edge_flip <- 0
    }
  }

  if (any(c("rf_mse", "rf_rsq") %in% metrics)) {
    time0 <- Sys.time()
    rfmse <- compute_rfmse(task, model)
    time1 <- Sys.time()
    summary_list$time_rf <- as.numeric(difftime(time1, time0, units = "sec"))
    summary_list$rf_mse <- rfmse$summary$rf_mse
    summary_list$rf_rsq <- rfmse$summary$rf_rsq
  }

  summary <- as_tibble(summary_list)

  list(coranking = coranking, summary = summary)
}
