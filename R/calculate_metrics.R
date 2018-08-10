#' Calculate the performance of a model with respect to a dataset
#'
#' @param dataset the original dataset
#' @param model the predicted model
#' @param metrics which metrics to evaluate:
#' \enumerate{
#'   \item Spearman correlation of geodesic distances: \code{"correlation"}
#'   \item Edge flip score: \code{"edge_flip"}
#'   \item RF MSE: \code{"rf_mse"}, \code{"rf_rsq"}
#'   \item Similarity in feature importance: \code{"featureimp_cor"}
#'   \item Overlap between branches: \code{"F1_branches"}
#'   \item Custom metric function. Format: \code{function(dataset, model) { 1 }}
#' }
#'
#' @importFrom igraph is_isomorphic_to graph_from_data_frame
#' @importFrom testthat expect_equal expect_true
#' @importFrom dynwrap is_wrapper_with_waypoint_cells compute_tented_geodesic_distances
#'
#' @export
calculate_metrics <- function(
  dataset,
  model,
  metrics = c("correlation", "edge_flip", "rf_mse", "rf_rsq", "lm_mse", "lm_rsq", "featureimp_cor", "F1_branches", "F1_milestones")
) {
  testthat::expect_true(dynwrap::is_wrapper_with_waypoint_cells(dataset))
  testthat::expect_true(is.null(model) || dynwrap::is_wrapper_with_waypoint_cells(model))

  if (!all(sapply(seq_along(metrics), function(i) !is.function(metrics[[i]]) || !is.null(names(metrics)[[i]])))) {
    stop("All custom metrics (functions) must be named!")
  }

  summary_list <- list()

  dataset <- dynwrap::simplify_trajectory(dataset)
  if (!is.null(model)) {
    model <- dynwrap::simplify_trajectory(model)
  }

  if (!is.null(model)) {
    testthat::expect_true(all(model$cell_ids %in% dataset$cell_ids))
    model$cell_ids <- dataset$cell_ids

    waypoints <- unique(c(dataset$waypoint_cells, model$waypoint_cells))

    # compute waypointed geodesic distances
    time0 <- Sys.time()
    dataset$geodesic_dist <- dynwrap::compute_tented_geodesic_distances(dataset, waypoints)
    model$geodesic_dist <- dynwrap::compute_tented_geodesic_distances(model, waypoints)
    time1 <- Sys.time()
    summary_list$time_waypointedgeodesic <- as.numeric(difftime(time1, time0, units = "sec"))
  }

  if (("correlation" %in% metrics)) {
    if (!is.null(model)) {
      dataset$geodesic_dist[is.infinite(dataset$geodesic_dist)] <- .Machine$double.xmax
      model$geodesic_dist[is.infinite(model$geodesic_dist)] <- .Machine$double.xmax

      # make sure the order of the cells is exactly equal
      testthat::expect_equal(rownames(dataset$geodesic_dist), rownames(model$geodesic_dist))
      testthat::expect_equal(colnames(dataset$geodesic_dist), colnames(model$geodesic_dist))

      # compute corrrelation
      time0 <- Sys.time()
      if (length(unique(c(model$geodesic_dist))) == 1 || length(unique(c(dataset$geodesic_dist))) == 1) {
        summary_list$correlation <- 0
      } else {
        summary_list$correlation <- cor(
          dataset$geodesic_dist %>% as.vector,
          model$geodesic_dist %>% as.vector,
          method = "spearman"
        ) %>% max(0)
      }

      time1 <- Sys.time()
      summary_list$time_correlation <- as.numeric(difftime(time1, time0, units = "sec"))

    } else {
      summary_list <- c(summary_list, list(correlation = 0))
    }
  }

  if ("edge_flip" %in% metrics) {
    if (!is.null(model)) {
      net1 <- model$milestone_network
      net2 <- dataset$milestone_network

      time0 <- Sys.time()
      summary_list$edge_flip <- calculate_edge_flip(net1, net2)
      time1 <- Sys.time()
      summary_list$time_edge_flip <- as.numeric(difftime(time1, time0, units = "sec"))
    } else {
      summary_list$edge_flip <- 0
    }
  }

  if (any(c("rf_mse", "rf_rsq", "rf_nmse", "lm_mse", "lm_rsq", "lm_nmse") %in% metrics)) {
    time0 <- Sys.time()
    position_predict <- compute_position_predict(dataset, model)
    time1 <- Sys.time()
    summary_list$time_pp <- as.numeric(difftime(time1, time0, units = "sec"))

    summary_list <- c(
      summary_list,
      position_predict$summary[intersect(metrics, names(position_predict$summary))]
    )
  }

  if ("featureimp_cor" %in% metrics) {
    time0 <- Sys.time()
    fimp <- compute_featureimp(dataset, model)
    time1 <- Sys.time()
    summary_list$time_featureimp <- as.numeric(difftime(time1, time0, units = "sec"))
    summary_list$featureimp_cor <- fimp$featureimp_cor
  }



  if (any(c("recovery_branches", "relevance_branches", "F1_branches", "recovery_milestones", "relevance_milestones", "F1_milestones") %in% metrics)) {
    dataset_simplified <- dynwrap::simplify_trajectory(dataset, allow_self_loops = TRUE)
    if (!is.null(model)) {
      model_simplified <- dynwrap::simplify_trajectory(model, allow_self_loops = TRUE)
    } else {
      model_simplified <- NULL
    }

    if (any(c("recovery_branches", "relevance_branches", "F1_branches") %in% metrics)) {
      time0 <- Sys.time()
      mapping_branches <- calculate_mapping_branches(dataset_simplified, model_simplified)
      time1 <- Sys.time()
      summary_list$time_mapping_branches <- as.numeric(difftime(time1, time0, units = "sec"))
      summary_list <- c(summary_list, mapping_branches)
    }

    if (any(c("recovery_milestones", "relevance_milestones", "F1_milestones") %in% metrics)) {
      time0 <- Sys.time()
      mapping_milestones <- calculate_mapping_milestones(dataset_simplified, model_simplified)
      time1 <- Sys.time()
      summary_list$time_mapping_milestones <- as.numeric(difftime(time1, time0, units = "sec"))
      summary_list <- c(summary_list, mapping_milestones)
    }
  }

  for (i in seq_along(metrics)) {
    f <- metrics[[i]]
    fn <- names(metrics)[[i]]
    if (is.function(f)) {

      if (!is.null(model)) {
        time0 <- Sys.time()
        output <- f(dataset, model)
        time1 <- Sys.time()
        summary_list[[paste0("time_", fn)]] <- as.numeric(difftime(time1, time0, units = "sec"))

        if (length(output) != 1) {
          stop("Metric ", sQuote(fn), " should return exactly 1 numeric score.")
        }
      } else {
        output <- 0
      }
      names(output) <- fn

      summary_list[names(output)] <- output
    }
  }

  summary <- as_tibble(summary_list)

  summary
}
