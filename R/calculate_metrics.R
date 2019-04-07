#' Calculate the similarity between two trajectory models.
#'
#' One use case for these metrics is to calculate the accuracy of a certain prediction compared to a reference trajectory. However, these metrics can also be used for other purposes, such as clustering of trajectories.
#'
#' Some metrics are asymmetric (see `dyneval::metrics$symmetric`), in which case the order of the dataset and model parameters matters.
#'
#' @param dataset The first trajectory, in most cases a gold standard trajectory
#' @param model The second trajectory, in most cases a predicted trajectory
#' @param metrics Which metrics to evaluate. Check `dyneval::metrics` for a list of possible metrics.
#'   Passing a custom metric function with format `function(dataset, model) { 1 }` is also supported. The name of this function within the list will be used as the name of the metric.
#' @param expression_source The expression data matrix, with features as columns.
#'   * If a matrix is provided, it is used as is.
#'   * If a character is provided, `dataset[[expression_source]]` should contain the matrix.
#'   * If a function is provided, that function will be called in order to obtain the expression (useful for lazy loading).
#'
#' @importFrom igraph is_isomorphic_to graph_from_data_frame
#' @importFrom testthat expect_equal expect_true
#' @importFrom dynwrap is_wrapper_with_waypoint_cells calculate_geodesic_distances
#'
#' @keywords metric
#'
#' @export
calculate_metrics <- function(
  dataset,
  model,
  metrics = dyneval::metrics$metric_id,
  expression_source = dataset$expression
) {
  # check if all function metrics are named
  if (!all(sapply(seq_along(metrics), function(i) !is.function(metrics[[i]]) || !is.null(names(metrics)[[i]])))) {
    stop("All custom metrics (functions) must be named!")
  }

  # check all character function metrics
  valid_metrics <- dyneval::metrics$metric_id
  character_metrics <- as.character(keep(metrics, is.character))
  if (!all(character_metrics %in% valid_metrics)) {
    stop("Invalid metrics: ", glue::glue_collapse(setdiff(character_metrics, valid_metrics), ", "))
  }

  summary_list <- list()

  dataset <- dynwrap::simplify_trajectory(dataset)
  if (!is.null(model)) {
    model <- dynwrap::simplify_trajectory(model)
  }

  if ("correlation" %in% metrics) {
    testthat::expect_true(dynwrap::is_wrapper_with_waypoint_cells(dataset))
    testthat::expect_true(is.null(model) || dynwrap::is_wrapper_with_waypoint_cells(model))

    # calculate geodesic distances
    if (!is.null(model)) {
      testthat::expect_true(all(model$cell_ids %in% dataset$cell_ids))
      model$cell_ids <- dataset$cell_ids

      waypoints <- unique(c(dataset$waypoint_cells, model$waypoint_cells))

      # compute waypointed geodesic distances
      time0 <- Sys.time()
      dataset$geodesic_dist <- dynwrap::calculate_geodesic_distances(dataset, waypoints)
      model$geodesic_dist <- dynwrap::calculate_geodesic_distances(model, waypoints)
      time1 <- Sys.time()
      summary_list$time_waypointedgeodesic <- as.numeric(difftime(time1, time0, units = "sec"))
    }


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

  if ("him" %in% metrics) {
    if (!is.null(model)) {
      net1 <- model$milestone_network
      net2 <- dataset$milestone_network

      time0 <- Sys.time()
      summary_list$him <- calculate_him(net1, net2)
      time1 <- Sys.time()
      summary_list$time_him <- as.numeric(difftime(time1, time0, units = "sec"))
    } else {
      summary_list$him <- 0
    }
  }

  if ("isomorphic" %in% metrics) {
    if (!is.null(model)) {
      graph1 <- model$milestone_network %>% igraph::graph_from_data_frame()
      graph2 <- dataset$milestone_network %>% igraph::graph_from_data_frame()

      time0 <- Sys.time()
      summary_list$isomorphic <- as.numeric(igraph::isomorphic(graph1, graph2))
      time1 <- Sys.time()
      summary_list$time_isomorphic <- as.numeric(difftime(time1, time0, units = "sec"))
    } else {
      summary_list$isomorphic <- 0
    }
  }

  if (any(c("rf_mse", "rf_rsq", "rf_nmse", "lm_mse", "lm_rsq", "lm_nmse") %in% metrics)) {
    time0 <- Sys.time()
    position_predict <- calculate_position_predict(dataset, model)
    time1 <- Sys.time()
    summary_list$time_pp <- as.numeric(difftime(time1, time0, units = "sec"))

    summary_list <- c(
      summary_list,
      position_predict$summary[intersect(metrics, names(position_predict$summary))]
    )
  }

  if (any(c("featureimp_cor", "featureimp_wcor") %in% metrics)) {
    time0 <- Sys.time()
    featureimp <- calculate_featureimp_cor(dataset, model, expression_source = expression_source)
    time1 <- Sys.time()
    summary_list$time_featureimp <- as.numeric(difftime(time1, time0, units = "sec"))
    summary_list$featureimp_cor <- featureimp$featureimp_cor
    summary_list$featureimp_wcor <- featureimp$featureimp_wcor
  }

  if (any(c("featureimp_ks", "featureimp_wilcox") %in% metrics)) {
    time0 <- Sys.time()
    featureimp <- calculate_featureimp_enrichment(dataset, model, expression_source = expression_source)
    time1 <- Sys.time()
    summary_list$time_featureimp_enrichment <- as.numeric(difftime(time1, time0, units = "sec"))
    summary_list$featureimp_ks <- featureimp$featureimp_ks
    summary_list$featureimp_wilcox <- featureimp$featureimp_wilcox
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
