#' Running an evaluation of a method on a set of datasets with a set of parameters
#'
#' @inheritParams calculate_metrics
#' @param datasets The datasets on which to evaluate.
#' @param method The method to evaluate.
#' @param parameters The parameters to evaluate with.
#' @param extra_metrics Extra metrics to calculate but not evaluate with.
#' @param output_model Whether or not the model will be outputted.
#' @param mc_cores The number of cores to use, allowing to parallellise the different datasets
#' @param verbose Whether or not to print extra information output.
#'
#' @export
#' @importFrom dynwrap infer_trajectories add_cell_waypoints
#' @importFrom parallel mclapply
#' @importFrom testthat expect_false expect_true
#' @importFrom readr write_rds
evaluate_ti_method <- function(
  datasets,
  method,
  parameters,
  metrics,
  extra_metrics = NULL,
  output_model = TRUE,
  mc_cores = 1,
  verbose = FALSE
) {
  testthat::expect_is(datasets, "tbl")
  testthat::expect_true(all(unlist(tmap(datasets, dynwrap::is_wrapper_with_waypoint_cells))))

  calc_metrics <- c(metrics, extra_metrics)
  calc_metrics <- calc_metrics[!duplicated(calc_metrics)]

  metric_names <- sapply(seq_along(metrics), function(i) {
    metric <- metrics[[i]]
    if (is.function(metric)) {
      names(metrics)[[i]]
    } else if (is.character(metric)) {
      metric
    } else {
      stop("Unexpected metric, check documentation.")
    }
  })

  method_outputs <- dynwrap::infer_trajectories(
    dataset = datasets,
    method = method,
    parameters = parameters,
    mc_cores = mc_cores,
    verbose = TRUE,
    capture_output = TRUE
  )

  # Calculate scores
  eval_outputs <- parallel::mclapply(seq_len(nrow(datasets)), mc.cores = mc_cores, function(i) {
    dataset <- dynutils::extract_row_to_list(datasets, i)

    # Fetch method outputs
    model <- method_outputs$model[[i]]

    if (!is.null(model)) {
      # Calculate geodesic distances
      time0 <- Sys.time()
      model <- model %>% dynwrap::add_cell_waypoints(num_cells_selected = length(dataset$waypoint_cells))
      time1 <- Sys.time()
      time_cellwaypoints <- as.numeric(difftime(time1, time0, units = "sec"))
      df_cellwaypoints <- data_frame(time_cellwaypoints)
    } else {
      df_cellwaypoints <- NULL
    }

    # Calculate metrics
    metrics_summary <- calculate_metrics(dataset, model, calc_metrics)

    # Create summary statistics
    summary <- bind_cols(
      method_outputs$summary[[i]],
      df_cellwaypoints,
      metrics_summary
    )

    # Return the output
    out <- list(summary = summary)
    if (output_model) {
      out$model <- model
    }
    out
  })

  summary <- map_dfr(eval_outputs, "summary")

  # Calculate the final score
  score <- summary %>%
    summarise_at(metric_names, funs(mean)) %>%
    as.matrix %>%
    as.vector %>%
    setNames(metric_names)

  # create output data structure
  out <- list(
    score = score,
    summary = summary
  )

  # add models if desired
  if (output_model) {
    out$models <- eval_outputs %>% map("model")
  }

  # return output
  out
}
