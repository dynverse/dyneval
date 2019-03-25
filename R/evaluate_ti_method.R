#' Running an evaluation of a method on a set of datasets with a set of parameters
#'
#' @inheritParams calculate_metrics
#' @inheritParams dynwrap::infer_trajectories
#' @param output_model Whether or not the model will be outputted.
#'
#' @keywords evaluation
#'
#' @importFrom dynwrap infer_trajectories add_cell_waypoints
#' @importFrom testthat expect_false expect_true
#'
#' @export
evaluate_ti_method <- function(
  dataset,
  method,
  parameters,
  metrics,
  give_priors = NULL,
  output_model = TRUE,
  seed = function() random_seed(),
  map_fun = map,
  verbose = FALSE
) {
  if (dynwrap::is_wrapper_with_trajectory(dataset)) {
    dataset <- list_as_tibble(list(dataset))
  }
  testthat::expect_true(all(mapdf_lgl(dataset, dynwrap::is_wrapper_with_waypoint_cells)))

  method_outputs <- dynwrap::infer_trajectories(
    dataset = dataset,
    method = method,
    parameters = parameters,
    give_priors = give_priors,
    seed = seed,
    map_fun = map_fun,
    verbose = verbose,
    return_verbose = TRUE
  )

  # Calculate scores
  eval_outputs <- map_fun(seq_len(nrow(dataset)), function(i) {
    dataseti <- dynutils::extract_row_to_list(dataset, i)

    # Fetch method outputs
    model <- method_outputs$model[[i]]

    if (!is.null(model)) {
      # Calculate geodesic distances
      time0 <- Sys.time()
      model <- model %>% dynwrap::add_cell_waypoints(num_cells_selected = length(dataseti$waypoint_cells))
      time1 <- Sys.time()
      time_cellwaypoints <- as.numeric(difftime(time1, time0, units = "sec"))
      df_cellwaypoints <- data_frame(time_cellwaypoints)
    } else {
      df_cellwaypoints <- NULL
    }

    # Calculate metrics
    metrics_summary <- calculate_metrics(dataseti, model, metrics)

    # Create summary statistics
    summary <- bind_cols(
      method_outputs$summary[[i]],
      df_cellwaypoints,
      metrics_summary
    )

    # if not interested in the created model,
    # do still return whether a model was created or not.
    if (!output_model && !is.null(model)) {
      model <- TRUE
    }

    # Return the output
    lst(summary, model)
  })

  # create output data structure
  list(
    summary = map_df(eval_outputs, "summary"),
    models = map(eval_outputs, "model")
  )
}
