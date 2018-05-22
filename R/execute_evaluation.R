#' Running an evaluation of a method on a set of tasks with a set of parameters
#'
#' @inheritParams calculate_metrics
#' @param tasks The tasks on which to evaluate.
#' @param method The method to evaluate.
#' @param parameters The parameters to evaluate with.
#' @param extra_metrics Extra metrics to calculate but not evaluate with.
#' @param output_model Whether or not the model will be outputted.
#' @param mc_cores The number of cores to use, allowing to parallellise the different tasks
#' @param verbose Whether or not to print extra information output.
#'
#' @export
#' @importFrom dynmethods execute_method
#' @importFrom parallel mclapply
#' @importFrom testthat expect_false expect_true
#' @importFrom readr write_rds
execute_evaluation <- function(
  tasks,
  method,
  parameters,
  metrics,
  extra_metrics = NULL,
  output_model = TRUE,
  mc_cores = 1,
  verbose = FALSE
) {
  testthat::expect_true("waypoint_cells" %in% colnames(tasks))
  testthat::expect_false(any(sapply(tasks$waypoint_cells, is.null)))

  calc_metrics <- unique(c(metrics, extra_metrics))

  method_outputs <- dynmethods::execute_method(
    tasks = tasks,
    method = method,
    parameters = parameters,
    mc_cores = mc_cores,
    verbose = verbose
  )

  # Calculate scores
  eval_outputs <- parallel::mclapply(seq_len(nrow(tasks)), mc.cores = mc_cores, function(i) {
    task <- dynutils::extract_row_to_list(tasks, i)

    # Fetch method outputs
    method_output <- method_outputs[[i]]

    if (any("try-error" %in% class(method_output))) {
      stop(method_output)
    }

    model <- method_output$model

    if (!is.null(model)) {
      # Calculate geodesic distances
      time0 <- Sys.time()
      model <- model %>% add_cell_waypoints(num_cells_selected = length(task$waypoint_cells))
      time1 <- Sys.time()
      time_cellwaypoints <- as.numeric(difftime(time1, time0, units = "sec"))
      df_cellwaypoints <- data_frame(time_cellwaypoints)
    } else {
      df_cellwaypoints <- NULL
    }

    # Calculate metrics
    metrics_summary <- calculate_metrics(task, model, calc_metrics)

    # Create summary statistics
    summary <- bind_cols(
      method_output$summary,
      df_cellwaypoints,
      metrics_summary
    )

    # Return the output
    lst(model, summary)
  })

  # Combine the different outputs in three lists/data frames
  models <- eval_outputs %>% map(~ .$model)
  summary <- eval_outputs %>% map_df(~ .$summary)

  # Calculate the final score
  score <- summary %>%
    summarise_at(metrics, funs(mean)) %>%
    as.matrix %>%
    as.vector %>%
    setNames(metrics)

  # Return extra information
  extras <- list(.summary = summary)

  # If output_model is true, add the model to the extras output
  if (output_model) {
    extras$.models <- models
  }

  # attach extras to score
  attr(score, "extras") <- extras

  # Return output
  score
}
