#' Running an evaluation of a method on a set of tasks with a set of parameters
#'
#' @inheritParams calculate_metrics
#' @param tasks The tasks on which to evaluate.
#' @param method The method to evaluate.
#' @param parameters The parameters to evaluate with.
#' @param timeout Kill execution after a given amount of time.
#' @param debug_timeout Setting debug to \code{TRUE} will avoid running the method in a separate R session
#'   using \code{\link[dynutils]{eval_with_timeout}} and run the method directly. Note that the timeout functionality
#'   will not work when \code{debug} is \code{TRUE}.#' @param metrics which metrics to use;
#'   see \code{\link{calculate_metrics}} for a list of which metrics are available.
#' @param output_model Whether or not the model will be outputted.
#'   If this is a character string, it will save the model in the requested folder.
#' @param error_score The aggregated score a method gets if it produces errors.
#'
#' @export
#' @importFrom dynmethods execute_method
execute_evaluation <- function(tasks, method, parameters, metrics, timeout, debug_timeout = FALSE, output_model = TRUE, error_score = 0) {
  method_outputs <- dynmethods::execute_method(tasks = tasks, method = method, parameters = parameters, timeout = timeout, debug_timeout = debug_timeout)

  # Calculate scores
  eval_outputs <- lapply(seq_len(nrow(tasks)), function(i) {
    task <- dynutils::extract_row_to_list(tasks, i)

    # Fetch method outputs
    method_output <- method_outputs[[i]]
    model <- method_output$model

    if (!is.null(model)) {
      # Calculate geodesic distances
      time0 <- Sys.time()
      model$geodesic_dist <- dynutils::compute_emlike_dist(model)
      time1 <- Sys.time()
      time_geodesic <- as.numeric(difftime(time1, time0, units = "sec"))

      # Calculate metrics
      metrics_output <- calculate_metrics(task, model, metrics)

      # Create summary statistics
      summary <- bind_cols(
        method_output$summary,
        data_frame(time_geodesic),
        metrics_output$summary
      )
    } else {
      summary <- method_output$summary %>%
        mutate_at(metrics, ~ error_score)
    }

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

  # If output_model is a boolean, it decides on whether
  # to add the model to the extras output
  if (is.logical(output_model) && output_model) {
    extras$.models <- models
  }

  # If output_model is a character, write the model
  # to the given destination
  if (is.character(output_model)) {
    if (!dir.exists(output_model)) {
      dir.create(output_model)
    }
    filename <- paste0(
      output_model,
      ifelse(grepl("/$", output_model), "", "/"),
      dynutils::random_time_string(), ".rds"
    )
    saveRDS(models, file = filename)
    extras$.models_file <- filename
  }

  # attach extras to score
  attr(score, "extras") <- extras

  # Return output
  score
}
