#' Used for wrapping an evaluation function around a TI method
#'
#' @importFrom smoof makeSingleObjectiveFunction makeMultiObjectiveFunction
#' @export
make_obj_fun <- function(method, noisy = F, suppress_output = T,
                         metrics = c("mean_R_nx", "auc_R_nx", "Q_global", "Q_local", "correlation", "isomorphic", "ged")) {
  # Use different makefunction if there are multiple metrics versus one
  if (length(metrics) > 1) {
    make_fun <- function(...) makeMultiObjectiveFunction(..., n.objectives = length(metrics))
  } else {
    make_fun <- makeSingleObjectiveFunction
  }

  # Wrap the method function in an evaluation function
  make_fun(
    name = "TItrain",
    vectorized = F,
    minimize = rep(F, length(metrics)),
    noisy = noisy,
    has.simple.signature = F,
    par.set = method$par_set,
    fn = function(x, tasks)
      execute_evaluation(
        tasks = tasks,
        method = method,
        parameters = x,
        metrics = metrics,
        suppress_output = suppress_output))
}

#' @export
impute_y_fun <- function(num_objectives) {
  function(x, y, opt.path, ...) {
    val <- rep(-1, num_objectives)
    attr(val, "extras") <- list(.summary = NA)
    val
  }
}

#' @importFrom future future value plan
execute_evaluation <- function(tasks, method, parameters,
                               metrics = c("mean_R_nx", "auc_R_nx", "Q_global", "Q_local", "correlation", "isomorphic", "ged"),
                               suppress_output = T) {
  # Run the method on each of the tasks
  method_futures <- future(
    {
      # Load required namespaces
      for (pack in method$package_required) {
        suppressMessages(do.call(requireNamespace, list(pack)))
      }

      # Disable seed setting
      dyneval:::my_assignin_namespace("set.seed", function(i) {}, ns = "base", envir = .BaseNamespaceEnv)

      # Run method on each task
      lapply(seq_len(nrow(tasks)), function(i) {
        task <- dyneval:::extract_row_to_list(tasks, i)
        dyneval:::run_method(task, method, parameters, suppress_output = suppress_output)
      })
    },
    globals = c("tasks", "method", "parameters", "suppress_output", "my_set_seed"),
    packages = c("dyneval", method$package_load),
    evaluator = plan("multisession")
  )
  method_outputs <- value(method_futures)

  # Calculate scores
  summary_outs <- lapply(seq_len(nrow(tasks)), function(i) {
    task <- extract_row_to_list(tasks, i)

    # Fetch method outputs
    method_output <- method_outputs[[i]]
    model <- method_output$model

    # Calculate geodesic distances
    time0 <- Sys.time()
    model$geodesic_dist <- compute_emlike_dist(model)
    time1 <- Sys.time()
    time_geodesic <- as.numeric(difftime(time1, time0, units = "sec"))

    # Calculate metrics
    metrics_output <- calculate_metrics(task, model, metrics)

    # Create summary statistics
    summary <- data.frame(
      method_name = method$name,
      method_short_name = method$short_name,
      task_id = task$id,
      method_output$summary,
      time_geodesic = time_geodesic,
      metrics_output$summary,
      stringsAsFactors = F,
      check.names = F
    )

    # Return the output
    lst(model, summary)
  })

  # Combine the different outputs in three lists/data frames
  models <- summary_outs %>% purrr::map(~ .$model)
  summary <- summary_outs %>% purrr::map_df(~ .$summary)

  # Calculate the final score
  score <- summary %>% summarise_at(metrics, funs(mean)) %>% as.matrix %>% as.vector %>% setNames(metrics)

  # Return extra information
  attr(score, "extras") <- list(.models = models, .summary = summary)

  # Return output
  score
}
