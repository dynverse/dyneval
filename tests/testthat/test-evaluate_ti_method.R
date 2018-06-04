context("Testing evaluate_ti_method")

custom_metric_1 <- function(task, model) {
  num_edges_task <- nrow(task$milestone_network)
  num_edges_model <- nrow(model$milestone_network)

  score <- 1 - abs(num_edges_task - num_edges_model) / num_edges_task

  ifelse(score < 0, 0, score)
}

custom_metric_2 <- function(task, model) {
  num_nodes_task <- length(task$milestone_ids)
  num_nodes_model <- length(model$milestone_ids)

  score <- 1 - abs(num_nodes_task - num_nodes_model) / num_nodes_task

  ifelse(score < 0, 0, score)
}

metrics <- list(
  "correlation",
  "edge_flip",
  "rf_mse",
  "featureimp_cor",
  num_edges = custom_metric_1,
  num_nodes = custom_metric_2
)

test_that(paste0("Testing evaluate_ti_method with random"), {
  out <- evaluate_ti_method(
    tasks = dyntoy::toy_tasks[5,],
    method = dynmethods::ti_random(),
    parameters = NULL,
    metrics = metrics,
    output_model = TRUE,
    extra_metrics = NULL,
    mc_cores = 2,
    verbose = TRUE
  )

  expect_is(out, "numeric")

  summary <- attr(out, "extras")$.summary
  models <- attr(out, "extras")$.models

  expect_null(summary$error[[1]])

  expect_is(summary$correlation, "numeric")

  expect_is(summary$edge_flip, "numeric")

  expect_is(summary$rf_mse, "numeric")

  expect_is(summary$featureimp_cor, "numeric")

  expect_is(summary$num_edges, "numeric")

  expect_is(summary$num_nodes, "numeric")

  expect_true(dynwrap::is_wrapper_with_trajectory(models[[1]]))
})



test_that(paste0("Testing evaluate_ti_method with error"), {
  out <- evaluate_ti_method(
    tasks = dyntoy::toy_tasks[5,],
    method = dynmethods::ti_error(),
    parameters = list(),
    metrics = metrics,
    output_model = TRUE,
    extra_metrics = NULL,
    mc_cores = 2,
    verbose = TRUE
  )

  expect_is(out, "numeric")
  expect_true(all(out == c(0, 0, 1, 0)))

  summary <- attr(out, "extras")$.summary
  models <- attr(out, "extras")$.models

  expect_true(!is.null(summary$error[[1]]))

  expect_is(summary$correlation, "numeric")

  expect_is(summary$edge_flip, "numeric")

  expect_is(summary$rf_mse, "numeric")

  expect_is(summary$featureimp_cor, "numeric")

  expect_true(is.null(models[[1]]))
})


test_that(paste0("Testing evaluate_ti_method with identity"), {
  out <- evaluate_ti_method(
    tasks = dyntoy::toy_tasks[5,],
    method = dynmethods::ti_identity(),
    parameters = list(),
    metrics = metrics,
    output_model = TRUE,
    extra_metrics = NULL,
    mc_cores = 2,
    verbose = TRUE
  )

  expect_is(out, "numeric")
  expect_true(all(out - c(1, 1, 0, 1) < .01))

  summary <- attr(out, "extras")$.summary
  models <- attr(out, "extras")$.models

  expect_null(summary$error[[1]])

  expect_is(summary$correlation, "numeric")

  expect_is(summary$edge_flip, "numeric")

  expect_is(summary$rf_mse, "numeric")

  expect_is(summary$featureimp_cor, "numeric")

  expect_true(dynwrap::is_wrapper_with_trajectory(models[[1]]))
})
