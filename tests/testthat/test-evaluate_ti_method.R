context("Testing evaluate_ti_method")

custom_metric_1 <- function(dataset, model) {
  num_edges_dataset <- nrow(dataset$milestone_network)
  num_edges_model <- nrow(model$milestone_network)

  score <- 1 - abs(num_edges_dataset - num_edges_model) / num_edges_dataset

  ifelse(score < 0, 0, score)
}

custom_metric_2 <- function(dataset, model) {
  num_nodes_dataset <- length(dataset$milestone_ids)
  num_nodes_model <- length(model$milestone_ids)

  score <- 1 - abs(num_nodes_dataset - num_nodes_model) / num_nodes_dataset

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
  tmp <- tempfile()

  sink(tmp)
  out <- evaluate_ti_method(
    datasets = dyntoy::toy_datasets[5,],
    method = dynwrap::ti_random(),
    parameters = NULL,
    metrics = metrics,
    output_model = TRUE,
    extra_metrics = NULL,
    mc_cores = 2,
    verbose = TRUE
  )
  sink()
  unlink(tmp)

  score <- out$score
  summary <- out$summary
  models <- out$models

  expect_is(score, "numeric")

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
    datasets = dyntoy::toy_datasets[5,],
    method = dynwrap::ti_error(),
    parameters = list(),
    metrics = metrics,
    output_model = TRUE,
    extra_metrics = NULL,
    mc_cores = 2,
    verbose = FALSE
  )

  score <- out$score
  summary <- out$summary
  models <- out$models

  expect_is(score, "numeric")
  expect_true(all(score == c(0, 0, 1, 0)))

  expect_true(!is.null(summary$error[[1]]))

  expect_is(summary$correlation, "numeric")

  expect_is(summary$edge_flip, "numeric")

  expect_is(summary$rf_mse, "numeric")

  expect_is(summary$featureimp_cor, "numeric")

  expect_true(is.null(models[[1]]))
})


test_that(paste0("Testing evaluate_ti_method with identity"), {
  out <- evaluate_ti_method(
    datasets = dyntoy::toy_datasets[5,],
    method = dynwrap::ti_identity(),
    parameters = list(),
    metrics = metrics,
    output_model = TRUE,
    extra_metrics = NULL,
    mc_cores = 2,
    verbose = FALSE
  )

  score <- out$score
  summary <- out$summary
  models <- out$models

  expect_is(score, "numeric")
  expect_true(all(score - c(1, 1, 0, 1) < .01))

  expect_null(summary$error[[1]])

  expect_is(summary$correlation, "numeric")

  expect_is(summary$edge_flip, "numeric")

  expect_is(summary$rf_mse, "numeric")

  expect_is(summary$featureimp_cor, "numeric")

  expect_true(dynwrap::is_wrapper_with_trajectory(models[[1]]))
})
