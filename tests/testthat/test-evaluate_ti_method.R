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
  "rf_nmse",
  "rf_rsq",
  "lm_mse",
  "lm_nmse",
  "lm_rsq",
  "featureimp_cor",
  "F1_branches",
  num_edges = custom_metric_1,
  num_nodes = custom_metric_2
)
metric_names <- ifelse(purrr::map_chr(metrics, is.character), metrics, names(metrics)) %>% unlist()

test_that(paste0("Testing evaluate_ti_method with random"), {
  tmp <- tempfile()

  sink(tmp)
  out <- evaluate_ti_method(
    dataset = dyntoy::toy_datasets[5, ],
    method = ti_random(),
    parameters = NULL,
    metrics = metrics,
    output_model = TRUE,
    map_fun = function(x, fun) parallel::mclapply(x, fun, mc.cores = 2),
    verbose = TRUE
  )
  sink()
  unlink(tmp)

  score <- as.list(out$summary)[metric_names] %>% unlist()
  summary <- out$summary
  models <- out$models

  expect_is(score, "numeric")

  expect_true(is.na(summary$error[[1]]))

  expect_is(summary$correlation, "numeric")

  expect_is(summary$edge_flip, "numeric")

  expect_is(summary$rf_mse, "numeric")

  expect_is(summary$featureimp_cor, "numeric")

  expect_is(summary$num_edges, "numeric")

  expect_is(summary$num_nodes, "numeric")

  expect_is(summary$F1_branches, "numeric")

  expect_true(dynwrap::is_wrapper_with_trajectory(models[[1]]))
})



test_that(paste0("Testing evaluate_ti_method with error"), {
  out <- evaluate_ti_method(
    dataset = dyntoy::toy_datasets[5, ],
    method = ti_error(),
    parameters = list(),
    metrics = metrics,
    output_model = TRUE,
    map_fun = map,
    verbose = FALSE
  )

  summary <- out$summary
  models <- out$models

  score <- as.list(out$summary)[metric_names] %>% unlist()

  expect_is(score, "numeric")

  expect_true(!is.na(summary$error[[1]]))

  expect_is(summary$correlation, "numeric")
  expect_is(summary$edge_flip, "numeric")
  expect_is(summary$rf_mse, "numeric")
  expect_is(summary$featureimp_cor, "numeric")

  expect_true(is.null(models[[1]]))
})


test_that(paste0("Testing evaluate_ti_method with identity"), {
  out <- evaluate_ti_method(
    dataset = dyntoy::toy_datasets[5, ],
    method = ti_identity(),
    parameters = list(),
    metrics = metrics,
    output_model = TRUE,
    map_fun = lapply,
    verbose = FALSE
  )

  summary <- out$summary
  models <- out$models

  score <- as.list(out$summary)[metric_names] %>% unlist()

  expect_is(score, "numeric")
  expect_true(all(score - c(1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1) < .01))

  expect_true(is.na(summary$error[[1]]))

  expect_is(summary$correlation, "numeric")

  expect_is(summary$edge_flip, "numeric")

  expect_is(summary$rf_mse, "numeric")

  expect_is(summary$featureimp_cor, "numeric")

  expect_true(dynwrap::is_wrapper_with_trajectory(models[[1]]))
})
