context("Testing evaluate_ti_method")

test_that(paste0("Testing evaluate_ti_method with random"), {
  out <- evaluate_ti_method(
    tasks = dyntoy::toy_tasks[5,],
    method = dynmethods::ti_random(),
    parameters = NULL,
    metrics = c("correlation", "edge_flip", "rf_mse"),
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

  expect_true(dynwrap::is_wrapper_with_trajectory(models[[1]]))
})



test_that(paste0("Testing evaluate_ti_method with error"), {
  out <- evaluate_ti_method(
    tasks = dyntoy::toy_tasks[5,],
    method = dynmethods::ti_error(),
    parameters = list(),
    metrics = c("correlation", "edge_flip", "rf_mse"),
    output_model = TRUE,
    extra_metrics = NULL,
    mc_cores = 2,
    verbose = TRUE
  )

  expect_is(out, "numeric")
  expect_true(all(out == c(0, 0, 1)))

  summary <- attr(out, "extras")$.summary
  models <- attr(out, "extras")$.models

  expect_true(!is.null(summary$error[[1]]))

  expect_is(summary$correlation, "numeric")

  expect_is(summary$edge_flip, "numeric")

  expect_is(summary$rf_mse, "numeric")

  expect_true(is.null(models[[1]]))
})


test_that(paste0("Testing evaluate_ti_method with identity"), {
  out <- evaluate_ti_method(
    tasks = dyntoy::toy_tasks[5,],
    method = dynmethods::ti_identity(),
    parameters = list(),
    metrics = c("correlation", "edge_flip", "rf_mse"),
    output_model = TRUE,
    extra_metrics = NULL,
    mc_cores = 2,
    verbose = TRUE
  )

  expect_is(out, "numeric")
  expect_true(all(out - c(1, 1, 0) < .001))

  summary <- attr(out, "extras")$.summary
  models <- attr(out, "extras")$.models

  expect_null(summary$error[[1]])

  expect_is(summary$correlation, "numeric")

  expect_is(summary$edge_flip, "numeric")

  expect_is(summary$rf_mse, "numeric")

  expect_true(dynwrap::is_wrapper_with_trajectory(models[[1]]))
})
