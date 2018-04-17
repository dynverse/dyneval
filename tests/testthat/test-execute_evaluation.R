context("Testing execute_evaluation")

test_that(paste0("Testing execute_evaluation"), {
  out <- execute_evaluation(
    tasks = dyntoy::toy_tasks[5,],
    method = dynmethods::description_random(),
    parameters = list(),
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
