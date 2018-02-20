context("Execute Evaluation")

test_that(paste0("Check execute_evaluation"), {
  out <- execute_evaluation(tasks = dyntoy::toy_tasks[5,], method = dynmethods::description_random(), parameters = list(), metrics = "correlation")
  summary <- attr(out, "extras")$.summary
  expect_null(summary$error[[1]])
})
