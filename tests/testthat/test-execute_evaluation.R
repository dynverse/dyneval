context("Execute Evaluation")

test_that(paste0("Check execute_evaluation"), {
  # out <- execute_evaluation(tasks = dyntoy::toy_tasks[5,], method = dynmethods::description_random(), parameters = list(), timeout = 0, metrics = "auc_R_nx")
  # summary <- attr(out, "extras")$.summary
  # expect_is(summary$error[[1]], "error")

  # out <- execute_evaluation(tasks = dyntoy::toy_tasks[5,], method = dynmethods::description_random(), parameters = list(), timeout = 60, metrics = "auc_R_nx")
  # summary <- attr(out, "extras")$.summary
  # expect_null(summary$error[[1]])
  #
  out <- execute_evaluation(tasks = dyntoy::toy_tasks[5,], method = dynmethods::description_random(), parameters = list(), metrics = "auc_R_nx")
  summary <- attr(out, "extras")$.summary
  expect_null(summary$error[[1]])
})
