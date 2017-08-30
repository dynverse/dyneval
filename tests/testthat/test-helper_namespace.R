context("Namespace helper")

test_that("Assign in namespace", {
  test_message <- "test set seed"
  orig_set_seed <- base::set.seed
  dyneval:::my_assignin_namespace("set.seed", function(i) test_message, ns = "base", envir = .BaseNamespaceEnv)

  expect_equal( set.seed(1), test_message )

  dyneval:::my_assignin_namespace("set.seed", orig_set_seed, ns = "base", envir = .BaseNamespaceEnv)
})
