context("Metric combinations")

test_that(paste0("calculate_combinations"), {
  testthat::expect_equal(calculate_arithmetic_mean(0.5, 0.6, 0.7), 0.6)
  testthat::expect_equal(calculate_harmonic_mean(0.5, 1, 1), 0.75)
  testthat::expect_equal(calculate_geometric_mean(0.001, 1, 1), 0.1)
})
