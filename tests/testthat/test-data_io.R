context("Data IO")

.datasets_location = "../dyngen/results/4/"

test_that("Loading datasets", {
  num_datasets <- 2
  datasets <- load_datasets(mc_cores = 1, num_datasets = num_datasets)

  expect_that( is_tibble(datasets), is_true() )
  expect_that( nrow(datasets), equals(num_datasets) )

  required_cols <- c("id", "cell_ids", "milestone_ids", "milestone_network", "milestone_percentages", "progressions", "counts", "geodesic_dist", "special_cells")
  expect_that( all(required_cols %in% colnames(datasets)), is_true() )
})
