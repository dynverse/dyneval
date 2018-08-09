context("Branch overlap")

# generate a simple dataset for comparison
dataset <- dyntoy::generate_dataset(model = "multifurcating", num_features = 2, allow_tented_progressions = FALSE, normalise = FALSE)

test_that(paste0("Branch overlap returns relevant results"), {
  # when exact gold standard -> all scores = 1
  prediction <- dataset

  results <- calculate_branch_overlap(dataset, prediction)
  testthat::expect_true(all(
    results$recovery == 1,
    results$relevance == 1,
    results$F1 == 1
  ))

  # when only one edge remains in the prediction -> relevance still 1
  prediction <- dataset
  prediction$progressions <- prediction$progressions %>% filter(from == first(from) & to == first(to))

  results <- calculate_branch_overlap(dataset, prediction)
  testthat::expect_true(all(
    results$recovery < 1,
    results$relevance == 1,
    results$F1 < 1
  ))

  # when all cells are removed -> all scores = 0
  prediction <- dataset
  prediction$progressions <- prediction$progressions %>% filter(FALSE)

  results <- calculate_branch_overlap(dataset, prediction)
  testthat::expect_true(all(
    results$recovery == 0,
    results$relevance == 0,
    results$F1 == 0
  ))

  # when all progressions are shuffled -> still scores = 1
  prediction <- dataset
  prediction$progressions <- prediction$progressions %>% mutate(percentage = ifelse(percentage > 0 & percentage < 1, runif(n()), percentage))

  results <- calculate_branch_overlap(dataset, prediction)
  testthat::expect_true(all(
    results$recovery_branches == 1,
    results$relevance_branches == 1,
    results$F1_branches == 1
  ))
})
