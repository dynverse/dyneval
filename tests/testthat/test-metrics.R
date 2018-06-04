context("Score metrics")

library(purrr)

test_that(paste0("Edge flip returns relevant results"), {
  topologies <- eval(formals(dyntoy:::generate_toy_milestone_network)$model)
  topologies <- topologies[topologies != "BA"]
  networks <- map(topologies, dyntoy:::generate_toy_milestone_network) %>% tibble(milestone_network = ., topology = topologies)
  design <- crossing(networks, networks)

  scores <- design %>% as.list() %>% pmap(function(milestone_network, milestone_network1, ...) {
    score <- dyneval:::calculate_edge_flip(milestone_network, milestone_network1)
    tibble(score = score)
  }) %>% bind_rows()
  design <- bind_cols(design, scores)

  expect_type(design$score, "double")
  expect_equal(design$score[design$topology == "simple_linear" & design$topology1 == "simple_linear"], 1)
  expect_equal(design$score[design$topology == "simple_linear" & design$topology1 == "linear"], 1)
  expect_lt(design$score[design$topology == "simple_linear" & design$topology1 == "bifurcating"], 1)
})
