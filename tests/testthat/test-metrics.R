context("Score metrics")
#
# test_that(paste0("Check lies network score"), {
#   net1 <- tibble(from=c(1), to=c(1), directed=TRUE, length=1)
#   net2 <- tibble(from=c(1), to=c(1), directed=TRUE, length=1)
#
#   expect_equal(dyneval:::calculate_lies_network_score(net1, net2), 1)
#
#   net1 <- tibble(from=c(1, 2, 2), to=c(2, 1, 1), directed=TRUE, length=1)
#   net2 <- tibble(from=c(1, 2, 2), to=c(2, 3, 4), directed=TRUE, length=1)
#
#   expect_equal(dyneval:::calculate_lies_network_score(net1, net2), 0.5)
#
#   net1 <- tibble(from=c(1), to=c(1), directed=TRUE, length=1)
#   net2 <- tibble(from=c(1), to=c(2), directed=TRUE, length=1)
#
#   expect_lt(dyneval:::calculate_lies_network_score(net1, net2), 1)
# })

test_that(paste0("Edge flip returns relevant results"), {
  library(tidyverse)
  topologies <- eval(formals(dyntoy:::generate_toy_milestone_network)$model)
  topologies <- topologies[topologies != "BA"]
  networks <- map(topologies, dyntoy:::generate_toy_milestone_network) %>% tibble(milestone_network = ., topology = topologies)
  design <- crossing(networks, networks)

  scores <- design %>% as.list() %>% pmap(function(milestone_network, milestone_network1, ...) {
    score <- dyneval:::calculate_edge_flip(milestone_network, milestone_network1)
    tibble(score=score)
  }) %>% bind_rows()
  design <- bind_cols(design, scores)

  expect_type(design$score, "double")
  expect_equal(design$score[design$topology == "simple_linear" & design$topology1 == "simple_linear"], 1)
  expect_equal(design$score[design$topology == "simple_linear" & design$topology1 == "linear"], 1)
  expect_lt(design$score[design$topology == "simple_linear" & design$topology1 == "bifurcating"], 1)
})



#
# library(ggplot2)
# design %>% ggplot() + geom_raster(aes(trajectory_type, trajectory_type1, fill=robbie_network_score))
#
#
# expand.grid(trajectory_types, trajectory_types) %>% map
#
# crossing(net1=networks, net2=networks)
