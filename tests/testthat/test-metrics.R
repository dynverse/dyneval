# context("Score metrics")
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

#
#
# library(tidyverse)
# trajectory_types <- eval(formals(dyntoy:::generate_toy_milestone_network)$trajectory_type)
# networks <- map(trajectory_types, dyntoy:::generate_toy_milestone_network) %>% tibble(milestone_network = ., trajectory_type = trajectory_types)
# design <- crossing(networks, networks)
# design <- design %>% bind_cols(pbapply::pblapply(cl=8, seq_len(nrow(design)), function(row_id) {
#   print(row_id)
#   row <- dynutils::extract_row_to_list(design, row_id)
#
#   tibble(
#     node_edit_score = dyneval:::calculate_node_edit_score(row$milestone_network, row$milestone_network1),
#     node_edge_edit_score = dyneval:::calculate_node_edge_edit_score(row$milestone_network, row$milestone_network1)
#   )
# }) %>% bind_rows())
#
# library(ggplot2)
# design %>% ggplot() + geom_raster(aes(trajectory_type, trajectory_type1, fill=robbie_network_score))
#
#
# expand.grid(trajectory_types, trajectory_types) %>% map
#
# crossing(net1=networks, net2=networks)
