
#' Compute the helenie network score (any resemblances to real-life persons is purely coincidental)
#'
#' @param net1 the first network to compare
#' @param net2 the second network to compare
calculate_hele_network_score <- function(net1, net2) {
  # net1 <- tibble(from=c(1, 2, 2), to=c(2, 1, 1), directed=TRUE, length=1)
  # net2 <- tibble(from=c(1, 2, 2), to=c(2, 3, 4), directed=TRUE, length=1)
  # net1 <- tibble(from=1, to=1, directed=TRUE, length=1)
  # net2 <- tibble(from=1, to=1, directed=TRUE, length=1)
  # net1 <- tibble(from=c(1), to=c(1), directed=TRUE, length=1)
  # net2 <- tibble(from=c(1), to=c(2), directed=TRUE, length=1)

  net1 <- dynutils::simplify_milestone_network(net1)
  net2 <- dynutils::simplify_milestone_network(net2)

  net1 <- net1 %>% igraph::graph_from_data_frame(directed=T) %>% igraph::make_line_graph()
  net2 <- net2 %>% igraph::graph_from_data_frame(directed=T) %>% igraph::make_line_graph()

  nodes1 <- net1 %>% igraph::V() %>% as.numeric()
  nodes2 <- net2 %>% igraph::V() %>% as.numeric()

  net1 <- net1 %>% igraph::as_data_frame() %>% mutate(length=1)
  net2 <- net2 %>% igraph::as_data_frame() %>% mutate(length=1)

  calculate_robbie_network_score(net1, net2, nodes1, nodes2, normalize_weights=FALSE)
}


