
#' Compute the lies network score (any resemblances to real-life persons is purely coincidental)
#'
#' @param net1 the first network to compare
#' @param net2 the second network to compare
calculate_lies_network_score <- function(
  net1,
  net2
) {
  nodes1=unique(c(net1$from, net1$to))
  nodes2=unique(c(net2$from, net2$to))
  if(length(nodes1) < length(nodes2)) {
    nodes3 <- nodes2
    net3 <- net2
    nodes2 <- nodes1
    net2 <- net1
    nodes1 <- nodes3
    net1 <- net3
  }

  if(any(nodes1 %in% nodes2)) {
    net1 <- net1 %>% mutate(from=paste0("__A", from), to=paste0("__A", to))
    net2 <- net2 %>% mutate(from=paste0("__B", from), to=paste0("__B", to))
    nodes1 <- unique(c(net1$from, net1$to))
    nodes2 <- unique(c(net2$from, net2$to))
  }

  net1$length <- net1$length / sum(net1$length)
  net2$length <- net2$length / sum(net2$length)

  adj1 <- get_adjacency(net1)
  adj2 <- get_adjacency(net2)

  nodes1 <- unique(c(net1$from, net1$to))
  nodes2 <- unique(c(net2$from, net2$to))

  map_grid <- expand.grid(node1=nodes1, node2=nodes2)
  map_size <- nrow(map_grid)

  number2binary = function(number, noBits) {
    binary_vector = rev(as.numeric(intToBits(number)))
    if(missing(noBits)) {
      return(binary_vector)
    } else {
      binary_vector[-(1:(length(binary_vector) - noBits))]
    }
  }

  if(map_size <= 11) {
    map_n <- 1:2^map_size
    combinations <- map(map_n, number2binary, map_size) %>% map(as.logical)
    scores <- map_dbl(combinations, score_map_lies, map_grid, net1, net2, nodes1, nodes2)
    score <- max(scores)
  } else {
    permutation2map <- function(perm) {
      (1:nrow(map_grid)) %in% ((0:(length(nodes1)-1)) * length(nodes1) + perm)
    }
    results <- optimize_robbie_network_score(net1, net2)
    perm <- results@solution
    initial <- apply(perm, 1, permutation2map) %>% t

    results <- GA::ga("binary", score_map_lies, map_grid, net1, net2, nodes1, nodes2, nBits = map_size, maxiter=50, monitor=FALSE, popSize=50, run=5, suggestions=initial, maxFitness=1)
    score <- results@fitnessValue
  }

  score
}


score_map_lies <- function(map, map_grid, net1, net2, nodes1, nodes2) {
  mapper <- map_grid[as.logical(map), ] %>% igraph::graph_from_data_frame(vertices=c(nodes1, nodes2)) %>% igraph::components() %>% .$membership

  net1_mapped <- net1 %>% mutate(from=as.character(mapper[from]), to=as.character(mapper[to]))
  net2_mapped <- net2 %>% mutate(from=as.character(mapper[from]), to=as.character(mapper[to]))

  adj1 <- get_adjacency(net1_mapped, unique(as.character(mapper)))
  adj2 <- get_adjacency(net2_mapped, unique(as.character(mapper)))

  adj_diff <- 1-sum(abs(adj1 - adj2))/((sum(adj1) + sum(adj2)))

  map_cost <- min(c(length(unique(mapper[nodes1]))/length(nodes1),length(unique(mapper[nodes2]))/length(nodes2)))

  edge_diff <- 1-(nrow(net1_mapped) - nrow(net2_mapped))/(max(nrow(net1_mapped), nrow(net2_mapped)))

  map_cost * adj_diff * edge_diff
}

