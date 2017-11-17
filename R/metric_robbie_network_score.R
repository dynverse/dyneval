

## Robie network score ----------------------------
score_map <- function(permutation, net, net_ref) {
  net_mapped <- net[permutation, permutation]

  1-sum(abs(net_mapped - net_ref))/(sum(net_mapped) + sum(net_ref))
}

complete_matrix <- function(mat, dim, fill=0) {
  mat <- rbind(mat, matrix(rep(0, ncol(mat) * (dim - nrow(mat))), ncol=ncol(mat)))
  mat <- cbind(mat, matrix(rep(0, nrow(mat) * (dim - ncol(mat))), nrow=nrow(mat)))
}

#' @importFrom reshape2 acast
get_adjacency <- function(net, nodes=unique(c(net$from, net$to))) {
  if(nrow(net) == 0) { # special case for circular
    newnet <- matrix(rep(0, length(nodes)))
    dimnames(newnet) <- list(nodes, nodes)
  } else {
    newnet <- net %>%
      mutate(from=factor(from, levels=nodes), to=factor(to, levels=nodes)) %>%
      reshape2::acast(from~to, value.var="length", fill=0, drop=FALSE, fun.aggregate=sum)
  }
  newnet + t(newnet)

}

#' Compute the robbie network score (any resemblances to real-life persons is purely coincidental)
#'
#' @param net1 the first network to compare
#' @param net2 the second network to compare
#' @param nodes1 nodes in `net1`, defaults to unique from and to
#' @param nodes2 nodes in `net2`, defaults to unique from and to
#' @param normalize_weights Whether to normalize the lengths of each network
#' @importFrom GA ga
calculate_robbie_network_score <- function(
  net1,
  net2,
  nodes1=unique(c(net1$from, net1$to)),
  nodes2=unique(c(net2$from, net2$to)),
  normalize_weights=TRUE
) {
  optimize_robbie_network_score(net1, net2, nodes1, nodes2, normalize_weights)@fitnessValue
}

optimize_robbie_network_score <- function(
  net1,
  net2,
  nodes1=unique(c(net1$from, net1$to)),
  nodes2=unique(c(net2$from, net2$to)),
  normalize_weights=TRUE
) {
  net <- get_adjacency(net1, nodes1)
  net_ref <- get_adjacency(net2, nodes2)
  net_ref <- complete_matrix(net_ref, nrow(net))

  if(nrow(net) == 1) {
    return(score_map(1, net, net_ref))
  }

  # normalize weights
  if(normalize_weights) {
    net <- net/sum(net)
    net_ref <- net_ref/sum(net_ref)
  }

  results <- GA::ga("permutation", score_map, net=net, net_ref=net_ref, min=rep(1, length(nodes1)), max=rep(length(nodes1), length(nodes1)), maxiter=100, monitor=FALSE, popSize=500, maxFitness = 1)

  results
}
