#' netdist scores
#'
#' @param net1 Network 1
#' @param net2 Network 2
#' @param simplify Whether or not to simplify the networks
#'
#' @examples
#' net1 <- dyntoy::generate_milestone_network("linear")
#' net2 <- dyntoy::generate_milestone_network("bifurcating")
#' calculate_him(net1, net2)
#'
#' net1 <- dyntoy::generate_milestone_network("cyclic")
#' net2 <- dyntoy::generate_milestone_network("cyclic")
#' calculate_him(net1, net2)
#'
#' @keywords metric
#'
#' @importFrom dynwrap simplify_igraph_network
#'
#' @export
calculate_him <- function(
  net1,
  net2,
  simplify = TRUE
) {
  requireNamespace("netdist", quietly = TRUE)

  # get the matched adjacencies
  adjacencies <- get_matched_adjacencies(
    net1,
    net2,
    simplify = simplify
  )

  # return 0 when the largest length of either graph is 0
  if (max(adjacencies[[2]]) == 0 || max(adjacencies[[1]]) == 0) {
    return(0)
  }

  netdist <- netdist::netdist(
    adjacencies[[1]] / sum(adjacencies[[1]]),
    adjacencies[[2]] / sum(adjacencies[[2]]),
    "HIM",
    n.cores = 1,
    ga = 0.1
  )["HIM"]

  netdist[netdist < 0] <- 0

  1 - netdist
}

