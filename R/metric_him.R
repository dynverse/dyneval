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
#' @importFrom dynwrap simplify_igraph_network
#'
#' @export
calculate_him <- function(
  net1,
  net2,
  simplify = TRUE
) {
  requireNamespace("netdist")

  # get the matched adjacencies
  adjacencies <- get_matched_adjacencies(
    net1,
    net2,
    simplify = simplify
  )

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

