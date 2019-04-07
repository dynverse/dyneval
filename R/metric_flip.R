
#' Edge flip score
#'
#' @param net1 Network 1
#' @param net2 Network 2
#' @param return Whether to return only the `score` or the full output (`all`)
#' @param simplify Whether or not to simplify the networks
#' @param limit_flips Maximal number of flips to check
#' @param limit_combinations Maximal number of combinations to check
#'
#' @keywords metric
#'
#' @examples
#' net1 <- dyntoy::generate_milestone_network("linear")
#' net2 <- dyntoy::generate_milestone_network("bifurcating")
#' calculate_edge_flip(net1, net2)
#'
#' net1 <- dyntoy::generate_milestone_network("cyclic")
#' net2 <- dyntoy::generate_milestone_network("diverging_with_loops")
#' calculate_edge_flip(net1, net2)
#'
#' @importFrom dynwrap simplify_igraph_network
#'
#' @export
calculate_edge_flip <- function(
  net1,
  net2,
  return = c("score", "all"),
  simplify = TRUE,
  limit_flips = 5,
  limit_combinations = choose(25, 4)
) {
  return <- match.arg(return, c("score", "all"))

  # get the matched adjacencies
  adjacencies <- get_matched_adjacencies(
    net1,
    net2,
    simplify = simplify
  )

  adj1 <- adjacencies[[1]] > 0
  adj2 <- adjacencies[[2]] > 0

  # calculate the mapping which nodes are connected to which edges
  edge_membership1 <- calculate_edge_membership(adj1)

  # substract the number of edges
  edge_difference <- sum(adj2[lower.tri(adj2, diag = FALSE)]) - sum(adj1[lower.tri(adj1, diag = FALSE)])

  # calculate the possible edges which can be added and removed to net1
  possible_edge_additions <- which(!adj1[lower.tri(adj1, diag = FALSE)])
  possible_edge_removes <- which(adj1[lower.tri(adj1, diag = FALSE)])

  graph2 <- adj2 %>% get_undirected_graph() # used later to calculate isomorphism
  sorted_degrees2 <- adj2 %>% rowSums() %>% sort(method = "radix") # used later to compare degree distributions

  # prepare for looping over the number of edges which can be flipped
  found <- FALSE
  n_flips <- abs(edge_difference) - 2

  # determine upper bound
  upper_bound <- sum(adj1[lower.tri(adj1, diag = FALSE)]) + sum(adj2[lower.tri(adj2, diag = FALSE)]) - 2
  if (upper_bound <= 0) {
    # this is only possible when one of the networks does not have any edges
    upper_bound <- 1
  }

  # now loop over the number of edge flips, starting with the minimal
  while (!found & n_flips <= upper_bound) {
    n_flips <- n_flips + 2
    # print(glue::glue("flips: {n_flips}"))

    # limit number of checked flips
    if (n_flips > limit_flips) {
      found <- TRUE
      n_flips <- upper_bound
    } else {
      # calculate the number of additions and removes
      n_additions <- (n_flips + edge_difference)/2
      n_removes <- n_additions - edge_difference

      if (n_additions < 0 | n_removes < 0) {stop("Edge additions and removes should be integer and higher than 0")}

      if (choose(length(possible_edge_additions), n_additions) > limit_combinations | choose(length(possible_edge_removes), n_removes) > limit_combinations ) {
        found <- TRUE
        n_flips <- upper_bound
      } else {
        # create the matrix which contains in the columns all possible flips, with in the rows the edge_id which will be flipped
        edge_additions <- combn_nice(possible_edge_additions, n_additions)
        edge_removes <- combn_nice(possible_edge_removes, n_removes)

        if (n_additions > 0 & n_removes > 0) {
          edge_flips <- rbind(
            matrix(rep(edge_additions, ncol(edge_removes)), nrow = nrow(edge_additions)),
            t(edge_removes) %>% rep(each = ncol(edge_additions)) %>% matrix(ncol = nrow(edge_removes)) %>% t
          )
        } else if (n_additions > 0) {
          edge_flips <- edge_additions
        } else {
          edge_flips <- edge_removes
        }

        # cut the edge_flips, avoiding huge memory consumption
        grouping <- ceiling(seq_len(ncol(edge_flips)) / 1000)
        ngroups <- max(grouping)
        group_id <- 0

        # loop over each group of edge_flips
        while(!found & group_id < ngroups) {
          group_id <- group_id + 1
          edge_flips_group <- edge_flips[, grouping == group_id, drop = F]

          # generate matrix with in the columns each flip and in the rows the vector format of the new adjacency of net1
          edge_flip_vectors1 <- generate_edge_flip_vectors(edge_flips_group, adj1)
          degree_vectors1 <- t(edge_flip_vectors1) %*% edge_membership1

          # now check several metrics of the new adjacency matrix, from fastest to slowest
          # after each check, the flips which are not OK are removed (in the selected object)

          selected <- seq_len(nrow(degree_vectors1))

          # quick max check
          degree_max_check <- check_degrees_max(degree_vectors1, sorted_degrees2)
          if (any(degree_max_check)) {
            selected <- selected[degree_max_check]

            # quick min check
            degree_min_check <- check_degrees_min(degree_vectors1[selected, , drop = F], sorted_degrees2)
            if (any(degree_min_check)) {
              selected <- selected[degree_min_check]

              # now check sorted
              degree_sorted_check <- check_degrees_sorted(degree_vectors1[selected, , drop = F], sorted_degrees2)

              if (any(degree_sorted_check)) {
                selected <- selected[degree_sorted_check]

                # now check isomorphic
                adj1s <- map(selected, ~edge_flips_group[, ., drop = F]) %>% map(flip_adj, adj1)

                isomorphic <- map_lgl(adj1s, function(adj1) {
                  graph1 <- adj1 %>% get_undirected_graph()

                  igraph::is_isomorphic_to(graph1, graph2)
                })

                if(any(isomorphic)) {
                  found <- TRUE
                  # print("found!")
                  newadj1 <- first(adj1s[isomorphic])
                }
              }
            }
          }
        }
      }
    }
  }

  # return output
  if (!found) {
    stop("Couldn't map, this shouldn't happen! Something went wrong!")
  }

  score <- 1-n_flips/upper_bound

  if (return == "score") {
    score
  } else if (return == "all") {
    list(score = score, newadj1 = newadj1, oldadj1 = adj1)
  }
}


get_adjacency_lengths <- function(net, nodes = sort(unique(c(net$from, net$to)))) {
  if(nrow(net) == 0) { # special case for circular
    newnet <- matrix(rep(0, length(nodes)))
    dimnames(newnet) <- list(nodes, nodes)
  } else {
    newnet <- net %>%
      mutate(from = factor(from, levels = nodes), to = factor(to, levels = nodes)) %>%
      reshape2::acast(from~to, value.var = "length", fill = 0, drop = FALSE, fun.aggregate = sum)
  }
  (newnet + t(newnet))
}

get_adjacency <- function(net) {
  get_adjacency_lengths(net) != 0
}

get_undirected_graph <- function(adj) {
  adj %>% igraph::graph_from_adjacency_matrix(mode = "upper", weighted = TRUE) %>% igraph::as.undirected()
}


# add extra rows and columns to matrix
complete_matrix <- function(mat, dim, fill = 0) {
  mat <- rbind(
    mat,
    matrix(
      rep(fill, ncol(mat) * (dim - nrow(mat))),
      ncol = ncol(mat),
      dimnames = list(sample.int(dim-nrow(mat)), colnames(mat))
    )
  )
  mat <- cbind(
    mat,
    matrix(
      rep(fill, nrow(mat) * (dim - ncol(mat))),
      nrow = nrow(mat),
      dimnames = list(rownames(mat), sample.int(dim-nrow(mat)))
    )
  )

  mat
}

# get the matched adjacency matrices between two networks
get_matched_adjacencies <- function(net1, net2, simplify = TRUE) {
  if (simplify) {
    directed1 <- any(net1$directed)
    directed2 <- any(net2$directed)
    net1 <- net1 %>%
      rename(weight = length) %>%
      filter(!(from == to & weight == 0)) %>% # remove self loop edges with length 0
      igraph::graph_from_data_frame(directed = F) %>%
      dynwrap::simplify_igraph_network() %>%
      igraph::as_data_frame() %>%
      rename(length = weight) %>%
      mutate(directed = directed1) %>%
      insert_two_nodes_into_selfloop() %>%
      change_single_edge_into_double() %>%
      insert_one_node_into_duplicate_edges()
    net2 <- net2 %>%
      rename(weight = length) %>%
      filter(!(from == to & weight == 0)) %>% # remove self loop edges with length 0
      igraph::graph_from_data_frame(directed = F) %>%
      dynwrap::simplify_igraph_network() %>%
      igraph::as_data_frame() %>%
      rename(length = weight) %>%
      mutate(directed = directed2) %>%
      insert_two_nodes_into_selfloop() %>%
      change_single_edge_into_double() %>%
      insert_one_node_into_duplicate_edges()
  }

  adj1 <- get_adjacency_lengths(net1)
  adj2 <- get_adjacency_lengths(net2)

  # make the adjacency matrices have the same dimensions
  if (nrow(adj1) > nrow(adj2)) {
    adj2 <- complete_matrix(adj2, nrow(adj1), fill = 0)
  } else {
    adj1 <- complete_matrix(adj1, nrow(adj2), fill = 0)
  }

  lst(adj1, adj2)
}


# what nodes are connecting to which edges
#' @importFrom purrr invoke
calculate_edge_membership <- function(adj) {
  adj_mapper <- adj
  tri_ids <- seq_len(sum(lower.tri(adj_mapper, diag = FALSE)))
  adj_mapper[lower.tri(adj_mapper, diag = FALSE)] <- tri_ids
  adj_mapper[upper.tri(adj_mapper, diag = TRUE)] <- 0
  purrr::invoke(cbind, map(
    seq_len(nrow(adj_mapper)),
    function(i) {
      as.integer(tri_ids %in% c(adj_mapper[i, ], adj_mapper[, i]))
    }
  ))
}

# flip edges (in adjacency format)
flip_adj <- function(i, adj) {
  adj[lower.tri(adj, diag = FALSE)][i] <- 1-adj[lower.tri(adj, diag = FALSE)][i]
  adj[upper.tri(adj, diag = TRUE)] <- 0
  adj + t(adj)
}

# flip edges (in edge vector format)
generate_edge_flip_vectors <- function(edge_flips, adj) {
  adjv <- adj[lower.tri(adj, diag = FALSE)]
  edge_flips %>% apply(2, function(x) {adjv[x] <- 1-adjv[x];adjv})
}

# check whether the maximal degree matches
check_degrees_max <- function(degree_vectors1, sorted_degrees2) {
  degree_vectors1 %>% apply(1, max) %>% {. == sorted_degrees2[length(sorted_degrees2)]}
}

# check whether the minimal degree matches
check_degrees_min <- function(degree_vectors1, sorted_degrees2) {
  degree_vectors1 %>% apply(1, min) %>% {. == sorted_degrees2[1]}
}

# check whether the degree distribution matches
check_degrees_sorted <- function(degree_vectors1, sorted_degrees2) {
  degree_vectors1 %>% apply(1, sort) %>% t %>% apply(1, function(x) all(x == sorted_degrees2))
}

# nice combn, which doesn't just use seq_len(x) if x is an integer -___-
#' @importFrom utils combn
combn_nice  <- function(x, m) {
  if(length(x) > 1 || m == 0) {
    utils::combn(x, m)
  } else {
    matrix(x, nrow = 1, ncol = 1)
  }
}

insert_two_nodes_into_selfloop <- function(df) {
  ix <- df$from == df$to
  new <- map_df(which(ix), function(i) {
    n <- df$from[[i]]
    l <- df$length[[i]]
    d <- df$directed[[i]]
    newn <- paste0(dynutils::random_time_string(), seq_len(2))
    data_frame(
      from = c(n, newn),
      to = c(newn, n),
      length = l/3,
      directed = d
    )
  })
  bind_rows(
    df[!ix,],
    new
  )
}

# df <- tibble(from = c("M1", "M2", "M2"), to = c("M2", "M1", "M1"), length = 1, directed = T)
insert_one_node_into_duplicate_edges <- function(df) {
  ix <- paste0(df$from, "#", df$to) %in% names(which(table(paste0(df$from, "#", df$to)) >= 2))
  new <- map_df(which(ix), function(i) {
    n <- df$from[[i]]
    t <- df$to[[i]]
    l <- df$length[[i]]
    d <- df$directed[[i]]
    newn <- paste0(dynutils::random_time_string(), seq_len(1))
    data_frame(
      from = c(n, newn),
      to = c(newn, t),
      length = l/2,
      directed = d
    )
  })
  bind_rows(
    df[!ix,],
    new
  )
}

change_single_edge_into_double <- function(df) {
  if (nrow(df) == 1 && df$from[[1]] != df$to[[1]]) {
    data_frame(
      from = c("a", "b"),
      to = c("b", "c"),
      length = df$length/2,
      directed = df$directed[[1]]
    )
  } else {
    df
  }
}
