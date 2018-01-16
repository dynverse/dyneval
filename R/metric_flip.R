get_adjacency_lengths <- function(net, nodes=sort(unique(c(net$from, net$to)))) {
  if(nrow(net) == 0) { # special case for circular
    newnet <- matrix(rep(0, length(nodes)))
    dimnames(newnet) <- list(nodes, nodes)
  } else {
    newnet <- net %>%
      mutate(from=factor(from, levels=nodes), to=factor(to, levels=nodes)) %>%
      reshape2::acast(from~to, value.var="length", fill=0, drop=FALSE, fun.aggregate=sum)
  }
  (newnet + t(newnet)) > 0
}

get_adjacency <- function(net) {
  get_adjacency_lengths(net) != 0
}

get_undirected_graph <- function(adj) {
  adj %>% igraph::graph_from_adjacency_matrix(mode="upper", weighted=TRUE) %>% igraph::as.undirected()
}

# check whether the degree distributions are equal
calculate_degree_equivalence <- function(adj1, sorted_degrees2) {
  degrees1 <- adj1 %>% rowSums()
  if(max(degrees1) == sorted_degrees2[length(sorted_degrees2)]) {
    all(sort(degrees1) == sorted_degrees2)
  }
  FALSE
}

# add extra rows and columns to matrix
complete_matrix <- function(mat, dim, fill=0) {
  mat <- rbind(mat, matrix(rep(fill, ncol(mat) * (dim - nrow(mat))), ncol=ncol(mat)))
  mat <- cbind(mat, matrix(rep(fill, nrow(mat) * (dim - ncol(mat))), nrow=nrow(mat)))
}

# what nodes are connecting to which edges
calculate_edge_membership <- function(adj) {
  adj_mapper <- adj
  tri_ids <- seq_len(sum(lower.tri(adj_mapper, diag = FALSE)))
  adj_mapper[lower.tri(adj_mapper, diag = FALSE)] <- tri_ids
  adj_mapper[upper.tri(adj_mapper, diag = TRUE)] <- 0
  map(seq_len(nrow(adj_mapper)), function(i) as.integer(tri_ids %in% c(adj_mapper[i, ], adj_mapper[, i]))) %>% invoke(cbind, .)
}

# flip edges (in adjacency format)
flip_adj <- function(i, adj) {
  adj[lower.tri(adj, diag=FALSE)][i] <- 1-adj[lower.tri(adj, diag=FALSE)][i]
  adj[upper.tri(adj, diag=TRUE)] <- 0
  adj + t(adj)
}

# flip edges (in edge vector format)
generate_edge_flip_vectors <- function(edge_flips, adj) {
  adjv <- adj[lower.tri(adj, diag=FALSE)]
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
combn_nice  <- function(x, m) {
  if(length(x) > 1 || m == 0) {
    combn(x, m)
  } else {
    matrix(x, nrow=1, ncol=1)
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

#' Edge flip score
#'
#' @param net1 Network 1
#' @param net2 Network 2
#' @param return Whether to return only the `score` or the full output (`all`)
#' @examples
#' net1 <- dyntoy:::generate_toy_milestone_network("linear_long")
#' net2 <- dyntoy:::generate_toy_milestone_network("trifurcating")
#' calculate_edge_flip(net1, net2)
calculate_edge_flip <- function(net1, net2, return=c("score", "all"), simplify=TRUE) {
  return <- match.arg(return, c("score", "all"))

  if (simplify) {
    directed1 <- any(net1$directed)
    directed2 <- any(net2$directed)
    net1 <- net1 %>%
      rename(weight = length) %>%
      igraph::graph_from_data_frame(directed = directed1) %>%
      dynutils::simplify_igraph_network() %>%
      igraph::as_data_frame() %>%
      rename(length = weight) %>%
      mutate(directed = directed1) %>%
      insert_two_nodes_into_selfloop()
    net2 <- net2 %>%
      rename(weight = length) %>%
      igraph::graph_from_data_frame(directed = directed2) %>%
      dynutils::simplify_igraph_network() %>%
      igraph::as_data_frame() %>%
      rename(length = weight) %>%
      mutate(directed = directed2) %>%
      insert_two_nodes_into_selfloop()

  }

  # calculate edge flip can't handle directed graphs
  net1$directed <- FALSE
  net2$directed <- FALSE

  adj1 <- get_adjacency(net1)
  adj2 <- get_adjacency(net2)

  if (nrow(adj1) > nrow(adj2)) {
    adj2 <- complete_matrix(adj2, nrow(adj1), fill = F)
  } else {
    adj1 <- complete_matrix(adj1, nrow(adj2), fill = F)
  }

  edge_membership1 <- calculate_edge_membership(adj1)

  edge_difference <- sum(adj2[lower.tri(adj2, diag = FALSE)]) - sum(adj1[lower.tri(adj1, diag = FALSE)])

  possible_edge_additions <- which(!adj1[lower.tri(adj1, diag = FALSE)])
  possible_edge_removes <- which(adj1[lower.tri(adj1, diag = FALSE)])

  graph2 <- adj2 %>% get_undirected_graph() # used later to calculate isomorphism
  sorted_degrees2 <- adj2 %>% rowSums() %>% sort(method="radix") # used later to compare degree distributions

  found <- FALSE
  n_flips <- abs(edge_difference) - 2

  max_flips <- sum(adj1[lower.tri(adj1, diag = FALSE)]) + sum(adj2[lower.tri(adj2, diag = FALSE)])

  baseline <- max_flips

  while (!found & n_flips <= max_flips) {
    n_flips <- n_flips + 2
    # print(glue::glue("flips: {n_flips}"))

    n_additions <- (n_flips + edge_difference)/2
    n_removes <- n_additions - edge_difference

    if (n_additions < 0 | n_removes < 0) {stop("Edge additions and removes should be integer and higher than 0")}

    edge_additions <- combn_nice(possible_edge_additions, n_additions)
    edge_removes <- combn_nice(possible_edge_removes, n_removes)

    if (n_additions > 0 & n_removes > 0) {
      edge_flips <- rbind(
        matrix(rep(edge_additions, ncol(edge_removes)), nrow=nrow(edge_additions)),
        t(edge_removes) %>% rep(each=ncol(edge_additions)) %>% matrix(ncol=nrow(edge_removes)) %>% t
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

    while(!found & group_id < ngroups) {
      group_id <- group_id + 1
      edge_flips_group <- edge_flips[, grouping == group_id, drop=F]

      edge_flip_vectors1 <- generate_edge_flip_vectors(edge_flips_group, adj1)
      degree_vectors1 <- t(edge_flip_vectors1) %*% edge_membership1

      # quick max check
      selected <- seq_len(nrow(degree_vectors1))

      degree_max_check <- check_degrees_max(degree_vectors1, sorted_degrees2)
      if (any(degree_max_check)) {
        selected <- selected[degree_max_check]

        # quick min check
        degree_min_check <- check_degrees_min(degree_vectors1[selected, , drop=F], sorted_degrees2)
        if (any(degree_min_check)) {
          selected <- selected[degree_min_check]

          # now check sorted
          degree_sorted_check <- check_degrees_sorted(degree_vectors1[selected, , drop=F], sorted_degrees2)

          if (any(degree_sorted_check)) {
            selected <- selected[degree_sorted_check]

            # now check isomorphic
            adj1s <- map(selected, ~edge_flips_group[, ., drop=F]) %>% map(flip_adj, adj1)

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

  if (!found) {
    stop("Couldn't map, this shouldn't happen! Something went wrong with the 'exact' algorithm!")
  } else {
    score <- 1-n_flips/baseline

    if (return == "score") {
      score
    } else if (return == "all") {
      list(score = score, newadj1 = newadj1, oldadj1 = adj1)
    }
  }
}

#' Plotting edge flips
#' @param oldadj Zouter needs to
#' @param newadj fill in these fields
#'
#' @importFrom tidygraph as_tbl_graph activate
#' @importFrom ggraph ggraph geom_edge_fan geom_edge_loop geom_node_label scale_edge_colour_manual
#' @importFrom cowplot theme_nothing
#'
#' @examples
#' net1 <- dyntoy:::generate_toy_milestone_network("linear_long")
#' net2 <- dyntoy:::generate_toy_milestone_network("trifurcating")
#' result <- calculate_edge_flip(net1, net2)
#' plot_edge_flips(result$newadj1, result$oldadj1)
plot_edge_flips <- function(oldadj, newadj) {
  # names are used for generating the network, make sure they are present and unique
  names <- seq_len(nrow(oldadj))
  dimnames(oldadj) <- list(names, names)
  dimnames(newadj) <- list(names, names)

  oldnet <- oldadj %>%
    reshape2::melt(varnames=c("from", "to"), value.name = "old") %>%
    mutate(old = old == 1) %>%
    filter(from >= to)

  newnet <- newadj %>%
    reshape2::melt(varnames=c("from", "to"), value.name = "new") %>%
    mutate(new = new == 1) %>%
    filter(from >= to)

  types <- tibble(old = c(F, F, T, T), new=c(F, T, F, T), type=c("irrelevant", "gained", "lost", "stayed"))
  net <- left_join(oldnet, newnet, by=c("from", "to")) %>%
    left_join(types, by=c("old", "new"))

  graph <- net %>% as_tbl_graph(directed=FALSE)

  graph <- graph %>%
    activate(edges) %>%
    mutate(weight = as.numeric(type != "irrelevant")) %>%
    arrange(weight)

  ggraph(graph, layout="fr") +
    geom_edge_fan(aes(color=type)) +
    geom_edge_loop(aes(colour = type)) +
    geom_node_label(aes(label=name)) +
    scale_edge_colour_manual(values=c(irrelevant = "#DDDDDD", lost="#FF4136", stayed="#0074D9", gained="#2ECC40")) +
    cowplot::theme_nothing()
}
