#' Used for wrapping an evaluation function around a TI method
#'
#' @inheritParams execute_evaluation
#' @param noisy whether or not the metric is noisy or not
#'
#' @importFrom smoof makeSingleObjectiveFunction makeMultiObjectiveFunction
#' @export
make_obj_fun <- function(method, metrics, noisy = FALSE) {
  # Use different makefunction if there are multiple metrics versus one
  if (length(metrics) > 1) {
    make_fun <- function(...) makeMultiObjectiveFunction(..., n.objectives = length(metrics))
  } else {
    make_fun <- makeSingleObjectiveFunction
  }

  # Wrap the method function in an evaluation function
  make_fun(
    name = "TItrain",
    vectorized = FALSE,
    minimize = rep(FALSE, length(metrics)),
    noisy = noisy,
    has.simple.signature = FALSE,
    par.set = method$par_set,
    fn = function(x, tasks, output_model, timeout) {
      execute_evaluation(
        tasks = tasks,
        method = method,
        parameters = x,
        metrics = metrics,
        timeout = timeout,
        output_model = output_model
      )
    }
  )
}

#' Running an evaluation of a method on a set of tasks with a set of parameters
#'
#' @inheritParams calculate_metrics
#' @param tasks The tasks on which to evaluate.
#' @param method The method to evaluate.
#' @param parameters The parameters to evaluate with.
#' @param timeout Kill execution after a given amount of time.
#' @param debug_timeout Setting debug to \code{TRUE} will avoid running the method in a separate R session
#'   using \code{\link[dynutils]{eval_with_timeout}} and run the method directly. Note that the timeout functionality
#'   will not work when \code{debug} is \code{TRUE}.#' @param metrics which metrics to use;
#'   see \code{\link{calculate_metrics}} for a list of which metrics are available.
#' @param output_model Whether or not the model will be outputted.
#'   If this is a character string, it will save the model in the requested folder.
#' @param error_score The aggregated score a method gets if it produces errors.
#'
#' @export
#' @importFrom netdist gdd net_emd
#' @importFrom dynmethods execute_method
execute_evaluation <- function(tasks, method, parameters, metrics, timeout, debug_timeout = FALSE, output_model = TRUE, error_score = 0) {
  method_outputs <- dynmethods::execute_method(tasks = tasks, method = method, parameters = parameters, timeout = timeout, debug_timeout = debug_timeout)

  # Calculate scores
  eval_outputs <- lapply(seq_len(nrow(tasks)), function(i) {
    task <- dynutils::extract_row_to_list(tasks, i)

    # Fetch method outputs
    method_output <- method_outputs[[i]]
    model <- method_output$model

    if (!is.null(model)) {
      # Calculate geodesic distances
      time0 <- Sys.time()
      model$geodesic_dist <- dynutils::compute_emlike_dist(model)
      time1 <- Sys.time()
      time_geodesic <- as.numeric(difftime(time1, time0, units = "sec"))

      # Calculate metrics
      metrics_output <- calculate_metrics(task, model, metrics)

      # Create summary statistics
      summary <- bind_cols(
        method_output$summary,
        data_frame(time_geodesic),
        metrics_output$summary
      )
    } else {
      summary <- method_output$summary %>%
        mutate_at(metrics, function(x) error_score)
    }

    # Return the output
    lst(model, summary)
  })

  # Combine the different outputs in three lists/data frames
  models <- eval_outputs %>% map(~ .$model)
  summary <- eval_outputs %>% map_df(~ .$summary)

  # Calculate the final score
  score <- summary %>%
    summarise_at(metrics, funs(mean)) %>%
    as.matrix %>%
    as.vector %>%
    setNames(metrics)

  # Return extra information
  extras <- list(.summary = summary)

  # If output_model is a boolean, it decides on whether
  # to add the model to the extras output
  if (is.logical(output_model) && output_model) {
    extras$.models <- models
  }

  # If output_model is a character, write the model
  # to the given destination
  if (is.character(output_model)) {
    if (!dir.exists(output_model)) {
      dir.create(output_model)
    }
    filename <- paste0(
      output_model,
      ifelse(grepl("/$", output_model), "", "/"),
      dynutils::random_time_string(), ".rds"
    )
    saveRDS(models, file = filename)
    extras$.models_file <- filename
  }

  # attach extras to score
  attr(score, "extras") <- extras

  # Return output
  score
}

#' Calculate the performance of a model with respect to a task
#'
#' @param task the original task
#' @param model the predicted model
#' @param metrics which metrics to evaluate:
#' \enumerate{
#'   \item Coranking of geodesic distances: \code{"mean_R_nx"}, \code{"auc_R_nx"}, \code{"Q_local"}, \code{"Q_global"}
#'   \item Spearman correlation of geodesic distances: \code{"correlation"}
#'   \item Isomorphism of the two networks: \code{"isomorphic"}
#'   \item GEDEVO Graph Edit Distance: \code{"ged"}
#'   \item Earth Mover's Distance on orbit counts: \code{"net_emd"}
#'   \item Genetic Algorithm for aligning small graphs: \code{"robbie_network_score"}
#'   \item Mantel test p-value: \code{"mantel_pval"}
#' }
#'
#' @importFrom igraph is_isomorphic_to graph_from_data_frame
#'
#' @export
calculate_metrics <- function(task, model, metrics) {
  summary_list <- list()

  if (any(c("mean_R_nx", "auc_R_nx", "Q_local", "Q_global", "correlation", "mantel_pval") %in% metrics)) {
    # compute coranking
    time0 <- Sys.time()
    coranking <- compute_coranking(task$geodesic_dist, model$geodesic_dist)
    summary_list <- c(summary_list, coranking$summary)
    time1 <- Sys.time()
    summary_list$time_coranking <- as.numeric(difftime(time1, time0, units = "sec"))

    # compute corrrelation
    time0 <- Sys.time()
    summary_list$correlation <- cor(task$geodesic_dist %>% as.vector, model$geodesic_dist %>% as.vector, method="spearman")
    time1 <- Sys.time()
    summary_list$time_correlation <- as.numeric(difftime(time1, time0, units = "sec"))

    # compute mantel pval
    time0 <- Sys.time()
    mantel <- vegan::mantel(task$geodesic_dist, model$geodesic_dist, permutations = 100)
    summary_list$mantel_pval <- -log10(mantel$signif)
    time1 <- Sys.time()
    summary_list$time_mantel <- as.numeric(difftime(time1, time0, units = "sec"))
  }

  net1 <- dynutils::simplify_milestone_network(model$milestone_network)
  net2 <- dynutils::simplify_milestone_network(task$milestone_network) %>% filter(to != "FILTERED_CELLS")

  # Compute the milestone network isomorphic
  if ("isomorphic" %in% metrics) {
    time0 <- Sys.time()

    summary_list$isomorphic <- (is_isomorphic_to(
      graph_from_data_frame(net1),
      graph_from_data_frame(net2)
    ))+0
    time1 <- Sys.time()
    summary_list$time_isomorphic <- as.numeric(difftime(time1, time0, units = "sec"))
  }

  # Compute the milestone network GED
  if ("ged" %in% metrics) {
    time0 <- Sys.time()
    summary_list$ged <- calculate_ged(net1, net2)
    time1 <- Sys.time()
    summary_list$time_ged <- as.numeric(difftime(time1, time0, units = "sec"))
  }

  # Compute Netdist EMD (see scripts/wouter/network_scores_tests.R)
  if ("net_emd" %in% metrics) {
    time0 <- Sys.time()
    gdd1 <- netdist::gdd(net1 %>% igraph::graph_from_data_frame())
    gdd2 <- netdist::gdd(net2 %>% igraph::graph_from_data_frame())
    summary_list$net_emd <- netdist::net_emd(gdd1, gdd2)
    time1 <- Sys.time()
    summary_list$time_net_emd <- 1-as.numeric(difftime(time1, time0, units = "sec"))
  }

  if ("robbie_network_score" %in% metrics) {
    time0 <- Sys.time()
    summary_list$robbie_network_score <- calculate_robbie_network_score(net1, net2)
    time1 <- Sys.time()
    summary_list$time_robbie_network_score <- 1-as.numeric(difftime(time1, time0, units = "sec"))
  }

  summary <- as_tibble(summary_list)

  list(coranking = coranking, summary = summary)
}

#' Compute the coranking matrix and
#'
#' @param gold_dist A data frame containing the pairwise distances in the original space
#' @param pred_dist A data frame containing the pairwise distances in the new space
#'
#' @importFrom coRanking coranking LCMC
#' @importFrom tibble lst
compute_coranking <- function(gold_dist, pred_dist) {
  gold_dist <- gold_dist + runif(length(gold_dist), 0, 1e-30)
  pred_dist <- pred_dist + runif(length(pred_dist), 0, 1e-30)
  gold_dist <- (gold_dist + t(gold_dist)) / 2
  pred_dist <- (pred_dist + t(pred_dist)) / 2

  diag(gold_dist) <- 0
  diag(pred_dist) <- 0

  Q <- coRanking::coranking(gold_dist, pred_dist, input = "dist")

  nQ <- nrow(Q)
  N <- nQ + 1

  # calculating Q_nx
  LCMC <- coRanking::LCMC(Q)
  Q_nx <- LCMC + seq_len(nQ) / nQ

  # calculating R_nx
  R_nx <- (nQ * Q_nx - seq_len(nQ)) / seq(nQ - 1, 0, -1)
  R_nx <- R_nx[-nQ]

  # calculating mean R_nx
  mean_R_nx <- mean(R_nx)

  # calculating AUC (space under R_nx curve)
  Ks <- seq_along(R_nx)
  auc_R_nx <- sum(R_nx / Ks) / sum(1 / Ks)

  # calculating Q_global and Q_local
  Kmax <- which.max(LCMC)
  ix <- seq(1, Kmax)
  Q_global <- mean(LCMC[-ix])
  Q_local <- mean(LCMC[ix])

  summary <- lst(mean_R_nx, auc_R_nx, Q_global, Q_local)

  lst(
    Q,
    LCMC,
    Q_nx,
    R_nx,
    summary
  )
}


#' Compute the graph edit distance using gedevo
#'
#' @param net1 the first network to compare
#' @param net2 the second network to compare
#'
#' @importFrom readr read_file
#' @importFrom utils write.table
#' @importFrom GEDEVO run_GEDEVO
calculate_ged <- function(net1, net2) {
  net1 <- net1 %>% mutate(dir="u") %>% select(from, dir, to)
  net2 <- net2 %>% mutate(dir="u") %>% select(from, dir, to)

  tempfolder <- tempfile()
  dir.create(tempfolder)
  write.table(net1, file.path(tempfolder, "net1.sif"), row.names = FALSE, col.names = FALSE)
  write.table(net2, file.path(tempfolder, "net2.sif"), row.names = FALSE, col.names = FALSE)

  cmd <- pritt("gedevo --groups a b --sif {tempfolder}/net1.sif a --sif {tempfolder}/net2.sif b --no-prematch --no-workfiles --save {tempfolder}/out --maxiter 100 --maxsecs 1 --maxsame 100")
  GEDEVO::run_GEDEVO(cmd)

  score <- readr::read_file(paste0(tempfolder, "/out.matching")) %>% gsub("^.*GED score:[ ]*([0-9\\.]*).*", "\\1", .) %>% as.numeric()

  1-score
}




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
