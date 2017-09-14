#' Used for wrapping an evaluation function around a TI method
#'
#' @param method the method to wrap an evaluation function around
#' @param noisy whether or not the metric is noisy or not
#' @param suppress_output whether or not to suppress the output
#' @param metrics which metrics to evaluate with
#'
#' @importFrom smoof makeSingleObjectiveFunction makeMultiObjectiveFunction
#' @export
make_obj_fun <- function(method, noisy = F, suppress_output = T,
                         metrics = c("mean_R_nx", "auc_R_nx", "Q_global", "Q_local", "correlation", "isomorphic", "robbie_network_score")) {
  # Use different makefunction if there are multiple metrics versus one
  if (length(metrics) > 1) {
    make_fun <- function(...) makeMultiObjectiveFunction(..., n.objectives = length(metrics))
  } else {
    make_fun <- makeSingleObjectiveFunction
  }

  # Wrap the method function in an evaluation function
  make_fun(
    name = "TItrain",
    vectorized = F,
    minimize = rep(F, length(metrics)),
    noisy = noisy,
    has.simple.signature = F,
    par.set = method$par_set,
    fn = function(x, tasks)
      execute_evaluation(
        tasks = tasks,
        method = method,
        parameters = x,
        metrics = metrics,
        suppress_output = suppress_output))
}

#' For returning a poor score when a method errors
#'
#' @param num_objectives the number of objectives used
#' @param error_score the score to return upon erroring
#'
#' @export
impute_y_fun <- function(num_objectives, error_score = -1) {
  function(x, y, opt.path, ...) {
    val <- rep(error_score, num_objectives)
    attr(val, "extras") <- list(.summary = NA)
    val
  }
}


#' Running an evaluation of a method on a set of tasks with a set of parameters
#'
#' @inheritParams execute_method
#' @param metrics which metrics to use
#'
#' @export
#' @importFrom dynutils override_setseed extract_row_to_list
#' @importFrom future future value plan
#' @importFrom netdist gdd net_emd
execute_evaluation <- function(
  tasks,
  method,
  parameters,
  metrics = c("mean_R_nx", "auc_R_nx", "Q_global", "Q_local", "correlation", "isomorphic", "robbie_network_score"),
  suppress_output = TRUE) {

  method_outputs <- execute_method(tasks, method, parameters, suppress_output)

  # Calculate scores
  summary_outs <- lapply(seq_len(nrow(tasks)), function(i) {
    task <- extract_row_to_list(tasks, i)

    # Fetch method outputs
    method_output <- method_outputs[[i]]
    model <- method_output$model

    # Calculate geodesic distances
    time0 <- Sys.time()
    model$geodesic_dist <- compute_emlike_dist(model)
    time1 <- Sys.time()
    time_geodesic <- as.numeric(difftime(time1, time0, units = "sec"))

    # Calculate metrics
    metrics_output <- calculate_metrics(task, model, metrics)

    # Create summary statistics
    summary <- data.frame(
      method_name = method$name,
      method_short_name = method$short_name,
      task_id = task$id,
      method_output$summary,
      time_geodesic = time_geodesic,
      metrics_output$summary,
      stringsAsFactors = F,
      check.names = F
    )

    # Return the output
    lst(model, summary)
  })

  # Combine the different outputs in three lists/data frames
  models <- summary_outs %>% purrr::map(~ .$model)
  summary <- summary_outs %>% purrr::map_df(~ .$summary)

  # Calculate the final score
  score <- summary %>% summarise_at(metrics, funs(mean)) %>% as.matrix %>% as.vector %>% setNames(metrics)

  # Return extra information
  attr(score, "extras") <- list(.models = models, .summary = summary)

  # Return output
  score
}

#' @importFrom igraph is_isomorphic_to graph_from_data_frame
calculate_metrics <- function(task, model, metrics) {
  summary <- data.frame(row.names = 1)

  # Compute coranking metrics
  if (any(c("mean_R_nx", "auc_R_nx", "Q_local", "Q_global") %in% metrics)) {
    time0 <- Sys.time()
    coranking <- compute_coranking(task$geodesic_dist, model$geodesic_dist)
    summary <- bind_cols(summary, coranking$summary)
    time1 <- Sys.time()
    summary$time_coranking <- as.numeric(difftime(time1, time0, units = "sec"))
  }

  # Compute the correlation of the geodesic distances
  if ("correlation" %in% metrics) {
    time0 <- Sys.time()
    summary$correlation <- cor(task$geodesic_dist %>% as.vector, model$geodesic_dist %>% as.vector, method="spearman")
    time1 <- Sys.time()
    summary$time_correlation <- as.numeric(difftime(time1, time0, units = "sec"))
  }

  # Compute the mantel test
  if ("mantel_pvalue" %in% metrics) {
    time0 <- Sys.time()
    mantel <- vegan::mantel(task$geodesic_dist, model$geodesic_dist, permutations = 1000, alternative="greater")
    summary$mantel_pval <- mantel$signif
    time1 <- Sys.time()
    summary$time_mantel <- as.numeric(difftime(time1, time0, units = "sec"))
  }

  net1 <- dynutils::simplify_network(model$milestone_network)
  net2 <- dynutils::simplify_network(task$milestone_network) %>% filter(to != "FILTERED_CELLS")

  # Compute the milestone network isomorphic
  if ("isomorphic" %in% metrics) {
    time0 <- Sys.time()

    summary$isomorphic <- (is_isomorphic_to(
      graph_from_data_frame(net1),
      graph_from_data_frame(net2)
    ))+0
    time1 <- Sys.time()
    summary$time_isomorphic <- as.numeric(difftime(time1, time0, units = "sec"))
  }

  # Compute the milestone network GED
  if ("ged" %in% metrics) {
    time0 <- Sys.time()
    summary$ged <- calculate_ged(net1, net2)
    time1 <- Sys.time()
    summary$time_ged <- as.numeric(difftime(time1, time0, units = "sec"))
  }

  # Compute Netdist EMD (see scripts/wouter/network_scores_tests.R)
  if ("net_emd" %in% metrics) {
    time0 <- Sys.time()
    gdd1 <- netdist::gdd(net1 %>% igraph::graph_from_data_frame())
    gdd2 <- netdist::gdd(net2 %>% igraph::graph_from_data_frame())
    summary$net_emd <- netdist::net_emd(gdd1, gdd2)
    time1 <- Sys.time()
    summary$time_net_emd <- 1-as.numeric(difftime(time1, time0, units = "sec"))
  }

  if ("robbie_network_score" %in% metrics) {
    time0 <- Sys.time()
    summary$robbie_network_score <- calculate_robbie_network_score(net1, net2)
    time1 <- Sys.time()
    summary$time_robbie_network_score <- 1-as.numeric(difftime(time1, time0, units = "sec"))
  }


  rownames(summary) <- NULL

  lst(coranking, summary)
}

#' Calculate Earth Mover's Distance between cells in a trajectory
#'
#' @param traj the trajectory
#'
#' @export
#'
#' @importFrom igraph graph_from_data_frame E distances
#' @importFrom transport transport
#' @import dplyr
#' @importFrom purrr %>% map map_df map_lgl
compute_emlike_dist <- function(traj) {
  cell_ids <- traj$cell_ids
  milestone_network <- traj$milestone_network
  milestone_ids <- traj$milestone_ids
  milestone_percentages <- traj$milestone_percentages

  # calculate the shortest path distances between milestones
  phantom_edges <- bind_rows(lapply(milestone_ids, function(sn) {
    sn_filt <- milestone_network %>% filter(from == sn)
    dis_vec <- setNames(sn_filt$length, sn_filt$to)
    phantom_edges <-
      expand.grid(from = sn_filt$to, to = sn_filt$to, stringsAsFactors = F) %>%
      filter(from < to) %>%
      left_join(milestone_network, by = c("from", "to")) %>%
      filter(is.na(length)) %>%
      mutate(length = dis_vec[from] + dis_vec[to])
    phantom_edges
  }))
  gr <- igraph::graph_from_data_frame(bind_rows(milestone_network %>% mutate(length = 2 * length), phantom_edges), directed = F, vertices = milestone_ids)
  milestone_distances <- igraph::distances(gr, weights = igraph::E(gr)$length, mode = "all")

  # transport percentages data
  milestone_percentages$milestone_id = factor(milestone_percentages$milestone_id, levels=milestone_ids) # make sure all milestones are included, even if none of the cells have a value for the milestone
  pct <- reshape2::acast(milestone_percentages, cell_id ~ milestone_id, value.var = "percentage", fill = 0, drop=FALSE)
  pct <- pct[cell_ids, milestone_ids]

  fromto_matrix <- matrix(0, nrow = length(milestone_ids), ncol = length(milestone_ids), dimnames = list(milestone_ids, milestone_ids))
  fromto2 <- milestone_network %>% reshape2::acast(from ~ to, value.var = "length", fun.aggregate = length)
  fromto_matrix[rownames(fromto2), colnames(fromto2)] <- fromto2
  diag(fromto_matrix) <- 1
  fromto_matrix[fromto_matrix > 0] <- 1

  froms <- setNames(lapply(rownames(pct), function(xi) {
    x <- pct[xi,]
    notzero <- x != 0
    wh <-
      if (sum(notzero) == 1) {
        which(notzero)
      } else {
        apply(fromto_matrix, 1, function(y) {
          all(!notzero | y)
        })
      }
    milestone_ids[wh]
  }), rownames(pct))

  closest <- bind_rows(lapply(milestone_ids, function(mid) {
    sample_node <- froms %>% map_lgl(~ mid %in% .)
    if (sum(sample_node) == 0) {
      NULL
    } else {
      milestones <- which(fromto_matrix[mid,] == 1)
      dist_milestones <- milestone_distances[milestones, milestones, drop = F]
      sample_pcts <- pct[sample_node, milestones, drop = F]
      closest_to_nodes <-
        sample_pcts %*% dist_milestones %>%
        reshape2::melt(varnames = c("from", "to"), value.name = "length") %>%
        mutate(from = as.character(from), to = as.character(to))

      closest_to_samples <-
        expand.grid(from = rownames(sample_pcts), to = rownames(sample_pcts), stringsAsFactors = F) %>%
        filter(from < to)
      closest_to_samples$length <- sapply(seq_len(nrow(closest_to_samples)), function(xi) {
        a <- sample_pcts[closest_to_samples$from[[xi]],]
        b <- sample_pcts[closest_to_samples$to[[xi]],]
        diff <- a - b
        which_from <- diff < 0
        num_from <- sum(which_from)
        which_to <- diff > 0
        num_to <- sum(which_to)

        if (num_from == 1) {
          sum(dist_milestones[which_from, which_to] * diff[which_to])
        } else if (num_to == 1) {
          -sum(dist_milestones[which_from, which_to] * diff[which_from])
        } else if (num_from == 0) {
          0
        } else {
          suppressWarnings({
            transport::transport(a, b, costm = dist_milestones) %>%
              mutate(dist = dist_milestones[cbind(from,to)], mult = mass * dist) %>%
              .$mult %>%
              sum
          })
        }
      })

      bind_rows(closest_to_nodes, closest_to_samples)
    }
  }))

  gr2 <- igraph::graph_from_data_frame(closest, directed = F, vertices = c(milestone_ids, cell_ids))
  gr2 %>% igraph::distances(v = cell_ids, to = cell_ids, weights = igraph::E(gr2)$length)
}

#' Compute the coranking matrix and
#'
#' @param gold_dist A data frame containing the pairwise distances in the original space
#' @param pred_dist A data frame containing the pairwise distances in the new space
#'
#' @export
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

  summary <- data_frame(mean_R_nx, auc_R_nx, Q_global, Q_local)

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
#' @importFrom glue glue
#' @importFrom utils write.table
#' @importFrom GEDEVO run_GEDEVO
#' @export
calculate_ged <- function(net1, net2) {
  net1 <- net1 %>% mutate(dir="u") %>% select(from, dir, to)
  net2 <- net2 %>% mutate(dir="u") %>% select(from, dir, to)

  tempfolder <- tempfile()
  dir.create(tempfolder)
  write.table(net1, file.path(tempfolder, "net1.sif"), row.names = FALSE, col.names = FALSE)
  write.table(net2, file.path(tempfolder, "net2.sif"), row.names = FALSE, col.names = FALSE)

  cmd <- glue::glue("gedevo --groups a b --sif {tempfolder}/net1.sif a --sif {tempfolder}/net2.sif b --no-prematch --no-workfiles --save {tempfolder}/out --maxiter 100 --maxsecs 1 --maxsame 100")
  GEDEVO::run_GEDEVO(cmd)

  score <- read_file(paste0(tempfolder, "/out.matching")) %>% gsub("^.*GED score:[ ]*([0-9\\.]*).*", "\\1", .) %>% as.numeric()

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

get_adjacency <- function(net, nodes=unique(c(net$from, net$to))) {
  newnet <- net %>%
    mutate(from=factor(from, levels=nodes), to=factor(to, levels=nodes)) %>%
    reshape2::acast(from~to, value.var="length", fill=0, drop=F)

  newnet[lower.tri(newnet)] = newnet[lower.tri(newnet)] + t(newnet)[lower.tri(newnet)] # make symmetric

  newnet
}

#' Compute the robbie network score (any resemblances to real-life persons is purely coincidental)
#'
#' @param net1 the first network to compare
#' @param net2 the second network to compare
#' @importFrom GA ga
#' @export
calculate_robbie_network_score <- function(net1, net2) {
  nodes1 <- unique(c(net1$from, net1$to))
  nodes2 <- unique(c(net2$from, net2$to))

  if(length(nodes1) < length(nodes2)) {
    nodes3 <- nodes2
    net3 <- net2
    nodes2 <- nodes1
    net2 <- net1
    nodes1 <- nodes3
    net1 <- net3
  }

  # if both networks have only one edge...
  if(nrow(net1) == 1) {
    # special cases: either the networks are a cycle, or contain one linear edge
    if(length(unique(c(net1$from, net1$to))) == length(unique(c(net2$from, net2$to)))) {
      return(1)
    } else {
      return(0)
    }
  }

  net <- get_adjacency(net1)
  net_ref <- get_adjacency(net2)
  net_ref <- complete_matrix(net_ref, nrow(net))

  # normalize weights
  net <- net/sum(net)
  net_ref <- net_ref/sum(net_ref)

  results <- GA::ga("permutation", score_map, net=net, net_ref=net_ref, min=rep(1, length(nodes1)), max=rep(length(nodes1), length(nodes1)), maxiter=10, monitor=FALSE)

  results@fitnessValue
}
