#' @export
make_obj_fun <- function(method) {
  makeSingleObjectiveFunction(
    name = "TItrain",
    vectorized = F,
    minimize = F,
    noisy = F,
    has.simple.signature = F,
    par.set = method$par_set,
    fn = function(x, tasks) {
      for (pack in method$package) do.call(library, list(pack))

      outs <- lapply(seq_len(nrow(tasks)), function(i) {
        arglist <- c(list(counts = tasks$counts[[i]]), x)
        capture.output({ model <- do.call(method$run_fun, arglist) })
        coranking <- compute_coranking(tasks$geodesic_dist[[i]], model$geodesic_dist)
        list(model = model, coranking = coranking)
      })
      df <- outs %>% map_df(~ .$coranking$summary) %>% mutate(task_name = tasks$name)

      score <- mean(df$max_lcmc)
      attr(score, "extras") <- list(.summary = df)
      score
    })
}

#' @export
impute_y_fun <- function(x, y, opt.path, ...) {
  val <- -1
  attr(val, "extras") <- list(.summary = NA)
  val
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
compute_emlike_dist <- function(traj) {
  state_network <- traj$state_network
  state_names <- traj$state_names
  state_percentages <- traj$state_percentages

  # calculate the shortest path distances between milestones
  phantom_edges <- bind_rows(lapply(state_names, function(sn) {
    sn_filt <- state_network %>% filter(from == sn)
    dis_vec <- setNames(sn_filt$length, sn_filt$to)
    phantom_edges <-
      expand.grid(from = sn_filt$to, to = sn_filt$to, stringsAsFactors = F) %>%
      filter(from < to) %>%
      mutate(length = dis_vec[from] + dis_vec[to])
    phantom_edges
  }))
  gr <- igraph::graph_from_data_frame(bind_rows(state_network %>% mutate(length = 2 * length), phantom_edges), directed = F, vertices = state_names)
  milestone_distances <- igraph::distances(gr, weights = igraph::E(gr)$length, mode = "all")

  # transport percentages data
  pct <- as.matrix(state_percentages[,-1])
  ids <- as.character(state_percentages$id)
  rownames(pct) <- ids

  fromto_matrix <- matrix(0, nrow = length(state_names), ncol = length(state_names), dimnames = list(state_names, state_names))
  fromto2 <- state_network %>% reshape2::acast(from~to, value.var = "length", fun.aggregate = length)
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
    state_names[wh]
  }), rownames(pct))

  closest <- bind_rows(lapply(state_names, function(state_name) {
    # cat("State ", state_name, "\n", sep="")
    sample_node <- froms %>% map_lgl(~ state_name %in% .)
    if (sum(sample_node) == 0) {
      NULL
    } else {
      milestones <- which(fromto_matrix[state_name,] == 1)
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
          transport::transport(a, b, costm = dist_milestones) %>%
            mutate(dist = dist_milestones[cbind(from,to)], mult = mass * dist) %>%
            .$mult %>%
            sum
        }
      })

      bind_rows(closest_to_nodes, closest_to_samples)
    }
  }))

  gr2 <- igraph::graph_from_data_frame(closest, directed = F, vertices = c(state_names, ids))
  gr2 %>% igraph::distances(v = ids, to = ids, weights = igraph::E(gr2)$length)
}

#' Plot the Earth Mover's distances in a heatmap
#'
#' @param traj the trajectory (less than 500 cells is recommended)
#' @param emdist the Earth Mover's distances as calculated by \code{\link{emdist}}
#' @param dimred the dimensionality reduction of the trajectory as produced by \code{\link{plotLearnerData.ti.default}}
#' @export
#'
#' @importFrom reshape2 acast
#' @importFrom pheatmap pheatmap
plot_emdist <- function(traj, dist, dimred = NULL, ...) {
  state_percentages <- traj$state_percentages
  pct <- as.data.frame(state_percentages[,-1])
  rownames(pct) <- state_percentages$id

  if (is.null(dimred)) {
    dimred <- plotLearnerData.ti.default(traj)
  }

  ann_colours <- setNames(lapply(dimred$space_states$colour, function(x) c("white", x)), dimred$space_states$id)

  pheatmap::pheatmap(
    dist,
    cluster_rows = F,
    cluster_cols = F,
    annotation_col = pct,
    annotation_row = pct,
    annotation_colors = ann_colours,
    legend = F,
    legend_labels = F,
    legend_breaks = F,
    border_color = NA,
    show_rownames = F,
    show_colnames = F,
    annotation_legend = F,
    annotation_names_row = F,
    annotation_names_col = F,
    fontsize = 20 / length(traj$state_names),
    ...)
}

#' Compute the coranking matrix and
#'
#' @param gold_dist A data frame containing the pairwise distances in the original space
#' @param pred_dist A data frame containing the pairwise distances in the new space
#'
#' @export
#'
#' @importFrom coRanking coranking LCMC
compute_coranking <- function(gold_dist, pred_dist) {
  gold_dist <- gold_dist + runif(length(gold_dist), 0, 1e-20)
  pred_dist <- pred_dist + runif(length(pred_dist), 0, 1e-20)
  diag(gold_dist) <- 0
  diag(pred_dist) <- 0

  corank <- coRanking::coranking(gold_dist, pred_dist, input = "dist")

  lcmc <- coRanking::LCMC(corank)

  summary <- data_frame(
    auc_lcmc = mean(lcmc),
    max_lcmc = max(lcmc)
  )

  list(
    corank = corank,
    lcmc = lcmc,
    summary = summary
  )
}

#' Convert percentages to a unique edge assignment
cal_branch_assignment <- function(percentages, network) {
  branch_assignment <- apply(percentages, 1, function(x) {
    positives <- names(which(x > 0))
    if(length(intersect(positives, network$from)) > 0) {
      intersect(positives, network$from)[[1]]
    } else {
      network %>% filter(to == positives[[1]]) %>% .$from %>% .[[1]]
    }
  }) %>% factor() %>% as.numeric()
}

#' Compute the F1rr based on two set of assignment labels
F1rr <- function(labels1, labels2) {
  overlaps <- map(unique(labels1), function(label1) {
    map_dbl(unique(labels2), function(label2) {
      sum((labels1 == label1) & (labels2 == label2)) / sum((labels1 == label1) | (labels2 == label2))
    })
  }) %>% invoke(cbind, .)
  # print(overlaps)
  2/(1/mean(apply(overlaps, 1, max)) + 1/mean(apply(overlaps, 2, max)))
}

#' Compute the F1rr score for overlap of branches
#'
#' Can currently not handle multiple linear edges without branching, which should be merged in the future
#'
#' @export
score_prediction_F1rr <- function(task, pred_output) {
  branch_assignment_true <- cal_branch_assignment(task$state_percentages %>% select(-id), task$state_network)
  branch_assignment_observed <- cal_branch_assignment(pred_output$state_percentages %>% select(-id), pred_output$state_network)
  F1rr(branch_assignment_true, branch_assignment_observed)
}
