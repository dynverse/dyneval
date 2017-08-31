#' Description for SLICER
#' @export
description_slicer <- function() {
  list(
    name = "SLICER",
    short_name = "SLICER",
    package_load = c("SLICER"),
    package_installed = c(),
    par_set = makeParamSet(
      makeIntegerParam(id = "kmin", lower = 2L, upper = 20L, default = 10),
      makeIntegerParam(id = "m", lower = 2L, upper = 20L, default = 2),
      makeNumericParam(id = "min_branch_len", lower = 0.5, upper = 20, default = 5),
      makeNumericParam(id = "min_representative_percentage", lower = 0.5, upper = 1, default = 0.8),
      makeNumericParam(id = "max_same_milestone_distance", lower = 0.1, upper = 10, default = 0.1)
    ),
    properties = c(),
    run_fun = run_slicer,
    plot_fun = plot_slicer
  )
}

run_slicer <- function(counts,
                      start_cell_id = sample(rownames(counts), 1),
                      kmin = 10,
                      m = 2,
                      min_branch_len = 5,
                      min_representative_percentage = 0.8,
                      max_same_milestone_distance = 0.1) {
  requireNamespace("SLICER")
  requireNamespace("lle")

  set.seed(1)

  expression <- log2(counts + 1)

  genes <- SLICER::select_genes(expression)
  expression_filtered <- expression[, genes]

  k <- SLICER::select_k(expression_filtered, kmin=kmin)
  traj_lle <- lle::lle(expression_filtered, m=m, k)$Y
  traj_graph <- SLICER::conn_knn_graph(traj_lle, k=k)
  ends <- SLICER::find_extreme_cells(traj_graph, traj_lle)

  #start_cell_id <- dataset$gs$cellinfo %>% filter(cell_id %in% rownames(expression_filtered)) %>% arrange(state_id, time) %>% pull(cell_id) %>% first
  start <- which(rownames(expression_filtered) == start_cell_id)
  #start <- 1
  cells_ordered <- SLICER::cell_order(traj_graph, start)
  branches <- SLICER::assign_branches(traj_graph, start, min_branch_len=min_branch_len) %>% factor %>% set_names(rownames(expression_filtered))

  # from SLICER we get the branch assignment, and the overall ordering of the cells from the start point
  progressions <- tibble(
    branch_id = branches,
    global_rank = match(rownames(expression), rownames(expression)[cells_ordered]),
    cell_id = rownames(expression)
  )
  progressions <- progressions %>%
    group_by(branch_id) %>%
    mutate(
      rank = rank(global_rank)-1,
      percentage = rank/max(rank)
    )

  # now extract the branch network
  # first create a temporary milestone network
  milestone_network <- progressions %>%
    group_by(branch_id) %>%
    summarise(length=n()) %>%
    mutate(
      from=paste0("M", seq_along(unique(progressions$branch_id))*2-1),
      to=paste0("M", seq_along(unique(progressions$branch_id))*2)
    )

  end_branches <- branches[rownames(expression_filtered)[ends]] # ignore these branches for making the connections
  start_branch <- branches[rownames(expression_filtered)[start]] # ignore these branches for making the connections
  connect_milestone_ids <- c(
    milestone_network %>% filter(!(branch_id %in% end_branches)) %>% pull(to),
    milestone_network %>% filter(!(branch_id %in% start_branch)) %>% pull(from)
  )

  milestone_ids <- unique(c(milestone_network$from, milestone_network$to))

  progressions <- progressions %>% left_join(milestone_network, by="branch_id")

  milestone_percentages <- convert_progressions_to_milestone_percentages(rownames(expression), milestone_ids, milestone_network, progressions)

  # we will now check which milestones are actually the "same" (ie. low distance to eachother)
  # use some milestone representatives, those which have a high percentage of the particular milestone
  # min_representative_percentage: between 0.5 and 1
  # to get at least one representative for each milestone, always get the one with the highest percentages
  milestone_representatives <- milestone_percentages %>%
    filter(percentage > (min_representative_percentage) | percentage == max(percentage)) %>%
    select(-percentage)

  # calculate distances between these representatives
  cell_distances <- map(seq_len(nrow(expression)), ~SLICER::process_distance(traj_graph, .)) %>%
    invoke(rbind, .) %>% set_rownames(rownames(expression)) %>% set_colnames(rownames(expression))
  representative_distances <- cell_distances[milestone_representatives$cell_id, milestone_representatives$cell_id] %>%
    reshape2::melt(varnames=c("cell_id_from", "cell_id_to"), value.name="distance") %>%
    mutate(cell_id_from=as.character(cell_id_from), cell_id_to=as.character(cell_id_to)) %>%
    left_join(milestone_representatives %>% rename(milestone_id_from=milestone_id), by=c("cell_id_from"="cell_id")) %>%
    left_join(milestone_representatives %>% rename(milestone_id_to=milestone_id), by=c("cell_id_to"="cell_id"))

  # group the representatives, according to close distance
  milestone_groups <- setNames(seq_along(milestone_ids), milestone_ids)

  close_representatives <- representative_distances %>%
    group_by(milestone_id_from, milestone_id_to) %>%
    summarise(distance=min(distance)) %>%
    filter(milestone_id_from != milestone_id_to) %>%
    filter((milestone_id_from %in% connect_milestone_ids) & (milestone_id_to %in% connect_milestone_ids)) %>%
    mutate(
      close =
        (distance < quantile(cell_distances, max_same_milestone_distance)) |
        (distance == min(distance))
      ) %>%
    filter(close)

  for (i in seq_len(nrow(close_representatives))) {
    milestone_groups[[close_representatives[i,]$milestone_id_from]] <- min(milestone_groups[[close_representatives[i,]$milestone_id_to]], milestone_groups[[close_representatives[i,]$milestone_id_from]])
  }

  # now map old milestones to new milestones
  milestone_mapper <- setNames(paste0("M", milestone_groups), names(milestone_groups))

  milestone_ids <- unique(milestone_mapper)
  milestone_percentages$milestone_id <- milestone_mapper[milestone_percentages$milestone_id]
  progressions$from <- milestone_mapper[progressions$from]
  progressions$to <- milestone_mapper[progressions$to]
  milestone_network$from <- milestone_mapper[milestone_network$from]
  milestone_network$to <- milestone_mapper[milestone_network$to]

  # tada!
  wrap_ti_prediction(
    ti_type = "linear",
    id = "SLICER",
    cell_ids = rownames(expression_filtered),
    milestone_ids = milestone_ids,
    milestone_network = milestone_network %>% select(from, to, length),
    milestone_percentages = milestone_percentages,
    dimred_samples = traj_lle,
    dimred_clust = branches
  )
}

#' @import ggplot2
plot_slicer <- function(ti_predictions) {
  dat <- as.data.frame(ti_predictions$dimred_samples) %>% mutate(branch = ti_predictions$dimred_clust)
  ggplot(dat) + geom_point(aes(V1, V2, color=branch))
}
