#' Description for SLICER
#' @export
description_slicer <- function() create_description(
  name = "SLICER",
  short_name = "SLICER",
  package_loaded = c("SLICER"),
  package_required = c("lle", "igraph"),
  par_set = makeParamSet(
    makeIntegerParam(id = "kmin", lower = 2L, upper = 20L, default = 10L),
    makeIntegerParam(id = "m", lower = 2L, upper = 20L, default = 2L),
    makeNumericParam(id = "min_branch_len", lower = 0.5, upper = 20, default = 5),
    makeNumericParam(id = "min_representative_percentage", lower = 0.5, upper = 1, default = 0.8),
    makeNumericParam(id = "max_same_milestone_distance", lower = 0.1, upper = 10, default = 0.1)
  ),
  properties = c(),
  run_fun = run_slicer,
  plot_fun = plot_slicer
)

run_slicer <- function(counts,
                      start_cell_id,
                      end_cell_ids = NULL,
                      kmin = 10,
                      m = 2,
                      min_branch_len = 5,
                      min_representative_percentage = 0.8,
                      max_same_milestone_distance = 0.1) {
  requireNamespace("SLICER")
  requireNamespace("lle")
  requireNamespace("igraph")

  expression <- log2(counts + 1)
  genes <- SLICER::select_genes(expression)
  expression_filtered <- expression[, genes]

  k <- SLICER::select_k(expression_filtered, kmin=kmin)
  traj_lle <- lle::lle(expression_filtered, m=m, k)$Y
  traj_graph <- SLICER::conn_knn_graph(traj_lle, k=k)

  if (is.null(end_cell_ids)) {
    ends <- SLICER::find_extreme_cells(traj_graph, traj_lle)
  } else {
    ends <- match(c(start_cell_id, end_cell_ids), rownames(counts))
  }

  start <- which(rownames(expression_filtered) == start_cell_id)
  cells_ordered <- SLICER::cell_order(traj_graph, start)
  branches <- SLICER::assign_branches(traj_graph, start, min_branch_len=min_branch_len) %>% factor %>% set_names(rownames(expression_filtered))

  # TODO: clean up code, add comments

  # from SLICER we get the branch assignment, and the overall ordering of the cells from the start point
  progressions <- tibble(
    branch_id = branches,
    global_rank = match(rownames(expression), rownames(expression)[cells_ordered]),
    cell_id = rownames(expression)
  ) %>%
    group_by(branch_id) %>%
    mutate(
      rank = rank(global_rank)-1,
      percentage = rank/max(rank)
    ) %>% ungroup()

  # now extract the branch network
  # first create a temporary milestone network
  milestone_network <- progressions %>%
    group_by(branch_id) %>%
    summarise(length=n()) %>%
    mutate(
      from=paste0("M", seq_along(unique(progressions$branch_id))*2-1),
      to=paste0("M", seq_along(unique(progressions$branch_id))*2),
      directed=FALSE
    )

  milestone_ids <- unique(c(milestone_network$from, milestone_network$to))

  progressions <- progressions %>% left_join(milestone_network, by="branch_id")

  milestone_percentages <- convert_progressions_to_milestone_percentages(rownames(expression), milestone_ids, milestone_network, progressions)

  representatives <- milestone_percentages %>% group_by(cell_id) %>% summarise(milestone_id=milestone_id[which.max(percentage)], percentage=max(percentage))

  conn <- traj_graph %>%
    igraph::as_data_frame() %>%
    mutate(
      from=rownames(counts)[from],
      to=rownames(counts)[to]
    ) %>%
    left_join(representatives %>% rename(milestone_id_from=milestone_id), by=c("from"="cell_id")) %>%
    left_join(representatives %>% rename(milestone_id_to=milestone_id), by=c("to"="cell_id"))

  branch_milestone_combinations <- c(paste0(milestone_network$from, milestone_network$to),paste0(milestone_network$to, milestone_network$from)) # disallowed merges of milestones

  weight_cutoff <- 0.5
  close_representatives <- conn %>%
    filter(milestone_id_from != milestone_id_to) %>%
    arrange(-weight) %>%
    filter(weight >= weight_cutoff) %>%
    filter(!(paste0(milestone_id_from, milestone_id_to) %in% branch_milestone_combinations))

  # group the representatives, according to close distance
  gr <- igraph::graph_from_data_frame(close_representatives %>% select(from = milestone_id_from, to = milestone_id_to), directed = FALSE, vertices = milestone_ids)
  milestone_groups <- igraph::components(gr)$membership

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
    milestone_network = milestone_network %>% select(from, to, length, directed),
    milestone_percentages = milestone_percentages %>% select(cell_id, milestone_id, percentage),
    dimred_samples = traj_lle,
    dimred_clust = branches
  )
}

#' @import ggplot2
plot_slicer <- function(ti_predictions) {
  dat <- as.data.frame(ti_predictions$dimred_samples) %>% mutate(branch = ti_predictions$dimred_clust)
  ggplot(dat) + geom_point(aes(V1, V2, color=branch))
}
