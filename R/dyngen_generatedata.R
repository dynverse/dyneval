#' Generate toy datasets with dyngen
#'
#' @importFrom dyngen generate_toy_milestone_network random_progressions_tented generate_expression
#' @importFrom dynutils extract_row_to_list list_as_tibble
#'
#' @export
generate_toy_datasets <- function() {
  settings <- expand.grid(ti_type = c("linear", "bifurcating", "cycle"), replicate = 1:5, stringsAsFactors = F)

  list_as_tibble(lapply(seq_len(nrow(settings)), function(rowi) {
    list2env(extract_row_to_list(settings, rowi), environment())

    milestone_network <- dyngen::generate_toy_milestone_network(ti_type)
    progressions <- dyngen::random_progressions_tented(milestone_network, ncells = 99)
    expression <- dyngen::generate_expression(milestone_network, progressions, ngenes = 101)
    counts <- dyngen::generate_counts(expression)
    colnames(counts) <- paste0("G", seq_len(ncol(counts)))

    cell_ids <- unique(progressions$cell_id)
    milestone_ids <- unique(c(milestone_network$from, milestone_network$to))

    # todo: replace with real special_cells
    special_cells <- list(
      start_cell_id = progressions %>% arrange(from, to, percentage) %>% pull(cell_id) %>% first
    )

    task <- wrap_ti_task_data(
      ti_type = "toy",
      id = paste0(ti_type, "_", replicate),
      counts = counts,
      cell_ids = cell_ids,
      milestone_ids = milestone_ids,
      milestone_network = milestone_network,
      progressions = progressions,
      special_cells = special_cells
    )

    task$cell_grouping = get_cell_grouping(task$milestone_percentages)
    task$geodesic_dist <- compute_emlike_dist(task)
    task
  }))
}
