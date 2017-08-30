#' @importFrom dyngen generate_toy_milestone_network random_progressions_tented
generate_toy_datasets <- function() {
  settings <- expand.grid(ti_type = c("linear", "bifurcating", "cycle"), replicate = 1:5, stringsAsFactors = F)

  list_as_tibble(lapply(seq_len(nrow(settings)), function(rowi) {
    list2env(extract_row_to_list(settings, rowi), environment())

    milestone_network <- dyngen::generate_toy_milestone_network(ti_type)
    progressions <- dyngen::random_progressions_tented(milestone_network)
    expression <- dyngen::generate_expression(milestone_network, progressions)
    counts <- round(expression * 100)

    cell_ids <- unique(progressions$cell_id)
    milestone_ids <- unique(c(milestone_network$from, milestone_network$to))

    # todo: replace with real special_cells
    special_cells <- list()

    # milestone_percentages <- convert_progressions_to_milestone_percentages(cell_ids, milestone_ids, milestone_network, progressions)

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
    task$geodesic_dist <- dyneval::compute_emlike_dist(task)
    task
  }))
}
