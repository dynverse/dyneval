
# sample some datasets
# task_sel <- dplyr::sample_n(load_datasets_info(), 8)

generate_toy_datasets <- function() {
  settings <- expand.grid(ti_type = c("linear", "bifurcating", "cycle"), replicate = 1:5, stringsAsFactors = F)

  list_as_tibble(lapply(seq_len(nrow(settings)), function(rowi) {
    list2env(extract_row_to_list(settings, 1), environment())

    milestone_network <- dyngen::generate_toy_milestone_network(ti_type)
    progressions <- dyngen::random_progressions_tented(milestone_network)

    cell_ids <- unique(progressions$cell_id)
    milestone_ids <- unique(c(milestone_network$from, milestone_network$to))

    # todo: replace with real counts
    counts <- matrix(runif(length(cell_ids) * 100), nrow = length(cell_ids), dimnames = list(cell_ids, paste0("G", seq_len(100))))

    # todo: replace with real special_cells
    special_cells <- list()

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

set.seed(1)
tasks <- generate_toy_datasets()
saveRDS(tasks, paste0(tempdir(), "/dyneval_test_datasets.rds"))
