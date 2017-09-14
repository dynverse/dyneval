#' Generate toy datasets with dyngen
#'
#' @param ti_types The types of TI to generate
#' @param num_replicates How many replicates of each TI type to generate
#' @param num_cells The number of cells in each dataset
#' @param num_genes The number of genes in each dataset
#'
#' @importFrom dyngen generate_toy_milestone_network random_progressions_tented generate_expression
#' @importFrom dynutils extract_row_to_list list_as_tibble
#'
#' @export
generate_toy_datasets <- function(ti_types = c("linear", "bifurcating", "cycle"), num_replicates = 3, num_cells = 99, num_genes = 101) {
  settings <- expand.grid(ti_type = ti_types, replicate = seq_len(num_replicates), stringsAsFactors = FALSE)

  list_as_tibble(lapply(seq_len(nrow(settings)), function(rowi) {
    list2env(extract_row_to_list(settings, rowi), environment())

    milestone_network <- dyngen::generate_toy_milestone_network(ti_type)
    progressions <- dyngen::random_progressions_tented(milestone_network, ncells = num_cells)
    expression <- dyngen::generate_expression(milestone_network, progressions, ngenes = num_genes)
    counts <- dyngen::generate_counts(expression)
    colnames(counts) <- paste0("G", seq_len(ncol(counts)))

    cell_ids <- unique(progressions$cell_id)
    milestone_ids <- unique(c(milestone_network$from, milestone_network$to))

    end_milestones <- milestone_ids[!(milestone_ids %in% milestone_network$from)]
    special_cells <- list(
      start_cell_id = progressions %>% arrange(from, to, percentage) %>% pull(cell_id) %>% first,
      end_cell_ids = progressions %>% filter(to %in% end_milestones) %>% group_by(to) %>% arrange(percentage) %>% summarise(cell_id=cell_id[which.max(percentage)]) %>% pull(cell_id)
    )

    task <- wrap_ti_task_data(
      ti_type = ti_type,
      id = paste0(ti_type, "_", replicate),
      counts = counts,
      cell_ids = cell_ids,
      milestone_ids = milestone_ids,
      milestone_network = milestone_network,
      progressions = progressions,
      special_cells = special_cells
    )

    task$type <- "ti_toy"
    task$cell_grouping <- get_cell_grouping(task$milestone_percentages)
    task$geodesic_dist <- compute_emlike_dist(task)

    task
  }))
}
