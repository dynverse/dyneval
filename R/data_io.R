
#' Load dyngen datasets
#'
#' \strong{Important!} A variable \code{.datasets_location} needs to be defined in the global environment for this function
#' to work correctly.
#'
#' @param mc_cores The number of cores to use while loading in the datasets (default: 1)
#' @param num_datasets The maximum number of datasets to load (default: Inf)
#'
#' @return A tibble of datasets
#' @export
#'
#' @examples
#' datasets <- load_datasets(mc_cores = 1, num_datasets = Inf)
load_datasets <- function(mc_cores = 1, num_datasets = Inf) {
  # read datasets meta information
  datasets_info <- readRDS(paste0(.datasets_location, "/datasets.rds"))

  # determine the number of datasets to read
  num_datasets <- min(nrow(datasets_info), num_datasets)

  # load the datasets one by one
  task_wrapped <- parallel::mclapply(seq_len(num_datasets), mc.cores = mc_cores, function(dataset_num) {
    dataset_id <- datasets_info$id[[dataset_num]]
    dataset <- dyngen::load_dataset(dataset_id)

    list2env(dataset, environment())

    cell_ids <- rownames(counts)
    cell_info <- cellinfo %>% slice(match(step_id, rownames(counts)))
    milestone_ids <- gs$milestone_percentages$milestone_id %>% unique %>% as.character
    milestone_network <- gs$milestone_network %>%
      mutate(
        from = as.character(from),
        to = as.character(to)
      )
    milestone_percentages <- cell_info %>%
      left_join(gs$milestone_percentages, by = c("step_id" = "cell_id")) %>%
      filter(percentage > 0) %>%
      mutate(milestone_id = as.character(milestone_id)) %>%
      select(cell_id, milestone_id, percentage)
    sample_info <- cell_info %>%
      select(cell_id, step, simulation_time = simulationtime)

    out <- wrap_ti_task_data(
      ti_type = model$modulenetname,
      id = dataset_id,
      cell_ids = cell_ids,
      milestone_ids = milestone_ids,
      milestone_network = milestone_network,
      milestone_percentages = milestone_percentages,
      counts = counts,
      sample_info = sample_info,
      task_ix = dataset_num,
      modulenet_id = model$modulenetname,
      platform_id = platform$platform_id,
      takesetting_type = dataset$takesetting$type,
      model_replicate = model$modelsetting$replicate,
      special_cells = special_cells
    )
    out$geodesic_dist <- compute_emlike_dist(out)
    out
  })
  task_wrapped %>%
    list_as_tibble %>%
    left_join(datasets_info, by = c("id" = "id"))
}
