#' Load dyngen datasets meta information
#'
#' \strong{Important!} A variable \code{.datasets_location} needs to be defined in the global environment for this function
#' to work correctly.
#'
#' @return The meta information of the datasets
#' @export
#'
#' @examples
#' \dontrun{
#' .dataset_location <- "path_to_dyngen_results"
#' load_datasets_info()
#' }
load_datasets_info <- function() {
  readRDS(paste0(.datasets_location, "/datasets.rds"))
}


#' Load dyngen datasets
#'
#' \strong{Important!} A variable \code{.datasets_location} needs to be defined in the global environment for this function
#' to work correctly.
#'
#' @param mc_cores The number of cores to use while loading in the datasets (default: 1)
#' @param datasets_info The meta information of the datasets to read
#'
#' @return A tibble of datasets
#' @export
#'
#' @importFrom dyngen load_dataset
#'
#' @examples
#' \dontrun{
#' .dataset_location <- "path_to_dyngen_results"
#' datasets <- load_datasets(mc_cores = 1, datasets_info = load_datasets_info())
#' }
load_datasets <- function(mc_cores = 1, datasets_info = load_datasets_info()) {
  # load the datasets one by one
  task_wrapped <- parallel::mclapply(seq_len(nrow(datasets_info)), mc.cores = mc_cores, function(dataset_num) {
    dataset_id <- datasets_info$id[[dataset_num]]
    dataset <- dyngen::load_dataset(dataset_id)

    list2env(dataset, environment())

    cell_ids <- rownames(counts)

    sample_info <- cellinfo %>%
      slice(match(cell_ids, step_id)) %>%
      select(cell_id, step, simulation_time = simulationtime)

    milestone_ids <- gs$milestone_percentages$milestone_id %>%
      unique %>%
      as.character

    milestone_network <- gs$milestone_network %>%
      mutate(
        from = as.character(from),
        to = as.character(to)
      )

    progression <- sample_info %>%
      select(cell_id) %>%
      left_join(gs$progression, by = "cell_id") %>%
      mutate(
        from = as.character(from),
        to = as.character(to)
      )

    out <- wrap_ti_task_data(
      ti_type = model$modulenetname,
      id = dataset_id,
      cell_ids = cell_ids,
      milestone_ids = milestone_ids,
      milestone_network = milestone_network,
      progression = progression,
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
