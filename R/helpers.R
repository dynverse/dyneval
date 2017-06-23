#' @export
to_tibble <- function(list_of_rows) {
  list_names <- names(list_of_rows[[1]])
  # list_types <- sapply(list_of_rows[[1]], typeof)
  # list_lengths <- sapply(list_of_rows[[1]], length)

  list_of_cols <- lapply(seq_along(list_names), function(x) {
    colname <- list_names[[x]]
    list <- lapply(list_of_rows, function(z) z[[colname]])
    if (typeof(list[[1]]) != "list" && all(sapply(list, length) == 1)) {
      unlist(list, recursive = F)
    } else {
      list
    }
  }) %>% setNames(list_names)

  list_of_cols %>% as.tibble()
}

#' @export
load_datasets <- function() {
  datasets_info <- readRDS(paste0(.datasets_location, "/datasets.rds"))

  task_wrapped <- lapply(seq_len(nrow(datasets_info)), function(dataset_num) {
    dataset_id <- datasets_info$id[[dataset_num]]
    dataset <- dyngen::load_dataset(dataset_id, contents = dyngen::contents_dataset(experiment = F))

    with(dataset, dyneval::wrap_ti_task_data(
      ti_type = model$modulenetname,
      name = info$id,
      counts = counts,
      state_names = gs$milestone_names,
      state_net = gs$milestone_net,
      state_percentages = gs$milestone_percentages %>% slice(match(rownames(counts), id))
    ))
  })
  task_tib <- to_tibble(task_wrapped)
  task_tib %>% left_join(datasets_info, by = c("name"="id"))
}
