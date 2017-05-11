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
