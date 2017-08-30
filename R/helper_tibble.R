#' Attempts to convert a list of lists to a tibble
#'
#' @param list_of_rows The list to be converted to a tibble
#'
#' @return A tibble with the same number of rows as there were elements in \code{list_of_rows}
#' @export
#'
#' @examples
#' l <- list(list(a = 1, b = log10), list(a = 2, b = sqrt))
#' tib <- list_as_tibble(l)
#' tib
list_as_tibble <- function(list_of_rows) {
  list_names <- names(list_of_rows[[1]])

  list_of_cols <- lapply(seq_along(list_names), function(x) {
    colname <- list_names[[x]]
    list <- lapply(list_of_rows, function(z) z[[colname]])
    if (typeof(list[[1]]) != "list" && all(sapply(list, length) == 1)) {
      unlist(list, recursive = F)
    } else {
      list
    }
  }) %>% setNames(list_names)

  list_of_cols %>% as_tibble()
}

#' Extracts one row from a tibble and converts it to a list
#'
#' @param tib the tibble
#' @param row_id the index of the row to be selected
#'
#' @return the corresponding row from the tibble as a list
#' @export
#'
#' @examples
#' l <- list(list(a = 1, b = log10), list(a = 2, b = sqrt))
#' tib <- list_as_tibble(l)
#'
#' extract_row_to_list(tib, 2)
extract_row_to_list <- function(tib, row_id) {
  tib[row_id, ] %>% as.list %>% map(function(x) {
    if (is.null(x) | !is.list(x)) {
      x
    } else {
      x[[1]]
    }
  })
}
