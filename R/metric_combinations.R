#' Calculate combinations of metrics
#'
#' @param ... Can be:
#'   - One numeric vector
#'   - A list containg numeric vectors
#'   - Numeric vectors given as separate inputs
#'
#' @rdname calculate_combinations
#'
#' @examples
#' calculate_harmonic_mean(0.1, 0.5, 0.9)
#' calculate_geometric_mean(0.1, 0.5, 0.9)
#' calculate_arithmetic_mean(0.1, 0.5, 0.9)
#'
#' calculate_harmonic_mean(c(0.1, 0.9), c(0.2, 1))
#' calculate_arithmetic_mean(c(0.1, 10), c(0.9, 20))
calculate_harmonic_mean <- function(...) {
  x <- process_combination_input(...)
  ncol(x) / rowSums(1/x)
}

#' @rdname calculate_combinations
calculate_geometric_mean <- function(...) {
  x <- process_combination_input(...)
  apply(x, 1, prod)^(1/ncol(x))
}

#' @rdname calculate_combinations
calculate_arithmetic_mean <- function(...) {
  x <- process_combination_input(...)
  rowSums(x)/ncol(x)
}

# Processes:
#   - process_combination_input(list(1, 2, 3))
#   - process_combination_input(c(1, 2, 3))
#   - process_combination_input(1, 2, 3)
# all to matrix(c(1, 2, 3))
process_combination_input <- function(...) {
  dots <- list(...)
  if (length(dots) > 1 && all(map_lgl(dots, is.numeric))) {
    do.call(cbind, dots)
  } else if (is.numeric(..1)) {
    do.call(cbind, as.list(..1))
  } else if (is.list(..1) && all(map_lgl(..1, is.numeric))) {
    do.call(cbind, ..1)
  } else if (is.matrix(..1) && is.numeric(..1)) {
    ..1
  } else {
    stop("Invalid input")
  }
}
