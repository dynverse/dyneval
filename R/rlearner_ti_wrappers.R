#' @export
get_task_identifier <- function(task) {
  task[c("type", "ti_type", "name")]
}

#' @export
wrap_ti_task_data <- function(ti_type, name, counts, ids, state_names, state_network, state_percentages, sample_info = NULL, feature_info = NULL, ...) {
  abstract_wrapper(
    "ti",
    ti_type,
    name,
    ids,
    state_names,
    state_network,
    state_percentages,
    counts = counts,
    sample_info = sample_info,
    feature_info = feature_info,
    ...
  )
}

#' @export
wrap_ti_prediction <- function(ti_type, name, ids, state_names, state_network, state_percentages, ...) {
  abstract_wrapper(
    "ti_pred",
    ti_type,
    name,
    ids,
    state_names,
    state_network,
    state_percentages,
    ...
  )
}

abstract_wrapper <- function(type, ti_type, name, ids, state_names, state_network, state_percentages, ...) {
  if (!is.data.frame(state_network) || ncol(state_network) != 3 || any(colnames(state_network) != c("from", "to", "length"))) {
    stop(sQuote("state_network"), " should be a data frame with exactly three columns named ", sQuote("from"),
         ", ", sQuote("to"), " and ", sQuote("length"), ".")
  }
  if (any(!state_network$from %in% state_names) || any(!state_network$to %in% state_names)) {
    stop("Not all states in ", sQuote("state_network"), " are in ", sQuote("state_names"), ".")
  }
  if (!is.data.frame(state_percentages) || ncol(state_percentages) != 3 || any(colnames(state_percentages) != c("id", "state", "percentage"))) {
    stop(sQuote("state_percentages"), " should be a data frame with exactly three columns named ", sQuote("id"),
         ", ", sQuote("state"), " and ", sQuote("percentage"), ".")
  }
  state_nam <- unique(state_percentages$state)
  if (!setequal(state_nam, state_names)) {
    stop("The set of all state names in ", sQuote("state_percentages"), " should be equal to ", sQuote("state_names"), ".")
  }

  ## create a separate state if some cells have been filtered out
  na_ids <- setdiff(ids, unique(state_percentages$id))
  if (length(na_ids) != 0) {
    state_percentages <- bind_rows(
      state_percentages,
      data_frame(id = na_ids, state = "FILTERED_CELLS", percentage = 1)
    )
    state_network <- dplyr::bind_rows(
      state_network,
      data_frame(from = state_names, to = "FILTERED_CELLS", length = max(state_network$length)*5)
    )
    state_names <- c(state_names, "FILTERED_CELLS")
  }

  l <- list(
    type = type,
    ti_type = ti_type,
    name = name,
    ids = ids,
    state_names = state_names,
    state_network = state_network,
    state_percentages = state_percentages,
    ...
  )
  class(l) <- paste0("dyneval::ti_wrapper")

  ## Precompute geodesic distances
  l$geodesic_dist <- compute_emlike_dist(l)

  l
}

is_ti_wrapper <- function(object) {
  class(object) == "dyneval::ti_wrapper"
}


