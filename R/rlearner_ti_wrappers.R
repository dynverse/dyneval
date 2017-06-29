#' @export
get_task_identifier <- function(task) {
  task[c("type", "ti_type", "name")]
}

#' @export
wrap_ti_task_data <- function(ti_type, name, counts, state_names, state_network, state_percentages, sample_info = NULL, feature_info = NULL, ...) {
  abstract_wrapper(
    "ti",
    ti_type,
    name,
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
wrap_ti_prediction <- function(ti_type, name, state_names, state_network, state_percentages, ...) {
  abstract_wrapper(
    "ti_pred",
    ti_type,
    name,
    state_names,
    state_network,
    state_percentages,
    ...
  )
}

abstract_wrapper <- function(type, ti_type, name, state_names, state_network, state_percentages, ...) {
  if (!is.data.frame(state_network) || any(colnames(state_network) != c("from", "to", "length"))) {
    stop(sQuote("state_network"), " should be a data frame with exactly three columns named ", sQuote("from"),
         ", ", sQuote("to"), " and ", sQuote("length"), ".")
  }
  if (any(!state_network$from %in% state_names) || any(!state_network$to %in% state_names)) {
    stop("Not all states in ", sQuote("state_network"), " are in ", sQuote("state_names"), ".")
  }
  if (!is.data.frame(state_percentages) || ncol(state_percentages) != length(state_names)+1) {
    stop(sQuote("state_network"), " should be a data frame with exactly N+1 columns, where N is the number of states as defined by ", sQuote("state_names"), ".")
  }
  if (!all(colnames(state_percentages) == c("id", state_names))) {
    stop("The column names of ", sQuote("state_network"), " should be c(\"id\", state_names).")
  }
  if (any(!state_names %in% colnames(state_percentages))) {
    stop("Not all states in ", sQuote("state_percentages"), " are in ", sQuote("state_names"), ".")
  }

  ## create a separate state if some cells have been filtered out
  state_pct_nas <- apply(state_percentages, 1, function(x) any(is.na(x)))
  if (any(state_pct_nas)) {
    state_percentages[state_pct_nas,-1] <- 0
    state_percentages[,"FILTERED_CELLS"] <- ifelse(state_pct_nas, 1, 0)
    state_network <- dplyr::bind_rows(
      state_network,
      data.frame(from = state_names, to = "FILTERED_CELLS", length = max(state_network$length)*5, stringsAsFactors = F)
    )
    state_names <- c(state_names, "FILTERED_CELLS")
  }

  l <- list(
    type = type,
    ti_type = ti_type,
    name = name,
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


