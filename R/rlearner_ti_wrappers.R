#' @export
get_task_identifier <- function(task) {
  task[c("type", "ti_type", "id")]
}

#' @export
wrap_ti_task_data <- function(
  ti_type,
  id,
  counts,
  cell_ids,
  milestone_ids,
  milestone_network,
  milestone_percentages = NULL,
  progressions = NULL,
  sample_info = NULL,
  feature_info = NULL,
  ...
) {
  abstract_wrapper(
    "ti",
    ti_type,
    id,
    cell_ids,
    milestone_ids,
    milestone_network,
    milestone_percentages,
    progressions,
    counts = counts,
    sample_info = sample_info,
    feature_info = feature_info,
    ...
  )
}

#' @export
wrap_ti_prediction <- function(
  ti_type,
  id,
  cell_ids,
  milestone_ids,
  milestone_network,
  milestone_percentages = NULL,
  progressions = NULL,
  ...
) {
  abstract_wrapper(
    "ti_pred",
    ti_type,
    id,
    cell_ids,
    milestone_ids,
    milestone_network,
    milestone_percentages,
    progressions,
    ...
  )
}

abstract_wrapper <- function(
  type,
  ti_type,
  id,
  cell_ids,
  milestone_ids,
  milestone_network,
  milestone_percentages = NULL,
  progressions = NULL,
  ...
) {
  if (!is.data.frame(milestone_network) || ncol(milestone_network) != 3 || any(colnames(milestone_network) != c("from", "to", "length"))) {
    stop(sQuote("milestone_network"), " should be a data frame with exactly three columns named ", sQuote("from"),
         ", ", sQuote("to"), " and ", sQuote("length"), ".")
  }
  if (any(!milestone_network$from %in% milestone_ids) || any(!milestone_network$to %in% milestone_ids)) {
    stop("Not all states in ", sQuote("milestone_network"), " are in ", sQuote("milestone_ids"), ".")
  }

  if (is.null(milestone_percentages) == is.null(progressions)) {
    stop("Exactly one of ", sQuote("milestone_percentages"), " or ", sQuote("progressions"), " must be defined, the other must be NULL.")
  }

  if (is.null(progressions)) {
    progressions <- convert_milestone_percentages_to_progressions(cell_ids, milestone_ids, milestone_network, milestone_percentages)
  } else if (is.null(milestone_percentages)) {
    milestone_percentages <- convert_progressions_to_milestone_percentages(cell_ids, milestone_ids, milestone_network, progressions)
  }

  if (!is.data.frame(milestone_percentages) || ncol(milestone_percentages) != 3 || any(colnames(milestone_percentages) != c("cell_id", "milestone_id", "percentage"))) {
    stop(sQuote("milestone_percentages"), " should be a data frame with exactly three columns named ", sQuote("cell_id"),
         ", ", sQuote("milestone_id"), " and ", sQuote("percentage"), ".")
  }
  if (!is.data.frame(progressions) || ncol(progressions) != 4 || any(colnames(progressions) != c("cell_id", "from", "to", "percentage"))) {
    stop(sQuote("progressions"), " should be a data frame with exactly four columns named ", sQuote("cell_id"),
         ", ", sQuote("from"), ", ", sQuote("to"), " and ", sQuote("percentage"), ".")
  }

  ## create a separate state if some cells have been filtered out
  na_ids <- setdiff(cell_ids, unique(milestone_percentages$cell_id))
  if (length(na_ids) != 0) {
    milestone_percentages <- bind_rows(
      milestone_percentages,
      data_frame(cell_id = na_ids, milestone_id = "FILTERED_CELLS", percentage = 1)
    )
    milestone_network <- dplyr::bind_rows(
      milestone_network,
      data_frame(from = milestone_ids, to = "FILTERED_CELLS", length = max(milestone_network$length)*5)
    )
    milestone_ids <- c(milestone_ids, "FILTERED_CELLS")
  }

  # create output structure
  out <- list(
    type = type,
    ti_type = ti_type,
    id = id,
    cell_ids = cell_ids,
    milestone_ids = milestone_ids,
    milestone_network = milestone_network,
    milestone_percentages = milestone_percentages,
    progressions = progressions,
    ...
  )
  class(out) <- paste0("dyneval::ti_wrapper")

  ## Precompute geodesic distances
  #out$geodesic_dist <- compute_emlike_dist(out)
  # do not precompute geodesic distances :)

  out
}

convert_milestone_percentages_to_progressions <- function(cell_ids, milestone_ids, milestone_network, milestone_percentages) {
  bind_rows(lapply(cell_ids, function(cid) {
    relevant_pct <- milestone_percentages %>% filter(cell_id == cid)

    if (nrow(relevant_pct) >= 2) {
      relevant_progr <- milestone_network %>%
        filter(from %in% relevant_pct$milestone_id & to %in% relevant_pct$milestone_id) %>%
        left_join(relevant_pct, by = c("to" = "milestone_id")) %>%
        select(cell_id, from, to, percentage)
    } else if (nrow(relevant_pct) == 1) {
      relevant_net <- milestone_network %>% filter(to %in% relevant_pct$milestone_id)
      if (nrow(relevant_net) == 0) {
        relevant_net <- milestone_network %>% filter(from %in% relevant_pct$milestone_id)
      }
      relevant_progr <- relevant_net %>%
        mutate(cell_id = cid, percentage = 1) %>%
        select(cell_id, from, to, percentage)
    } else {
      relevant_progr <- NULL
    }
    relevant_progr
  }))
}

convert_progressions_to_milestone_percentages <- function(cell_ids, milestone_ids, milestone_network, progressions) {
  check_froms <- progressions %>% group_by(cell_id) %>% summarise(n = length(unique(from)))
  if (any(check_froms$n > 1)) {
    stop("In ", sQuote("progressions"), ", cells should only have 1 unique from milestone.")
  }

  bind_rows(lapply(cell_ids, function(cid) {
    relevant_prg <- progressions %>% filter(cell_id == cid)

    # there should be only 1 unique from, according to an earlier check.
    from <- unique(relevant_prg$from)

    bind_rows(
      data_frame(cell_id = cid, milestone_id = from, percentage = 1 - sum(relevant_prg$percentage)),
      data_frame(cell_id = cid, milestone_id = relevant_prg$to, percentage = relevant_prg$percentage)
    )
  })) %>% filter(percentage > 0)
}

is_ti_wrapper <- function(object) {
  class(object) == "dyneval::ti_wrapper"
}


