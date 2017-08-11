#' @import tidyverse
generate_linear <- function(ncells = 100) {
  milestone_network <- tibble(from="M1", to="M2", length=1)
  milestone_ids <- unique(c(milestone_network$from, milestone_network$to))
  cell_ids <- paste0("C", seq_len(ncells))

  progressions <- tibble(cell_id = cell_ids) %>%
    mutate(
      from="M1",
      to="M2",
      percentage=runif(nrow(.))
    )

  task <- dyneval::wrap_ti_prediction(
    "linear",
    "linear",
    cell_ids,
    milestone_ids,
    milestone_network,
    progressions = progressions
  )

  task
}

#' @import tidyverse
generate_bifurcating <- function(ncells = 100) {
  milestone_network <- tibble(from=c("M1", "M2", "M2"), to=c("M2", "M3", "M4"), length=1)
  milestone_ids <- unique(c(milestone_network$from, milestone_network$to))
  cell_ids <- paste0("C", seq_len(ncells))

  progressions <- tibble(cell_id = cell_ids) %>%
    bind_cols(milestone_network[sample(seq_len(nrow(milestone_network)), length(cell_ids), replace=TRUE), ]) %>%
    mutate(percentage = map_dbl(length, ~runif(1, 0, .)))

  task <- dyneval::wrap_ti_prediction(
    "bifurcating",
    "bifurcating",
    cell_ids,
    milestone_ids,
    milestone_network,
    progressions = progressions %>% select(cell_id, from, to, percentage)
  )

  task
}
