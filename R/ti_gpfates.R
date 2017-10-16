#' Description for GPfates
#' @export
description_gpfates <- function() create_description(
  name = "GPfates",
  short_name = "GPfates",
  package_loaded = c(),
  package_required = c("GPfates"),
  par_set = makeParamSet(
    makeNumericParam(id = "log_expression_cutoff", lower = 0.5, upper = 5, default = 2),
    makeNumericParam(id = "min_cells_expression_cutoff", lower = 0, upper = 20, default = 2),
    makeIntegerParam(id = "nfates", lower=1L, upper=20L, default=1L),
    makeIntegerParam(id = "ndims", lower=1L, upper=5L, default=2L)
  ),
  properties = c(),
  run_fun = run_gpfates,
  plot_fun = plot_gpfates
)

## TODO: give simulationtime as prior
#' @importFrom readr read_csv
#' @importFrom utils write.table
run_gpfates <- function(
  counts,
  nfates = 1,
  ndims = 2,
  log_expression_cutoff = 2,
  min_cells_expression_cutoff = 2
) {
  requireNamespace("GPfates")

  gp_out <- GPfates::GPfates(counts, nfates, ndims, log_expression_cutoff, min_cells_expression_cutoff)

  pseudotime <- gp_out$pseudotime
  phi <- gp_out$phi
  dr <- gp_out$dr

  pseudotime <- pseudotime %>% mutate(time = (time-min(time))/(max(time) - min(time)))

  # first get percentages of all final milestones, by getting their phi values, and multiplying by the pseudotime
  milestone_percentages <- phi %>%
    gather(milestone_id, percentage, -cell_id) %>%
    left_join(pseudotime, by="cell_id") %>%
    mutate(percentage = percentage*time)

  # now add the starting milestone, just 1-pseudotime
  milestone_percentages <-
    bind_rows(
      milestone_percentages,
      pseudotime %>% mutate(milestone_id="M0", percentage= 1-time)
    ) %>%
    select(-time)

  # create milestone network
  milestone_network <- data_frame(
    from = "M0",
    to = paste0("M", seq_len(nfates)),
    length = 1,
    directed = TRUE
  )
  milestone_ids <- paste0("M", seq(0, nfates))

  wrap_ti_prediction(
    ti_type = "GPfates",
    id = "GPfates",
    cell_ids = rownames(counts),
    milestone_ids = milestone_ids,
    milestone_network = milestone_network,
    milestone_percentages = milestone_percentages,
    dimred_samples = dr,
    pseudotime = pseudotime
  )
}

## TODO extract OMGP
plot_gpfates <- function(ti_predictions) {
  sample_df <- data.frame(
    ti_predictions$dimred_samples
  ) %>% left_join(ti_predictions$pseudotime, by = "cell_id")
  ggplot() +
    geom_point(aes(Comp1, Comp2, colour = time), sample_df) +
    coord_equal() +
    viridis::scale_color_viridis()
}
