#' Description for SCOUP
#' @export
description_scoup <- function() create_description(
  name = "SCOUP",
  short_name = "SCOUP",
  package_required = c(),
  package_loaded = c(),
  par_set = makeParamSet(
    makeIntegerParam(id = "ndim", lower = 2L, default = 2L, upper = 20L),
    makeIntegerParam(id = "nbranch", lower = 1L, default = 1L, upper = 20L),
    makeNumericParam(id = "max_ite1", lower = log(2), default = log(1000), upper = log(5000), trafo = function(x) round(exp(x))),
    makeNumericParam(id = "max_ite2", lower = log(2), default = log(10000), upper = log(50000), trafo = function(x) round(exp(x))),
    makeNumericParam(id = "alpha_min", lower = log(.001), default = log(.1), upper = log(10), trafo = exp),
    makeNumericParam(id = "alpha_max", lower = log(1), default = log(100), upper = log(10000), trafo = exp),
    makeNumericParam(id = "t_min", lower = log(.00001), default = log(.001), upper = log(1), trafo = exp),
    makeNumericParam(id = "t_max", lower = log(.1), default = log(2), upper = log(100), trafo = exp),
    makeNumericParam(id = "sigma_squared_min", lower = log(.001), default = log(.1), upper = log(10), trafo = exp),
    makeNumericParam(id = "thresh", lower = log(.01), default = log(.01), upper = log(10), trafo = exp)
  ),
  properties = c(),
  run_fun = run_scoup,
  plot_fun = plot_scoup
)

#' @importFrom utils read.table write.table
#' @importFrom stats var
run_scoup <- function(
  counts,
  cell_grouping,
  start_cell_id,
  ndim = 3,
  nbranch = 1,
  max_ite1 = 1000,
  max_ite2 = 10000,
  alpha_min = .1,
  alpha_max = 100,
  t_min = .001,
  t_max = 2,
  sigma_squared_min = .1,
  thresh = .01
) {
  requireNamespace("SCOUP")

  # figure out indices of starting population
  # from the cell_grouping and the start_cell_id
  start_ix <- cell_grouping %>%
    filter(cell_id == start_cell_id) %>%
    select(group_id) %>%
    left_join(cell_grouping, by = "group_id") %>%
    .$cell_id

  # log transform counts
  expr <- log2(counts + 1)

  # run SP and SCOUP
  model <- SCOUP::run_SCOUP(
    expr = expr,
    start_ix = start_ix,
    ndim = ndim,
    nbranch = nbranch,
    max_ite1 = max_ite1,
    max_ite2 = max_ite2,
    alpha_min = alpha_min,
    alpha_max = alpha_max,
    t_min = t_min,
    t_max = t_max,
    sigma_squared_min = sigma_squared_min,
    thresh = thresh
  ) %>%
    rownames_to_column("cell_id")

  maxtime <- max(model$time)

  # for visualisation purposes, run pca
  space <- stats::prcomp(expr)$x[,1:2]

  # create progressions
  progressions <- model %>%
    gather("to", "percentage", -cell_id, -time) %>%
    mutate(percentage = percentage * time / maxtime) %>%
    mutate(from = "M0") %>%
    select(cell_id, from, to, percentage)

  # create milestone ids
  milestone_ids <- c("M0", paste0("M", seq_len(nbranch)))

  # create milestone network
  milestone_network <- data_frame(
    from = milestone_ids[[1]],
    to = milestone_ids[-1],
    length = 1,
    directed = TRUE
  )

  # return output
  wrap_ti_prediction(
    ti_type = "multiforcating",
    id = "SCOUP",
    cell_ids = rownames(counts),
    milestone_ids = milestone_ids,
    milestone_network = milestone_network,
    progressions = progressions,
    model = model,
    space = space
  )
}

#' @importFrom cowplot theme_cowplot
plot_scoup <- function(prediction, type = c("dimred", "percentages")) {
  type <- match.arg(type)

  switch(
    type,
    dimred = {
      # find most likely milestone of each cell
      celltypes <- prediction$milestone_percentages %>%
        group_by(cell_id) %>%
        arrange(desc(percentage)) %>%
        slice(1) %>%
        ungroup()

      # combine data
      space_df <- prediction$space %>%
        as.data.frame %>%
        rownames_to_column(var = "cell_id") %>%
        left_join(celltypes, by = "cell_id")

      # make plot
      ggplot(space_df) +
        geom_point(aes(PC1, PC2, colour = milestone_id)) +
        cowplot::theme_cowplot() +
        labs(colour = "Milestone")
    },
    percentages = {
      ggplot(prediction$model %>% gather("endstate", "percentage", -cell_id, -time)) +
        geom_point(aes(time, percentage, color=endstate)) +
        facet_wrap(~endstate, ncol = 1) +
        cowplot::theme_cowplot()
    }
  )
}
