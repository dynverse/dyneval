#' Description for Wishbone
#' @export
description_wishbone <- function() create_description(
  name = "Wishbone",
  short_name = "Wishbone",
  package_loaded = c(),
  package_required = c("jsonlite", "Wishbone"),
  par_set = makeParamSet(
    makeIntegerParam(id = "knn", lower=2L, upper=100L, default=10L),
    makeIntegerParam(id = "n_diffusion_components", lower=2L, upper=20L, default=10L),
    makeIntegerParam(id = "n_pca_components", lower=2L, upper=30L, default=15L),
    makeLogicalParam(id = "branch", default = TRUE),
    makeIntegerParam(id = "k", lower=2L, upper=100L, default=15L),
    makeIntegerParam(id = "num_waypoints", lower=2L, upper=500L, default=250L),
    makeLogicalParam(id = "normalize", default = TRUE),
    makeNumericParam(id = "epsilon", lower=0.1, upper=10, default=1)
  ),
  properties = c(),
  run_fun = run_wishbone,
  plot_fun = plot_wishbone
)

run_wishbone <- function(
  counts,
  start_cell_id,
  knn = 10,
  n_diffusion_components = 2,
  n_pca_components = 15,
  markers="~",
  branch=TRUE,
  k=15,
  num_waypoints=50,
  normalize=TRUE,
  epsilon=1
) {
  if(is.null(start_cell_id)) stop("Give start cell id")

  requireNamespace("Wishbone")

  wb_out <- Wishbone::Wishbone(
    counts = counts,
    start_cell_id = start_cell_id,
    knn = knn,
    n_diffusion_components = n_diffusion_components,
    n_pca_components = n_pca_components,
    markers = markers,
    branch = branch,
    k = k,
    num_waypoints = num_waypoints,
    normalize = normalize,
    epsilon = epsilon
  )

  branch_assignment <- wb_out$branch_assignment
  trajectory <- wb_out$trajectory
  space <- wb_out$space

  model <- left_join(branch_assignment, trajectory, by="cell_id")

  if (branch) {
    milestone_network <- tibble(from=c("M1", "M2", "M2"), to=c("M2", "M3", "M4"), branch=c(1, 2, 3))
  } else {
    milestone_network <- tibble(from=c("M1"), to=c("M2"), branch=c(1))
  }
  milestone_network <- milestone_network %>% mutate(directed=TRUE)

  progressions <- left_join(model, milestone_network, by="branch")

  # get lengths of milestone network
  milestone_network <- progressions %>% group_by(branch) %>% summarise(length=max(time) - min(time)) %>% left_join(milestone_network, by="branch")

  # now scale the times between 0 and 1 => percentages
  progressions <- progressions %>%
    group_by(branch) %>%
    mutate(percentage=(time - min(time))/(max(time) - min(time))) %>%
    ungroup()

  milestone_ids <- unique(c(milestone_network$from, milestone_network$to))

  wrap_ti_prediction(
    ti_type = "linear",
    id = "Wishbone",
    cell_ids = rownames(counts),
    milestone_ids = milestone_ids,
    milestone_network = milestone_network %>% select(from, to, length, directed),
    progressions = progressions %>% select(cell_id, from, to, percentage),
    space = space,
    model=model
  )
}

plot_wishbone <- function(ti_predictions) {
  ggplot(left_join(ti_predictions$model, ti_predictions$space, by="cell_id")) + geom_point(aes(Comp0, Comp1, color=branch))
}
