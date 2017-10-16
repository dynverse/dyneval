#' Description for DPT
#' @export
description_dpt <- function() create_description(
  name = "DPT",
  short_name = "DPT",
  package_loaded = c("destiny"),
  package_required = c("dynutils", "reshape2"),
  par_set = makeParamSet(
    makeDiscreteParam(id = "sigma", default = "local", values = c("local", "global")),
    makeDiscreteParam(id = "distance", default = "euclidean", values = c("euclidean", "cosine", "rankcor")),
    makeIntegerParam(id = "n_eigs", lower = 3L, upper = 100L, default = 20L),
    makeLogicalParam(id = "density_norm", default = TRUE),
    makeIntegerParam(id = "n_local_lower", lower = 2L, upper = 20L, default = 5L),
    makeIntegerParam(id = "n_local_upper", lower = 2L, upper = 20L, default = 7L),
    makeNumericParam(id = "w_width", lower = -4, upper = 0, default = log(.1), trafo = exp),
    forbidden = quote(n_local_lower > n_local_upper)
  ),
  properties = c(),
  run_fun = run_dpt,
  plot_fun = plot_dpt
)

#' @importFrom reshape2 melt
run_dpt <- function(counts,
                    start_cell_id = NULL,
                    sigma = "local",
                    distance = "euclidean",
                    n_eigs = 20,
                    density_norm = TRUE,
                    n_local_lower = 5,
                    n_local_upper = 7,
                    w_width = .1) {
  requireNamespace("destiny")

  # create n_local vector
  n_local <- seq(n_local_lower, n_local_upper, by = 1)

  expr <- log2(counts+1)

  dm <- destiny::DiffusionMap(
    expr,
    sigma = sigma,
    distance = distance,
    n_eigs = n_eigs,
    density_norm = density_norm,
    n_local = n_local
  )

  dpt_params <- lst(dm, w_width)

  if (!is.null(start_cell_id)) {
    dpt_params$tips <- which(rownames(counts) %in% start_cell_id)
  }

  dpt <- do.call(destiny::DPT, dpt_params)

  tips <- destiny::tips(dpt)
  milestone_ids <- paste0("DPT", tips)

  milestone_percentages <- map_df(milestone_ids, function(mid) {
    d <- dpt[[mid]]
    dsc <- dynutils::scale_minmax(d)
    data_frame(cell_id = rownames(expr), milestone_id = mid, percentage = 1 - dsc)
  }) %>%
    group_by(cell_id) %>%
    mutate(percentage = percentage / sum(percentage)) %>%
    ungroup()

  milestone_network <- crossing(from = milestone_ids, to = milestone_ids) %>%
    filter(from < to) %>%
    rowwise() %>%
    mutate(
      length = dpt[[from]][[setNames(tips, milestone_ids)[[to]]]],
      directed = FALSE
    ) %>%
    ungroup()

  wrap_ti_prediction(
    ti_type = "triangle",
    id = "DPT",
    cell_ids = rownames(expr),
    milestone_ids = milestone_ids,
    milestone_network = milestone_network,
    milestone_percentages = milestone_percentages,
    dpt = dpt
  )
}

plot_dpt <- function(prediction) {
  requireNamespace("destiny")
  destiny::plot.DPT(prediction$dpt)
}

