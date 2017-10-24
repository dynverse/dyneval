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

  # transform counts
  expr <- log2(counts+1)

  # run diffusion maps
  dm <- destiny::DiffusionMap(
    expr,
    sigma = sigma,
    distance = distance,
    n_eigs = n_eigs,
    density_norm = density_norm,
    n_local = n_local
  )

  # run DPT
  dpt_params <- lst(dm, w_width)
  if (!is.null(start_cell_id)) {
    dpt_params$tips <- which(rownames(counts) %in% start_cell_id)
  }
  dpt <- do.call(destiny::DPT, dpt_params)

  # find DPT tips
  tips <- destiny::tips(dpt)
  milestone_ids <- paste0("DPT", tips)

  # transform to percentages
  milestone_percentages <- map_df(milestone_ids, function(mid) {
    data_frame(cell_id = rownames(expr), milestone_id = mid, percentage = 1 - dynutils::scale_minmax(dpt[[mid]]))
  }) %>%
    group_by(cell_id) %>%
    mutate(percentage = percentage / sum(percentage)) %>%
    ungroup()

  # generate milestone network
  milestone_network <- crossing(from = milestone_ids, to = milestone_ids) %>%
    filter(from < to) %>%
    rowwise() %>%
    mutate(
      length = dpt[[from]][[setNames(tips, milestone_ids)[[to]]]],
      directed = FALSE
    ) %>%
    ungroup()

  # extract dimred for visualisation
  space <- dpt@dm@eigenvectors[,1:3] %>%
    as.data.frame %>%
    mutate(
      is_tip = dpt@tips[,1],
      branch = dpt@branch[,1],
      col_lab = ifelse(is_tip, "Tip", ifelse(is.na(branch), "Unassigned", paste0("Branch ", branch)))
    )

  # return output
  wrap_ti_prediction(
    ti_type = "triangle",
    id = "DPT",
    cell_ids = rownames(expr),
    milestone_ids = milestone_ids,
    milestone_network = milestone_network,
    milestone_percentages = milestone_percentages,
    space = space
  )
}

plot_dpt <- function(prediction) {
  # based on destiny::plot.DPT(prediction$dpt, col_by = "branch")

  palette <- c("#8DD3C7", "#FFED6F", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#BC80BD", "#FCCDE5", "gray85", "#CCEBC5", "#FFFFB3")
  ann_cols <- c(
    setNames(palette, paste0("Branch ", seq_along(palette))),
    Unassigned = "lightgray",
    Tip = "red"
  )

  g <- ggplot(prediction$space) +
    geom_point(aes(DC1, DC2, colour = col_lab), size = 2) +
    scale_colour_manual(values = ann_cols) +
    labs(colour = "Branch") +
    theme(legend.position = c(0.9, 0.1))
  process_dyneval_plot(g, prediction$id)
}

