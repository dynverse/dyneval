#' Description for DPT
#' @export
description_dpt <- function() create_description(
  name = "DPT",
  short_name = "DPT",
  package_loaded = c("destiny"),
  package_required = c(),
  par_set = makeParamSet(
    makeDiscreteParam(id = "sigma", default = "local", values = c("local", "global")),
    makeDiscreteParam(id = "distance", default = "euclidean", values = c("euclidean", "cosine", "rankcor")),
    makeIntegerParam(id = "n_eigs", lower = 3, upper = 100, default = 20),
    makeLogicalParam(id = "density_norm", default = T),
    makeIntegerParam(id = "n_local_lower", lower = 2L, upper = 20L, default = 5),
    makeIntegerParam(id = "n_local_upper", lower = 2L, upper = 20L, default = 7),
    makeNumericParam(id = "w_width", lower = -4, upper = 0, default = log(.1), trafo = exp),
    forbidden = quote(n_local_lower > n_local_upper)# | (sigma == "global" && distance %in% c("cosine", "rankcor")))
  ),
  properties = c(),
  run_fun = run_dpt,
  plot_fun = plot_dpt
)

run_dpt <- function(counts,
                    sigma = "local",
                    distance = "euclidean",
                    n_eigs = 20,
                    density_norm = T,
                    n_local_lower = 5,
                    n_local_upper = 7,
                    w_width = .1) {
  requireNamespace("destiny")

  n_local <- seq(n_local_lower, n_local_upper, by = 1)

  expr <- log2(counts+1)

  dm <- destiny::DiffusionMap(expr, sigma = sigma, distance = distance, n_eigs = n_eigs, density_norm = density_norm, n_local = n_local)
  dpt <- destiny::DPT(dm, w_width = w_width)

  tips <- destiny::tips(dpt)
  milestone_ids <- paste0("DPT", tips)

  milestone_percentages <- bind_rows(lapply(milestone_ids, function(x) data_frame(cell_id = rownames(expr), milestone_id = x, percentage = dpt[[x]]))) %>%
    group_by(cell_id) %>%
    mutate(percentage = percentage / sum(percentage)) %>%
    ungroup()

  milestone_network <- bind_rows(lapply(seq_along(milestone_ids), function(i) {
    from <- milestone_ids[[i]]
    index <- tips[[i]]
    data_frame(from = from, to = milestone_ids, length = dpt[[from]][tips])
  })) %>% filter(row_number() %in% c(2, 3)) %>% mutate(directed=FALSE)

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

