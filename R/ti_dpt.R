#' @import ParamHelpers
#' @export
description_dpt <- function() {
  list(
    name = "DPT",
    short_name = "DPT",
    package = c("destiny", "tidyverse"),
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
}

#' @import tidyverse
#' @import destiny
#'
#' @export
run_dpt <- function(counts,
                    sigma = "local",
                    distance = "euclidean",
                    n_eigs = 20,
                    density_norm = T,
                    n_local_lower = 5,
                    n_local_upper = 7,
                    w_width = .1) {
  n_local <- seq(n_local_lower, n_local_upper, by = 1)

  expr <- log2(counts+1)

  dm <- destiny::DiffusionMap(expr, sigma = sigma, distance = distance, n_eigs = n_eigs, density_norm = density_norm, n_local = n_local)
  dpt <- destiny::DPT(dm, w_width = w_width)

  tips <- tips(dpt)
  state_names <- paste0("DPT", tips)

  state_percentages <- bind_rows(lapply(state_names, function(x) data_frame(id = rownames(expr), state = x, percentage = dpt[[x]]))) %>%
    group_by(id) %>%
    mutate(percentage = percentage / sum(percentage)) %>%
    ungroup()

  state_network <- bind_rows(lapply(seq_along(state_names), function(i) {
    state <- state_names[[i]]
    index <- tips[[i]]
    data_frame(from = state, to = state_names, length = dpt[[state]][tips])
  })) %>% filter(from < to)

  ids <- rownames(expr)

  wrap_ti_prediction(
    ti_type = "triangle",
    name = "DPT",
    ids = ids,
    state_names = state_names,
    state_network = state_network,
    state_percentages = state_percentages
  )
}

#' @export
plot_dpt <- function(ti_predictions) {
  #qplot(percent_rank(ti_predictions$state_percentages[,1]), ti_predictions$state_percentages[,1], colour = data$sample_info$group.name)
  # todo
}

