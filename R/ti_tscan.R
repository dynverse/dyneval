#' @import ParamHelpers
#' @export
description_tscan <- function() {
  list(
    name = "TSCAN",
    short_name = "TSCAN",
    package = c("TSCAN", "mclust", "igraph"),
    par_set = makeParamSet(
      # makeNumericParam(id = "minexpr_value", lower = 0, upper = 5, default = 1),
      # makeNumericParam(id = "minexpr_percent", lower = 0, upper = .5, default = .5),
      # makeNumericParam(id = "cvcutoff", lower = 0, upper = 5, default = 1),
      makeIntegerParam(id = "exprmclust_clusternum_lower", lower = 2L, upper = 20L, default = 2),
      makeIntegerParam(id = "exprmclust_clusternum_upper", lower = 2L, upper = 20L, default = 9),
      makeDiscreteParam(id = "modelNames", default = "VVV", values = mclust::mclust.options("emModelNames")),
      forbidden = quote(exprmclust_clusternum_lower > exprmclust_clusternum_upper)
    ),
    properties = c(),
    run_fun = run_tscan,
    plot_fun = plot_tscan
  )
}

#' @import tidyverse
#' @import TSCAN
#'
#' @export
run_tscan <- function(counts,
                      minexpr_percent = 1,
                      minexpr_value = .5,
                      cvcutoff = 1,
                      exprmclust_clusternum_lower = 2,
                      exprmclust_clusternum_upper = 9,
                      modelNames = "VVV") {
  expr <- t(as.matrix(counts))
  # cds_1 <- TSCAN::preprocess(
  #   expr,
  #   takelog = T,
  #   logbase = 2,
  #   clusternum = preprocess_clusternum,
  #   pseudocount = pseudocount,
  #   minexpr_value = minexpr_value,
  #   minexpr_percent = minexpr_percent,
  #   cvcutoff = cvcutoff)
  cds_1 <- TSCAN::preprocess(
    expr,
    takelog = T,
    logbase = 2,
    pseudocount = 1,
    clusternum = NULL,
    minexpr_value = 0,
    minexpr_percent = 0,
    cvcutoff = 0)

  exprmclust_clusternum <- seq(exprmclust_clusternum_lower, exprmclust_clusternum_upper, 1)
  cds_2 <- TSCAN::exprmclust(
    cds_1,
    clusternum = exprmclust_clusternum,
    modelNames = modelNames,
    reduce = T)

  cds_3 <- TSCAN::TSCANorder(cds_2)

  ids <- rownames(counts)
  state_names <- paste0("state_", c(head(cds_3, 1), tail(cds_3, 1)))
  state_network <- data_frame(from = state_names[[1]], to = state_names[[2]], length = 1)

  pseudotime <- setNames(percent_rank(match(rownames(counts), cds_3)), rownames(counts))
  state_percentages <- bind_rows(
    tibble::data_frame(id = rownames(counts), state = state_names[[1]], percentage = 1 - pseudotime),
    tibble::data_frame(id = rownames(counts), state = state_names[[2]], percentage = pseudotime)
  )

  wrap_ti_prediction(
    ti_type = "linear",
    name = "TSCAN",
    ids = ids,
    state_names = state_names,
    state_network = state_network,
    state_percentages = state_percentages,
    dimred_samples = cds_2$pcareduceres,
    dimred_clust = cds_2$clusterid,
    clust_centers = cds_2$clucenter,
    cds = cds_3
  )
}

#' @export
plot_tscan <- function(ti_predictions) {
  qplot(percent_rank(ti_predictions$state_percentages[,1]), ti_predictions$state_percentages[,1], colour = data$sample_info$group.name)
}

