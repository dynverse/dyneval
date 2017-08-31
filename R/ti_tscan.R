#' @export
description_tscan <- function() {
  modelNames_values <- c("EII", "VII", "EEI", "VEI", "EVI", "VVI", "EEE", "EVE", "VEE", "VVE", "EEV", "VEV", "EVV", "VVV")
  list(
    name = "TSCAN",
    short_name = "TSCAN",
    package_load = c("mclust", "igraph", "ggplot2"),
    package_installed = c("TSCAN"),
    par_set = makeParamSet(
      # makeNumericParam(id = "minexpr_value", lower = 0, upper = 5, default = 1),
      # makeNumericParam(id = "minexpr_percent", lower = 0, upper = .5, default = .5),
      # makeNumericParam(id = "cvcutoff", lower = 0, upper = 5, default = 1),
      makeIntegerParam(id = "exprmclust_clusternum_lower", lower = 2L, upper = 20L, default = 2),
      makeIntegerParam(id = "exprmclust_clusternum_upper", lower = 2L, upper = 20L, default = 9),
      makeDiscreteParam(id = "modelNames", default = "VVV", values = modelNames_values),
      forbidden = quote(exprmclust_clusternum_lower > exprmclust_clusternum_upper)
    ),
    properties = c(),
    run_fun = run_tscan,
    plot_fun = plot_tscan
  )
}

run_tscan <- function(counts,
                      minexpr_percent = 1,
                      minexpr_value = .5,
                      cvcutoff = 1,
                      exprmclust_clusternum_lower = 2,
                      exprmclust_clusternum_upper = 9,
                      modelNames = "VVV") {
  requireNamespace("TSCAN")
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

  milestone_ids <- paste0("milestone_", c(head(cds_3, 1), tail(cds_3, 1)))
  milestone_network <- data_frame(from = milestone_ids[[1]], to = milestone_ids[[2]], length = 1)

  pseudotime <- setNames(percent_rank(match(rownames(counts), cds_3)), rownames(counts))
  milestone_percentages <- bind_rows(
    tibble::data_frame(cell_id = rownames(counts), milestone_id = milestone_ids[[1]], percentage = 1 - pseudotime),
    tibble::data_frame(cell_id = rownames(counts), milestone_id = milestone_ids[[2]], percentage = pseudotime)
  )

  milestone_percentages <- milestone_percentages %>% filter(!is.na(percentage))

  wrap_ti_prediction(
    ti_type = "linear",
    id = "TSCAN",
    cell_ids = rownames(counts),
    milestone_ids = milestone_ids,
    milestone_network = milestone_network,
    milestone_percentages = milestone_percentages,
    dimred_samples = cds_2$pcareduceres,
    dimred_clust = cds_2$clusterid,
    clust_centers = cds_2$clucenter,
    cds = cds_3
  )
}

plot_tscan <- function(ti_predictions) {
  qplot(percent_rank(ti_predictions$milestone_percentages[,1]), ti_predictions$milestone_percentages[,1], colour = data$sample_info$group.name)
}

