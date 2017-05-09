#' @import ParamHelpers
#' @export
makeRLearner.ti.tscan <- function() {
  makeRLearnerTI(
    cl = "ti.tscan",
    package = c("TSCAN"),
    par.set = makeParamSet(
      makeIntegerLearnerParam(id = "preprocess_clusternum", lower = NULL, upper=1000 ,default = NULL),
      makeLogicalLearnerParam(id = "preprocess_takelog", default = TRUE),
      makeIntegerLearnerParam(id = "preprocess_logbase", lower = 2, upper=10 ,default = 2),
      makeIntegerLearnerParam(id = "preprocess_pseudocount", lower = 1, upper=5 ,default = 1),
      makeIntegerLearnerParam(id = "preprocess_minexpr_value", lower = 1, upper=5 ,default = 1),
      makeIntegerLearnerParam(id = "preprocess_minexpr_percent", lower = 1, upper=5 ,default = .5),
      makeIntegerLearnerParam(id = "preprocess_cvcutoff", lower = 1, upper=5 ,default = 1),
      makeIntegerLearnerParam(id = "clustering_min_clusternum", lower = 2, upper=5 ,default = 2),
      makeIntegerLearnerParam(id = "clustering_max_clusternum", lower = 7, upper=15 ,default = 9),
      makeLogicalLearnerParam(id = "clustering_reduce", default = TRUE)
      #makeNumericVectorLearnerParam(id = "ordering_MSTorder", default = TRUE),
    ),
    properties = c("dimred", "pseudotime"),
    name = "TSCAN",
    short.name = "TSCAN"
  )
}

#' @import tidyverse
#' @import TSCAN
#'
#' @export
trainLearner.ti.tscan <- function(.task, num_dimensions,preprocess_clusternum, preprocess_takelog,
                                  preprocess_logbase, preprocess_pseudocount,preprocess_minexpr_value,
                                  preprocess_minexpr_percent, preprocess_cvcutoff, clustering_min_clusternum,
                                  clustering_max_clusternum, clustering_reduce ) {
  # subsetting will not work yet, but the function is already provided
  data <- get_task_data(.task)

  expression <- data$expression

  # cds_1 <- TSCAN::preprocess(t(as.matrix(expression)),clusternum = preprocess_clusternum,
  #                            takelog = F, logbase = preprocess_logbase,
  #                            pseudocount =preprocess_pseudocount, minexpr_value = preprocess_minexpr_value,
  #                            minexpr_percent = preprocess_minexpr_percent, cvcutoff = preprocess_cvcutoff)

  clusternums <- clustering_min_clusternum:clustering_max_clusternum
  cds_2 <- TSCAN::exprmclust(t(as.matrix(expression)), clusternum = clusternums, reduce = clustering_reduce)

  cds_3 <- TSCAN::TSCANorder(cds_2)
  # problem: not all cells are in cds_3!

  state_names <- paste0("state_", c(head(cds_3, 1), tail(cds_3, 1)))
  state_network <- data_frame(from = state_names[[1]], to = state_names[[2]], length = 1)

  pseudotime <- setNames(percent_rank(match(rownames(expression), cds_3)), rownames(expression))
  state_percentages <- data.frame(id = names(pseudotime), from = 1-pseudotime, to = pseudotime)
  colnames(state_percentages) <- c("id", state_names)

  wrap_ti_prediction(
    ti_type = "linear",
    name = "TSCAN",
    state_names = state_names,
    state_network = state_network,
    state_percentages = state_percentages,
    task_id = get_task_identifier(.task),
    dimred_samples = cds_2$pcareduceres,
    dimred_clust = cds_2$clusterid,
    clust_centers = cds_2$clucenter,
    cds = cds_3
  )
}

#' @export
plotLearner.ti.TSCAN <- function(ti_predictions) {
  qplot(percent_rank(ti_predictions$state_percentages[,1]), ti_predictions$state_percentages[,1], colour = data$sample_info$group.name)
}

