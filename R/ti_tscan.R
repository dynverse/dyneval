#' @import ParamHelpers
#' @importFrom mclust mclust.options
#' @export
makeRLearner.ti.tscan <- function() {
  makeRLearnerTI(
    cl = "ti.tscan",
    package = c("TSCAN"),
    par.set = makeParamSet(
      makeIntegerLearnerParam(id = "pseudocount", lower = 1, upper = 5, default = 1),
      makeNumericLearnerParam(id = "minexpr_value", lower = 0, upper = 20, default = 1),
      makeNumericLearnerParam(id = "minexpr_percent", lower = 0, upper = 1, default = .5),
      makeNumericLearnerParam(id = "cvcutoff", lower = 0, upper = 20, default = 1),
      makeIntegerVectorLearnerParam(id = "clusternum", len = as.integer(NA), lower = 2, upper = 5, default = 2),
      makeDiscreteLearnerParam(id = "modelNames", default = "VVV", values = mclust::mclust.options("emModelNames"))
    ),
    #properties = c("dimred", "pseudotime"),
    properties = c("tibble"),#, "dimred", "dimred_traj", "pseudotime"), # todo: implement other outputs
    name = "TSCAN",
    short.name = "TSCAN"
  )
}

#' @import tidyverse
#' @import TSCAN
#'
#' @export
trainLearner.ti.tscan <- function(.learner, .task, .subset, ...) {
  NULL
}

#' @export
predictLearner.ti.tscan <- function(.learner, .model, .newdata, ...) {
  outs <- lapply(.newdata$expression, run.ti.tscan, ...)
  outs_tib <- to_tibble(outs)
  outs_tib
}

run.ti.tscan <- function(count, num_dimensions, pseudocount, minexpr_percent,
                         minexpr_value, cvcutoff, clusternum, modelNames) {
  cds_1 <- TSCAN::preprocess(
    t(as.matrix(count)), takelog = T, logbase = 2,
    clusternum = NULL , pseudocount = pseudocount, # clusternum needs to be added
    minexpr_value = preprocess_minexpr_value, minexpr_percent = minexpr_value,
    cvcutoff = cvcutoff)

  cds_2 <- TSCAN::exprmclust(cds_1, clusternum = clusternum, modelNames = modelNames, reduce = T)

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

