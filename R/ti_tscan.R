#' @export
makeRLearner.ti.tscan <- function() {
  makeRLearnerTI(
    cl = "ti.tscan",
    package = c("TSCAN"),
    par.set = makeParamSet(
      makeIntegerLearningParam(id = "preprocess_clusternum", lower = NULL, upper=1000 ,default = NULL),
      makeLogicalLearningParam(id = "preprocess_takelog", default = TRUE),
      makeIntegerLearningParam(id = "preprocess_logbase", lower = 2, upper=10 ,default = 2),
      makeIntegerLearningParam(id = "preprocess_pseudocount", lower = 1, upper=5 ,default = 1),
      makeIntegerLearningParam(id = "preprocess_minexpr_value", lower = 1, upper=5 ,default = 1),
      makeIntegerLearningParam(id = "preprocess_minexpr_percent", lower = 1, upper=5 ,default = .5),
      makeIntegerLearningParam(id = "preprocess_cvcutoff", lower = 1, upper=5 ,default = 1),
      makeIntegerLearningParam(id = "clustering_min_clusternum", lower = 2, upper=5 ,default = 2),
      makeIntegerLearningParam(id = "clustering_max_clusternum", lower = 7, upper=15 ,default = 9),
      makeLogicalLearningParam(id = "clustering_reduce", default = TRUE)
      #makeNumericVectorLearningParam(id = "ordering_MSTorder", default = TRUE),
    ),
    properties = c("linear", "dimred_samples", "pseudotime"),
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

  cds_2 <- TSCAN::exprmclust(t(as.matrix(expression)), clusternum= 2:9,reduce=clustering_reduce)

  cds_3 <- TSCAN::TSCANorder(cds_2)
  len <- length(cds_3)

  state_network <- data_frame(from = paste0("state_", cds_3[1]), to = paste0("state_", cds_3[len]), length = 1)
  state_percentages <- matrix(rep(0,len*2),ncol=2)
  state_percentages[,1] <- 1-(1:len/len)
  state_percentages[,2] <- 1:len/len
  state_names <- paste0("state_", c(cds_3[1],cds_3[len]))
  colnames(state_percentages) <- state_names
  rownames(state_percentages) <- cds_3
  state_percentages <- state_percentages[rownames(expression),]
  state_percentages_df <- data.frame(id = rownames(state_percentages), state_percentages, row.names = NULL, stringsAsFactors = F)

  wrap_ti_prediction(
    ti_type = "linear",
    name = "TSCAN",
    state_names = state_names,
    state_network = state_network,
    state_percentages = state_percentages_df,
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

