#' Description for TSCAN
#' @export
description_tscan <- function() create_description(
  name = "TSCAN",
  short_name = "TSCAN",
  package_loaded = c("mclust", "igraph", "ggplot2"),
  package_required = c("TSCAN"),
  par_set = makeParamSet(
    makeNumericParam(id = "minexpr_percent", lower=0, upper=1, default=0),
    makeNumericParam(id = "minexpr_value", lower=0, upper=10, default=0),
    makeNumericParam(id = "cvcutoff", lower=0, upper=5, default=0),
    makeIntegerParam(id = "exprmclust_clusternum_lower", lower = 2L, upper = 20L, default = 2L),
    makeIntegerParam(id = "exprmclust_clusternum_upper", lower = 2L, upper = 20L, default = 9L, requires = expression(exprmclust_clusternum_lower <= exprmclust_clusternum_upper)),
    makeDiscreteParam(id = "modelNames", default = "VVV", values = c("EII", "VII", "EEI", "VEI", "EVI", "VVI", "EEE", "EVE", "VEE", "VVE", "EEV", "VEV", "EVV", "VVV"))
  ),
  properties = c(),
  run_fun = run_tscan,
  plot_fun = plot_tscan
)

run_tscan <- function(counts,
                      minexpr_percent = 0,
                      minexpr_value = 0,
                      cvcutoff = 0,
                      exprmclust_clusternum_lower = 2,
                      exprmclust_clusternum_upper = 9,
                      modelNames = "VVV"
  ) {
  requireNamespace("TSCAN")
  expr <- t(as.matrix(counts))
  # cds_1 <- TSCAN::preprocess(
  #   expr,
  #   takelog = TRUE,
  #   logbase = 2,
  #   clusternum = preprocess_clusternum,
  #   pseudocount = pseudocount,
  #   minexpr_value = minexpr_value,
  #   minexpr_percent = minexpr_percent,
  #   cvcutoff = cvcutoff)
  cds_1 <- TSCAN::preprocess(
    expr,
    takelog = TRUE,
    logbase = 2,
    pseudocount = 1,
    clusternum = NULL,
    minexpr_value = minexpr_value,
    minexpr_percent = minexpr_percent,
    cvcutoff = cvcutoff
  )

  exprmclust_clusternum <- seq(exprmclust_clusternum_lower, exprmclust_clusternum_upper, 1)
  cds_2 <- TSCAN::exprmclust(
    cds_1,
    clusternum = 2,
    modelNames = modelNames,
    reduce = T
  )

  cds_3 <- TSCAN::TSCANorder(cds_2)

  milestone_network <- cds_2$MSTtree %>% igraph::as_data_frame() %>% select(-weight) %>% mutate(directed=TRUE)
  end_milestones <- setdiff(milestone_network$to, milestone_network$from)
  milestone_network <- bind_rows(
    milestone_network,
    tibble(
      from=end_milestones,
      to=as.character(seq(max(as.numeric(milestone_network$to))+1, length.out=length(end_milestones)))
    )
  )

  model <- tibble(cluster_id = cds_2$clusterid, cell_id = names(cds_2$clusterid)) %>%
    mutate(rank = match(cell_id, cds_3), cluster_id=as.character(cluster_id)) %>%
    group_by(cluster_id) %>%
    mutate(percentage=order(rank)/n()) %>%
    ungroup()
  progressions <- model %>%
    left_join(milestone_network, by=c("cluster_id"="from")) %>%
    group_by(cell_id) %>%
    mutate(percentage=percentage/n()) %>%
    rename(from=cluster_id)

  milestone_network <- progressions %>% group_by(from, to) %>% summarise(length=n(), directed=TRUE) %>% ungroup()

  wrap_ti_prediction(
    ti_type = "branching",
    id = "TSCAN",
    cell_ids = rownames(counts),
    milestone_ids = unique(c(milestone_network$from, milestone_network$to)),
    milestone_network = milestone_network,
    progressions = progressions %>% select(cell_id, from, to, percentage),
    dimred_samples = cds_2$pcareduceres,
    dimred_clust = cds_2$clusterid,
    clust_centers = cds_2$clucenter,
    cds = cds_3
  )
}

plot_tscan <- function(prediction) {
  plotdata_cells <- prediction$dimred_samples %>% as.data.frame() %>% tibble::rownames_to_column("cell_id") %>% mutate(cluster_id=factor(prediction$dimred_clust)) %>% mutate(rank=match(cell_id, prediction$cds))
  cluster_order <- plotdata_cells %>% group_by(cluster_id) %>% summarise(rank=mean(rank)) %>% arrange(rank) %>% pull(cluster_id)
  plotdata_clusters <- prediction$clust_centers[cluster_order, ] %>% as.data.frame()

  ggplot() +
    geom_point(aes(PC1, PC2, color=cluster_id), data=plotdata_cells) +
    geom_path(aes(V1, V2), data=plotdata_clusters)
}

