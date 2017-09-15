#' Description for Mpath
#' @export
description_mpath <- function() create_description(
  name = "Mpath",
  short_name = "Mpath",
  package_loaded = c("Mpath"),
  package_required = c(),
  par_set = makeParamSet(
    makeDiscreteParam(id = "distMethod", default = "euclidean", values = c("pearson", "kendall", "spearman", "euclidean")),
    makeDiscreteParam(id = "method", default = "diversity_size", values = c("kmeans", "diversity", "size", "diversity_size")),
    makeIntegerParam(id = "numcluster", lower = 3L, default = 11L, upper = 30L),
    makeNumericParam(id = "diversity_cut", lower = .1, default = .6, upper = 1),
    makeNumericParam(id = "size_cut", lower = .01, default = .05, upper = 1)
  ),
  properties = c(),
  run_fun = run_mpath,
  plot_fun = plot_mpath
)

run_mpath <- function(counts, cell_grouping,
                      distMethod = "euclidean", method = "kmeans",
                      numcluster = 11, diversity_cut = .6, size_cut = .05) {
  # function to save a data.frame in a temporary directory and return the file's location
  fakeFile <- function(x) {
    loc <- tempfile()
    write.table(x, file=loc, sep="\t")
    loc
  }

  sampleInfo <- cell_grouping %>% rename(GroupID=group_id)

  landmark_cluster <- landmark_designation(
    fakeFile(t(counts)),
    "output_disabled",
    fakeFile(sampleInfo),
    distMethod = distMethod,
    method = method,
    numcluster = numcluster,
    diversity_cut = diversity_cut,
    size_cut = size_cut,
    saveRes = FALSE)

  network <- build_network(t(counts), landmark_cluster, distMethod = distMethod)
  trimmed_network <- trim_net(network)

  ordering <- nbor_order(t(counts), landmark_cluster, unique(landmark_cluster$landmark_cluster))

  trimmed_network[upper.tri(trimmed_network)] = 0
  attr(trimmed_network, "class") = "matrix"

  mpath_network <- trimmed_network %>%
    as.matrix() %>%
    reshape2::melt(varnames = c("from", "to")) %>%
    filter(value == 1) %>%
    select(-value) %>%
    mutate(from = as.character(from), to = as.character(to))

  milestone_network <- mpath_network

  # add milestones for landmarks with only outgoing edges
  beginning_milestones <- unique(c(milestone_network$from, milestone_network$to)) %>%
    keep(~!(. %in% milestone_network$from))

  milestone_network <- bind_rows(
    milestone_network,
    tibble(
      from = beginning_milestones,
      to = paste0("extra_", seq_along(beginning_milestones))
    )
  )

  # add ordering and from and to from milestone_network
  progressions1 <- landmark_cluster %>%
    rename(cell_id=cell) %>%
    mutate(global_rank=match(cell_id, ordering)) %>%
    filter(!is.na(global_rank)) %>%
    rename(from=landmark_cluster) %>%
    left_join(milestone_network, by=c("from"))

  # calculate time based on ordering
  progressions <- progressions1 %>%
    group_by(from, to) %>%
    mutate(percentage=(global_rank - min(global_rank))/(max(global_rank) - min(global_rank))) %>%
    ungroup() %>%
    group_by(cell_id) %>%
    mutate(percentage=percentage/n()) %>%
    ungroup()

  milestone_network <- progressions %>% group_by(from, to) %>% summarise(length=n()) %>%
    right_join(milestone_network, by=c("from", "to")) %>%
    ungroup() %>%
    mutate(directed=FALSE)

  wrap_ti_prediction(
    ti_type = "tree",
    id = "Mpath",
    cell_ids = rownames(counts),
    milestone_ids = unique(c(milestone_network$from, milestone_network$to)),
    milestone_network = milestone_network %>% select(from, to, length, directed),
    progressions = progressions %>% select(cell_id, from, to, percentage) %>% mutate(cell_id=as.character(cell_id)),
    landmark_cluster = landmark_cluster,
    mpath_network = mpath_network,
    cell_grouping = cell_grouping
  )
}

plot_mpath <- function(prediction) {
  requireNamespace("igraph")

  pie_sizes <- prediction$landmark_cluster %>% left_join(prediction$cell_grouping, by=c("cell"="cell_id")) %>%
    mutate(group_id=factor(group_id)) %>%
    group_by(landmark_cluster) %>%
    summarise(counts=list(as.numeric(table(group_id)))) %>%
    {set_names(.$counts, .$landmark_cluster)}

  g <- prediction$mpath_network %>% igraph::graph_from_data_frame()
  pie_sizes <- pie_sizes[names(igraph::V(g))]
  g %>% igraph::plot.igraph(vertex.shape="pie",vertex.pie=pie_sizes)
}
