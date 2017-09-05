#' Description for Mpath
#' @export
description_mpath <- function() {
  list(
    name = "Mpath",
    short_name = "Mpath",
    package_load = c("Mpath"),
    package_installed = c(),
    par_set = makeParamSet(
    ),
    properties = c(),
    run_fun = run_mpath,
    plot_fun = plot_mpath
  )
}
run_mpath <- function(counts, cell_grouping, numcluster=15, method="diversity") {
  oldwd <- getwd() # Mpath randomly generates plots and output files in the working directory, so change the wd here
  setwd(tempdir())
  fakeFile <- function(x) {
    loc <- tempfile()
    write.table(x, file=loc, sep="\t")
    loc
  }

  sampleInfo <- cell_grouping %>% rename(GroupID=group_id)

  #find_optimal_cluster_number(fakeFile(t(counts)), fakeFile(sampleInfo))

  # never use method = kmeans, dark monsters from the abiss reside there!
  landmark_cluster <- landmark_designation(fakeFile(t(counts)), "wat", fakeFile(sampleInfo), numcluster=numcluster, method = method, saveRes=FALSE)

  network <- build_network(t(counts), landmark_cluster)
  trimmed_network <- trim_net(network)

  ordering <- nbor_order(t(counts), landmark_cluster, unique(landmark_cluster$landmark_cluster))

  trimmed_network[upper.tri(trimmed_network)] = 0
  attr(trimmed_network, "class") = "matrix"
  milestone_network <- trimmed_network %>% as.matrix() %>% reshape2::melt(varnames=c("from", "to")) %>%
    filter(value == 1) %>% select(-value) %>% mutate(from=as.character(from), to=as.character(to))

  # milestone_network %>% igraph::graph_from_data_frame() %>% plot
  #
  # milestone_network %>% {unique(.$from, .$to) %in% .$to}

  # add milestones for landmarks with only outgoing edges
  beginning_milestones <- unique(c(milestone_network$from, milestone_network$to)) %>% keep(~!(. %in% milestone_network$from))
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
    ungroup()

  setwd(oldwd)

  wrap_ti_prediction(
    ti_type = "tree",
    id = "Mpath",
    cell_ids = rownames(counts),
    milestone_ids = unique(c(milestone_network$from, milestone_network$to)),
    milestone_network = milestone_network %>% select(from, to, length),
    progressions = progressions %>% select(cell_id, from, to, percentage) %>% mutate(cell_id=as.character(cell_id))
  )
}


plot_mpath <- function(prediction) {

}
