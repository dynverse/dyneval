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

#' @importFrom utils write.table
run_mpath <- function(counts,
                      cell_grouping,
                      distMethod = "euclidean",
                      method = "kmeans",
                      numcluster = 11,
                      diversity_cut = .6,
                      size_cut = .05) {
  requireNamespace("igraph")

  # write data to files
  tmp <- tempfile(pattern = "mpath")
  counts_file <- paste0(tmp, "_counts.tsv")
  sampleinfo_file <- paste0(tmp, "_sampleinfo.tsv")

  utils::write.table(t(counts), counts_file, sep = "\t")
  utils::write.table(cell_grouping %>% rename(GroupID = group_id), sampleinfo_file, sep = "\t")

  # designate landmarks
  tryCatch({
    landmark_cluster <- Mpath::landmark_designation(
      rpkmFile = counts_file,
      baseName = NULL,
      sampleFile = sampleinfo_file,
      distMethod = distMethod,
      method = method,
      numcluster = numcluster,
      diversity_cut = diversity_cut,
      size_cut = size_cut,
      saveRes = FALSE
    )
  },
  finally = {
    file.remove(counts_file)
    file.remove(sampleinfo_file)
  })

  milestone_ids <- unique(landmark_cluster$landmark_cluster)

  # build network
  network <- Mpath::build_network(
    exprs = t(counts),
    baseName = NULL,
    landmark_cluster = landmark_cluster,
    distMethod = distMethod,
    writeRes = FALSE
  )

  # trim network
  trimmed_network <- Mpath::trim_net(
    nb12 = network,
    writeRes = FALSE
  )

  # create final milestone network
  class(trimmed_network) <- NULL
  milestone_network <- trimmed_network %>%
    reshape2::melt(varnames = c("from", "to"), value.name = "length") %>%
    mutate_if(is.factor, as.character) %>%
    filter(length > 0, from < to) %>%
    mutate(directed = FALSE)

  # find an edge for each cell to sit on
  connections <- bind_rows(
    milestone_network %>% mutate(landmark_cluster = from, percentage = 0),
    milestone_network %>% mutate(landmark_cluster = to, percentage = 1)
  )
  progressions <-
    landmark_cluster %>%
    left_join(connections, by = "landmark_cluster") %>%
    select(cell_id = cell, from, to, percentage) %>%
    na.omit %>%
    group_by(cell_id) %>%
    arrange(desc(percentage)) %>%
    slice(1) %>%
    ungroup()

  # The nbor ordering can't really be used in this
  # gr <- milestone_network %>%
  #   igraph::graph_from_data_frame(directed = FALSE, vertices = milestone_ids)
  #
  # # randomly select start milestone
  # start_node <- igraph::degree(gr) %>% keep(~. == 1) %>% names() %>% sample(1)
  #
  # # create ordering of landmarks
  # lm_order <- igraph::distances(gr, v = start_node)[1,] %>% sort %>% names
  #
  # # create cell ordering
  # ordering <- Mpath::nbor_order(
  #   t(counts),
  #   landmark_cluster,
  #   lm_order = lm_order,
  #   writeRes = FALSE
  # )

  # return output
  wrap_ti_prediction(
    ti_type = "tree",
    id = "Mpath",
    cell_ids = rownames(counts),
    milestone_ids = milestone_ids,
    milestone_network = milestone_network,
    progressions = progressions,
    landmark_cluster = landmark_cluster,
    cell_grouping = cell_grouping
  )
}

# TODO: migrate to ggplot!
plot_mpath <- function(prediction) {
  requireNamespace("igraph")

  # milestone net as igraph
  g <- prediction$milestone_network %>%
    filter(to != "FILTERED_CELLS") %>%
    igraph::graph_from_data_frame(directed = F)

  # calculate sizes of pie chunks
  pie_sizes <- prediction$landmark_cluster %>%
    left_join(prediction$cell_grouping, by = c("cell"="cell_id")) %>%
    mutate(group_id = factor(group_id)) %>%
    group_by(landmark_cluster) %>%
    summarise(counts = list(as.numeric(table(group_id)))) %>%
    {set_names(.$counts, .$landmark_cluster)}

  pie_sizes <- pie_sizes[names(igraph::V(g))]

  # plot graph
  igraph::plot.igraph(g, vertex.shape = "pie",vertex.pie = pie_sizes)
}
