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
#' @importFrom stats na.omit
run_mpath <- function(counts,
                      cell_grouping,
                      distMethod = "euclidean",
                      method = "kmeans",
                      numcluster = 11,
                      diversity_cut = .6,
                      size_cut = .05) {
  requireNamespace("igraph")

  sample_info <- cell_grouping %>% rename(GroupID = group_id) %>% as.data.frame

  landmark_cluster <- Mpath::landmark_designation(
    rpkmFile = t(counts),
    baseName = NULL,
    sampleFile = sample_info,
    distMethod = distMethod,
    method = method,
    numcluster = numcluster,
    diversity_cut = diversity_cut,
    size_cut = size_cut,
    saveRes = FALSE
  )

  milestone_ids <- unique(landmark_cluster$landmark_cluster)

  # catch situation where mpath only detects 1 landmark
  if (length(milestone_ids) == 1) {
    cell_ids <- rownames(counts)
    milestone_network <- data_frame(
      from = milestone_ids,
      to = milestone_ids,
      length = 1,
      directed = FALSE
    )
    progressions <- data_frame(
      cell_id = cell_ids,
      from = milestone_ids,
      to = milestone_ids,
      percentage = 1
    )
  } else {
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
      stats::na.omit %>%
      group_by(cell_id) %>%
      arrange(desc(percentage)) %>%
      slice(1) %>%
      ungroup()
  }

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

plot_mpath <- function(prediction) {
  requireNamespace("igraph")

  # milestone net as igraph
  gr <- igraph::graph_from_data_frame(
    prediction$milestone_network %>% filter(to != "FILTERED_CELLS"),
    directed = FALSE,
    vertices = prediction$milestone_ids
  )

  # collect info on cells
  cell_ids <- prediction$cell_ids
  labels <- prediction$cell_grouping %>% slice(match(cell_ids, cell_id)) %>% .$group_id
  clustering <- prediction$progressions %>%
    slice(match(cell_ids, cell_id)) %>%
    {with(., ifelse(percentage == 0, from, to))}

  # make plot
  make_piegraph_plot(gr, clustering, labels)
}

#' @importFrom cowplot theme_nothing
make_pie_plot <- function(x_count, annotation_colours, add_label = FALSE) {
  count_df <- data_frame(
    celltype = factor(names(x_count), levels = names(x_count)),
    count = x_count,
    percent = x_count / sum(x_count),
    right = cumsum(percent),
    left = right - percent
  )
  g <- ggplot() +
    geom_rect(aes(xmin = 0, xmax = 1, ymin = 0, ymax = 1), fill = "white") +
    geom_rect(aes(xmin = left, xmax = right, ymin = 0, ymax = 1, fill = celltype), count_df) +
    geom_hline(yintercept = c(0, 1)) +
    scale_fill_manual(values = annotation_colours) +
    xlim(0, 1) +
    cowplot::theme_nothing() +
    coord_polar()
  if (add_label) {
    g +
      geom_text(aes((left + right)/2, 1.5, label = celltype, colour = celltype), count_df) +
      ylim(0, 2) +
      scale_colour_manual(values = annotation_colours)
  } else {
    g + ylim(0, 1)
  }
}

#' @importFrom cowplot theme_nothing
make_legend_plot <- function(annotation_colours) {
  count <- setNames(rep(1, length(annotation_colours)), names(annotation_colours))
  make_pie_plot(
    count,
    annotation_colours,
    add_label = TRUE
  ) + labs(title = "Legend")
}

#' @importFrom RColorBrewer brewer.pal
#' @importFrom cowplot theme_nothing plot_grid
#' @importFrom ggimage geom_subview
make_piegraph_plot <- function(gr, clustering, labels) {
  requireNamespace("igraph")

  # plot de graph
  lay <- gr %>%
    igraph::layout_with_kk() %>%
    dynutils::scale_quantile(0)
  dimnames(lay) <- list(
    igraph::V(gr)$name,
    c("X", "Y")
  )
  lay_df <- lay %>% as.data.frame %>% rownames_to_column()

  mst_df <- igraph::as_data_frame(gr)
  colnames(mst_df) <- c("i", "j", "dist")

  # kijken hoeveel van welke category er in iedere node zit
  categories <- if (is.factor(labels)) levels(labels) else sort(unique(labels))
  clusters <- names(igraph::V(gr))
  node_counts <-
    sapply(categories, function(ca) {
      sapply(clusters, function(cl) {
        sum(labels == ca & clustering == cl)
      })
    })
  num_labels <- rowSums(node_counts)

  # stel een kleurenschema op
  annotation_colours <- setNames(RColorBrewer::brewer.pal(ncol(node_counts), "Set2"), colnames(node_counts))

  # attach positions of mst_df
  mst_df_with_pos <- data.frame(
    mst_df,
    i = lay[mst_df[,1],,drop=F],
    j = lay[mst_df[,2],,drop=F]
  )

  # maak een plot van de lijntjes
  max_size <- .075
  p <- ggplot(data.frame(lay), aes(X, Y)) +
    geom_segment(aes(x = i.X, xend = j.X, y = i.Y, yend = j.Y), data.frame(mst_df_with_pos)) +
    cowplot::theme_nothing() +
    coord_equal()

  # maak subplots voor elk van de nodes
  subplots <- lapply(seq_len(nrow(node_counts)), function(i) {
    x <- lay[i,"X"]
    y <- lay[i,"Y"]
    x_count <- node_counts[i,]
    subplot <- make_pie_plot(x_count, annotation_colours)
    area <- num_labels[[i]] / max(num_labels) * pi/4*max_size^2
    width <- sqrt(area) * 4 / pi
    ggimage::geom_subview(subplot, x, y, width, width)
  })

  # maak een plot van de legende
  legends <- make_legend_plot(annotation_colours)

  # combineer alle plots
  combined_plot <- Reduce("+", c(list(p), subplots)) + coord_equal()

  combined_plot +
    ggimage::geom_subview(legends, .15, .15, .2, .2) +
    geom_text(aes(X, Y, label = rowname), lay_df)
}


