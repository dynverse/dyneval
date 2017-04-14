library(cowplot)
library(tidyverse)
library(dyneval)
library(monocle)

make_node_graph <- function(linegraph, edge_names) {
  # add 2 nodes to each edge, call them <edgename>_start and <edgename>_end
  edges_str <- sort(c(paste0(edge_names, "_start"), paste0(edge_names, "_end")))

  # load the linegraph in igraph, and determine which nodes are equivalent to eachother
  gr <- linegraph %>%
    transmute(from_str = paste0(from, "_end"), to_str = paste0(to, "_start")) %>%
    igraph::graph_from_data_frame(vertices = edges_str)

  edge_to_node <- igraph::components(gr)$membership

  # give each node a unique name
  edge_to_node <- setNames(paste0("node_", edge_to_node), names(edge_to_node))

  # translate each edge to its 'from -> to' notation
  node_graph <- bind_rows(lapply(edge_names, function(i) {
    data_frame(from = edge_to_node[paste0(i, "_start")], to = edge_to_node[paste0(i, "_end")], edge_id = i)
  }))

  list(
    node_graph = node_graph,
    node_names = sort(unique(edge_to_node))
  )
}

dyngen_dataset_to_task <- function(dataset, name) {
  linegraph <- dataset$gs$piecestatenet
  #edge_names <- sort(unique(with(linegraph, c(from, to))))
  edge_names <- seq_along(dataset$gs$piecestates)

  node_graph <- make_node_graph(linegraph, edge_names)

  state_names <- node_graph$node_names

  state_network <- dataset$gs$cellinfo %>%
    group_by(piecestateid) %>%
    summarise(length = max(progression) - min(progression)) %>%
    mutate(length = length / min(length)) %>%
    mutate(piecestateid = as.integer(piecestateid)) %>%
    left_join(node_graph$node_graph, by = c("piecestateid"="edge_id")) %>%
    select(from, to, length)

  # expression <- log2(dataset$counts + 1)
  expression <- dataset$expression

  state_percentages <- dataset$gs$cellinfo %>%
    mutate(piecestateid = as.integer(piecestateid)) %>%
    left_join(node_graph$node_graph, by = c("piecestateid"="edge_id")) %>%
    as_data_frame() %>%
    group_by(piecestateid) %>%
    mutate(progression_sc = (progression - min(progression)) / (max(progression) - min(progression))) %>%
    ungroup() %>%
    select(cell, from, to, progression_sc) %>%
    gather(direct, node, from, to) %>%
    mutate(value = ifelse(direct == "from", 1 - progression_sc, progression_sc)) %>%
    select(-progression_sc, -direct) %>%
    spread(node, value, fill = 0) %>%
    rename(id = cell) %>%
    slice(match(rownames(expression), id))

  wrap_ti_task_data(
    ti_type = "manyfurcating",
    name = name,
    expression = expression,
    state_names = state_names,
    state_network = state_network,
    state_percentages = state_percentages
  )
}

datasets <- readRDS("../dyngen/results/datasets.rds")

# dataset_num <- 10
for (dataset_num in seq_along(datasets)) {
  task <- dyngen_dataset_to_task(datasets[[dataset_num]], paste0("dataset_", dataset_num))
  cat("Processing ", task$name, "\n", sep="")

  ## compare with scorpius
  cat("  Running SCORPIUS\n")
  pred_output1 <- trainLearner.ti.scorpius(task, .subset = NULL, num_dimensions = 3, num_clusters = 4)
  cat("  Running Monocle\n")
  pred_output2 <- trainLearner.ti.monocle(task, .subset = NULL)

  ## calculate EM distances
  cat("  Computing EM distances: gold\n")
  task_emdist <- compute_em_dist(task)
  cat("  Computing EM distances: SCORPIUS\n")
  pred_emdist1 <- compute_em_dist(pred_output1)
  cat("  Computing EM distances: Monocle\n")
  pred_emdist2 <- compute_em_dist(pred_output2)

  cat("  Computing coranking: SCORPIUS\n")
  corank1 <- compute_coranking(task_emdist, pred_emdist1)
  cat("  Computing coranking: Monocle\n")
  corank2 <- compute_coranking(task_emdist, pred_emdist2)

  cat("  Generating plots\n")
  scores <- data_frame(
    method = c("SCORPIUS", "monocle"),
    auc_lcmc = c(corank1$auc_lcmc, corank2$auc_lcmc),
    cor = cor(task_emdist %>% as.vector, cbind(pred_emdist1 %>% as.vector, pred_emdist2 %>% as.vector))[1,])

  # row 1
  plsz <- 500
  png("scratch/tmp/plot1.png", plsz, plsz)
  plotLearner.ti.default(task) + labs(title = "Gold")
  dev.off()
  png("scratch/tmp/plot2.png", plsz, plsz)
  plotLearner.ti.default(pred_output1) + labs(title = "SCORPIUS prediction")
  dev.off()
  png("scratch/tmp/plot3.png", plsz, plsz)
  plotLearner.ti.default(pred_output2) + labs(title = "monocle prediction")
  dev.off()


  plotdata_task <- plotLearnerData.ti.default(task)
  plotdata_output1 <- plotLearnerData.ti.default(pred_output1)
  plotdata_output2 <- plotLearnerData.ti.default(pred_output2)

  plotdata_overlay_output1 <- plotdata_output1
  plotdata_overlay_output1$space_samples$colour <- plotdata_task$space_samples$colour
  plotdata_overlay_output2 <- plotdata_output2
  plotdata_overlay_output2$space_samples$colour <- plotdata_task$space_samples$colour

  pl2 <- function(object, insert_phantom_edges = T) {
    if (dyneval:::is_ti_wrapper(object)) {
      dimred_object <- dyneval:::plotLearnerData.ti.default(object, insert_phantom_edges = insert_phantom_edges)
    } else if (dyneval:::is_ti_dimred_wrapper(object)) {
      dimred_object <- object
    } else {
      stop(sQuote("object"), " is not an object as generated by plotLearnerData.ti.default, wrap_ti_task_data or wrap_ti_prediction.")
    }
    with(dimred_object, {
      ggplot() +
        geom_segment(aes(x = from.Comp1, xend = to.Comp1, y = from.Comp2, yend = to.Comp2), space_lines,
                     size = 10, colour = "#444444", arrow = arrow(length = unit(.5, "cm"), type="closed")) +
        geom_point(aes(Comp1, Comp2, colour = colour), space_samples, size = 3, alpha = .5) +
        geom_text(aes(Comp1, Comp2, label = id), space_states, nudge_y = .05) +
        scale_colour_identity() +
        scale_fill_identity() +
        scale_x_continuous(limits = c(-.55, .55)) +
        scale_y_continuous(limits = c(-.55, .55)) +
        coord_equal() +
        labs(x = "Component 1", y = "Component 2")
    })
  }

  # row 1b
  png("scratch/tmp/plot1b.png", plsz, plsz)
  plot.new()
  dev.off()
  png("scratch/tmp/plot2b.png", plsz, plsz)
  pl2(plotdata_overlay_output1) + labs(title = "SCORPIUS prediction")
  dev.off()
  png("scratch/tmp/plot3b.png", plsz, plsz)
  pl2(plotdata_overlay_output2) + labs(title = "monocle prediction")
  dev.off()


  # row 2
  png("scratch/tmp/plot4.png", plsz, plsz)
  plot.new()
  dev.off()
  png("scratch/tmp/plot5.png", plsz, plsz)
  plotLearner.ti.scorpius(pred_output1) + labs(title = "SCORPIUS own plot") + theme(legend.position = "none")
  dev.off()
  png("scratch/tmp/plot6.png", plsz, plsz)
  plotLearner.ti.monocle(pred_output2) + coord_equal() + labs(title = "monocle own plot")
  dev.off()

  # row 3
  png("scratch/tmp/plot7.png", plsz, plsz)
  plot_emdist(task, task_emdist)
  dev.off()
  png("scratch/tmp/plot8.png", plsz, plsz)
  plot_emdist(pred_output1, pred_emdist1)
  dev.off()
  png("scratch/tmp/plot9.png", plsz, plsz)
  plot_emdist(pred_output2, pred_emdist2)
  dev.off()

  # row 4
  png("scratch/tmp/plot10.png", plsz, plsz)
  ggplot(scores) + geom_bar(aes(method, auc_lcmc), stat = "identity")
  dev.off()
  png("scratch/tmp/plot11.png", plsz, plsz)
  coRanking::imageplot(corank1$corank)
  dev.off()
  png("scratch/tmp/plot12.png", plsz, plsz)
  coRanking::imageplot(corank2$corank)
  dev.off()

  system(paste0(
    "cd scratch/tmp\n",
    "bash << HERE\n",
    "convert plot[123].png +append result_a.png\n",
    "convert plot[123]b.png +append result_b.png\n",
    "convert plot[456].png +append result_c.png\n",
    "convert plot[789].png +append result_d.png\n",
    "convert plot1[012].png +append result_e.png\n",
    "convert result_[abcde].png -append result.png\n",
    "HERE\n"))

  system(paste0(
    "cd scratch/tmp\n",
    "bash << HERE\n",
    "convert plot[123].png -append result_a.png\n",
    "convert plot[123]a.png -append result_b.png\n",
    "convert plot[456].png -append result_c.png\n",
    "convert plot[789].png -append result_d.png\n",
    "convert plot1[012].png -append result_e.png\n",
    "convert result_[abcde].png +append result2.png\n",
    "HERE\n"))

  system("rm scratch/tmp/*")

  file.copy("scratch/tmp/result.png", paste0("scratch/output/", task$name, "_plot.png"))
  file.copy("scratch/tmp/result2.png", paste0("scratch/output/", task$name, "_plot2.png"))
  save(task, pred_output1, pred_output2, task_emdist, pred_emdist1, pred_emdist2, corank1, corank2, scores, file = paste0("scratch/output/", task$name, "_data.RData"))
}
