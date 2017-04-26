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

  expression <- t(SCORPIUS::quant.scale(t(log2(dataset$counts + 1))))
  # expression <- dataset$expression

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
output_root_folder <- "scratch/output_dyngentest/"

# dataset_num <- 10
for (dataset_num in seq_along(datasets)) {
  task_name <- paste0("dataset_", dataset_num)
  cat("Processing ", task_name, "\n", sep="")
  dataset <- datasets[[dataset_num]]
  data_dir <- paste0(output_root_folder, task_name)
  dir.create(data_dir)
  data_file <- paste0(data_dir, "/data.RData")

  if (file.exists(data_file)) {
    cat("  Skipping; dataset already executed\n")
  } else if (any(dataset$gs$piecestatenet$from == dataset$gs$piecestatenet$to)) {
    cat("  Skipping; self-edges are not allowed\n")
  } else {
    task <- dyngen_dataset_to_task(datasets[[dataset_num]], task_name)

    ## compare with scorpius
    cat("  Running SCORPIUS\n")
    pred_output1 <- trainLearner.ti.scorpius(task, .subset = NULL, num_dimensions = 3, num_clusters = 4)
    cat("  Running Monocle\n")
    pred_output2 <- trainLearner.ti.monocle(task, .subset = NULL, num_dimensions = 3)

    ## calculate EM distances
    cat("  Computing EM distances: gold\n")
    task_emdist <- compute_emlike_dist(task)
    cat("  Computing EM distances: SCORPIUS\n")
    pred_emdist1 <- compute_emlike_dist(pred_output1)
    cat("  Computing EM distances: Monocle\n")
    pred_emdist2 <- compute_emlike_dist(pred_output2)

    cat("  Computing coranking: SCORPIUS\n")
    corank1 <- compute_coranking(task_emdist, pred_emdist1)
    cat("  Computing coranking: Monocle\n")
    corank2 <- compute_coranking(task_emdist, pred_emdist2)

    cat("  Saving data\n")
    scores <- data.frame(
      method = c("SCORPIUS", "monocle"),
      bind_rows(corank1$summary, corank2$summary),
      cor = cor(task_emdist %>% as.vector, cbind(pred_emdist1 %>% as.vector, pred_emdist2 %>% as.vector))[1,])

    save(task, pred_output1, pred_output2, task_emdist, pred_emdist1, pred_emdist2, corank1, corank2, scores, file = data_file)
  }
}

for (dataset_num in seq_along(datasets)) {
  task_name <- paste0("dataset_", dataset_num)
  cat("Processing ", task_name, "\n", sep="")
  data_dir <- paste0(output_root_folder, task_name)
  data_file <- paste0(data_dir, "/data.RData")
  if (file.exists(data_file)) {
    load(data_file)

    plotdata_task <- plotLearnerData.ti.default(task)
    plotdata_output1 <- plotLearnerData.ti.default(pred_output1)
    plotdata_output2 <- plotLearnerData.ti.default(pred_output2)

    # row 1
    plsz <- 1200
    png(paste0(data_dir, "/plot1.png"), plsz, plsz, res = 300)
    print(plotLearner.ti.default(plotdata_task) + labs(title = "Gold"))
    dev.off()
    png(paste0(data_dir, "/plot2.png"), plsz, plsz, res = 300)
    print(plotLearner.ti.default(plotdata_output1) + labs(title = "SCORPIUS prediction"))
    dev.off()
    png(paste0(data_dir, "/plot3.png"), plsz, plsz, res = 300)
    print(plotLearner.ti.default(plotdata_output2) + labs(title = "monocle prediction"))
    dev.off()

    # row 2
    png(paste0(data_dir, "/plot4.png"), plsz, plsz, res = 300)
    plot.new()
    dev.off()
    png(paste0(data_dir, "/plot5.png"), plsz, plsz, res = 300)
    print(plotLearner.ti.combined(plotdata_task, plotdata_output1) + labs(title = "SCORPIUS prediction"))
    dev.off()
    png(paste0(data_dir, "/plot6.png"), plsz, plsz, res = 300)
    print(plotLearner.ti.combined(plotdata_task, plotdata_output2) + labs(title = "monocle prediction"))
    dev.off()

    # row 3
    png(paste0(data_dir, "/plot7.png"), plsz, plsz, res = 300)
    plot.new()
    dev.off()
    png(paste0(data_dir, "/plot8.png"), plsz, plsz, res = 300)
    print(plotLearner.ti.scorpius(pred_output1) + labs(title = "SCORPIUS own plot") + theme(legend.position = "none"))
    dev.off()
    png(paste0(data_dir, "/plot9.png"), plsz, plsz, res = 300)
    print(plotLearner.ti.monocle(pred_output2) + coord_equal() + labs(title = "monocle own plot") + theme(legend.position = "none"))
    dev.off()

    # row 4
    plot_emdist(task, task_emdist, plotdata_task, filename = paste0(data_dir, "/plot10.png"), width = plsz / 300, height = plsz / 300)
    plot_emdist(pred_output1, pred_emdist1, plotdata_output1, filename = paste0(data_dir, "/plot11.png"), width = plsz / 300, height = plsz / 300)
    plot_emdist(pred_output2, pred_emdist2, plotdata_output2, filename = paste0(data_dir, "/plot12.png"), width = plsz / 300, height = plsz / 300)

    # row 5
    png(paste0(data_dir, "/plot13.png"), plsz, plsz, res = 300)
    print(ggplot(scores) + geom_bar(aes(method, max_lcmc), stat = "identity"))
    dev.off()
    png(paste0(data_dir, "/plot14.png"), plsz, plsz, res = 300)
    coRanking::imageplot(corank1$corank)
    dev.off()
    png(paste0(data_dir, "/plot15.png"), plsz, plsz, res = 300)
    coRanking::imageplot(corank2$corank)
    dev.off()

    # row 6
    png(paste0(data_dir, "/plot16.png"), plsz, plsz, res = 300)
    print(ggplot(scores) + geom_bar(aes(method, cor), stat = "identity"))
    dev.off()
    png(paste0(data_dir, "/plot17.png"), plsz, plsz, res = 300)
    (ggplot() + geom_point(aes(as.vector(task_emdist), as.vector(pred_emdist1)), alpha = .1, size = .5) + labs(x = "original EM dist", y = "SCORPIUS EM dist")) %>% print
    dev.off()
    png(paste0(data_dir, "/plot18.png"), plsz, plsz, res = 300)
    (ggplot() + geom_point(aes(as.vector(task_emdist), as.vector(pred_emdist2)), alpha = .1, size = .5) + labs(x = "original EM dist", y = "monocle EM dist")) %>% print
    dev.off()

    system(paste0(
      "cd ", data_dir, "\n",
      "bash << HERE\n",
      "convert plot[123].png -append result_a.png\n",
      "convert plot[456].png -append result_b.png\n",
      "convert plot[789].png -append result_c.png\n",
      "convert plot1[012].png -append result_d.png\n",
      "convert plot1[345].png -append result_e.png\n",
      "convert plot1[678].png -append result_f.png\n",
      "convert result_[abcdef].png +append result.png\n",
      "HERE\n"))

    file.copy(paste0(data_dir, "/result.png"), paste0(data_dir, "_plot.png"), overwrite = T)
  }
}


dat_evals <- lapply(seq_along(datasets), function(dataset_num) {
  task_name <- paste0("dataset_", dataset_num)
  data_dir <- paste0(output_root_folder, task_name)
  data_file <- paste0(data_dir, "/data.RData")
  if (!file.exists(data_file)) {
    NULL
  } else {
    vals <- load(data_file)
    environment() %>% as.list() %>% .[vals]
  }
})

scores <- dat_evals %>% map_df(~ data.frame(dataset = .$task$name, ti_type = .$task$ti_type, .$scores))
ggplot(scores) +
  geom_path(aes(max_lcmc, cor, group = dataset)) +
  geom_point(aes(max_lcmc, cor, colour = method), size = 3)
#
# de <- dat_evals[[1]]
# qplot(as.vector(de$task_emdist), as.vector(de$pred_emdist1))
# qplot(as.vector(de$task_emdist), as.vector(de$pred_emdist2))
#
#
# dataset <- datasets[[10]]
#
# # expression_c <- t(SCORPIUS::quant.scale(t(log2(dataset$counts + 1))))
# # expression_r <- dataset$expression[rownames(expression_c),]
# #
# # pheatmap::pheatmap(t(SCORPIUS::quant.scale(expression_c)), cluster_cols = F)
# # pheatmap::pheatmap(t(SCORPIUS::quant.scale(expression_r)), cluster_cols = F)
#
# task <- dyngen_dataset_to_task(dataset, "freiu")
# plotdata <- plotLearnerData.ti.default(task)
# plotLearner.ti.default(plotdata)
# space_states <- plotdata$space_states
#
# expr <- task$expression
# ddr <- DDRTree::DDRTree(SCORPIUS::quant.scale(t(expr)), dimensions = 20, sigma = .001, maxIter = 20, lambda = NULL, param.gamma = 10, tol = .001)
# sample_coord <- data.frame(t(ddr$Z))
# traj_coord <- data.frame(t(ddr$Y))
# colnames(sample_coord) <- colnames(traj_coord) <- paste0("Comp", seq_len(ncol(sample_coord)))
# sample_coord$type <- task$state_names[apply(task$state_percentages[,-1], 1, which.max)]
# tree_links <- ddr$stree %>% as.matrix %>% reshape2::melt(varnames=c("from", "to"), value.name = "value") %>% filter(value != 0)
# tree_links_coord <- tree_links %>% left_join(data.frame(from = seq_len(nrow(traj_coord)), from=traj_coord)) %>% left_join(data.frame(to = seq_len(nrow(traj_coord)), to=traj_coord))
#
# ggplot() +
#   geom_segment(aes(x = from.Comp1, xend = to.Comp1, y = from.Comp2, yend = to.Comp2), tree_links_coord) +
#   geom_point(aes(Comp1, Comp2, colour = type), sample_coord) +
#   scale_colour_manual(values = setNames(space_states$colour, space_states$id)) +
#   theme(panel.background = element_rect(fill = "#777777"))
# #
