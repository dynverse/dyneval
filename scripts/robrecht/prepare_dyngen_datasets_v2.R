library(cowplot)
library(tidyverse)
library(dyneval)
library(monocle)

datasets_info <- readRDS("../dyngen/results/datasets.rds")

output_root_folder <- "results/output_dyngentest/"
.datasets_location <- "../dyngen/results"

# dataset_num <- 10
for (dataset_num in seq_len(nrow(datasets_info))) {
  dataset_id <- datasets_info$id[[dataset_num]]
  cat("Processing ", dataset_num, "/", nrow(datasets_info), ": ", dataset_id, "\n", sep="")

  data_dir <- paste0(output_root_folder, dataset_id)
  dir.create(data_dir)
  data_file <- paste0(data_dir, "/data.RData")

  if (file.exists(data_file)) {
    cat("  Skipping; dataset already executed\n")
  } else {
    dataset <- dyngen::load_dataset(dataset_id, contents = dyngen::contents_dataset(experiment = F))

    task <- with(dataset, dyneval::wrap_ti_task_data(
      ti_type = model$modulenetname,
      name = info$id,
      expression = log2(counts+1),
      state_names = gs$milestone_names,
      state_net = gs$milestone_net,
      state_percentages = gs$milestone_percentages %>% slice(match(rownames(counts), id))
    ))

    ## compare with scorpius
    cat("  Running SCORPIUS\n")
    pred_output1 <- trainLearner.ti.scorpius(task, .subset = NULL, num_dimensions = 3, num_clusters = 4)
    cat("  Running Monocle\n")
    pred_output2 <- trainLearner.ti.monocle(task, .subset = NULL, num_dimensions = 3)
    # cat("  Running TSCAN\n")
    # pred_output3 <- trainLearner.ti.tscan(task, num_dimensions = 2, preprocess_clusternum = NULL,
    #                                       preprocess_takelog = T, preprocess_logbase = 2,
    #                                       preprocess_pseudocount = 1, preprocess_minexpr_value = 1,
    #                                       preprocess_minexpr_percent = 1, preprocess_cvcutoff = 1,
    #                                       clustering_min_clusternum = 2, clustering_max_clusternum = 9,
    #                                       clustering_reduce = T)

    ## calculate EM distances
    cat("  Computing EM distances: gold\n")
    task_emdist <- compute_emlike_dist(task)
    cat("  Computing EM distances: SCORPIUS\n")
    pred_emdist1 <- compute_emlike_dist(pred_output1)
    cat("  Computing EM distances: Monocle\n")
    pred_emdist2 <- compute_emlike_dist(pred_output2)
    # cat("  Computing EM distances: TSCAN\n")
    # pred_emdist3 <- compute_emlike_dist(pred_output3)

    cat("  Computing coranking: SCORPIUS\n")
    corank1 <- compute_coranking(task_emdist, pred_emdist1)
    cat("  Computing coranking: Monocle\n")
    corank2 <- compute_coranking(task_emdist, pred_emdist2)
    # cat("  Computing coranking: TSCAN\n")
    # corank3 <- compute_coranking(task_emdist, pred_emdist3)

    cat("  Saving data\n")
    # scores <- data.frame(
    #   method = c("SCORPIUS", "monocle", "TSCAN"),
    #   bind_rows(corank1$summary, corank2$summary, corank3$summary),
    #   cor = cor(task_emdist %>% as.vector, cbind(pred_emdist1 %>% as.vector, pred_emdist2 %>% as.vector, pred_emdist3 %>% as.vector))[1,])
    #
    # save(task, pred_output1, pred_output2, pred_output3, task_emdist, pred_emdist1, pred_emdist2, pred_emdist3, corank1, corank2, corank3, scores,
    #      file = data_file)

    scores <- data.frame(
      method = c("SCORPIUS", "monocle"),
      bind_rows(corank1$summary, corank2$summary),
      cor = cor(task_emdist %>% as.vector, cbind(pred_emdist1 %>% as.vector, pred_emdist2 %>% as.vector))[1,])

    save(task, pred_output1, pred_output2,  task_emdist, pred_emdist1, pred_emdist2, corank1, corank2, scores,
         file = data_file)
  }
}

for (dataset_num in seq_len(nrow(datasets_info))) {
  dataset_id <- datasets_info$id[[dataset_num]]
  cat("Processing ", dataset_num, "/", nrow(datasets_info), ": ", dataset_id, "\n", sep="")

  data_dir <- paste0(output_root_folder, dataset_id)
  dir.create(data_dir)
  data_file <- paste0(data_dir, "/data.RData")

  if (file.exists(data_file)) {
    load(data_file)

    plotdata_task <- plotLearnerData.ti.default(task)
    plotdata_output1 <- plotLearnerData.ti.default(pred_output1)
    plotdata_output2 <- plotLearnerData.ti.default(pred_output2)
    plotdata_output3 <- plotLearnerData.ti.default(pred_output3)

    # row 1
    plsz <- 1200
    png(paste0(data_dir, "/plot1a.png"), plsz, plsz, res = 300)
    print(plotLearner.ti.default(plotdata_task) + labs(title = "Gold"))
    dev.off()
    png(paste0(data_dir, "/plot1b.png"), plsz, plsz, res = 300)
    print(plotLearner.ti.default(plotdata_output1) + labs(title = "SCORPIUS prediction"))
    dev.off()
    png(paste0(data_dir, "/plot1c.png"), plsz, plsz, res = 300)
    print(plotLearner.ti.default(plotdata_output2) + labs(title = "monocle prediction"))
    dev.off()
    # png(paste0(data_dir, "/plot1d.png"), plsz, plsz, res = 300)
    # print(plotLearner.ti.default(plotdata_output3) + labs(title = "TSCAN prediction"))
    # dev.off()

    # row 2
    png(paste0(data_dir, "/plot2a.png"), plsz, plsz, res = 300)
    plot.new()
    dev.off()
    png(paste0(data_dir, "/plot2b.png"), plsz, plsz, res = 300)
    print(plotLearner.ti.combined(plotdata_task, plotdata_output1) + labs(title = "SCORPIUS prediction"))
    dev.off()
    png(paste0(data_dir, "/plot2c.png"), plsz, plsz, res = 300)
    print(plotLearner.ti.combined(plotdata_task, plotdata_output2) + labs(title = "monocle prediction"))
    dev.off()
    # png(paste0(data_dir, "/plot2d.png"), plsz, plsz, res = 300)
    # print(plotLearner.ti.combined(plotdata_task, plotdata_output3) + labs(title = "TSCAN prediction"))
    # dev.off()

    # row 3
    png(paste0(data_dir, "/plot3a.png"), plsz, plsz, res = 300)
    plot.new()
    dev.off()
    png(paste0(data_dir, "/plot3b.png"), plsz, plsz, res = 300)
    print(plotLearner.ti.scorpius(pred_output1) + labs(title = "SCORPIUS own plot") + theme(legend.position = "none"))
    dev.off()
    png(paste0(data_dir, "/plot3c.png"), plsz, plsz, res = 300)
    print(plotLearner.ti.monocle(pred_output2) + coord_equal() + labs(title = "monocle own plot") + theme(legend.position = "none"))
    dev.off()
    # png(paste0(data_dir, "/plot3d.png"), plsz, plsz, res = 300)
    # plot.new()
    # dev.off()

    # row 4
    plot_emdist(task, task_emdist, plotdata_task, filename = paste0(data_dir, "/plot4a.png"), width = plsz / 300, height = plsz / 300)
    plot_emdist(pred_output1, pred_emdist1, plotdata_output1, filename = paste0(data_dir, "/plot4b.png"), width = plsz / 300, height = plsz / 300)
    plot_emdist(pred_output2, pred_emdist2, plotdata_output2, filename = paste0(data_dir, "/plot4c.png"), width = plsz / 300, height = plsz / 300)
    # plot_emdist(pred_output3, pred_emdist3, plotdata_output3, filename = paste0(data_dir, "/plot4d.png"), width = plsz / 300, height = plsz / 300)

    # row 5
    png(paste0(data_dir, "/plot5a.png"), plsz, plsz, res = 300)
    print(ggplot(scores) + geom_bar(aes(method, max_lcmc), stat = "identity"))
    dev.off()
    png(paste0(data_dir, "/plot5b.png"), plsz, plsz, res = 300)
    coRanking::imageplot(corank1$corank)
    dev.off()
    png(paste0(data_dir, "/plot5c.png"), plsz, plsz, res = 300)
    coRanking::imageplot(corank2$corank)
    dev.off()
    # png(paste0(data_dir, "/plot5d.png"), plsz, plsz, res = 300)
    # coRanking::imageplot(corank3$corank)
    # dev.off()

    # row 6
    png(paste0(data_dir, "/plot6a.png"), plsz, plsz, res = 300)
    print(ggplot(scores) + geom_bar(aes(method, cor), stat = "identity"))
    dev.off()
    png(paste0(data_dir, "/plot6b.png"), plsz, plsz, res = 300)
    (ggplot() + geom_point(aes(as.vector(task_emdist), as.vector(pred_emdist1)), alpha = .1, size = .5) + labs(x = "original EM dist", y = "SCORPIUS EM dist")) %>% print
    dev.off()
    png(paste0(data_dir, "/plot6c.png"), plsz, plsz, res = 300)
    (ggplot() + geom_point(aes(as.vector(task_emdist), as.vector(pred_emdist2)), alpha = .1, size = .5) + labs(x = "original EM dist", y = "monocle EM dist")) %>% print
    dev.off()
    # png(paste0(data_dir, "/plot6d.png"), plsz, plsz, res = 300)
    # (ggplot() + geom_point(aes(as.vector(task_emdist), as.vector(pred_emdist3)), alpha = .1, size = .5) + labs(x = "original EM dist", y = "TSCAN EM dist")) %>% print
    # dev.off()

    system(paste0(
      "cd ", data_dir, "\n",
      "bash << HERE\n",
      "convert plot1[abcd].png -append result_a.png\n",
      "convert plot2[abcd].png -append result_b.png\n",
      "convert plot3[abcd].png -append result_c.png\n",
      "convert plot4[abcd].png -append result_d.png\n",
      "convert plot5[abcd].png -append result_e.png\n",
      "convert plot6[abcd].png -append result_f.png\n",
      "convert result_[abcdef].png +append result.png\n",
      "HERE\n"))

    file.copy(paste0(data_dir, "/result.png"), paste0(data_dir, "_plot.png"), overwrite = T)
  }
}

dat_evals <- lapply(seq_len(nrow(datasets_info)), function(dataset_num) {
  dataset_id <- datasets_info$id[[dataset_num]]
  data_dir <- paste0(output_root_folder, dataset_id)
  dir.create(data_dir)
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
  geom_path(aes(max_lcmc, cor, group = dataset, colour = ti_type)) +
  geom_point(aes(max_lcmc, cor, shape = method), size = 3)
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
