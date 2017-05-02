library(dyneval)
library(tidyverse)
library(cowplot)
set.seed(1)

task <- read_task_data("ti", "linear", "ginhoux")
task_emdist <- compute_em_dist(task)

pred_output <- trainLearner.ti.scorpius(task, .subset = NULL, num_dimensions = 3, num_clusters = 4)
pred_emdist <- compute_em_dist(pred_output)

plot_grid(
  plotLearner.ti.default(task),
  plotLearner.ti.default(pred_output),
  plotLearner.ti.scorpius(pred_output),
  nrow = 1
)


corank <- compute_coranking(task_emdist, pred_emdist)

corank$lcmc
corank$auc_lcmc

plot(corank$lcmc)
coRanking::imageplot(corank$corank)
pheatmap::pheatmap(corank$corank, cluster_rows = F, cluster_cols = F, color = colorRampPalette(c("white", "black"))(100), show_rownames = F, show_colnames = F)
