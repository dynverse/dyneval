library(dyneval)
library(tidyverse)
library(cowplot)
set.seed(1)

task <- read_task_data("ti", "linear", "ginhoux")

plotLearner.ti.default(task)

scorpius_output <- trainLearner.ti.scorpius(task, .subset = NULL, num_dimensions = 3, num_clusters = 4)

plotLearner.ti.default(scorpius_output)
plotLearner.ti.scorpius(scorpius_output)

scorpius_output$state_names
scorpius_output$state_network
scorpius_output$state_percentages

plot_grid(
  plotLearner.ti.default(task),
  plotLearner.ti.default(scorpius_output),
  plotLearner.ti.scorpius(scorpius_output),
  nrow = 1
)

plot_emdist2 <- function(traj, emdist) {
  state_percentages <- traj$state_percentages
  pct <- as.data.frame(state_percentages[,-1])
  rownames(pct) <- state_percentages$id
  dist_mat <- reshape2::acast(emdist, from~to, value.var = "dist")
  pheatmap::pheatmap(dist_mat, annotation_col = pct, annotation_row = pct)
}

task_emdist <- emdist(task)
plot_emdist2(task, task_emdist)
scorpius_emdist <- emdist(scorpius_output)
plot_emdist2(scorpius_output, scorpius_emdist)

joined <- task_emdist %>% rename(gold = dist) %>%
  left_join(scorpius_emdist %>% rename(pred = dist), by = c("from", "to")) %>%
  as_data_frame

cor(joined$gold, joined$pred)
qplot(joined$gold, joined$pred)
