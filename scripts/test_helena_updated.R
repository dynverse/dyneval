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

# todo:
#  - compare task$state_network and task$state_percentages
#  - to output$state_network and output$state_percentages
