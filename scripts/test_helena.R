library(dyneval)
library(tidyverse)

task <- dyneval:::read_ti_task_data("linear_ginhoux")

output <- dyneval:::trainLearner.ti.scorpius(task, .subset = NULL, num_dimensions = 3, num_clusters = 4)

output$state_network
output$state_percentages

# todo:
#  - compare task$state_network and task$state_percentages
#  - to output$state_network and output$state_percentages
