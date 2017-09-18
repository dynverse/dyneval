library(dyneval)
library(dyntoy)
library(tidyverse)

out_dir <- "~/Workspace/dynverse/dynresults/cluster_toy/"
dir.create(out_dir, recursive = T)

# easy test
tasks <- generate_toy_datasets(num_replicates = 2)
task_group <- rep("group", nrow(tasks))
task_fold <- gsub(".*_", "", tasks$id) %>% as.integer()
methods <- get_descriptions(as_tibble = T)

benchmark_suite_submit(
  tasks,
  task_group,
  task_fold,
  out_dir = out_dir,
  methods = methods,
  metrics = c("auc_R_nx", "robbie_network_score"),
  timeout = 600,
  memory = "16G",
  num_cores = 2,
  num_iterations = 5,
  num_init_params = 16
)

benchmark_suite_retrieve(out_dir)
