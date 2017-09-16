library(dyneval)
library(dyntoy)
library(tidyverse)

out_dir <- "~/Workspace/dynverse/dynresults/cluster_tasks/"
dir.create(out_dir, recursive = T)

# load datasets
# .datasets_location = "~/Workspace/dynverse/dynresults/dyngen_v4/4/"
# tasks <- load_datasets(8) # this function takes way too long due to the geodesic distances being calculated
# saveRDS(tasks, paste0(out_dir, "tasks.rds"))
tasks <- readRDS(paste0(out_dir, "tasks.rds")) %>%
  filter(
    platform_id == "fluidigm_c1",
    takesetting_type == "snapshot",
    ti_type == "consecutive_bifurcating",
    model_replicate %in% c(1,2))


task_group <- rep("group", nrow(tasks))
task_fold <- tasks$model_replicate
methods <- get_descriptions(as_tibble = T)

benchmark_suite_submit(
  tasks,
  task_group,
  task_fold,
  out_dir = out_dir,
  methods = methods,
  metrics = c("auc_R_nx", "robbie_network_score"),
  timeout = 1200,
  memory = "16G",
  num_cores = 2,
  num_iterations = 10,
  num_init_params = 1100
)

benchmark_suite_retrieve(out_dir)
