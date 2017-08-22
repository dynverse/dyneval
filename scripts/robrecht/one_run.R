library(cowplot)
library(tidyverse)
library(dyneval)
library(mlrMBO)
library(parallelMap)

output_root_folder <- "results/output_dyngentest/"

## load datasets
.datasets_location = "../dyngen/results/4/" # needs to be defined, to let dyngen know where the datasets are
tasks <- load_datasets(8, ndatasets = 10) # this function takes way too long due to the geodesic distances being calculated
# saveRDS(tasks, paste0(.datasets_location, "tasks.rds"))
#tasks <- readRDS(paste0(.datasets_location, "tasks.rds")) %>% mutate(group = paste0(platform_id, "_", takesetting_type)) %>% group_by(group, ti_type) %>% mutate(subtask_ix = seq_len(n())) %>% ungroup()

## take one task
counts <- tasks$counts[[1]]

## choose a method
# method <- description_scuba()
# method <- description_waterfall()
#method <- description_monocle_ddrtree()
method <- description_scorpius()
method <- description_slingshot()
method_fun <- make_obj_fun(method)

## apply default params
default_params <- generateDesignOfDefaults(method$par_set)
params <- generateDesign(n = 20, par.set = method$par_set)

#eval_out <- method_fun(default_params[1,] %>% as.list, tasks = tasks)
eval_out <- method_fun(list(), tasks = tasks[1:min(nrow(tasks), 10),])
eval_extras <- attr(eval_out,"extras")
attr(eval_out,"extras") <- NULL
eval_out
eval_extras$.summary

## apply it manually
method_out <- method$run_fun(counts)

corank_out <- compute_coranking(tasks$geodesic_dist[[1]], method_out$geodesic_dist)
