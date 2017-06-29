library(cowplot)
library(tidyverse)
library(dyneval)
library(mlrMBO)
library(parallelMap)

output_root_folder <- "results/output_dyngentest/"

## load datasets
.datasets_location = "../dyngen/results" # needs to be defined, to let dyngen know where the datasets are
tasks <- load_datasets() %>% mutate(dataset_i = seq_len(n())) %>% group_by(ti_type) %>% mutate(subdataset_i = seq_len(n())) %>% ungroup

## choose a method
method <- description_scorpius()
method_fun <- make_obj_fun(method)

## apply it
method_out <- method$run_fun(tasks$counts[[1]])

corank_out <- compute_coranking(tasks$geodesic_dist[[1]], method_out$geodesic_dist)
