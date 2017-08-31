library(cowplot)
library(tidyverse)
library(dyneval)
library(mlrMBO)
library(parallelMap)

tasks <- generate_toy_datasets()
counts <- tasks$counts[[1]]

## choose a method
# method <- description_scuba()
# method <- description_waterfall()
#method <- description_monocle_ddrtree()
# method <- description_scorpius()
# method <- description_slingshot()
method <- description_scuba()
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
