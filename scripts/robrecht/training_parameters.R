library(cowplot)
library(tidyverse)
library(dyneval)
library(mlrMBO)
library(parallelMap)

output_root_folder <- "results/output_dyngentest/"

## load datasets
.datasets_location = "../dyngen/results" # needs to be defined, to let dyngen know where the datasets are
tasks <- load_datasets()

## choose a method
# method <- description_scorpius()
method <- description_monocle_ddrtree()

## MBO settings
control <- makeMBOControl(propose.points = 8) %>% setMBOControlTermination(iters = 2L)
design <- generateDesign(n = 8, par.set = method$par_set)

## start parameter optimisation
parallelStartMulticore(cpus = 8, show.info = TRUE)
tune_out <- mbo(make_obj_fun(method), design = design, control = control, show.info = T, more.args = list(tasks = tasks %>% filter(ti_type == "linear")))
parallelStop()

## look at results
tune_out

## look at all of the parameters that have been tried
tune_out$opt.path$env$path
