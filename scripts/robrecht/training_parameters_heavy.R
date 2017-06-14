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
method <- description_scorpius()
# method <- description_monocle_ddrtree()

## MBO settings
control <- makeMBOControl(propose.points = 8) %>% setMBOControlTermination(iters = 20L)
design <- generateDesign(n = 40, par.set = method$par_set)

## start parameter optimisation
parallelStartMulticore(cpus = 8, show.info = TRUE)
tune_out <- mbo(make_obj_fun(method), design = design, control = control, show.info = T, more.args = list(tasks = tasks %>% filter(ti_type == "linear")))
parallelStop()

## look at results
tune_out

## look at all of the parameters that have been tried
tune_out$opt.path$env$path

path <- tune_out$opt.path$env$path %>% mutate(i = seq_len(n())) %>% as_data_frame()

ggplot(path) + geom_point(aes(i, y))
ggplot(path) + geom_point(aes(i, y, colour = num_dimensions)) + scale_colour_distiller(palette = "RdBu")
ggplot(path) + geom_point(aes(i, y, colour = norm_method))
ggplot(path) + geom_point(aes(i, y, colour = maxIter)) + scale_colour_distiller(palette = "RdBu")
ggplot(path) + geom_point(aes(i, y, colour = sigma)) + scale_colour_distiller(palette = "RdBu")
ggplot(path) + geom_point(aes(i, y, colour = lambda)) + scale_colour_distiller(palette = "RdBu")
ggplot(path) + geom_point(aes(i, y, colour = ncenter)) + scale_colour_distiller(palette = "RdBu")
ggplot(path) + geom_point(aes(i, y, colour = param.gamma)) + scale_colour_distiller(palette = "RdBu")
ggplot(path) + geom_point(aes(i, y, colour = tol)) + scale_colour_distiller(palette = "RdBu")
ggplot(path) + geom_point(aes(i, y, colour = factor(auto_param_selection)))


path
