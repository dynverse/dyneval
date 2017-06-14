library(cowplot)
library(tidyverse)
library(dyneval)
library(mlrMBO)
library(parallelMap)
library(PRISM)

output_root_folder <- "results/output_dyngentest/"

## load datasets
.datasets_location = "../dyngen/results" # needs to be defined, to let dyngen know where the datasets are
tasks <- load_datasets()
tasks <- tasks %>% group_by(ti_type) %>% mutate(i = seq_len(n())) %>% ungroup

## choose a method
# method <- description_scorpius()
# method <- description_monocle_ddrtree()
method <- description_celltree_maptpx()
obj_fun <- make_obj_fun(method)

## select datasets
select_tasks <- tasks %>% arrange(ti_type, i)

## MBO settings
control_train <- makeMBOControl(
  propose.points = 24,
  impute.y.fun = impute_y_fun
) %>% setMBOControlTermination(iters = 20L)
control_eval <- makeMBOControl(
  propose.points = 2,
  impute.y.fun = impute_y_fun
) %>% setMBOControlTermination(iters = 1L)
design <- generateDesign(n = 100, par.set = method$par_set)

## Fold setting
test_set <- 1

## Run MBO
qsub_handle = qsub_run(qsub_config = override_qsub_config(wait = F), FUN = function() {
  library(dplyr)
  library(dyneval)
  library(mlrMBO)
  library(parallelMap)

  ## start parameter optimisation
  parallelStartMulticore(cpus = 24, show.info = TRUE)
  tune_train <- mbo(obj_fun, design = design, control = control_train, show.info = T, more.args = list(tasks = select_tasks %>% filter(!i %in% test_set)))
  tune_test <- mbo(obj_fun, design = tune_train$opt.path$env$path %>% select(-y), control = control_eval, show.info = T, more.args = list(tasks = select_tasks %>% filter(i %in% test_set)))
  parallelStop()

  list(tune_train = tune_train, tune_test = tune_test)
})




opt_path_test <- tune_test$opt.path$env$path
opt_path$y_test <- head(opt_path_test$y, nrow(opt_path))



individual_scores <- bind_rows(lapply(seq_len(nrow(opt_path)), function(i) {
  if (opt_path$y[[i]] != 0) {
    tune_out$opt.path$env$extra[[i]]$.summary %>% mutate(point = i, fold_type = "train")
  } else {
    NULL
  }
})) %>% left_join(tasks %>% select(type, ti_type, name, experimentid, platformname, version, i), by = c("task_name"="name"))
grouped_scores <- individual_scores %>% group_by(ti_type, point, fold_type) %>% summarise(auc_lcmc = mean(auc_lcmc), max_lcmc = mean(max_lcmc))
ggplot(grouped_scores) + geom_point(aes(point, max_lcmc, colour = ti_type)) + scale_colour_brewer(palette = "Dark2")
ggplot(grouped_scores) + geom_point(aes(point, max_lcmc, colour = ti_type)) + scale_colour_brewer(palette = "Dark2") + facet_wrap(~ti_type, nrow = 1) + geom_vline(xintercept = which(diff(opt_path$iteration) != 0)+.5)
ggplot(grouped_scores) + geom_point(aes(point, max_lcmc, colour = ti_type)) + scale_colour_brewer(palette = "Dark2") + facet_wrap(~ti_type, ncol = 1) + geom_vline(xintercept = which(diff(opt_path$iteration) != 0)+.5)

ggplot(opt_path) + geom_point(aes(i, y)) + geom_vline(xintercept = which(diff(opt_path$iteration) != 0)+.5)
ggplot(opt_path) + geom_point(aes(i, y, colour = num_topics_lower)) + geom_vline(xintercept = which(diff(opt_path$iteration) != 0)+.5) + scale_colour_distiller(palette="RdBu")
ggplot(opt_path) + geom_point(aes(i, y, colour = num_topics_upper)) + geom_vline(xintercept = which(diff(opt_path$iteration) != 0)+.5) + scale_colour_distiller(palette="RdBu")
ggplot(opt_path) + geom_point(aes(i, y, colour = sd_filter)) + geom_vline(xintercept = which(diff(opt_path$iteration) != 0)+.5) + scale_colour_distiller(palette="RdBu")
ggplot(opt_path) + geom_point(aes(i, y, colour = tot_iter)) + geom_vline(xintercept = which(diff(opt_path$iteration) != 0)+.5) + scale_colour_distiller(palette="RdBu")
ggplot(opt_path) + geom_point(aes(i, y, colour = tolerance)) + geom_vline(xintercept = which(diff(opt_path$iteration) != 0)+.5) + scale_colour_distiller(palette="RdBu")
ggplot(opt_path) + geom_point(aes(i, y, colour = width_scale_factor)) + geom_vline(xintercept = which(diff(opt_path$iteration) != 0)+.5) + scale_colour_distiller(palette="RdBu")
ggplot(opt_path) + geom_point(aes(i, y, colour = outlier_tolerance_factor)) + geom_vline(xintercept = which(diff(opt_path$iteration) != 0)+.5) + scale_colour_distiller(palette="RdBu")





individual_scores_test <- bind_rows(lapply(seq_len(nrow(opt_path)), function(i) {
  if (opt_path_test$y[[i]] != 0) {
    tune_test$opt.path$env$extra[[i]]$.summary %>% mutate(point = i, fold_type = "test")
  } else {
    NULL
  }
})) %>% left_join(tasks %>% select(type, ti_type, name, experimentid, platformname, version, i), by = c("task_name"="name"))
grouped_scores_test <- individual_scores_test %>% group_by(ti_type, point, fold_type) %>% summarise(auc_lcmc = mean(auc_lcmc), max_lcmc = mean(max_lcmc))

ggplot(opt_path) + geom_point(aes(i, y, colour = "train")) + geom_point(aes(i, y_test, colour = "test")) + geom_vline(xintercept = which(diff(opt_path$iteration) != 0)+.5)
ggplot(grouped_scores_test) + geom_point(aes(point, max_lcmc, colour = ti_type)) + scale_colour_brewer(palette = "Dark2") + facet_wrap(~ti_type, ncol = 1) + geom_vline(xintercept = which(diff(opt_path$iteration) != 0)+.5)

ggplot(bind_rows(grouped_scores, grouped_scores_test)) + geom_point(aes(point, max_lcmc, colour = fold_type)) + scale_colour_brewer(palette = "Dark2") + facet_wrap(~ti_type, ncol = 1) + geom_vline(xintercept = which(diff(opt_path$iteration) != 0)+.5)
ggplot(bind_rows(grouped_scores, grouped_scores_test)) + geom_point(aes(point, max_lcmc, colour = fold_type)) + scale_colour_brewer(palette = "Dark2") + facet_grid(ti_type~fold_type) + geom_vline(xintercept = which(diff(opt_path$iteration) != 0)+.5)
