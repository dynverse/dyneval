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
num_cores <- 24
control_train <- makeMBOControl(propose.points = num_cores, impute.y.fun = impute_y_fun) %>% setMBOControlTermination(iters = 20L)
control_eval <- makeMBOControl(propose.points = 2, impute.y.fun = impute_y_fun) %>% setMBOControlTermination(iters = 1L)

# ## Run MBO
# qsub_handle <- qsub_lapply(
#   X = 1:4,
#   qsub_environment = c("method", "obj_fun", "select_tasks", "num_cores", "control_train", "control_eval", "design", "test_set"),
#   qsub_config = override_qsub_config(wait = F, num_cores = num_cores, memory = "10G", name = "dyneval", remove_tmp_folder = F, max_wall_time = "24:00:00"),
#   FUN = function(test_set) {
#   library(dplyr)
#   library(purrr)
#   library(dyneval)
#   library(mlrMBO)
#   library(parallelMap)
#
#   design <- bind_rows(
#     generateDesignOfDefaults(method$par_set),
#     generateDesign(n = num_cores * 4, par.set = method$par_set)
#   )
#   ## start parameter optimisation
#   parallelStartMulticore(cpus = num_cores, show.info = TRUE)
#   tune_train <- mbo(obj_fun, design = design, control = control_train, show.info = T, more.args = list(tasks = select_tasks %>% filter(!i %in% test_set)))
#   tune_test <- mbo(obj_fun, design = tune_train$opt.path$env$path %>% select(-y), control = control_eval, show.info = T, more.args = list(tasks = select_tasks %>% filter(i %in% test_set)))
#   parallelStop()
#
#   list(design = design, tune_train = tune_train, tune_test = tune_test)
# })
#
# save.image("results/temp.RData")
load("results/temp.RData")

# default_design <- generateDesignOfDefaults(method$par_set)
# many_dd <- bind_rows(lapply(1:4, function(x) default_design))
#
# parallelStartMulticore(cpus = 8, show.info = TRUE)
# default_train <- mbo(obj_fun, design = many_dd, control = control_eval, show.info = T, more.args = list(tasks = select_tasks %>% filter(!i %in% test_set)))
# default_test <- mbo(obj_fun, design = many_dd, control = control_eval, show.info = T, more.args = list(tasks = select_tasks %>% filter(i %in% test_set)))
# parallelStop()

output <- qsub_retrieve(qsub_handle)

gathered <- lapply(seq_len(4), function(repeat_i) {
  tune_train <- output[[repeat_i]]$tune_train
  tune_test <- output[[repeat_i]]$tune_test

  opt_path <- tune_train$opt.path$env$path %>%
    mutate(iteration = tune_train$opt.path$env$dob, time = tune_train$opt.path$env$exec.time) %>%
    rename(y_train = y) %>%
    left_join(tune_test$opt.path$env$path %>% rename(y_test = y), by = setdiff(colnames(tune_train$opt.path$env$path), "y")) %>%
    mutate(repeat_i = repeat_i, run_i = seq_len(n()))

  individual_scores <- bind_rows(lapply(seq_len(nrow(opt_path)), function(i) {
    iteration <- opt_path$iteration[[i]]
    bind_rows(
      if (opt_path$y_train[[i]] != 0) {
        tune_train$opt.path$env$extra[[i]]$.summary %>% mutate(repeat_i = repeat_i, run_i = i, iteration, fold_type = "train")
      } else {
        NULL
      },
      if (opt_path$y_test[[i]] != 0) {
        tune_test$opt.path$env$extra[[i]]$.summary %>% mutate(repeat_i = repeat_i, run_i = i, iteration, fold_type = "test")
      } else {
        NULL
      }
    )
  })) %>% left_join(tasks %>% select(type, ti_type, name, experimentid, platformname, version, i), by = c("task_name"="name"))
  grouped_scores <- individual_scores %>% group_by(iteration, repeat_i, ti_type, run_i, fold_type) %>% summarise(auc_lcmc = mean(auc_lcmc), max_lcmc = mean(max_lcmc)) %>% ungroup()

  list(opt_path = opt_path, individual_scores = individual_scores, grouped_scores = grouped_scores)
})

opt_path <- gathered %>% map_df(~ .$opt_path)
individual_scores <- gathered %>% map_df(~ .$individual_scores)
grouped_scores <- gathered %>% map_df(~ .$grouped_scores)

png("rplot.png", 1000, 400)
ggplot(opt_path %>% filter(y_train > 0, y_test > 0)) + geom_point(aes(y_train, y_test, colour = iteration)) + scale_colour_distiller(palette = "RdBu") + facet_wrap(~repeat_i, nrow = 1) + coord_equal()
dev.off()
png("rplot2.png", 1600, 1200)
ggplot(grouped_scores %>% select(-auc_lcmc) %>% spread(fold_type, max_lcmc)) + geom_point(aes(train, test, colour = iteration)) + scale_colour_distiller(palette = "RdBu") + facet_grid(repeat_i~ti_type) + coord_equal()
dev.off()

best_opt_path <- opt_path %>% arrange(desc(y_train)) %>% slice(1)
best_param <- best_opt_path %>% select(-y_train:-iteration)
best_individual_scores <- individual_scores %>% filter(run_i == best_opt_path$run_i, repeat_i == best_opt_path$repeat_i)
best_grouped_scores <- grouped_scores %>% filter(run_i == best_opt_path$run_i, repeat_i == best_opt_path$repeat_i)

# ggplot(best_individual_scores) + geom_bar(aes(factor(i), max_lcmc, fill = ti_type), stat = "identity") + facet_wrap(~ti_type)
png("rplot3.png", 1000, 300)
ggplot(best_grouped_scores) + geom_bar(aes(ti_type, max_lcmc, fill = factor(fold_type, levels = c("train", "test"))), stat = "identity", position = "dodge") + labs(fill = "Fold type")
dev.off()

g <- ggplot() +
  geom_bar(aes(ti_type, max_lcmc, group = factor(fold_type, levels = c("train", "test")), colour = "best"), fill = NA, stat = "identity", position = "dodge", best_grouped_scores) +
  geom_boxplot(aes(ti_type, max_lcmc, fill = factor(fold_type, levels = c("train", "test"))), position = "dodge", grouped_scores) +
  labs(fill = "Fold type") +
  scale_fill_brewer(palette = "Set2") +
  scale_colour_brewer(palette = "Set1")
g

png("rplot4.png", 1000, 600)
g
dev.off()

ggplot(individual_scores) + geom_point(aes(run_i, max_lcmc, colour = fold_type)) + facet_grid(i ~ ti_type)

# ggplot(grouped_scores) + geom_point(aes(point, max_lcmc, colour = ti_type)) + scale_colour_brewer(palette = "Dark2")
# ggplot(grouped_scores) + geom_point(aes(point, max_lcmc, colour = ti_type)) + scale_colour_brewer(palette = "Dark2") + facet_wrap(~ti_type, nrow = 1) + geom_vline(xintercept = which(diff(opt_path$iteration) != 0)+.5)
# ggplot(grouped_scores) + geom_point(aes(point, max_lcmc, colour = ti_type)) + scale_colour_brewer(palette = "Dark2") + facet_wrap(~ti_type, ncol = 1) + geom_vline(xintercept = which(diff(opt_path$iteration) != 0)+.5)
#
ggplot(opt_path) + geom_point(aes(i, y_train)) + geom_vline(xintercept = which(diff(opt_path$iteration) != 0)+.5)
ggplot(opt_path) + geom_point(aes(i, y_test)) + geom_vline(xintercept = which(diff(opt_path$iteration) != 0)+.5)
ggplot(opt_path) + geom_point(aes(y_train, y_test))
ggplot(opt_path %>% filter(y_train > 0)) + geom_point(aes(y_train, y_test, colour = iteration)) + scale_colour_distiller(palette = "RdBu")
ggplot(grouped_scores %>% select(-auc_lcmc) %>% spread(fold_type, max_lcmc)) + geom_point(aes(train, test, colour = iteration)) + scale_colour_distiller(palette = "RdBu") + facet_wrap(~ti_type)
# ggplot(opt_path) + geom_point(aes(i, y, colour = num_topics_lower)) + geom_vline(xintercept = which(diff(opt_path$iteration) != 0)+.5) + scale_colour_distiller(palette="RdBu")
# ggplot(opt_path) + geom_point(aes(i, y, colour = num_topics_upper)) + geom_vline(xintercept = which(diff(opt_path$iteration) != 0)+.5) + scale_colour_distiller(palette="RdBu")
# ggplot(opt_path) + geom_point(aes(i, y, colour = sd_filter)) + geom_vline(xintercept = which(diff(opt_path$iteration) != 0)+.5) + scale_colour_distiller(palette="RdBu")
# ggplot(opt_path) + geom_point(aes(i, y, colour = tot_iter)) + geom_vline(xintercept = which(diff(opt_path$iteration) != 0)+.5) + scale_colour_distiller(palette="RdBu")
# ggplot(opt_path) + geom_point(aes(i, y, colour = tolerance)) + geom_vline(xintercept = which(diff(opt_path$iteration) != 0)+.5) + scale_colour_distiller(palette="RdBu")
# ggplot(opt_path) + geom_point(aes(i, y, colour = width_scale_factor)) + geom_vline(xintercept = which(diff(opt_path$iteration) != 0)+.5) + scale_colour_distiller(palette="RdBu")
# ggplot(opt_path) + geom_point(aes(i, y, colour = outlier_tolerance_factor)) + geom_vline(xintercept = which(diff(opt_path$iteration) != 0)+.5) + scale_colour_distiller(palette="RdBu")

ggplot(opt_path) + geom_point(aes(i, y, colour = "train")) + geom_point(aes(i, y_test, colour = "test")) + geom_vline(xintercept = which(diff(opt_path$iteration) != 0)+.5)
ggplot(grouped_scores_test) + geom_point(aes(point, max_lcmc, colour = ti_type)) + scale_colour_brewer(palette = "Dark2") + facet_wrap(~ti_type, ncol = 1) + geom_vline(xintercept = which(diff(opt_path$iteration) != 0)+.5)

ggplot(grouped_scores) + geom_point(aes(point, max_lcmc, colour = fold_type)) + scale_colour_brewer(palette = "Dark2") + facet_wrap(~ti_type, ncol = 1) + geom_vline(xintercept = which(diff(opt_path$iteration) != 0)+.5)
ggplot(grouped_scores) + geom_point(aes(point, max_lcmc, colour = fold_type)) + scale_colour_brewer(palette = "Dark2") + facet_grid(ti_type~fold_type) + geom_vline(xintercept = which(diff(opt_path$iteration) != 0)+.5)
