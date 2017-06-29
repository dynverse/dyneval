library(cowplot)
library(tidyverse)
library(dyneval)
library(mlrMBO)
library(parallelMap)
library(PRISM)

output_root_folder <- "results/output_dyngen_paramtraincv/"
dir.create(output_root_folder)

## load datasets
.datasets_location = "../dyngen/results" # needs to be defined, to let dyngen know where the datasets are
tasks <- load_datasets() %>% mutate(dataset_i = seq_len(n())) %>% group_by(ti_type) %>% mutate(subdataset_i = seq_len(n())) %>% ungroup

## choose a method
method <- description_scorpius()
# method <- description_monocle_ddrtree()
# method <- description_celltree_maptpx()
# method <- description_celltree_gibbs()
# method <- description_celltree_vem()

## select datasets
select_tasks <- tasks %>% arrange(ti_type, subdataset_i)

## MBO settings
num_cores <- 6

## set up evaluation
metrics <- c("Q_global", "Q_local", "correlation")
obj_fun <- make_obj_fun(method, metrics = metrics)
impute_fun <- impute_y_fun(obj_fun)

## MBO settings
control_train <- makeMBOControl(
  n.objectives = length(metrics),
  propose.points = num_cores,
  impute.y.fun = impute_fun) %>%
  setMBOControlTermination(iters = 20L) %>%
  setMBOControlInfill(makeMBOInfillCritDIB())
control_test <- control_train
control_test$iters <- 1
control_test$propose.points <- 1

grid <- expand.grid(repeat_i = seq_len(4), fold_i = seq_len(4))

## Run MBO
qsub_handle <- qsub_lapply(
  X = seq_len(nrow(grid)),
  qsub_environment = c("method", "select_tasks", "num_cores", "metrics", "obj_fun", "impute_fun", "control_train", "control_test", "design", "grid"),
  qsub_config = override_qsub_config(wait = F, num_cores = num_cores, memory = "20G", name = "dyneval", remove_tmp_folder = F, max_wall_time = "24:00:00"),
  FUN = function(grid_i) {

  library(dplyr)
  library(purrr)
  library(dyneval)
  library(mlrMBO)
  library(parallelMap)

  repeat_i <- grid$repeat_i[[grid_i]]
  fold_i <- grid$fold_i[[grid_i]]

  design <- bind_rows(
    generateDesignOfDefaults(method$par_set),
    generateDesign(n = 100, par.set = method$par_set)
  )
  ## start parameter optimisation
  parallelStartMulticore(cpus = num_cores, show.info = TRUE)
  tune_train <- mbo(
    obj_fun,
    design = design,
    control = control_train,
    show.info = T,
    more.args = list(tasks = select_tasks %>% filter(!subdataset_i %in% fold_i))
  )
  tune_test <- mbo(
    obj_fun,
    design = tune_train$opt.path$env$path %>% dplyr::select(-starts_with("y_"), -one_of("y")),
    control = control_test,
    show.info = T,
    more.args = list(tasks = select_tasks %>% filter(subdataset_i %in% fold_i))
  )
  parallelStop()

  list(design = design, tune_train = tune_train, tune_test = tune_test)
})

save.image(paste0(output_root_folder, "temp.RData"))
load(paste0(output_root_folder, "temp.RData"))

output <- qsub_retrieve(qsub_handle)
#
# save(output, file = paste0(output_root_folder, "temp_output.RData"))
load(paste0(output_root_folder, "temp_output.RData"))


gathered <- lapply(seq_len(nrow(grid)), function(grid_i) {
  fold_i <- grid$fold_i[[grid_i]]
  repeat_i <- grid$repeat_i[[grid_i]]
  design <- output[[grid_i]]$design
  tune_train <- output[[grid_i]]$tune_train
  tune_test <- output[[grid_i]]$tune_test

  scores <- tune_train$opt.path$env$path %>%
    rename(y_train = y) %>%
    left_join(tune_test$opt.path$env$path %>% rename(y_test = y), by = setdiff(colnames(tune_train$opt.path$env$path), "y")) %>%
    mutate(
      iteration_i = tune_train$opt.path$env$dob,
      time = tune_train$opt.path$env$exec.time,
      grid_i, repeat_i, fold_i, param_i = seq_len(n()))

  individual_scores <- bind_rows(lapply(seq_len(nrow(scores)), function(param_i) {
    iteration_i <- scores$iteration[[param_i]]
    bind_rows(
      if (scores$y_train[[param_i]] >= 0) {
        tune_train$opt.path$env$extra[[param_i]]$.summary %>% mutate(grid_i, repeat_i, fold_i, param_i, iteration_i, fold_type = "train")
      } else {
        NULL
      },
      if (scores$y_test[[param_i]] >= 0) {
        tune_test$opt.path$env$extra[[param_i]]$.summary %>% mutate(grid_i, repeat_i, fold_i, param_i, iteration_i, fold_type = "test")
      } else {
        NULL
      }
    )
  })) %>% left_join(tasks %>% select(type, ti_type, name, experimentid, platformname, version, dataset_i = i), by = c("task_name" = "name"))
  grouped_scores <- individual_scores %>% group_by(grid_i, repeat_i, fold_i, ti_type, iteration_i, param_i, fold_type) %>% summarise(auc_lcmc = mean(auc_lcmc), max_lcmc = mean(max_lcmc)) %>% ungroup()

  list(scores = scores, individual_scores = individual_scores, grouped_scores = grouped_scores)
})

scores <- gathered %>% map_df(~ .$scores) %>% mutate(grid_str = paste0("Repeat ", repeat_i, ", fold ", fold_i))
individual_scores <- gathered %>% map_df(~ .$individual_scores) %>% mutate(grid_str = paste0("Repeat ", repeat_i, ", fold ", fold_i))
grouped_scores <- gathered %>% map_df(~ .$grouped_scores) %>% mutate(grid_str = paste0("Repeat ", repeat_i, ", fold ", fold_i))

best_scores <- scores %>% group_by(fold_i) %>%  arrange(desc(y_train)) %>% slice(1)
best_param <- best_scores %>% select(-y_train:-param_i)
best_individual_scores <- individual_scores %>% filter(param_i == best_scores$param_i, grid_i == best_scores$grid_i)
best_grouped_scores <- grouped_scores %>% filter(param_i == best_scores$param_i, grid_i == best_scores$grid_i)

pdf(paste0(output_root_folder, "1_celltree_overallscore.pdf"), 10, 4)
# png(paste0(output_root_folder, "1_celltree_overallscore.png"), 1000, 400)
ggplot(scores %>% filter(y_train >= 0, y_test >= 0)) +
  geom_vline(aes(xintercept = y_train), scores %>% filter(param_i == 1)) +
  geom_hline(aes(yintercept = y_test), scores %>% filter(param_i == 1)) +
  geom_vline(aes(xintercept = y_train), best_scores, colour = "red") +
  geom_hline(aes(yintercept = y_test), best_scores, colour = "red") +
  geom_point(aes(y_train, y_test, colour = iteration_i)) +
  scale_colour_distiller(palette = "RdBu") +
  facet_wrap(~grid_str, nrow = 1) +
  coord_equal()
dev.off()

pdf(paste0(output_root_folder, "2_celltree_groupedscore.pdf"), 16, 12)
# png(paste0(output_root_folder, "2_celltree_groupedscore.png"), 1600, 1200)
ggplot(grouped_scores %>% select(-auc_lcmc) %>% spread(fold_type, max_lcmc)) +
  geom_vline(aes(xintercept = train), grouped_scores %>% select(-auc_lcmc) %>% spread(fold_type, max_lcmc) %>% filter(param_i == 1)) +
  geom_hline(aes(yintercept = test), grouped_scores %>% select(-auc_lcmc) %>% spread(fold_type, max_lcmc) %>% filter(param_i == 1)) +
  geom_vline(aes(xintercept = train), best_grouped_scores %>% select(-auc_lcmc) %>% spread(fold_type, max_lcmc), colour = "red") +
  geom_hline(aes(yintercept = test), best_grouped_scores %>% select(-auc_lcmc) %>% spread(fold_type, max_lcmc), colour = "red") +
  geom_point(aes(train, test, colour = iteration_i)) +
  scale_colour_distiller(palette = "RdBu") +
  facet_grid(grid_str~ti_type) +
  coord_equal()
dev.off()



# ggplot(best_individual_scores) + geom_bar(aes(factor(i), max_lcmc, fill = ti_type), stat = "identity") + facet_wrap(~ti_type)
# png("rplot3.png", 1000, 300)
# ggplot(best_grouped_scores) + geom_bar(aes(ti_type, max_lcmc, fill = factor(fold_type, levels = c("train", "test"))), stat = "identity", position = "dodge") + labs(fill = "Fold type")
# dev.off()
ggplot() +
  geom_bar(aes(ti_type, max_lcmc, group = factor(fold_type, levels = c("train", "test")), colour = "best"), fill = NA, stat = "identity", position = "dodge", best_grouped_scores) +
  geom_boxplot(aes(ti_type, max_lcmc, fill = factor(fold_type, levels = c("train", "test"))), position = "dodge", grouped_scores) +
  labs(fill = "Fold type") +
  scale_fill_brewer(palette = "Set2") +
  scale_colour_brewer(palette = "Set1")


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
# ggplot(grouped_scores) + geom_point(aes(point, max_lcmc, colour = ti_type)) + scale_colour_brewer(palette = "Dark2") + facet_wrap(~ti_type, nrow = 1) + geom_vline(xintercept = which(diff(scores$iteration) != 0)+.5)
# ggplot(grouped_scores) + geom_point(aes(point, max_lcmc, colour = ti_type)) + scale_colour_brewer(palette = "Dark2") + facet_wrap(~ti_type, ncol = 1) + geom_vline(xintercept = which(diff(scores$iteration) != 0)+.5)
#
ggplot(scores) + geom_point(aes(i, y_train)) + geom_vline(xintercept = which(diff(scores$iteration) != 0)+.5)
ggplot(scores) + geom_point(aes(i, y_test)) + geom_vline(xintercept = which(diff(scores$iteration) != 0)+.5)
ggplot(scores) + geom_point(aes(y_train, y_test))
ggplot(scores %>% filter(y_train > 0)) + geom_point(aes(y_train, y_test, colour = iteration)) + scale_colour_distiller(palette = "RdBu")
ggplot(grouped_scores %>% select(-auc_lcmc) %>% spread(fold_type, max_lcmc)) + geom_point(aes(train, test, colour = iteration)) + scale_colour_distiller(palette = "RdBu") + facet_wrap(~ti_type)
# ggplot(scores) + geom_point(aes(i, y, colour = num_topics_lower)) + geom_vline(xintercept = which(diff(scores$iteration) != 0)+.5) + scale_colour_distiller(palette="RdBu")
# ggplot(scores) + geom_point(aes(i, y, colour = num_topics_upper)) + geom_vline(xintercept = which(diff(scores$iteration) != 0)+.5) + scale_colour_distiller(palette="RdBu")
# ggplot(scores) + geom_point(aes(i, y, colour = sd_filter)) + geom_vline(xintercept = which(diff(scores$iteration) != 0)+.5) + scale_colour_distiller(palette="RdBu")
# ggplot(scores) + geom_point(aes(i, y, colour = tot_iter)) + geom_vline(xintercept = which(diff(scores$iteration) != 0)+.5) + scale_colour_distiller(palette="RdBu")
# ggplot(scores) + geom_point(aes(i, y, colour = tolerance)) + geom_vline(xintercept = which(diff(scores$iteration) != 0)+.5) + scale_colour_distiller(palette="RdBu")
# ggplot(scores) + geom_point(aes(i, y, colour = width_scale_factor)) + geom_vline(xintercept = which(diff(scores$iteration) != 0)+.5) + scale_colour_distiller(palette="RdBu")
# ggplot(scores) + geom_point(aes(i, y, colour = outlier_tolerance_factor)) + geom_vline(xintercept = which(diff(scores$iteration) != 0)+.5) + scale_colour_distiller(palette="RdBu")

ggplot(scores) + geom_point(aes(i, y, colour = "train")) + geom_point(aes(i, y_test, colour = "test")) + geom_vline(xintercept = which(diff(scores$iteration) != 0)+.5)
ggplot(grouped_scores_test) + geom_point(aes(point, max_lcmc, colour = ti_type)) + scale_colour_brewer(palette = "Dark2") + facet_wrap(~ti_type, ncol = 1) + geom_vline(xintercept = which(diff(scores$iteration) != 0)+.5)

ggplot(grouped_scores) + geom_point(aes(point, max_lcmc, colour = fold_type)) + scale_colour_brewer(palette = "Dark2") + facet_wrap(~ti_type, ncol = 1) + geom_vline(xintercept = which(diff(scores$iteration) != 0)+.5)
ggplot(grouped_scores) + geom_point(aes(point, max_lcmc, colour = fold_type)) + scale_colour_brewer(palette = "Dark2") + facet_grid(ti_type~fold_type) + geom_vline(xintercept = which(diff(scores$iteration) != 0)+.5)
