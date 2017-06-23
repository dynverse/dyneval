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
# method <- description_scorpius()
# method <- description_ddrtree()
method <- description_celltree_maptpx()
method_fun <- make_obj_fun(method)

## MBO settings
control_train <- makeMBOControl(propose.points = 8, impute.y.fun = impute_y_fun) %>% setMBOControlTermination(iters = 10L)
control_test <- makeMBOControl(propose.points = 1, impute.y.fun = impute_y_fun) %>% setMBOControlTermination(iters = 1L)
design <- generateDesign(n = 8, par.set = method$par_set)

## start parameter optimisation
parallelStartMulticore(cpus = 8, show.info = TRUE)
train_out <- mbo(
  fun = method_fun,
  design = design,
  control = control_train,
  show.info = T,
  more.args = list(tasks = tasks %>% filter(ti_type == "linear", subdataset_i != 4)) # use training datasets
)
test_out <- mbo(
  fun = method_fun,
  design = train_out$opt.path$env$path %>% dplyr::select(-y), # evaluate the parameters that have been evaluated in the training
  control = control_test,
  show.info = T,
  more.args = list(tasks = tasks %>% filter(ti_type == "linear", subdataset_i == 4)) # use test datasets
)
parallelStop()

## look at results
train_out

## combine the scores for the different parameters that have been tried, for train and test
train_path <- train_out$opt.path$env$path %>% dplyr::rename(train = y)
test_path <- test_out$opt.path$env$path %>% dplyr::rename(test = y)

eval_summ <- train_path %>%
  left_join(test_path, by = setdiff(colnames(train_path), "train")) %>%
  mutate(
    iteration_i = train_out$opt.path$env$dob,
    param_i = seq_len(n()),
    time_train = train_out$opt.path$env$exec.time,
    time_test = test_out$opt.path$env$exec.time[param_i]
  )
eval_summ_gath <- eval_summ %>% gather(fold_type, y, train, test)

## collect the scores per dataset individually
eval_ind <- bind_rows(lapply(seq_len(nrow(eval_summ)), function(param_i) {
  iteration_i <- eval_summ$iteration[[param_i]]
  bind_rows(
    if (eval_summ$y_train[[param_i]] >= 0) { # did this execution finish correctly?
      train_out$opt.path$env$extra[[param_i]]$.summary %>% mutate(param_i, iteration_i, fold_type = "train")
    } else {
      NULL
    },
    if (eval_summ$y_test[[param_i]] >= 0) { # did this execution finish correctly?
      test_out$opt.path$env$extra[[param_i]]$.summary %>% mutate(param_i, iteration_i, fold_type = "test")
    } else {
      NULL
    }
  )
})) %>% left_join(tasks %>% dplyr::select(type, ti_type, name, experimentid, platformname, version, dataset_i, subdataset_i), by = c("task_name" = "name"))

## group them together per ti_type
eval_grp <- eval_ind %>% group_by(ti_type, iteration_i, param_i, fold_type) %>% summarise(auc_lcmc = mean(auc_lcmc), max_lcmc = mean(max_lcmc)) %>% ungroup()

## Make several plots
ggplot(eval_summ_gath %>% filter(y >= 0)) + geom_point(aes(param_i, y, colour = fold_type)) + facet_wrap(~fold_type, ncol = 1)
ggplot(eval_summ %>% filter(train >= 0, test >= 0)) + geom_point(aes(train, test))
