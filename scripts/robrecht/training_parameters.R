library(tidyverse)
library(dyneval)
library(mlrMBO)
library(parallelMap)

output_root_folder <- "results/output_dyngentest/"

## load datasets
.datasets_location = "../dyngen/results" # needs to be defined, to let dyngen know where the datasets are
tasks <- load_datasets(mc_cores = 8)

## choose a method
method <- description_scorpius()
# method <- description_monocle_ddrtree()
# method <- description_celltree_maptpx()
# method <- description_celltree_gibbs()
# method <- description_celltree_vem()
# method <- description_tscan()
# method <- description_stemid()
# method <- description_scuba()
# method <- description_dpt()
# method <- description_embeddr()

## MBO settings
num_cores <- 8

## set up evaluation
metrics <- c("Q_global", "Q_local", "correlation")
obj_fun <- make_obj_fun(method, metrics = metrics, suppress_output = T)
impute_fun <- impute_y_fun(length(metrics))

# try first
try <- obj_fun(list(), tasks = tasks)
try_extra <- attr(try, "extras")
attr(try, "extras") <- NULL
try


## MBO settings
control_train <- makeMBOControl(
  n.objectives = length(metrics),
  propose.points = num_cores,
  impute.y.fun = impute_fun) %>%
  setMBOControlTermination(iters = 2) %>%
  setMBOControlInfill(makeMBOInfillCritDIB())
control_test <- control_train
control_test$iters <- 1
control_test$propose.points <- 1
design <- bind_rows(
  generateDesignOfDefaults(method$par_set),
  generateDesign(n = 4 * num_cores, par.set = method$par_set)
)

## start parameter optimisation
parallelStartMulticore(cpus = 8, show.info = TRUE)
train_out <- mbo(
  fun = obj_fun,
  design = design,
  control = control_train,
  show.info = T,
  more.args = list(tasks = tasks %>% filter(ti_type == "linear", subdataset_i != 4)) # use training datasets
)
test_out <- mbo(
  fun = obj_fun,
  design = train_out$opt.path$env$path %>% dplyr::select(-starts_with("y_"), -one_of("y")), # evaluate the parameters that have been evaluated in the training
  control = control_test,
  show.info = T,
  more.args = list(tasks = tasks %>% filter(ti_type == "linear", subdataset_i == 4)) # use test datasets
)
parallelStop()

library(cowplot)

## look at results
train_out

## combine the scores for the different parameters that have been tried, for train and test
eval_summ_gath <- bind_rows(
  data.frame(type = "train", train_out$opt.path$env$path) %>% mutate(
    param_i = seq_len(n()),
    time = train_out$opt.path$env$exec.time
  ),
  data.frame(type = "test", test_out$opt.path$env$path) %>% mutate(
    param_i = seq_len(n()),
    time = test_out$opt.path$env$exec.time
  )
) %>% filter(param_i <= nrow(train_out$opt.path$env$path)) %>%
  dplyr::as_data_frame()

eval_summ <- eval_summ_gath %>% gather(eval_metric, score, y_1:y_3, time) %>%
  mutate(comb = paste0(type, "_", eval_metric)) %>%
  select(-type, -eval_metric) %>%
  spread(comb, score) %>%
  mutate(iteration_i = train_out$opt.path$env$dob[param_i]) %>%
  arrange(param_i)

## collect the scores per dataset individually
eval_ind <- bind_rows(lapply(seq_len(nrow(eval_summ)), function(param_i) {
  iteration_i <- eval_summ$iteration_i[[param_i]]
  bind_rows(
    if (eval_summ$train_y_1[[param_i]] >= 0) { # did this execution finish correctly?
      train_out$opt.path$env$extra[[param_i]]$.summary %>% mutate(param_i, iteration_i, fold_type = "train")
    } else {
      NULL
    },
    if (eval_summ$test_y_1[[param_i]] >= 0) { # did this execution finish correctly?
      test_out$opt.path$env$extra[[param_i]]$.summary %>% mutate(param_i, iteration_i, fold_type = "test")
    } else {
      NULL
    }
  )
})) %>% left_join(tasks %>% dplyr::select(type, ti_type, name, experimentid, platformname, version, dataset_i, subdataset_i), by = c("task_name" = "name"))

## group them together per ti_type
eval_grp <- eval_ind %>% group_by(ti_type, iteration_i, param_i, fold_type) %>% summarise_at(metrics, mean) %>% ungroup()

## Make several plots
ggplot(eval_summ_gath %>% filter(y_1 >= 0)) + geom_point(aes(param_i, y_1, colour = type))
ggplot(eval_summ %>% filter(train_y_1 >= 0, test_y_1 >= 0)) + geom_point(aes(train_y_1, test_y_1))
ggplot(eval_grp) + geom_point(aes(param_i, Q_global, colour = fold_type)) + facet_wrap(~ti_type)
