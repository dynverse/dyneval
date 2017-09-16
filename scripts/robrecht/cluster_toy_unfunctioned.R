library(dyneval)
library(dyntoy)
library(tidyverse)

# easy test
tasks <- generate_toy_datasets(num_replicates = 2)
task_group <- rep("group", nrow(tasks))
task_fold <- gsub(".*_", "", tasks$id) %>% as.integer()
methods <- get_descriptions(as_tibble = T)

out_dir = out_dir
methods = methods
metrics = c("auc_R_nx", "robbie_network_score")
timeout = 600
memory = "16G"
num_cores = 8
num_iterations = 5
num_init_params = 16
max_wall_time = "72:00:00"
num_repeats = 1

## MBO settings
control_train <- mlrMBO::makeMBOControl(
  n.objectives = length(metrics),
  propose.points = 1,
  impute.y.fun = impute_y_fun(length(metrics))) %>%
  mlrMBO::setMBOControlTermination(iters = num_iterations) %>%
  mlrMBO::setMBOControlInfill(mlrMBO::makeMBOInfillCritDIB())
control_test <- control_train
control_test$iters <- 1
control_test$propose.points <- 1
learner <- mlr::makeLearner(
  "regr.randomForest",
  se.method = "jackknife",
  predict.type = "se",
  keep.inbag = TRUE)

## Grid settings
grid <- expand.grid(
  fold_i = sort(unique(task_fold)),
  group_sel = sort(unique(task_group)),
  repeat_i = seq_len(num_repeats),
  stringsAsFactors = F
)

## Run MBO
methodi <- 2
method <- dynutils::extract_row_to_list(methods, methodi)

cat("Submitting ", method$name, "\n", sep="")

# create an objective function
obj_fun <- make_obj_fun(method = method, metrics = metrics, timeout = timeout)

# generate initial parameters
design <- bind_rows(
  ParamHelpers::generateDesignOfDefaults(method$par_set),
  ParamHelpers::generateDesign(n = num_init_params, par.set = method$par_set)
)

# cluster params
library(mlrMBO)
library(parallelMap)

# submit to the cluster
grid_i <- 1
fold_i <- grid[grid_i,]$fold_i
group_sel <- grid[grid_i,]$group_sel
repeat_i <- grid[grid_i,]$repeat_i

## start parameter optimisation
parallelStartMulticore(cpus = num_cores, show.info = TRUE)
tune_train <- mbo(
  obj_fun,
  learner = learner,
  design = design,
  control = control_train,
  show.info = TRUE,
  more.args = list(tasks = tasks[task_group == group_sel & task_fold != fold_i,])
)
tune_test <- mbo(
  obj_fun,
  learner = learner,
  design = tune_train$opt.path$env$path %>% select(-starts_with("y_"), -one_of("y")),
  control = control_test,
  show.info = TRUE,
  more.args = list(tasks = tasks[task_group == group_sel & task_fold == fold_i,])
)
parallelStop()

